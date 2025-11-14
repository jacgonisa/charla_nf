#!/usr/bin/env python3
"""
Extract large indels (>threshold bp) from BAM files and extract their sequences.

This script:
1. Reads indel position TSV files from non_hybrid_reads_analysis
2. Filters for indels larger than the specified threshold (default 50bp)
3. Extracts the indel sequences from the original reads
4. Creates a catalog of large indels with metadata
5. Outputs indel sequences in FASTA format for downstream analysis
"""

import os
import sys
import pysam
import csv
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd


def load_fasta_sequences(fasta_file):
    """Load all sequences from FASTA file into a dictionary."""
    print(f"\n[INFO] Loading sequences from {fasta_file}...")
    sequences = {}
    count = 0
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)
        count += 1
        if count % 10000 == 0:
            print(f"  Loaded {count} sequences so far...")
    print(f"[SUCCESS] Loaded {len(sequences)} total sequences")
    return sequences


def extract_large_indels_from_tsv(tsv_file, threshold=50):
    """
    Extract large indels from the indel positions TSV file.
    Returns a list of indel records with metadata.
    """
    large_indels = []

    if not os.path.exists(tsv_file):
        print(f"  [WARNING] TSV file not found: {tsv_file}")
        return large_indels

    total_indels = 0
    with open(tsv_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            total_indels += 1
            size = int(row['Size'])
            if size >= threshold:
                large_indels.append({
                    'read_id': row['Read_ID'],
                    'type': row['Type'],
                    'size': size,
                    'read_pos': int(row['Read_Pos']),
                    'ref_pos': int(row['Ref_Pos']),
                    'source_file': os.path.basename(tsv_file)
                })

    print(f"  [INFO] Processed {total_indels} total indels, found {len(large_indels)} large indels (>={threshold}bp)")
    return large_indels


def extract_indel_sequence_from_bam(bam_file, read_id, read_pos, indel_type, indel_size):
    """
    Extract the actual indel sequence from the BAM file.
    For insertions, extract the inserted sequence.
    For deletions, we'll mark it as a deletion (no sequence to extract from read).
    """
    try:
        bamfile = pysam.AlignmentFile(bam_file, "rb")
    except Exception as e:
        print(f"  [ERROR] Opening BAM file {bam_file}: {e}")
        return None

    for read in bamfile.fetch(until_eof=True):
        if read.query_name == read_id:
            if indel_type == "Insertion":
                # Extract the inserted sequence
                start_pos = read_pos
                end_pos = read_pos + indel_size
                if end_pos <= len(read.query_sequence):
                    sequence = read.query_sequence[start_pos:end_pos]
                    bamfile.close()
                    return sequence
                else:
                    print(f"  [WARNING] Insertion coordinates out of bounds for read {read_id}")
            else:  # Deletion
                # For deletions, we don't have the deleted sequence in the read
                # We could potentially extract it from the reference, but for now
                # we'll just mark it as a deletion
                bamfile.close()
                return None

    bamfile.close()
    return None


def process_nonhybrid_directory(nonhybrid_dir, fasta_sequences, threshold=50):
    """
    Process all indel TSV files in the non_hybrid directory.
    Extract large indels and their sequences.
    """
    mappings_dir = os.path.join(nonhybrid_dir, "mappings")

    if not os.path.exists(mappings_dir):
        print(f"[ERROR] Mappings directory not found: {mappings_dir}")
        return [], []

    print(f"\n[INFO] Searching for indel TSV files in: {mappings_dir}")

    large_indels = []
    indel_sequences = []

    # Find all indel position TSV files
    all_tsv_files = list(Path(mappings_dir).glob("*_indel_positions.tsv"))

    # Filter to only analyze:
    # 1. mm_ont and wm_mapont modes
    # 2. Centromeric haplotypes only (CC* and LC*, not CA* or LA*)
    tsv_files = []
    for f in all_tsv_files:
        # Check for mapping mode
        if not ('mm_ont' in f.name or 'wm_mapont' in f.name):
            continue

        # Extract haplotype code (e.g., "CC1_mm_ont_indel_positions.tsv" -> "CC1")
        basename = f.stem.replace("_indel_positions", "")
        haplotype = basename.split('_')[0]

        # Check if it's a centromeric haplotype (second character is 'C')
        if len(haplotype) >= 2 and haplotype[1] == 'C':
            tsv_files.append(f)

    print(f"[INFO] Found {len(all_tsv_files)} total indel position TSV files")
    print(f"[INFO] Filtering for: mm_ont and wm_mapont modes only")
    print(f"[INFO] Filtering for: Centromeric haplotypes only (CC*, LC*)")
    print(f"[INFO] Analyzing {len(tsv_files)} files after filtering")

    if len(tsv_files) == 0:
        print(f"[WARNING] No TSV files to process!")
        return [], []

    for file_idx, tsv_file in enumerate(tsv_files, 1):
        print(f"\n[PROGRESS] Processing file {file_idx}/{len(tsv_files)}: {tsv_file.name}")

        # Extract haplotype from filename (e.g., "CA1_against_own_indel_positions.tsv" -> "CA1")
        basename = tsv_file.stem.replace("_indel_positions", "")
        haplotype = basename.split('_')[0]
        print(f"  [INFO] Haplotype: {haplotype}")

        # Find corresponding BAM file
        bam_file = tsv_file.with_suffix('.bam').with_name(basename + '.bam')

        if not bam_file.exists():
            print(f"  [WARNING] BAM file not found: {bam_file}")
            print(f"  [WARNING] Skipping this file...")
            continue

        print(f"  [INFO] BAM file found: {bam_file.name}")

        # Extract large indels from TSV
        indels = extract_large_indels_from_tsv(str(tsv_file), threshold)

        if len(indels) == 0:
            print(f"  [INFO] No large indels found in this file, moving to next...")
            continue

        # For each large indel, extract the sequence
        insertions_count = 0
        deletions_count = 0
        sequences_extracted = 0

        print(f"  [INFO] Extracting sequences for {len(indels)} large indels...")

        for idx, indel in enumerate(indels):
            read_id = indel['read_id']

            # Try to get sequence from BAM first (for insertions)
            if indel['type'] == "Insertion":
                insertions_count += 1
                sequence = extract_indel_sequence_from_bam(
                    str(bam_file),
                    read_id,
                    indel['read_pos'],
                    indel['type'],
                    indel['size']
                )
                if sequence:
                    sequences_extracted += 1
            else:
                deletions_count += 1
                sequence = None

            # Create unique indel ID
            indel_id = f"{haplotype}_{indel['type'][:3]}_{indel['size']}bp_{idx+1}"

            # Add metadata to indel record
            indel_record = {
                'indel_id': indel_id,
                'haplotype': haplotype,
                'read_id': read_id,
                'type': indel['type'],
                'size': indel['size'],
                'read_pos': indel['read_pos'],
                'ref_pos': indel['ref_pos'],
                'multiple_of_178': (indel['size'] % 178 == 0),
                'closest_178_multiple': round(indel['size'] / 178, 2),
                'source_bam': os.path.basename(bam_file)
            }
            large_indels.append(indel_record)

            # Create sequence record if we have a sequence
            if sequence:
                seq_record = SeqRecord(
                    Seq(sequence),
                    id=indel_id,
                    description=f"type={indel['type']} size={indel['size']}bp haplotype={haplotype} read={read_id}"
                )
                indel_sequences.append(seq_record)

        print(f"  [SUMMARY] {haplotype}: {insertions_count} insertions, {deletions_count} deletions")
        print(f"  [SUMMARY] {haplotype}: {sequences_extracted} sequences extracted")

    print(f"\n[SUCCESS] Finished processing all files!")
    return large_indels, indel_sequences


def main():
    if len(sys.argv) != 5:
        print("Usage: python3 13.1-extract_large_indels.py <nonhybrid_dir> <fasta_file> <output_dir> <threshold>")
        sys.exit(1)

    nonhybrid_dir = sys.argv[1]
    fasta_file = sys.argv[2]
    output_dir = sys.argv[3]
    threshold = int(sys.argv[4])

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    print("\n" + "="*80)
    print("=== EXTRACTING LARGE INDELS ===")
    print("="*80)
    print(f"[CONFIG] Threshold: {threshold} bp")
    print(f"[CONFIG] Input directory: {nonhybrid_dir}")
    print(f"[CONFIG] Output directory: {output_dir}")
    print(f"[CONFIG] FASTA file: {fasta_file}")
    print("="*80)

    # Load FASTA sequences
    fasta_sequences = load_fasta_sequences(fasta_file)

    # Process non_hybrid directory
    large_indels, indel_sequences = process_nonhybrid_directory(
        nonhybrid_dir,
        fasta_sequences,
        threshold
    )

    # Save catalog
    print("\n" + "="*80)
    print("=== SAVING RESULTS ===")
    print("="*80)

    catalog_file = os.path.join(output_dir, "large_indels_catalog.tsv")
    if large_indels:
        print(f"[INFO] Saving catalog with {len(large_indels)} large indels...")
        df = pd.DataFrame(large_indels)
        df.to_csv(catalog_file, sep='\t', index=False)
        print(f"[SUCCESS] Catalog saved to: {catalog_file}")

        # Print summary statistics
        print("\n" + "="*80)
        print("=== SUMMARY STATISTICS ===")
        print("="*80)
        print(f"Total large indels: {len(large_indels)}")
        insertions = sum(1 for x in large_indels if x['type'] == 'Insertion')
        deletions = sum(1 for x in large_indels if x['type'] == 'Deletion')
        print(f"  Insertions: {insertions} ({insertions/len(large_indels)*100:.1f}%)")
        print(f"  Deletions: {deletions} ({deletions/len(large_indels)*100:.1f}%)")

        exact_multiples = sum(1 for x in large_indels if x['multiple_of_178'])
        print(f"\nExact 178bp multiples: {exact_multiples} ({exact_multiples/len(large_indels)*100:.1f}%)")

        sizes = [x['size'] for x in large_indels]
        print(f"\nSize distribution:")
        print(f"  Min: {min(sizes)} bp")
        print(f"  Max: {max(sizes)} bp")
        print(f"  Mean: {sum(sizes)/len(sizes):.1f} bp")
        print(f"  Median: {sorted(sizes)[len(sizes)//2]} bp")

        # Show size distribution by multiples of 178
        print(f"\nIndels by 178bp multiples:")
        for mult in range(1, 6):
            count = sum(1 for s in sizes if abs(s - mult*178) <= 20)
            if count > 0:
                print(f"  ~{mult}×178bp ({mult*178}bp ±20): {count} indels")
    else:
        print("[WARNING] No large indels found.")
        # Create empty catalog
        pd.DataFrame(columns=['indel_id', 'haplotype', 'read_id', 'type', 'size',
                             'read_pos', 'ref_pos', 'multiple_of_178',
                             'closest_178_multiple', 'source_bam']).to_csv(
                                 catalog_file, sep='\t', index=False)

    # Save sequences
    sequences_file = os.path.join(output_dir, "indel_sequences.fa")
    print(f"\n[INFO] Saving indel sequences...")
    if indel_sequences:
        SeqIO.write(indel_sequences, sequences_file, "fasta")
        print(f"[SUCCESS] Saved {len(indel_sequences)} indel sequences to: {sequences_file}")
    else:
        print("[INFO] No indel sequences to save (deletions have no sequence to extract).")
        # Create empty file
        open(sequences_file, 'w').close()

    print("\n" + "="*80)
    print("=== EXTRACTION COMPLETE ===")
    print("="*80)


if __name__ == "__main__":
    main()
