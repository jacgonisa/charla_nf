#!/usr/bin/env python3
"""
Analyze indel sequences for 178bp satellite block structure.

This script:
1. Analyzes if indels are multiples of 178bp
2. Breaks indel sequences into 178bp blocks
3. Performs MSA (Multiple Sequence Alignment) on blocks
4. Checks if blocks align to known centromeric satellite repeats
5. Identifies if indels are intact satellite blocks or broken/chimeric
6. Generates reports and visualizations
"""

import os
import sys
import subprocess
from pathlib import Path
import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns


SATELLITE_LENGTH = 178  # Arabidopsis centromeric satellite repeat length


def split_sequence_into_blocks(sequence, block_size=178):
    """
    Split a sequence into blocks of specified size.
    Returns list of (block_index, block_sequence) tuples.
    """
    blocks = []
    seq_len = len(sequence)
    num_complete_blocks = seq_len // block_size

    for i in range(num_complete_blocks):
        start = i * block_size
        end = start + block_size
        blocks.append((i, sequence[start:end]))

    # Handle remainder
    remainder_start = num_complete_blocks * block_size
    if remainder_start < seq_len:
        remainder = sequence[remainder_start:]
        if len(remainder) >= block_size * 0.5:  # Only keep if at least 50% of block size
            blocks.append((num_complete_blocks, remainder))

    return blocks


def analyze_indel_structure(indel_sequences_file, catalog_file, output_dir):
    """
    Analyze structure of each indel sequence.
    """
    # Load catalog
    catalog_df = pd.read_csv(catalog_file, sep='\t')

    # Load indel sequences
    indel_seqs = {record.id: str(record.seq) for record in SeqIO.parse(indel_sequences_file, 'fasta')}

    print(f"Analyzing {len(indel_seqs)} indel sequences...")

    analysis_results = []
    all_blocks = []

    for indel_id, sequence in indel_seqs.items():
        # Get indel info from catalog
        indel_info = catalog_df[catalog_df['indel_id'] == indel_id]
        if indel_info.empty:
            continue

        indel_size = indel_info.iloc[0]['size']
        indel_type = indel_info.iloc[0]['type']

        # Split into 178bp blocks
        blocks = split_sequence_into_blocks(sequence, SATELLITE_LENGTH)

        # Analyze block structure
        num_blocks = len(blocks)
        complete_blocks = sum(1 for _, block in blocks if len(block) == SATELLITE_LENGTH)
        partial_blocks = num_blocks - complete_blocks

        # Calculate expected blocks
        expected_blocks = indel_size / SATELLITE_LENGTH

        # Check if size is multiple of 178
        is_exact_multiple = (indel_size % SATELLITE_LENGTH == 0)
        deviation_from_multiple = indel_size % SATELLITE_LENGTH

        # Store block sequences for MSA
        for block_idx, block_seq in blocks:
            if len(block_seq) >= SATELLITE_LENGTH * 0.8:  # At least 80% of expected length
                all_blocks.append(SeqRecord(
                    Seq(block_seq),
                    id=f"{indel_id}_block{block_idx}",
                    description=f"indel={indel_id} block={block_idx} length={len(block_seq)}"
                ))

        # Analyze sequence composition
        gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence) if len(sequence) > 0 else 0

        # Calculate sequence complexity (simple entropy-like measure)
        base_counts = Counter(sequence)
        complexity = len(base_counts) / 4.0  # Normalized by 4 possible bases

        analysis_results.append({
            'indel_id': indel_id,
            'indel_size': indel_size,
            'indel_type': indel_type,
            'is_exact_multiple': is_exact_multiple,
            'expected_blocks': round(expected_blocks, 2),
            'actual_blocks': num_blocks,
            'complete_blocks': complete_blocks,
            'partial_blocks': partial_blocks,
            'deviation_bp': deviation_from_multiple,
            'gc_content': round(gc_content, 4),
            'sequence_complexity': round(complexity, 4),
            'block_completeness': round(complete_blocks / expected_blocks, 2) if expected_blocks > 0 else 0
        })

    return analysis_results, all_blocks


def run_msa_clustal(sequences, output_file, threads=4):
    """
    Run multiple sequence alignment using Clustal Omega.
    Returns True if successful, False otherwise.
    """
    if len(sequences) < 2:
        print("Not enough sequences for MSA (need at least 2)")
        return False

    # Write sequences to temporary file
    temp_fasta = output_file + ".temp.fa"
    SeqIO.write(sequences, temp_fasta, "fasta")

    # Run Clustal Omega
    cmd = [
        'clustalo',
        '-i', temp_fasta,
        '-o', output_file,
        '--outfmt=clustal',
        '--threads', str(threads),
        '--force'
    ]

    print(f"Running MSA with Clustal Omega: {len(sequences)} sequences...")

    try:
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        print(f"MSA complete: {output_file}")
        os.remove(temp_fasta)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running Clustal Omega: {e}")
        print(f"stderr: {e.stderr}")
        return False
    except FileNotFoundError:
        print("Warning: Clustal Omega not found. Trying MAFFT...")
        return run_msa_mafft(sequences, output_file, threads)


def run_msa_mafft(sequences, output_file, threads=4):
    """
    Run multiple sequence alignment using MAFFT as fallback.
    """
    if len(sequences) < 2:
        print("Not enough sequences for MSA (need at least 2)")
        return False

    # Write sequences to temporary file
    temp_fasta = output_file + ".temp.fa"
    SeqIO.write(sequences, temp_fasta, "fasta")

    # Run MAFFT
    cmd = [
        'mafft',
        '--auto',
        '--thread', str(threads),
        temp_fasta
    ]

    print(f"Running MSA with MAFFT: {len(sequences)} sequences...")

    try:
        with open(output_file, 'w') as outfile:
            result = subprocess.run(
                cmd,
                stdout=outfile,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
        print(f"MSA complete: {output_file}")
        os.remove(temp_fasta)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running MAFFT: {e}")
        print(f"stderr: {e.stderr}")
        return False
    except FileNotFoundError:
        print("Error: Neither Clustal Omega nor MAFFT found. Skipping MSA.")
        return False


def analyze_alignment_quality(alignment_file):
    """
    Analyze quality and conservation of MSA alignment.
    """
    if not os.path.exists(alignment_file):
        return None

    try:
        # Try reading as clustal format first
        try:
            alignment = AlignIO.read(alignment_file, "clustal")
        except:
            # Try FASTA format if clustal fails
            alignment = AlignIO.read(alignment_file, "fasta")

        # Calculate alignment statistics
        num_sequences = len(alignment)
        alignment_length = alignment.get_alignment_length()

        # Calculate conservation score (simple identity-based)
        conservation_scores = []
        for i in range(alignment_length):
            column = alignment[:, i]
            most_common = Counter(column).most_common(1)[0]
            conservation = most_common[1] / num_sequences
            conservation_scores.append(conservation)

        avg_conservation = np.mean(conservation_scores)
        highly_conserved_positions = sum(1 for x in conservation_scores if x >= 0.8)

        return {
            'num_sequences': num_sequences,
            'alignment_length': alignment_length,
            'avg_conservation': round(avg_conservation, 4),
            'highly_conserved_positions': highly_conserved_positions,
            'conservation_percentage': round(highly_conserved_positions / alignment_length, 4) if alignment_length > 0 else 0
        }
    except Exception as e:
        print(f"Error analyzing alignment: {e}")
        return None


def create_satellite_analysis_plots(analysis_df, output_dir):
    """
    Create comprehensive plots for satellite structure analysis.
    """
    plt.style.use('default')
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))

    # Plot 1: Size distribution with 178bp multiples
    ax = axes[0, 0]
    sizes = analysis_df['indel_size'].values
    ax.hist(sizes, bins=30, alpha=0.7, color='steelblue', edgecolor='black')

    # Add vertical lines for 178bp multiples
    max_size = max(sizes) if len(sizes) > 0 else 178
    for i in range(1, int(max_size / 178) + 2):
        ax.axvline(i * 178, color='red', linestyle='--', alpha=0.5, linewidth=1)
        ax.text(i * 178, ax.get_ylim()[1] * 0.95, f'{i}×178',
               rotation=90, va='top', ha='right', fontsize=8)

    ax.set_xlabel('Indel Size (bp)', fontweight='bold')
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_title('Indel Size Distribution with 178bp Multiples', fontweight='bold')

    # Plot 2: Exact multiples vs non-multiples
    ax = axes[0, 1]
    exact_counts = analysis_df['is_exact_multiple'].value_counts()
    colors = ['#2ecc71' if x else '#e74c3c' for x in exact_counts.index]
    ax.bar(['Exact\nMultiple', 'Not Exact\nMultiple'], exact_counts.values,
           color=colors, alpha=0.7, edgecolor='black')
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_title('Exact 178bp Multiples', fontweight='bold')

    # Plot 3: Block completeness
    ax = axes[0, 2]
    block_completeness = analysis_df['block_completeness'].values
    ax.hist(block_completeness, bins=20, alpha=0.7, color='purple', edgecolor='black')
    ax.axvline(1.0, color='red', linestyle='--', linewidth=2, label='100% complete')
    ax.set_xlabel('Block Completeness Ratio', fontweight='bold')
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_title('Satellite Block Completeness', fontweight='bold')
    ax.legend()

    # Plot 4: Expected vs actual blocks
    ax = axes[1, 0]
    ax.scatter(analysis_df['expected_blocks'], analysis_df['actual_blocks'],
              alpha=0.6, s=50, color='steelblue')
    max_val = max(analysis_df['expected_blocks'].max(), analysis_df['actual_blocks'].max())
    ax.plot([0, max_val], [0, max_val], 'r--', linewidth=2, label='Perfect agreement')
    ax.set_xlabel('Expected Blocks', fontweight='bold')
    ax.set_ylabel('Actual Blocks', fontweight='bold')
    ax.set_title('Expected vs Actual Block Count', fontweight='bold')
    ax.legend()

    # Plot 5: GC content distribution
    ax = axes[1, 1]
    gc_content = analysis_df['gc_content'].values
    ax.hist(gc_content, bins=20, alpha=0.7, color='green', edgecolor='black')
    ax.axvline(0.5, color='red', linestyle='--', linewidth=2, label='50% GC')
    ax.set_xlabel('GC Content', fontweight='bold')
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_title('GC Content Distribution', fontweight='bold')
    ax.legend()

    # Plot 6: Deviation from 178bp multiples
    ax = axes[1, 2]
    deviations = analysis_df['deviation_bp'].values
    ax.hist(deviations, bins=30, alpha=0.7, color='orange', edgecolor='black')
    ax.set_xlabel('Deviation from 178bp Multiple (bp)', fontweight='bold')
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_title('Size Deviation from Exact Multiples', fontweight='bold')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'satellite_structure_analysis.png'),
               dpi=300, bbox_inches='tight')
    plt.close()

    print("Generated satellite structure analysis plots")


def main():
    if len(sys.argv) != 6:
        print("Usage: python3 13.3-analyze_satellite_structure.py <indel_sequences.fa> <catalog.tsv> <msa_output_dir> <output_analysis.tsv> <threads>")
        sys.exit(1)

    indel_sequences_file = sys.argv[1]
    catalog_file = sys.argv[2]
    msa_output_dir = sys.argv[3]
    output_analysis_file = sys.argv[4]
    threads = int(sys.argv[5])

    print("=== Analyzing Satellite Structure ===")
    print(f"Indel sequences: {indel_sequences_file}")
    print(f"Catalog: {catalog_file}")
    print(f"MSA output directory: {msa_output_dir}")

    # Check if files exist
    if not os.path.exists(indel_sequences_file) or os.path.getsize(indel_sequences_file) == 0:
        print("No indel sequences to analyze. Creating empty results.")
        pd.DataFrame(columns=['indel_id', 'indel_size', 'is_exact_multiple',
                             'expected_blocks', 'actual_blocks', 'alignment_quality']).to_csv(
                                 output_analysis_file, sep='\t', index=False)
        return

    # Analyze indel structure
    print("\n=== Step 1: Analyzing indel structure ===")
    analysis_results, all_blocks = analyze_indel_structure(
        indel_sequences_file,
        catalog_file,
        msa_output_dir
    )

    # Save blocks for MSA (will be done by MAFFT module in Nextflow)
    print("\n=== Step 2: Saving 178bp blocks for MSA ===")
    blocks_file = os.path.join(msa_output_dir, "satellite_blocks.fa")

    if len(all_blocks) >= 2:
        SeqIO.write(all_blocks, blocks_file, "fasta")
        print(f"[SUCCESS] Saved {len(all_blocks)} blocks to: {blocks_file}")
        print(f"[INFO] MSA will be performed by MAFFT module in Nextflow pipeline")

        # Initialize alignment fields (will be updated after MAFFT runs)
        for result in analysis_results:
            result['alignment_included'] = False
            result['avg_block_conservation'] = 0
    else:
        print(f"[WARNING] Not enough blocks for MSA (found {len(all_blocks)}, need at least 2)")
        # Create empty file
        open(blocks_file, 'w').close()

        for result in analysis_results:
            result['alignment_included'] = False
            result['avg_block_conservation'] = 0

    # Save analysis results
    analysis_df = pd.DataFrame(analysis_results)
    analysis_df.to_csv(output_analysis_file, sep='\t', index=False)
    print(f"\n=== Step 3: Saving results ===")
    print(f"Saved satellite analysis to: {output_analysis_file}")

    # Generate summary statistics
    print("\n=== Summary Statistics ===")
    print(f"Total indels analyzed: {len(analysis_results)}")
    exact_multiples = sum(1 for x in analysis_results if x['is_exact_multiple'])
    print(f"Exact 178bp multiples: {exact_multiples} ({exact_multiples/len(analysis_results)*100:.1f}%)")

    complete_satellites = sum(1 for x in analysis_results if x['block_completeness'] >= 0.95)
    print(f"Nearly complete satellites (≥95% blocks): {complete_satellites}")

    avg_completeness = np.mean([x['block_completeness'] for x in analysis_results])
    print(f"Average block completeness: {avg_completeness:.2%}")

    # Create plots
    print("\n=== Step 4: Creating plots ===")
    create_satellite_analysis_plots(analysis_df, os.path.dirname(output_analysis_file))

    print("\n=== Analysis complete ===")


if __name__ == "__main__":
    main()
