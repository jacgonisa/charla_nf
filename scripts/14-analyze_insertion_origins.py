#!/usr/bin/env python3
"""
Analyze Insertion Origins Using Cenhapmer Profiling

This script extracts large insertion sequences (>50bp) from BAM files and profiles
them against cenhapmer databases to determine their genomic origin. This can reveal:

1. Homologous recombination-mediated expansions
2. Ectopic insertions from different chromosomes/regions
3. Centromeric satellite repeat expansions
4. Evidence of recombination between non-homologous regions

The key insight: If an insertion in chromosome 1 contains k-mers specific to
chromosome 2 centromere, this suggests recombination-mediated duplication!

Usage:
    python 14-analyze_insertion_origins.py \\
        --bam-dir <mapping_directory> \\
        --cenhapmer-dir <cenhapmer_database_dir> \\
        --kmer-size <k> \\
        --output-prefix <output_prefix>

Example:
    python scripts/14-analyze_insertion_origins.py \\
        --bam-dir results/non_hybrid_mapping/sample_nonhybrid/mappings \\
        --cenhapmer-dir 03-cenhapmers/k31 \\
        --kmer-size 31 \\
        --output-prefix results/insertion_origins/analysis
"""

import os
import sys
import argparse
import pysam
import subprocess
import tempfile
from pathlib import Path
from collections import defaultdict, Counter
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
import warnings
warnings.filterwarnings('ignore')

# Publication-ready plot settings
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans']
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


def parse_haplotype(haplotype_code):
    """Parse haplotype code (e.g., CA1 -> Col, ARMS, Chr1)"""
    if len(haplotype_code) != 3:
        return None

    background = 'Col-0' if haplotype_code[0] == 'C' else 'Ler-0'
    region = 'ARMS' if haplotype_code[1] == 'A' else 'CEN'
    chromosome = haplotype_code[2]

    return {
        'background': background,
        'region': region,
        'chromosome': chromosome,
        'code': haplotype_code
    }


def extract_insertions_from_bam(bam_path, min_size=50):
    """
    Extract insertion sequences from BAM file.

    Returns:
        List of dictionaries with insertion details:
        - read_id: Read identifier
        - ref_chr: Reference chromosome
        - ref_pos: Reference position
        - size: Insertion size
        - sequence: Inserted sequence
        - context: 50bp surrounding context
    """

    if not os.path.exists(bam_path):
        print(f"Warning: BAM file not found: {bam_path}")
        return []

    try:
        bamfile = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        print(f"Error opening BAM file {bam_path}: {e}")
        return []

    insertions = []

    for read in bamfile.fetch(until_eof=True):
        if read.is_unmapped or read.cigartuples is None:
            continue

        query_pos = 0
        ref_pos = read.reference_start
        ref_chr = read.reference_name

        for op, length in read.cigartuples:
            if op == 0:  # Match
                query_pos += length
                ref_pos += length
            elif op == 1 and length >= min_size:  # Insertion >= min_size
                # Extract inserted sequence
                inserted_seq = read.query_sequence[query_pos:query_pos + length]

                # Get surrounding context (50bp before and after)
                context_start = max(0, query_pos - 50)
                context_end = min(len(read.query_sequence), query_pos + length + 50)
                context_seq = read.query_sequence[context_start:context_end]

                insertions.append({
                    'read_id': read.query_name,
                    'ref_chr': ref_chr if ref_chr else 'unknown',
                    'ref_pos': ref_pos,
                    'size': length,
                    'sequence': inserted_seq,
                    'context': context_seq,
                    'query_start': query_pos,
                    'query_end': query_pos + length
                })

                query_pos += length
            elif op == 2:  # Deletion
                ref_pos += length
            elif op in [4, 5]:  # Soft/hard clip
                query_pos += length
            elif op == 3:  # Skip
                ref_pos += length

    bamfile.close()

    return insertions


def profile_sequence_with_cenhapmers(sequence, cenhapmer_dir, kmer_size, temp_dir):
    """
    Profile a sequence against all cenhapmer databases.

    Returns:
        Dictionary mapping haplotype codes to k-mer counts
    """

    # Write sequence to temporary FASTA
    temp_fasta = os.path.join(temp_dir, f"seq_{np.random.randint(10000)}.fa")
    with open(temp_fasta, 'w') as f:
        f.write(f">insertion\n{sequence}\n")

    # Find all cenhapmer databases
    cenhapmer_files = list(Path(cenhapmer_dir).glob("*.kmc_pre"))

    if not cenhapmer_files:
        print(f"Warning: No cenhapmer databases found in {cenhapmer_dir}")
        return {}

    profiles = {}

    for kmc_pre in cenhapmer_files:
        kmc_prefix = str(kmc_pre).replace('.kmc_pre', '')
        basename = os.path.basename(kmc_prefix)

        # Parse haplotype from filename
        # Format: unique_Parent_Region_ChrN_kSize or Parent_RegionChr
        if basename.startswith('unique_'):
            parts = basename.split('_')
            if len(parts) >= 4:
                parent = parts[1]
                region_type = parts[2]
                chr_part = parts[3]
                chr_num = chr_part.replace('Chr', '')

                # Create haplotype code
                parent_code = 'C' if 'Col' in parent else 'L'
                region_code = 'A' if region_type == 'ARMS' else 'C'
                haplotype = f"{parent_code}{region_code}{chr_num}"
        else:
            parts = basename.split('_')
            if len(parts) >= 2:
                parent = parts[0]
                region_chr = parts[1]
                haplotype = region_chr[:3] if len(region_chr) >= 3 else 'UNK'
            else:
                continue

        # Count k-mers using kmc_tools
        try:
            cmd = f"kmc_tools transform {kmc_prefix} dump /dev/stdout | grep -c '^'"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=10)

            # Alternative: use intersect to count shared k-mers
            temp_kmc_prefix = os.path.join(temp_dir, f"temp_kmc_{np.random.randint(10000)}")

            # Generate k-mers from insertion sequence
            kmc_cmd = f"kmc -k{kmer_size} -m1 -fm {temp_fasta} {temp_kmc_prefix} {temp_dir} 2>/dev/null"
            subprocess.run(kmc_cmd, shell=True, capture_output=True, timeout=10)

            # Intersect with cenhapmer database
            intersect_prefix = os.path.join(temp_dir, f"intersect_{np.random.randint(10000)}")
            intersect_cmd = f"kmc_tools simple {temp_kmc_prefix} {kmc_prefix} intersect {intersect_prefix} 2>/dev/null"
            subprocess.run(intersect_cmd, shell=True, capture_output=True, timeout=10)

            # Count intersected k-mers
            count_cmd = f"kmc_tools transform {intersect_prefix} dump /dev/stdout 2>/dev/null | wc -l"
            result = subprocess.run(count_cmd, shell=True, capture_output=True, text=True, timeout=10)

            count = int(result.stdout.strip()) if result.stdout.strip() else 0
            profiles[haplotype] = count

            # Cleanup temp KMC files
            for ext in ['.kmc_pre', '.kmc_suf']:
                for prefix in [temp_kmc_prefix, intersect_prefix]:
                    try:
                        os.remove(f"{prefix}{ext}")
                    except:
                        pass

        except Exception as e:
            print(f"Warning: Error profiling against {haplotype}: {e}")
            profiles[haplotype] = 0

    # Cleanup temp fasta
    try:
        os.remove(temp_fasta)
    except:
        pass

    return profiles


def analyze_bam_insertions(bam_path, cenhapmer_dir, kmer_size, min_size=50):
    """
    Analyze all insertions from a BAM file and profile their origins.
    """

    print(f"\nAnalyzing insertions from: {os.path.basename(bam_path)}")

    # Extract insertions
    insertions = extract_insertions_from_bam(bam_path, min_size=min_size)

    if not insertions:
        print(f"  No insertions >={min_size}bp found")
        return []

    print(f"  Found {len(insertions)} insertions >={min_size}bp")

    # Create temp directory for KMC operations
    with tempfile.TemporaryDirectory() as temp_dir:
        results = []

        for i, insertion in enumerate(insertions, 1):
            if i % 10 == 0:
                print(f"  Profiling insertion {i}/{len(insertions)}...")

            # Profile insertion sequence
            profile = profile_sequence_with_cenhapmers(
                insertion['sequence'],
                cenhapmer_dir,
                kmer_size,
                temp_dir
            )

            # Determine likely origin (highest k-mer count)
            if profile:
                sorted_profile = sorted(profile.items(), key=lambda x: x[1], reverse=True)
                top_hit = sorted_profile[0] if sorted_profile[0][1] > 0 else ('UNKNOWN', 0)

                # Calculate enrichment scores
                total_kmers = sum(profile.values())
                enrichments = {k: (v / total_kmers * 100) if total_kmers > 0 else 0
                              for k, v in profile.items()}
            else:
                top_hit = ('UNKNOWN', 0)
                enrichments = {}

            results.append({
                'read_id': insertion['read_id'],
                'ref_chr': insertion['ref_chr'],
                'ref_pos': insertion['ref_pos'],
                'size': insertion['size'],
                'top_origin': top_hit[0],
                'top_origin_kmers': top_hit[1],
                'profile': profile,
                'enrichments': enrichments,
                'sequence': insertion['sequence']
            })

    return results


def detect_ectopic_insertions(results, mapped_haplotype):
    """
    Detect insertions where the origin differs from the mapped location.

    These represent potential ectopic recombination events or expansions.
    """

    ectopic = []

    hap_info = parse_haplotype(mapped_haplotype)
    if not hap_info:
        return ectopic

    expected_code = hap_info['code']

    for result in results:
        origin = result['top_origin']

        if origin == 'UNKNOWN' or result['top_origin_kmers'] == 0:
            continue

        # Check if origin differs from expected location
        if origin != expected_code:
            origin_info = parse_haplotype(origin)

            if origin_info:
                ectopic.append({
                    'insertion': result,
                    'mapped_location': hap_info,
                    'origin_location': origin_info,
                    'is_different_chr': origin_info['chromosome'] != hap_info['chromosome'],
                    'is_different_region': origin_info['region'] != hap_info['region'],
                    'is_different_background': origin_info['background'] != hap_info['background']
                })

    return ectopic


def create_origin_sankey_plot(all_results, output_prefix):
    """
    Create Sankey-like flow diagram showing insertion origins.
    """

    fig, ax = plt.subplots(figsize=(16, 12))

    # Aggregate data
    flow_data = defaultdict(lambda: defaultdict(int))

    for bam_name, results in all_results.items():
        haplotype = bam_name.split('_')[0]
        hap_info = parse_haplotype(haplotype)

        if not hap_info:
            continue

        mapped_loc = f"{hap_info['background']}\nChr{hap_info['chromosome']}\n{hap_info['region']}"

        for result in results:
            origin = result['top_origin']
            if origin == 'UNKNOWN' or result['top_origin_kmers'] == 0:
                continue

            origin_info = parse_haplotype(origin)
            if origin_info:
                origin_loc = f"{origin_info['background']}\nChr{origin_info['chromosome']}\n{origin_info['region']}"
                flow_data[mapped_loc][origin_loc] += 1

    # Create matrix for heatmap
    all_locations = sorted(set(list(flow_data.keys()) +
                              [dest for dests in flow_data.values() for dest in dests.keys()]))

    matrix = np.zeros((len(all_locations), len(all_locations)))

    for i, source in enumerate(all_locations):
        for j, dest in enumerate(all_locations):
            if source in flow_data and dest in flow_data[source]:
                matrix[i, j] = flow_data[source][dest]

    # Plot heatmap
    sns.heatmap(matrix, annot=True, fmt='.0f', cmap='YlOrRd', ax=ax,
               xticklabels=all_locations, yticklabels=all_locations,
               cbar_kws={'label': 'Number of Insertions'},
               linewidths=1, linecolor='white')

    ax.set_xlabel('Insertion Origin (by cenhapmer profile)', fontweight='bold', fontsize=12)
    ax.set_ylabel('Mapped Location (where insertion was found)', fontweight='bold', fontsize=12)
    ax.set_title('Insertion Origin Analysis: Evidence of Ectopic Recombination',
                fontweight='bold', fontsize=14, pad=20)

    # Highlight diagonal (same location = expansion, not ectopic)
    for i in range(len(all_locations)):
        rect = Rectangle((i, i), 1, 1, fill=False, edgecolor='blue', linewidth=3)
        ax.add_patch(rect)

    plt.tight_layout()
    plt.savefig(f"{output_prefix}_origin_heatmap.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_origin_heatmap.pdf", bbox_inches='tight')
    plt.close()

    print(f"Saved origin heatmap: {output_prefix}_origin_heatmap.png/pdf")


def create_ectopic_summary_plots(ectopic_events, output_prefix):
    """
    Create comprehensive plots for ectopic insertion events.
    """

    if not ectopic_events:
        print("No ectopic events detected")
        return

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Panel A: Counts by type
    ax = axes[0, 0]
    event_types = {
        'Different Chr': sum(1 for e in ectopic_events if e['is_different_chr']),
        'Different Region': sum(1 for e in ectopic_events if e['is_different_region'] and not e['is_different_chr']),
        'Different Background': sum(1 for e in ectopic_events if e['is_different_background']),
        'Same Location': len(ectopic_events) - sum(1 for e in ectopic_events
                                                    if e['is_different_chr'] or e['is_different_region'])
    }

    colors_map = ['#e74c3c', '#f39c12', '#3498db', '#27ae60']
    ax.bar(event_types.keys(), event_types.values(), color=colors_map, alpha=0.7, edgecolor='black', linewidth=1.5)
    ax.set_ylabel('Number of Insertions', fontweight='bold')
    ax.set_title('A. Ectopic Insertion Types', fontweight='bold', loc='left')
    ax.grid(axis='y', alpha=0.3)
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')

    # Panel B: Size distribution of ectopic insertions
    ax = axes[0, 1]
    sizes = [e['insertion']['size'] for e in ectopic_events]
    ax.hist(sizes, bins=30, color='#3498db', alpha=0.7, edgecolor='black', linewidth=1)
    ax.set_xlabel('Insertion Size (bp)', fontweight='bold')
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_title('B. Size Distribution of Ectopic Insertions', fontweight='bold', loc='left')
    ax.axvline(np.median(sizes), color='red', linestyle='--', linewidth=2, label=f'Median: {np.median(sizes):.0f} bp')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)

    # Panel C: Chromosome origin matrix
    ax = axes[1, 0]
    chr_flows = defaultdict(lambda: defaultdict(int))
    for e in ectopic_events:
        mapped_chr = e['mapped_location']['chromosome']
        origin_chr = e['origin_location']['chromosome']
        chr_flows[mapped_chr][origin_chr] += 1

    all_chrs = sorted(set([e['mapped_location']['chromosome'] for e in ectopic_events] +
                         [e['origin_location']['chromosome'] for e in ectopic_events]))

    chr_matrix = np.zeros((len(all_chrs), len(all_chrs)))
    for i, mapped_chr in enumerate(all_chrs):
        for j, origin_chr in enumerate(all_chrs):
            chr_matrix[i, j] = chr_flows[mapped_chr][origin_chr]

    sns.heatmap(chr_matrix, annot=True, fmt='.0f', cmap='YlOrRd', ax=ax,
               xticklabels=[f'Chr{c}' for c in all_chrs],
               yticklabels=[f'Chr{c}' for c in all_chrs],
               cbar_kws={'label': 'Count'}, linewidths=1, linecolor='white')
    ax.set_xlabel('Origin Chromosome', fontweight='bold')
    ax.set_ylabel('Mapped Chromosome', fontweight='bold')
    ax.set_title('C. Inter-Chromosomal Insertion Origins', fontweight='bold', loc='left')

    # Panel D: Region type transitions
    ax = axes[1, 1]
    region_flows = defaultdict(int)
    for e in ectopic_events:
        mapped_reg = e['mapped_location']['region']
        origin_reg = e['origin_location']['region']
        region_flows[f"{origin_reg}→{mapped_reg}"] += 1

    if region_flows:
        ax.bar(region_flows.keys(), region_flows.values(), color='#9b59b6', alpha=0.7,
              edgecolor='black', linewidth=1.5)
        ax.set_ylabel('Number of Insertions', fontweight='bold')
        ax.set_title('D. Region Type Transitions', fontweight='bold', loc='left')
        ax.grid(axis='y', alpha=0.3)
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')

    plt.suptitle(f'Ectopic Insertion Analysis: Evidence of Homologous Recombination\n'
                f'Total Ectopic Events: {len(ectopic_events)}',
                fontweight='bold', fontsize=14, y=0.995)

    plt.tight_layout()
    plt.savefig(f"{output_prefix}_ectopic_summary.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_ectopic_summary.pdf", bbox_inches='tight')
    plt.close()

    print(f"Saved ectopic summary: {output_prefix}_ectopic_summary.png/pdf")


def save_results_table(all_results, ectopic_events, output_prefix):
    """Save detailed results to TSV files."""

    # All insertions with profiles
    all_insertions = []
    for bam_name, results in all_results.items():
        haplotype = bam_name.split('_')[0]
        for result in results:
            all_insertions.append({
                'haplotype': haplotype,
                'read_id': result['read_id'],
                'ref_chr': result['ref_chr'],
                'ref_pos': result['ref_pos'],
                'size': result['size'],
                'origin': result['top_origin'],
                'origin_kmers': result['top_origin_kmers'],
                'is_ectopic': any(e['insertion']['read_id'] == result['read_id']
                                for e in ectopic_events)
            })

    df_all = pd.DataFrame(all_insertions)
    df_all.to_csv(f"{output_prefix}_all_insertions.tsv", sep='\t', index=False)
    print(f"Saved all insertions: {output_prefix}_all_insertions.tsv")

    # Ectopic events only
    if ectopic_events:
        ectopic_data = []
        for e in ectopic_events:
            ectopic_data.append({
                'read_id': e['insertion']['read_id'],
                'size': e['insertion']['size'],
                'mapped_background': e['mapped_location']['background'],
                'mapped_chr': e['mapped_location']['chromosome'],
                'mapped_region': e['mapped_location']['region'],
                'origin_background': e['origin_location']['background'],
                'origin_chr': e['origin_location']['chromosome'],
                'origin_region': e['origin_location']['region'],
                'origin_kmers': e['insertion']['top_origin_kmers'],
                'is_inter_chr': e['is_different_chr'],
                'is_inter_region': e['is_different_region']
            })

        df_ectopic = pd.DataFrame(ectopic_data)
        df_ectopic.to_csv(f"{output_prefix}_ectopic_events.tsv", sep='\t', index=False)
        print(f"Saved ectopic events: {output_prefix}_ectopic_events.tsv")


def main():
    parser = argparse.ArgumentParser(
        description='Analyze insertion origins using cenhapmer profiling',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument('--bam-dir', required=True,
                       help='Directory containing BAM files from non-hybrid mapping')
    parser.add_argument('--cenhapmer-dir', required=True,
                       help='Directory containing cenhapmer KMC databases')
    parser.add_argument('--kmer-size', type=int, default=31,
                       help='K-mer size used for cenhapmers (default: 31)')
    parser.add_argument('--min-size', type=int, default=50,
                       help='Minimum insertion size to analyze (default: 50bp)')
    parser.add_argument('--mapping-mode', default='wm_mapont',
                       help='Mapping mode to analyze (default: wm_mapont for winnowmap)')
    parser.add_argument('--output-prefix', required=True,
                       help='Output file prefix')

    args = parser.parse_args()

    # Create output directory
    output_dir = os.path.dirname(args.output_prefix)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    print("="*70)
    print("INSERTION ORIGIN ANALYSIS - Detecting Ectopic Recombination")
    print("="*70)
    print(f"\nParameters:")
    print(f"  BAM directory: {args.bam_dir}")
    print(f"  Cenhapmer DB: {args.cenhapmer_dir}")
    print(f"  K-mer size: {args.kmer_size}")
    print(f"  Min insertion size: {args.min_size}bp")

    # Find BAM files matching the specified mapping mode
    if args.mapping_mode:
        bam_pattern = f"*_{args.mapping_mode}.bam"
        bam_files = list(Path(args.bam_dir).glob(bam_pattern))

        if not bam_files:
            print(f"\nWarning: No BAM files matching pattern '{bam_pattern}' found.")
            print(f"Looking for any BAM files instead...")
            bam_files = list(Path(args.bam_dir).glob("*.bam"))
    else:
        bam_files = list(Path(args.bam_dir).glob("*.bam"))

    if not bam_files:
        print(f"\nError: No BAM files found in {args.bam_dir}")
        sys.exit(1)

    print(f"\nFound {len(bam_files)} BAM files to analyze (mapping mode: {args.mapping_mode})")

    # Show which files we're analyzing
    print(f"\nBAM files to process:")
    for bam in sorted(bam_files)[:5]:  # Show first 5
        print(f"  - {bam.name}")
    if len(bam_files) > 5:
        print(f"  ... and {len(bam_files) - 5} more")

    # Analyze each BAM file
    all_results = {}
    all_ectopic = []

    for bam_file in bam_files:
        basename = bam_file.stem
        haplotype = basename.split('_')[0]

        results = analyze_bam_insertions(
            str(bam_file),
            args.cenhapmer_dir,
            args.kmer_size,
            min_size=args.min_size
        )

        if results:
            all_results[basename] = results

            # Detect ectopic insertions
            ectopic = detect_ectopic_insertions(results, haplotype)
            all_ectopic.extend(ectopic)

            if ectopic:
                print(f"  → Found {len(ectopic)} ectopic insertions!")

    # Generate visualizations
    print(f"\n{''*70}")
    print("Generating visualizations...")

    if all_results:
        create_origin_sankey_plot(all_results, args.output_prefix)

    if all_ectopic:
        create_ectopic_summary_plots(all_ectopic, args.output_prefix)

        print(f"\n{'='*70}")
        print(f"ECTOPIC RECOMBINATION EVENTS DETECTED: {len(all_ectopic)}")
        print(f"{'='*70}")

        # Summary statistics
        inter_chr = sum(1 for e in all_ectopic if e['is_different_chr'])
        inter_region = sum(1 for e in all_ectopic if e['is_different_region'])

        print(f"\nEvent Types:")
        print(f"  Inter-chromosomal: {inter_chr}")
        print(f"  Inter-regional: {inter_region}")
        print(f"\nThis provides direct evidence of:")
        print(f"  1. Homologous recombination-mediated expansions")
        print(f"  2. Ectopic recombination between non-homologous regions")
        print(f"  3. Potential centromeric satellite repeat mobility")
    else:
        print("\nNo ectopic events detected (all insertions match their mapped location)")

    # Save results
    save_results_table(all_results, all_ectopic, args.output_prefix)

    print(f"\n{'='*70}")
    print("ANALYSIS COMPLETE")
    print(f"{'='*70}")
    print(f"\nResults saved with prefix: {args.output_prefix}")


if __name__ == '__main__':
    main()
