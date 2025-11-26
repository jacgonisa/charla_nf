#!/usr/bin/env python3
"""
Analyze cenhapmer (marker) database balance and distribution.

This script diagnoses potential biases in the marker database that could affect
crossover detection:
1. Unequal marker counts between parents (Col vs Ler)
2. Variable marker density across genome regions (CEN vs ARMS)
3. Chromosome-specific imbalances

Outputs:
- Balance metrics (ratios, coefficients of variation)
- Visualization plots
- Normalization recommendations
"""

import os
import sys
import argparse
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict


def count_kmers_in_database(kmc_prefix):
    """
    Count total k-mers in a KMC database using kmc_tools.

    Parameters
    ----------
    kmc_prefix : str
        Path to KMC database (without .kmc_pre/.kmc_suf extension)

    Returns
    -------
    int
        Total number of unique k-mers
    """
    try:
        # Use kmc_tools info to get stats
        cmd = f"kmc_tools info {kmc_prefix}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)

        # Parse output for total unique k-mers
        for line in result.stdout.split('\n'):
            if 'Total number of k-mers' in line or 'Total no. of unique k-mers' in line:
                # Extract number
                parts = line.split(':')
                if len(parts) > 1:
                    count = int(parts[1].strip())
                    return count

        # Fallback: count lines in dump
        cmd_dump = f"kmc_tools transform {kmc_prefix} dump /dev/stdout | wc -l"
        result = subprocess.run(cmd_dump, shell=True, capture_output=True, text=True, check=True)
        return int(result.stdout.strip())

    except subprocess.CalledProcessError as e:
        print(f"Warning: Could not count k-mers in {kmc_prefix}: {e}", file=sys.stderr)
        return 0
    except Exception as e:
        print(f"Error processing {kmc_prefix}: {e}", file=sys.stderr)
        return 0


def analyze_marker_database(cenhapmer_dir, output_prefix):
    """
    Analyze marker balance across the cenhapmer database.

    Expected directory structure (supports two filename formats):

    Format 1 (compact):
    cenhapmer_dir/
    ├── Col-0_CA1.kmc_pre/suf  # Col chromosome 1 arms
    ├── Col-0_CC1.kmc_pre/suf  # Col chromosome 1 centromere
    ├── Ler-0_LA1.kmc_pre/suf  # Ler chromosome 1 arms
    ├── Ler-0_LC1.kmc_pre/suf  # Ler chromosome 1 centromere
    └── ...

    Format 2 (verbose):
    cenhapmer_dir/
    ├── unique_Col-0_ARMS_Chr1_k31.kmc_pre/suf
    ├── unique_Col-0_CEN_Chr1_k31.kmc_pre/suf
    ├── unique_Ler-0_ARMS_Chr1_k31.kmc_pre/suf
    ├── unique_Ler-0_CEN_Chr1_k31.kmc_pre/suf
    └── ...
    """

    print(f"Analyzing cenhapmer databases in: {cenhapmer_dir}")

    # Find all KMC databases
    kmc_files = []
    for root, dirs, files in os.walk(cenhapmer_dir):
        for f in files:
            if f.endswith('.kmc_pre'):
                kmc_prefix = os.path.join(root, f.replace('.kmc_pre', ''))
                kmc_files.append(kmc_prefix)

    if not kmc_files:
        print(f"Error: No KMC databases found in {cenhapmer_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(kmc_files)} KMC databases")

    # Parse filenames and count k-mers
    results = []
    for kmc_file in sorted(kmc_files):
        basename = os.path.basename(kmc_file)

        # Parse filename - support two formats:
        # Format 1: Parent_RegionChr (e.g., Col-0_CA1, Ler-0_LC2)
        # Format 2: unique_Parent_RegionType_ChrN_kSize (e.g., unique_Col-0_ARMS_Chr1_k31)

        parts = basename.split('_')

        # Detect format based on number of parts
        if len(parts) >= 5 and parts[0] == 'unique':
            # Format 2: unique_Parent_RegionType_ChrN_kSize
            parent = parts[1]  # Col-0 or Ler-0
            region_type_str = parts[2]  # ARMS or CEN
            chr_part = parts[3]  # Chr1, Chr2, etc.

            # Extract chromosome number
            chr_num = chr_part.replace('Chr', '')  # 1, 2, 3, 4, 5

            # Determine region type and code
            if region_type_str == 'ARMS':
                region_type = 'ARMS'
                # Create region code: CA (Col ARMS), LA (Ler ARMS)
                if 'Col' in parent:
                    region_code = 'CA' + chr_num
                elif 'Ler' in parent:
                    region_code = 'LA' + chr_num
                else:
                    region_code = 'XA' + chr_num  # Unknown parent
            elif region_type_str == 'CEN':
                region_type = 'CEN'
                # Create region code: CC (Col CEN), LC (Ler CEN)
                if 'Col' in parent:
                    region_code = 'CC' + chr_num
                elif 'Ler' in parent:
                    region_code = 'LC' + chr_num
                else:
                    region_code = 'XC' + chr_num  # Unknown parent
            else:
                region_type = 'UNKNOWN'
                region_code = 'XX' + chr_num

        elif len(parts) >= 2:
            # Format 1: Parent_RegionChr (e.g., Col-0_CA1, Ler-0_LC2)
            parent = parts[0]  # Col-0 or Ler-0
            region_chr = parts[1]  # CA1, CC1, LA1, LC1, etc.

            # Extract region (C=CEN, A=ARMS) and chromosome number
            if len(region_chr) < 3:
                print(f"Warning: Cannot parse region/chr from {region_chr}", file=sys.stderr)
                continue

            region_code = region_chr[:2]  # CA, CC, LA, LC
            chr_num = region_chr[2:]      # 1, 2, 3, 4, 5

            # Determine region type
            if region_code[1] == 'A':
                region_type = 'ARMS'
            elif region_code[1] == 'C':
                region_type = 'CEN'
            else:
                region_type = 'UNKNOWN'
        else:
            print(f"Warning: Unexpected filename format: {basename}", file=sys.stderr)
            continue

        print(f"Counting k-mers in {basename}...", end=' ')
        kmer_count = count_kmers_in_database(kmc_file)
        print(f"{kmer_count:,} k-mers")

        results.append({
            'database': basename,
            'parent': parent,
            'region_code': region_code,
            'region_type': region_type,
            'chromosome': chr_num,
            'kmer_count': kmer_count,
            'kmc_path': kmc_file
        })

    # Create DataFrame
    df = pd.DataFrame(results)

    # Save raw counts
    output_tsv = f"{output_prefix}_marker_counts.tsv"
    df.to_csv(output_tsv, sep='\t', index=False)
    print(f"\nSaved marker counts to: {output_tsv}")

    return df


def calculate_balance_metrics(df, output_prefix):
    """
    Calculate balance metrics and identify potential biases.
    """

    print("\n" + "="*70)
    print("MARKER BALANCE ANALYSIS")
    print("="*70)

    # Overall parent balance
    parent_totals = df.groupby('parent')['kmer_count'].sum()
    print("\n1. OVERALL PARENT BALANCE")
    print("-" * 40)
    for parent, count in parent_totals.items():
        print(f"  {parent}: {count:,} markers")

    if len(parent_totals) == 2:
        ratio = parent_totals.max() / parent_totals.min()
        imbalance_pct = (ratio - 1) * 100
        print(f"\n  Imbalance ratio: {ratio:.2f}:1")
        print(f"  Imbalance: {imbalance_pct:.1f}%")

        if ratio > 1.2:
            print(f"  ⚠️  WARNING: >20% imbalance between parents!")
        elif ratio > 1.1:
            print(f"  ⚡ CAUTION: >10% imbalance detected")
        else:
            print(f"  ✓ Parents are well-balanced")

    # Region-specific balance
    print("\n2. REGION-SPECIFIC BALANCE (CEN vs ARMS)")
    print("-" * 40)
    region_totals = df.groupby(['parent', 'region_type'])['kmer_count'].sum().unstack(fill_value=0)
    print(region_totals)

    if 'ARMS' in region_totals.columns and 'CEN' in region_totals.columns:
        for parent in region_totals.index:
            arms_count = region_totals.loc[parent, 'ARMS']
            cen_count = region_totals.loc[parent, 'CEN']
            ratio = arms_count / cen_count if cen_count > 0 else float('inf')
            print(f"\n  {parent} ARMS:CEN ratio: {ratio:.2f}:1")

    # Chromosome-specific balance
    print("\n3. CHROMOSOME-SPECIFIC BALANCE")
    print("-" * 40)
    chr_totals = df.groupby(['parent', 'chromosome'])['kmer_count'].sum().unstack(fill_value=0)
    print(chr_totals)

    # Calculate coefficient of variation per parent
    for parent in chr_totals.index:
        counts = chr_totals.loc[parent].values
        cv = np.std(counts) / np.mean(counts) * 100
        print(f"\n  {parent} CV across chromosomes: {cv:.1f}%")
        if cv > 30:
            print(f"    ⚠️  High variability across chromosomes!")

    # Save summary metrics
    summary_file = f"{output_prefix}_balance_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("MARKER BALANCE SUMMARY\n")
        f.write("=" * 70 + "\n\n")

        f.write("Parent Totals:\n")
        for parent, count in parent_totals.items():
            f.write(f"  {parent}: {count:,}\n")

        if len(parent_totals) == 2:
            f.write(f"\nImbalance Ratio: {ratio:.2f}:1\n")
            f.write(f"Imbalance Percentage: {imbalance_pct:.1f}%\n")

        f.write("\n\nRegion Totals:\n")
        f.write(region_totals.to_string())

        f.write("\n\nChromosome Totals:\n")
        f.write(chr_totals.to_string())

    print(f"\nSaved balance summary to: {summary_file}")

    return parent_totals, region_totals, chr_totals


def plot_marker_distribution(df, output_prefix):
    """
    Create visualizations of marker distribution.
    """

    print("\nGenerating visualizations...")

    # Set style
    sns.set_style("whitegrid")

    # 1. Overall parent comparison
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # 1a. Total markers per parent
    ax = axes[0, 0]
    parent_totals = df.groupby('parent')['kmer_count'].sum()
    colors = ['#3498db', '#e74c3c'] if len(parent_totals) == 2 else None
    parent_totals.plot(kind='bar', ax=ax, color=colors, edgecolor='black', alpha=0.7)
    ax.set_title('Total Markers per Parent', fontsize=13, fontweight='bold')
    ax.set_xlabel('Parent', fontsize=11)
    ax.set_ylabel('Total Marker Count', fontsize=11)
    ax.tick_params(axis='x', rotation=0)
    for i, v in enumerate(parent_totals):
        ax.text(i, v, f'{v:,.0f}', ha='center', va='bottom', fontweight='bold')

    # 1b. Markers per region type
    ax = axes[0, 1]
    region_data = df.groupby(['parent', 'region_type'])['kmer_count'].sum().unstack(fill_value=0)
    region_data.plot(kind='bar', ax=ax, color=['#2ecc71', '#f39c12'], edgecolor='black', alpha=0.7)
    ax.set_title('Markers by Region Type', fontsize=13, fontweight='bold')
    ax.set_xlabel('Parent', fontsize=11)
    ax.set_ylabel('Marker Count', fontsize=11)
    ax.legend(title='Region', fontsize=9)
    ax.tick_params(axis='x', rotation=0)

    # 1c. Markers per chromosome
    ax = axes[1, 0]
    chr_data = df.groupby(['parent', 'chromosome'])['kmer_count'].sum().unstack(fill_value=0)
    chr_data.T.plot(kind='bar', ax=ax, color=['#3498db', '#e74c3c'], edgecolor='black', alpha=0.7, width=0.8)
    ax.set_title('Markers per Chromosome', fontsize=13, fontweight='bold')
    ax.set_xlabel('Chromosome', fontsize=11)
    ax.set_ylabel('Marker Count', fontsize=11)
    ax.legend(title='Parent', fontsize=9)
    ax.tick_params(axis='x', rotation=0)

    # 1d. Heatmap of markers by parent x region x chromosome
    ax = axes[1, 1]
    pivot = df.pivot_table(values='kmer_count', index=['parent', 'region_type'],
                           columns='chromosome', fill_value=0)
    sns.heatmap(pivot, annot=True, fmt='.0f', cmap='YlOrRd', ax=ax,
                cbar_kws={'label': 'Marker Count'}, linewidths=0.5)
    ax.set_title('Marker Heatmap: Parent x Region x Chr', fontsize=13, fontweight='bold')
    ax.set_xlabel('Chromosome', fontsize=11)
    ax.set_ylabel('Parent / Region', fontsize=11)

    plt.tight_layout()
    plot_file = f"{output_prefix}_marker_distribution.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved distribution plot: {plot_file}")

    # 2. Detailed per-database plot
    fig, ax = plt.subplots(figsize=(16, 8))

    # Sort by parent, region, chromosome
    df_sorted = df.sort_values(['parent', 'region_type', 'chromosome'])
    x_labels = df_sorted['database'].tolist()
    y_values = df_sorted['kmer_count'].tolist()

    # Color by parent
    colors_list = ['#3498db' if 'Col' in db else '#e74c3c' for db in x_labels]

    bars = ax.bar(range(len(x_labels)), y_values, color=colors_list, edgecolor='black', alpha=0.7)
    ax.set_xticks(range(len(x_labels)))
    ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=8)
    ax.set_xlabel('Database', fontsize=12, fontweight='bold')
    ax.set_ylabel('Marker Count', fontsize=12, fontweight='bold')
    ax.set_title('Marker Count per Database', fontsize=14, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)

    # Add value labels on bars
    for i, (bar, val) in enumerate(zip(bars, y_values)):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{val:,.0f}', ha='center', va='bottom', fontsize=7, rotation=90)

    plt.tight_layout()
    plot_file2 = f"{output_prefix}_per_database.png"
    plt.savefig(plot_file2, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved per-database plot: {plot_file2}")


def generate_normalization_recommendations(df, parent_totals, output_prefix):
    """
    Generate recommendations for normalization strategies.
    """

    rec_file = f"{output_prefix}_normalization_recommendations.txt"

    with open(rec_file, 'w') as f:
        f.write("MARKER BALANCE NORMALIZATION RECOMMENDATIONS\n")
        f.write("=" * 70 + "\n\n")

        # Check imbalance level
        if len(parent_totals) == 2:
            ratio = parent_totals.max() / parent_totals.min()
            imbalance_pct = (ratio - 1) * 100

            f.write(f"Overall Imbalance: {imbalance_pct:.1f}% ({ratio:.2f}:1 ratio)\n\n")

            if ratio < 1.1:
                f.write("STATUS: ✓ Well-balanced database\n")
                f.write("  - Imbalance <10% is acceptable for most analyses\n")
                f.write("  - No normalization required\n\n")

            elif ratio < 1.2:
                f.write("STATUS: ⚡ Moderate imbalance (10-20%)\n")
                f.write("  - Consider normalization for high-precision analyses\n\n")
                f.write("RECOMMENDED STRATEGIES:\n")
                f.write("  1. Weighted scoring (EASIEST):\n")
                f.write("     - Multiply k-mer counts by inverse frequency\n")
                f.write("     - Weight = 1 / (marker_count_for_region / total_markers)\n\n")
                f.write("  2. Add to confidence scoring:\n")
                f.write("     - Use marker balance as additional component\n")
                f.write("     - Penalize reads in over-represented regions\n\n")

            else:
                f.write("STATUS: ⚠️  HIGH IMBALANCE (>20%)\n")
                f.write("  - Normalization STRONGLY RECOMMENDED\n\n")
                f.write("RECOMMENDED STRATEGIES (in order of preference):\n\n")
                f.write("  1. REGENERATE CENHAPMER DATABASE (BEST):\n")
                f.write("     - Use stricter filtering for over-represented parent\n")
                f.write("     - Or more lenient filtering for under-represented parent\n")
                f.write("     - Aim for <10% imbalance\n\n")
                f.write("  2. DOWN-SAMPLE over-represented regions:\n")
                f.write("     - Randomly remove markers from abundant regions\n")
                f.write("     - Match marker density to under-represented parent\n\n")
                f.write("  3. NORMALIZE k-mer counts:\n")
                f.write("     - Use probability scores instead of raw counts\n")
                f.write("     - Score = (count_parent_A / markers_parent_A) vs\n")
                f.write("               (count_parent_B / markers_parent_B)\n\n")
                f.write("  4. REGION-SPECIFIC thresholds:\n")
                f.write("     - Adjust min_token_count based on marker density\n")
                f.write("     - Higher threshold for marker-rich regions\n\n")

        # Chromosome-specific recommendations
        f.write("\nCHROMOSOME-SPECIFIC CONSIDERATIONS:\n")
        f.write("-" * 40 + "\n")
        chr_totals = df.groupby('chromosome')['kmer_count'].sum().sort_values(ascending=False)
        cv = np.std(chr_totals) / np.mean(chr_totals) * 100

        f.write(f"Coefficient of variation across chromosomes: {cv:.1f}%\n\n")

        if cv > 30:
            f.write("⚠️  High variability detected!\n")
            f.write("  - Consider chromosome-specific normalization\n")
            f.write("  - Flag crossovers in low-marker chromosomes\n")
        else:
            f.write("✓ Chromosomes have similar marker densities\n")

        f.write("\n\nREGION-SPECIFIC CONSIDERATIONS (CEN vs ARMS):\n")
        f.write("-" * 40 + "\n")
        region_totals = df.groupby('region_type')['kmer_count'].sum()
        if 'ARMS' in region_totals and 'CEN' in region_totals:
            region_ratio = region_totals['ARMS'] / region_totals['CEN']
            f.write(f"ARMS:CEN marker ratio: {region_ratio:.2f}:1\n\n")

            if region_ratio > 2:
                f.write("⚠️  ARMS have much higher marker density than CEN\n")
                f.write("  - Expected: centromeres are repetitive, fewer unique markers\n")
                f.write("  - Impact: Crossovers in CEN may be less confident\n")
                f.write("  - Mitigation: Use region-aware confidence scoring\n\n")

    print(f"\nSaved normalization recommendations to: {rec_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Analyze cenhapmer marker database balance',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze marker balance
  python 16-analyze_marker_balance.py \\
      --cenhapmer-dir /path/to/cenhapmers/k31 \\
      --output-prefix results/marker_analysis/k31_balance

  # The script will generate:
  #   - k31_balance_marker_counts.tsv (raw counts)
  #   - k31_balance_balance_summary.txt (metrics)
  #   - k31_balance_marker_distribution.png (plots)
  #   - k31_balance_normalization_recommendations.txt (guidance)
        """
    )

    parser.add_argument('--cenhapmer-dir', required=True,
                       help='Directory containing cenhapmer KMC databases')
    parser.add_argument('--output-prefix', required=True,
                       help='Output file prefix for results')

    args = parser.parse_args()

    # Create output directory
    output_dir = os.path.dirname(args.output_prefix)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # 1. Analyze marker database
    df = analyze_marker_database(args.cenhapmer_dir, args.output_prefix)

    # 2. Calculate balance metrics
    parent_totals, region_totals, chr_totals = calculate_balance_metrics(df, args.output_prefix)

    # 3. Generate visualizations
    plot_marker_distribution(df, args.output_prefix)

    # 4. Generate normalization recommendations
    generate_normalization_recommendations(df, parent_totals, args.output_prefix)

    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
    print(f"\nAll results saved with prefix: {args.output_prefix}")


if __name__ == '__main__':
    main()
