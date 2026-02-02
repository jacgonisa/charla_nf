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
from scipy import stats
import warnings
warnings.filterwarnings('ignore', category=UserWarning)

# Publication-ready plot settings
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans']
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['xtick.labelsize'] = 9
plt.rcParams['ytick.labelsize'] = 9
plt.rcParams['legend.fontsize'] = 9
plt.rcParams['figure.titlesize'] = 13
plt.rcParams['pdf.fonttype'] = 42  # TrueType fonts for PDF
plt.rcParams['ps.fonttype'] = 42


def load_bed_file_sizes(bed_file):
    """
    Load genomic region sizes from BED file.

    BED format: Chr_Region_Parent_Arm  Start  End
    Example: Chr1_ARMS_Col_1  0  14841147

    Returns dict: {(parent, region_type, chr): size_in_bp}
    """
    region_sizes = {}

    if not os.path.exists(bed_file):
        return region_sizes

    with open(bed_file) as f:
        for line in f:
            if not line.strip():
                continue

            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue

            name = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            size_bp = end - start

            # Parse name: Chr1_ARMS_Col_1 or Chr1_CEN_Col
            name_parts = name.split('_')
            if len(name_parts) >= 3:
                chr_part = name_parts[0]  # Chr1, Chr2, etc.
                region_part = name_parts[1]  # ARMS or CEN
                parent_part = name_parts[2]  # Col or Ler

                chr_num = chr_part.replace('Chr', '')

                # Normalize parent name
                if parent_part == 'Col':
                    parent = 'Col-0'
                elif parent_part == 'Ler':
                    parent = 'Ler-0'
                else:
                    parent = parent_part

                region_type = region_part  # ARMS or CEN

                # Sum sizes for multiple arms (Chr1_ARMS_Col_1 + Chr1_ARMS_Col_2)
                key = (parent, region_type, chr_num)
                if key in region_sizes:
                    region_sizes[key] += size_bp
                else:
                    region_sizes[key] = size_bp

    return region_sizes


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


def analyze_marker_database(cenhapmer_dir, output_prefix, col_bed=None, ler_bed=None):
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

    # Load genomic region sizes from BED files if provided
    region_sizes = {}
    if col_bed and os.path.exists(col_bed):
        print(f"\nLoading Col-0 region sizes from: {col_bed}")
        region_sizes.update(load_bed_file_sizes(col_bed))
    if ler_bed and os.path.exists(ler_bed):
        print(f"Loading Ler-0 region sizes from: {ler_bed}")
        region_sizes.update(load_bed_file_sizes(ler_bed))

    # Calculate marker density (markers per Mb) if region sizes available
    if region_sizes:
        print(f"\nCalculating marker density (markers per Mb)...")
        df['region_size_bp'] = df.apply(
            lambda row: region_sizes.get((row['parent'], row['region_type'], row['chromosome']), 0),
            axis=1
        )
        df['region_size_mb'] = df['region_size_bp'] / 1_000_000
        df['marker_density'] = df.apply(
            lambda row: row['kmer_count'] / row['region_size_mb'] if row['region_size_mb'] > 0 else 0,
            axis=1
        )
        print(f"Added marker density column (markers/Mb)")
    else:
        print(f"\nWarning: No BED files provided, marker density not calculated")
        df['region_size_bp'] = 0
        df['region_size_mb'] = 0.0
        df['marker_density'] = 0.0

    # Save raw counts with density
    output_tsv = f"{output_prefix}_marker_counts.tsv"
    df.to_csv(output_tsv, sep='\t', index=False)
    print(f"Saved marker counts to: {output_tsv}")

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

    # ========== DENSITY-BASED ANALYSIS (NORMALIZED BY GENOMIC SIZE) ==========
    if 'marker_density' in df.columns and df['marker_density'].sum() > 0:
        print("\n" + "="*70)
        print("MARKER DENSITY ANALYSIS (Normalized by genomic region size)")
        print("="*70)

        # Overall parent density balance
        parent_density = df.groupby('parent', group_keys=False).apply(
            lambda x: x['kmer_count'].sum() / x['region_size_mb'].sum() if x['region_size_mb'].sum() > 0 else 0,
            include_groups=False
        )
        print("\n4. OVERALL PARENT DENSITY (markers/Mb)")
        print("-" * 40)
        for parent, density in parent_density.items():
            print(f"  {parent}: {density:,.1f} markers/Mb")

        if len(parent_density) == 2:
            ratio = parent_density.max() / parent_density.min()
            imbalance_pct = (ratio - 1) * 100
            print(f"\n  Density imbalance ratio: {ratio:.2f}:1")
            print(f"  Density imbalance: {imbalance_pct:.1f}%")

            if ratio > 1.2:
                print(f"  ⚠️  WARNING: >20% density imbalance between parents!")
            elif ratio > 1.1:
                print(f"  ⚡ CAUTION: >10% density imbalance detected")
            else:
                print(f"  ✓ Parents have well-balanced marker density")

        # Region-specific density
        print("\n5. REGION-SPECIFIC DENSITY (markers/Mb)")
        print("-" * 40)
        region_density = df.groupby(['parent', 'region_type'], group_keys=False).apply(
            lambda x: x['kmer_count'].sum() / x['region_size_mb'].sum() if x['region_size_mb'].sum() > 0 else 0,
            include_groups=False
        ).unstack(fill_value=0)
        print(region_density)

        if 'ARMS' in region_density.columns and 'CEN' in region_density.columns:
            for parent in region_density.index:
                arms_density = region_density.loc[parent, 'ARMS']
                cen_density = region_density.loc[parent, 'CEN']
                ratio = arms_density / cen_density if cen_density > 0 else float('inf')
                print(f"\n  {parent} ARMS:CEN density ratio: {ratio:.2f}:1")
                print(f"    (This is the TRUE marker density difference)")

        # Chromosome-specific density
        print("\n6. CHROMOSOME-SPECIFIC DENSITY (markers/Mb)")
        print("-" * 40)
        chr_density = df.groupby(['parent', 'chromosome'], group_keys=False).apply(
            lambda x: x['kmer_count'].sum() / x['region_size_mb'].sum() if x['region_size_mb'].sum() > 0 else 0,
            include_groups=False
        ).unstack(fill_value=0)
        print(chr_density)

        # Calculate coefficient of variation per parent (DENSITY)
        for parent in chr_density.index:
            densities = chr_density.loc[parent].values
            cv = np.std(densities) / np.mean(densities) * 100 if np.mean(densities) > 0 else 0
            print(f"\n  {parent} CV across chromosomes (density): {cv:.1f}%")
            if cv > 30:
                print(f"    ⚠️  High density variability across chromosomes!")
    else:
        print("\n⚠️  Marker density not calculated (provide --col-bed and --ler-bed)")

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


def format_number(num):
    """Format large numbers for plot labels (e.g., 1.2M, 450K)"""
    if num >= 1_000_000:
        return f'{num/1_000_000:.1f}M'
    elif num >= 1_000:
        return f'{num/1_000:.0f}K'
    else:
        return f'{num:.0f}'


def add_significance_bar(ax, x1, x2, y, p_value, height_offset=0):
    """Add significance bar and p-value annotation to plot"""
    y_max = y + height_offset
    h = y_max * 0.02

    # Draw significance bar
    ax.plot([x1, x1, x2, x2], [y_max, y_max + h, y_max + h, y_max],
            lw=1.5, c='black')

    # Add p-value text
    if p_value < 0.001:
        sig_text = '***'
    elif p_value < 0.01:
        sig_text = '**'
    elif p_value < 0.05:
        sig_text = '*'
    else:
        sig_text = 'ns'

    ax.text((x1 + x2) / 2, y_max + h, sig_text,
            ha='center', va='bottom', fontsize=11, fontweight='bold')


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


def plot_marker_density_publication(df, output_prefix, format='png'):
    """
    Create publication-ready marker density visualizations with statistical annotations.

    Parameters
    ----------
    df : pd.DataFrame
        Marker count data with density information
    output_prefix : str
        Output file prefix
    format : str
        Output format ('png', 'pdf', 'svg', or 'all')
    """

    if 'marker_density' not in df.columns or df['marker_density'].sum() == 0:
        print("\nSkipping density plots (no BED files provided)")
        return

    print("\nGenerating publication-ready density visualizations...")

    # Define color palette
    colors = {
        'Col-0': '#2C7BB6',  # Blue
        'Ler-0': '#D7191C',  # Red
        'ARMS': '#66C2A5',   # Teal
        'CEN': '#FC8D62'     # Orange
    }

    # ========== FIGURE 1: Comprehensive 4-panel density analysis ==========
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)

    # Panel A: Marker density by parent
    ax1 = fig.add_subplot(gs[0, 0])
    parent_density = df.groupby('parent', group_keys=False).apply(
        lambda x: x['kmer_count'].sum() / x['region_size_mb'].sum()
        if x['region_size_mb'].sum() > 0 else 0,
        include_groups=False
    )

    bars = ax1.bar(range(len(parent_density)), parent_density.values,
                   color=[colors.get(p, '#888888') for p in parent_density.index],
                   edgecolor='black', linewidth=1.5, alpha=0.8)

    ax1.set_xticks(range(len(parent_density)))
    ax1.set_xticklabels(parent_density.index, fontweight='bold')
    ax1.set_ylabel('Marker Density (markers/Mb)', fontweight='bold')
    ax1.set_title('A. Overall Marker Density by Parent', fontweight='bold', loc='left')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # Add value labels
    for i, (bar, val) in enumerate(zip(bars, parent_density.values)):
        ax1.text(bar.get_x() + bar.get_width()/2, val,
                f'{val:,.0f}', ha='center', va='bottom',
                fontweight='bold', fontsize=10)

    # Add statistical test if two parents
    if len(parent_density) == 2:
        parent1_data = df[df['parent'] == parent_density.index[0]]['marker_density'].values
        parent2_data = df[df['parent'] == parent_density.index[1]]['marker_density'].values
        t_stat, p_val = stats.ttest_ind(parent1_data, parent2_data)

        y_max = max(parent_density.values)
        add_significance_bar(ax1, 0, 1, y_max, p_val, height_offset=y_max*0.05)

    # Panel B: Marker density by region type
    ax2 = fig.add_subplot(gs[0, 1])
    region_density = df.groupby(['parent', 'region_type'], group_keys=False).apply(
        lambda x: x['kmer_count'].sum() / x['region_size_mb'].sum()
        if x['region_size_mb'].sum() > 0 else 0,
        include_groups=False
    ).unstack(fill_value=0)

    x = np.arange(len(region_density))
    width = 0.35

    if 'ARMS' in region_density.columns and 'CEN' in region_density.columns:
        bars1 = ax2.bar(x - width/2, region_density['ARMS'], width,
                       label='ARMS', color=colors['ARMS'],
                       edgecolor='black', linewidth=1.5, alpha=0.8)
        bars2 = ax2.bar(x + width/2, region_density['CEN'], width,
                       label='CEN', color=colors['CEN'],
                       edgecolor='black', linewidth=1.5, alpha=0.8)

        # Add value labels
        for bars in [bars1, bars2]:
            for bar in bars:
                height = bar.get_height()
                ax2.text(bar.get_x() + bar.get_width()/2, height,
                        f'{height:,.0f}', ha='center', va='bottom',
                        fontsize=8, fontweight='bold')

    ax2.set_xticks(x)
    ax2.set_xticklabels(region_density.index, fontweight='bold')
    ax2.set_ylabel('Marker Density (markers/Mb)', fontweight='bold')
    ax2.set_title('B. Marker Density by Region Type', fontweight='bold', loc='left')
    ax2.legend(frameon=False, loc='upper right')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # Panel C: Marker density across chromosomes
    ax3 = fig.add_subplot(gs[1, 0])
    chr_density = df.groupby(['parent', 'chromosome'], group_keys=False).apply(
        lambda x: x['kmer_count'].sum() / x['region_size_mb'].sum()
        if x['region_size_mb'].sum() > 0 else 0,
        include_groups=False
    ).unstack(fill_value=0)

    chr_order = sorted(chr_density.columns, key=lambda x: int(x))
    x_pos = np.arange(len(chr_order))
    width = 0.35

    parents = chr_density.index.tolist()
    if len(parents) == 2:
        bars1 = ax3.bar(x_pos - width/2, chr_density.loc[parents[0], chr_order], width,
                       label=parents[0], color=colors.get(parents[0], '#888'),
                       edgecolor='black', linewidth=1.5, alpha=0.8)
        bars2 = ax3.bar(x_pos + width/2, chr_density.loc[parents[1], chr_order], width,
                       label=parents[1], color=colors.get(parents[1], '#888'),
                       edgecolor='black', linewidth=1.5, alpha=0.8)

    ax3.set_xticks(x_pos)
    ax3.set_xticklabels([f'Chr{c}' for c in chr_order])
    ax3.set_xlabel('Chromosome', fontweight='bold')
    ax3.set_ylabel('Marker Density (markers/Mb)', fontweight='bold')
    ax3.set_title('C. Chromosome-Specific Marker Density', fontweight='bold', loc='left')
    ax3.legend(frameon=False, loc='upper right')
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)

    # Panel D: Violin plot of marker density distribution
    ax4 = fig.add_subplot(gs[1, 1])

    # Prepare data for violin plot
    plot_data = []
    for parent in df['parent'].unique():
        for region in df['region_type'].unique():
            subset = df[(df['parent'] == parent) & (df['region_type'] == region)]
            if len(subset) > 0:
                for density in subset['marker_density']:
                    plot_data.append({
                        'Parent': parent,
                        'Region': region,
                        'Density': density,
                        'Group': f'{parent}\n{region}'
                    })

    plot_df = pd.DataFrame(plot_data)

    # Create violin plot
    groups = plot_df['Group'].unique()
    positions = np.arange(len(groups))

    parts = ax4.violinplot([plot_df[plot_df['Group'] == g]['Density'].values
                            for g in groups],
                           positions=positions, widths=0.7,
                           showmeans=True, showmedians=True)

    # Color the violins
    for i, pc in enumerate(parts['bodies']):
        group = groups[i]
        if 'Col' in group:
            color = colors['Col-0']
        elif 'Ler' in group:
            color = colors['Ler-0']
        else:
            color = '#888888'
        pc.set_facecolor(color)
        pc.set_alpha(0.6)
        pc.set_edgecolor('black')
        pc.set_linewidth(1.5)

    ax4.set_xticks(positions)
    ax4.set_xticklabels(groups, fontsize=8)
    ax4.set_ylabel('Marker Density (markers/Mb)', fontweight='bold')
    ax4.set_title('D. Marker Density Distribution', fontweight='bold', loc='left')
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)
    ax4.grid(axis='y', alpha=0.3, linestyle='--')

    # Save figure in requested format(s)
    formats_to_save = ['png', 'pdf', 'svg'] if format == 'all' else [format]

    for fmt in formats_to_save:
        output_file = f"{output_prefix}_density_analysis_publication.{fmt}"
        dpi = 300 if fmt == 'png' else None
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight', format=fmt)
        print(f"  Saved: {output_file}")

    plt.close()

    # ========== FIGURE 2: Genome-wide density landscape ==========
    fig, axes = plt.subplots(len(df['parent'].unique()), 1,
                             figsize=(16, 4 * len(df['parent'].unique())),
                             sharex=True)

    if len(df['parent'].unique()) == 1:
        axes = [axes]

    for idx, parent in enumerate(sorted(df['parent'].unique())):
        ax = axes[idx]
        parent_df = df[df['parent'] == parent].copy()
        parent_df = parent_df.sort_values(['chromosome', 'region_type'])

        x_pos = 0
        x_labels = []
        x_ticks = []
        chr_boundaries = []

        for chr_num in sorted(parent_df['chromosome'].unique(), key=int):
            chr_df = parent_df[parent_df['chromosome'] == chr_num]

            for _, row in chr_df.iterrows():
                width = row['region_size_mb']
                density = row['marker_density']

                # Color by region type
                color = colors.get(row['region_type'], '#888888')

                ax.barh(0, width, left=x_pos, height=0.8,
                       color=color, edgecolor='black', linewidth=1,
                       alpha=0.8)

                # Add density text
                if width > 1:  # Only label if wide enough
                    ax.text(x_pos + width/2, 0, f'{density:.0f}',
                           ha='center', va='center', fontweight='bold',
                           fontsize=8)

                x_labels.append(f"{row['region_type']}")
                x_ticks.append(x_pos + width/2)
                x_pos += width

            chr_boundaries.append(x_pos)

        # Add chromosome boundaries
        for boundary in chr_boundaries[:-1]:
            ax.axvline(boundary, color='black', linewidth=2, linestyle='--')

        ax.set_xlim(0, x_pos)
        ax.set_ylim(-0.5, 0.5)
        ax.set_yticks([0])
        ax.set_yticklabels([parent], fontweight='bold', fontsize=12)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.tick_params(left=False)

        if idx == len(df['parent'].unique()) - 1:
            ax.set_xlabel('Genomic Position (Mb)', fontweight='bold')

        # Add chromosome labels at the top
        if idx == 0:
            chr_centers = []
            x_pos = 0
            for chr_num in sorted(parent_df['chromosome'].unique(), key=int):
                chr_df = parent_df[parent_df['chromosome'] == chr_num]
                chr_width = chr_df['region_size_mb'].sum()
                chr_centers.append((x_pos + chr_width/2, chr_num))
                x_pos += chr_width

            for x, chr_num in chr_centers:
                ax.text(x, 0.6, f'Chr{chr_num}', ha='center', va='bottom',
                       fontweight='bold', fontsize=10)

    plt.suptitle('Genome-wide Marker Density Landscape',
                fontweight='bold', fontsize=14, y=0.995)

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=colors['ARMS'], edgecolor='black',
                            label='Chromosome Arms', alpha=0.8),
                      Patch(facecolor=colors['CEN'], edgecolor='black',
                            label='Centromere', alpha=0.8)]
    axes[0].legend(handles=legend_elements, loc='upper right',
                  frameon=False, fontsize=10)

    # Save figure
    for fmt in formats_to_save:
        output_file = f"{output_prefix}_genome_landscape.{fmt}"
        dpi = 300 if fmt == 'png' else None
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight', format=fmt)
        print(f"  Saved: {output_file}")

    plt.close()

    # ========== FIGURE 3: Statistical comparison heatmap ==========
    fig, ax = plt.subplots(figsize=(10, 8))

    # Create pivot table for heatmap
    pivot = df.pivot_table(values='marker_density',
                          index=['parent', 'region_type'],
                          columns='chromosome',
                          fill_value=0)

    # Create heatmap with annotations
    sns.heatmap(pivot, annot=True, fmt='.0f', cmap='YlOrRd', ax=ax,
               cbar_kws={'label': 'Marker Density (markers/Mb)'},
               linewidths=1, linecolor='white',
               annot_kws={'fontsize': 9, 'fontweight': 'bold'})

    ax.set_xlabel('Chromosome', fontweight='bold')
    ax.set_ylabel('Parent / Region', fontweight='bold')
    ax.set_title('Marker Density Heatmap: Parent × Region × Chromosome',
                fontweight='bold', pad=20)

    # Improve tick labels
    ax.set_xticklabels([f'Chr{c}' for c in pivot.columns], rotation=0)

    plt.tight_layout()

    # Save figure
    for fmt in formats_to_save:
        output_file = f"{output_prefix}_density_heatmap.{fmt}"
        dpi = 300 if fmt == 'png' else None
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight', format=fmt)
        print(f"  Saved: {output_file}")

    plt.close()


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
  # Basic analysis (counts and summary plots)
  python 16-analyze_marker_balance.py \\
      --cenhapmer-dir /path/to/cenhapmers/k31 \\
      --output-prefix results/marker_analysis/k31_balance

  # Full analysis with publication-ready figures (requires BED files)
  python 16-analyze_marker_balance.py \\
      --cenhapmer-dir /path/to/cenhapmers/k31 \\
      --output-prefix results/marker_analysis/k31_balance \\
      --col-bed index/Col-0.bed \\
      --ler-bed index/Ler-0.bed \\
      --publication \\
      --format all

  # Generate outputs:
  #   - *_marker_counts.tsv (raw counts with density)
  #   - *_balance_summary.txt (statistical metrics)
  #   - *_marker_distribution.png (basic plots)
  #   - *_normalization_recommendations.txt (guidance)
  #
  # Publication-ready figures (with --publication):
  #   - *_density_analysis_publication.[png/pdf/svg] (4-panel comprehensive)
  #   - *_genome_landscape.[png/pdf/svg] (genome-wide density map)
  #   - *_density_heatmap.[png/pdf/svg] (statistical heatmap)
        """
    )

    parser.add_argument('--cenhapmer-dir', required=True,
                       help='Directory containing cenhapmer KMC databases')
    parser.add_argument('--output-prefix', required=True,
                       help='Output file prefix for results')
    parser.add_argument('--col-bed', required=False,
                       help='Col-0 BED file with genomic region coordinates (for density calculation)')
    parser.add_argument('--ler-bed', required=False,
                       help='Ler-0 BED file with genomic region coordinates (for density calculation)')
    parser.add_argument('--format', default='png', choices=['png', 'pdf', 'svg', 'all'],
                       help='Output format for plots (default: png)')
    parser.add_argument('--publication', action='store_true',
                       help='Generate publication-ready figures with statistical annotations')

    args = parser.parse_args()

    # Create output directory
    output_dir = os.path.dirname(args.output_prefix)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # 1. Analyze marker database
    df = analyze_marker_database(args.cenhapmer_dir, args.output_prefix,
                                  args.col_bed, args.ler_bed)

    # 2. Calculate balance metrics
    parent_totals, region_totals, chr_totals = calculate_balance_metrics(df, args.output_prefix)

    # 3. Generate visualizations
    plot_marker_distribution(df, args.output_prefix)

    # 4. Generate publication-ready density visualizations if requested
    if args.publication or (args.col_bed and args.ler_bed):
        plot_marker_density_publication(df, args.output_prefix, format=args.format)

    # 5. Generate normalization recommendations
    generate_normalization_recommendations(df, parent_totals, args.output_prefix)

    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
    print(f"\nAll results saved with prefix: {args.output_prefix}")

    if args.publication or (args.col_bed and args.ler_bed):
        print("\nPublication-ready figures generated:")
        formats = ['png', 'pdf', 'svg'] if args.format == 'all' else [args.format]
        for fmt in formats:
            print(f"  - {args.output_prefix}_density_analysis_publication.{fmt}")
            print(f"  - {args.output_prefix}_genome_landscape.{fmt}")
            print(f"  - {args.output_prefix}_density_heatmap.{fmt}")


if __name__ == '__main__':
    main()
