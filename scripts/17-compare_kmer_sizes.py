#!/usr/bin/env python3
"""
Compare marker balance and density across different k-mer sizes.

This script analyzes cenhapmer databases for multiple k-mer sizes and generates
comparative visualizations to help select the optimal k-mer size for CHARLA analysis.

Key comparisons:
1. Total marker counts vs k-mer size
2. Marker density vs k-mer size
3. Parent balance vs k-mer size
4. ARMS vs CEN ratios across k-mer sizes
5. Chromosome-specific patterns

Outputs:
- Comparative plots showing trends across k-mer sizes
- Statistical summaries
- Recommendations for k-mer size selection
"""

import os
import sys
import argparse
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
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


def count_kmers_in_database(kmc_prefix):
    """Count total k-mers in a KMC database."""
    try:
        cmd = f"kmc_tools info {kmc_prefix}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)

        for line in result.stdout.split('\n'):
            if 'Total number of k-mers' in line or 'Total no. of unique k-mers' in line:
                parts = line.split(':')
                if len(parts) > 1:
                    return int(parts[1].strip())

        # Fallback: count lines in dump
        cmd_dump = f"kmc_tools transform {kmc_prefix} dump /dev/stdout | wc -l"
        result = subprocess.run(cmd_dump, shell=True, capture_output=True, text=True, check=True)
        return int(result.stdout.strip())

    except Exception as e:
        print(f"Warning: Could not count k-mers in {kmc_prefix}: {e}", file=sys.stderr)
        return 0


def load_bed_file_sizes(bed_file):
    """Load genomic region sizes from BED file."""
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
                chr_part = name_parts[0]
                region_part = name_parts[1]
                parent_part = name_parts[2]

                chr_num = chr_part.replace('Chr', '')

                if parent_part == 'Col':
                    parent = 'Col-0'
                elif parent_part == 'Ler':
                    parent = 'Ler-0'
                else:
                    parent = parent_part

                region_type = region_part

                key = (parent, region_type, chr_num)
                region_sizes[key] = region_sizes.get(key, 0) + size_bp

    return region_sizes


def analyze_single_kmer_size(kmer_size, cenhapmer_base_dir, col_bed, ler_bed):
    """Analyze marker database for a single k-mer size."""

    kmer_dir = os.path.join(cenhapmer_base_dir, f'k{kmer_size}')

    if not os.path.exists(kmer_dir):
        print(f"Warning: Directory not found: {kmer_dir}", file=sys.stderr)
        return None

    print(f"\nAnalyzing k={kmer_size}...")

    # Find all KMC databases
    kmc_files = []
    for root, dirs, files in os.walk(kmer_dir):
        for f in files:
            if f.endswith('.kmc_pre'):
                kmc_prefix = os.path.join(root, f.replace('.kmc_pre', ''))
                kmc_files.append(kmc_prefix)

    if not kmc_files:
        print(f"  No KMC databases found in {kmer_dir}", file=sys.stderr)
        return None

    print(f"  Found {len(kmc_files)} databases")

    # Parse filenames and count k-mers
    results = []
    for kmc_file in sorted(kmc_files):
        basename = os.path.basename(kmc_file)
        parts = basename.split('_')

        # Format: unique_Parent_RegionType_ChrN_kSize
        if len(parts) >= 5 and parts[0] == 'unique':
            parent = parts[1]
            region_type_str = parts[2]
            chr_part = parts[3]
            chr_num = chr_part.replace('Chr', '')
            region_type = region_type_str
        else:
            continue

        kmer_count = count_kmers_in_database(kmc_file)

        results.append({
            'kmer_size': kmer_size,
            'parent': parent,
            'region_type': region_type,
            'chromosome': chr_num,
            'kmer_count': kmer_count
        })

    if not results:
        return None

    df = pd.DataFrame(results)

    # Load genomic region sizes
    region_sizes = {}
    if col_bed and os.path.exists(col_bed):
        region_sizes.update(load_bed_file_sizes(col_bed))
    if ler_bed and os.path.exists(ler_bed):
        region_sizes.update(load_bed_file_sizes(ler_bed))

    # Calculate marker density if region sizes available
    if region_sizes:
        df['region_size_bp'] = df.apply(
            lambda row: region_sizes.get((row['parent'], row['region_type'], row['chromosome']), 0),
            axis=1
        )
        df['region_size_mb'] = df['region_size_bp'] / 1_000_000
        df['marker_density'] = df.apply(
            lambda row: row['kmer_count'] / row['region_size_mb'] if row['region_size_mb'] > 0 else 0,
            axis=1
        )
    else:
        df['region_size_bp'] = 0
        df['region_size_mb'] = 0.0
        df['marker_density'] = 0.0

    print(f"  Total markers: {df['kmer_count'].sum():,}")

    return df


def plot_kmer_size_comparison(all_data, output_prefix, format='png'):
    """Create comprehensive k-mer size comparison plots."""

    print("\nGenerating k-mer size comparison plots...")

    # Define colors
    colors = {
        'Col-0': '#2C7BB6',
        'Ler-0': '#D7191C',
        'ARMS': '#66C2A5',
        'CEN': '#FC8D62'
    }

    kmer_sizes = sorted(all_data['kmer_size'].unique())

    # ========== FIGURE 1: Comprehensive 6-panel comparison ==========
    fig = plt.figure(figsize=(18, 12))
    gs = fig.add_gridspec(3, 2, hspace=0.35, wspace=0.25)

    # Panel A: Total markers vs k-mer size (by parent)
    ax1 = fig.add_subplot(gs[0, 0])
    parent_totals = all_data.groupby(['kmer_size', 'parent'])['kmer_count'].sum().unstack()

    x = np.arange(len(kmer_sizes))
    width = 0.35

    for i, parent in enumerate(parent_totals.columns):
        offset = (i - 0.5) * width
        bars = ax1.bar(x + offset, parent_totals[parent], width,
                      label=parent, color=colors.get(parent, '#888'),
                      edgecolor='black', linewidth=1.5, alpha=0.8)

        # Add value labels
        for bar in bars:
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2, height,
                    f'{height/1e6:.1f}M', ha='center', va='bottom',
                    fontsize=8, fontweight='bold')

    ax1.set_xticks(x)
    ax1.set_xticklabels([f'k={k}' for k in kmer_sizes])
    ax1.set_xlabel('K-mer Size', fontweight='bold')
    ax1.set_ylabel('Total Markers', fontweight='bold')
    ax1.set_title('A. Total Marker Count vs K-mer Size', fontweight='bold', loc='left')
    ax1.legend(frameon=False, loc='upper left')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.grid(axis='y', alpha=0.3, linestyle='--')

    # Panel B: Parent balance ratio
    ax2 = fig.add_subplot(gs[0, 1])
    balance_data = []
    for k in kmer_sizes:
        k_data = all_data[all_data['kmer_size'] == k]
        parent_totals_k = k_data.groupby('parent')['kmer_count'].sum()
        if len(parent_totals_k) == 2:
            ratio = parent_totals_k.max() / parent_totals_k.min()
            imbalance_pct = (ratio - 1) * 100
        else:
            imbalance_pct = 0
        balance_data.append(imbalance_pct)

    bars = ax2.bar(x, balance_data, color='#95a5a6', edgecolor='black', linewidth=1.5, alpha=0.8)

    # Color bars based on threshold
    for i, (bar, val) in enumerate(zip(bars, balance_data)):
        if val < 10:
            bar.set_color('#27ae60')  # Green - good
        elif val < 20:
            bar.set_color('#f39c12')  # Orange - caution
        else:
            bar.set_color('#e74c3c')  # Red - warning

        ax2.text(bar.get_x() + bar.get_width()/2, val,
                f'{val:.1f}%', ha='center', va='bottom',
                fontsize=9, fontweight='bold')

    ax2.axhline(y=10, color='orange', linestyle='--', linewidth=2, label='10% threshold')
    ax2.axhline(y=20, color='red', linestyle='--', linewidth=2, label='20% threshold')

    ax2.set_xticks(x)
    ax2.set_xticklabels([f'k={k}' for k in kmer_sizes])
    ax2.set_xlabel('K-mer Size', fontweight='bold')
    ax2.set_ylabel('Parent Imbalance (%)', fontweight='bold')
    ax2.set_title('B. Parent Balance vs K-mer Size', fontweight='bold', loc='left')
    ax2.legend(frameon=False, loc='upper left', fontsize=8)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.grid(axis='y', alpha=0.3, linestyle='--')

    # Panel C: Marker density vs k-mer size (if density available)
    ax3 = fig.add_subplot(gs[1, 0])
    if all_data['marker_density'].sum() > 0:
        density_totals = all_data.groupby(['kmer_size', 'parent'], group_keys=False).apply(
            lambda x: x['kmer_count'].sum() / x['region_size_mb'].sum() if x['region_size_mb'].sum() > 0 else 0,
            include_groups=False
        ).unstack()

        for i, parent in enumerate(density_totals.columns):
            offset = (i - 0.5) * width
            ax3.bar(x + offset, density_totals[parent], width,
                   label=parent, color=colors.get(parent, '#888'),
                   edgecolor='black', linewidth=1.5, alpha=0.8)

        ax3.set_xticks(x)
        ax3.set_xticklabels([f'k={k}' for k in kmer_sizes])
        ax3.set_xlabel('K-mer Size', fontweight='bold')
        ax3.set_ylabel('Marker Density (markers/Mb)', fontweight='bold')
        ax3.set_title('C. Marker Density vs K-mer Size', fontweight='bold', loc='left')
        ax3.legend(frameon=False, loc='upper left')
        ax3.spines['top'].set_visible(False)
        ax3.spines['right'].set_visible(False)
        ax3.grid(axis='y', alpha=0.3, linestyle='--')
    else:
        ax3.text(0.5, 0.5, 'Density data not available\n(provide BED files)',
                ha='center', va='center', transform=ax3.transAxes,
                fontsize=12, style='italic', color='gray')
        ax3.set_xticks([])
        ax3.set_yticks([])
        ax3.spines['top'].set_visible(False)
        ax3.spines['right'].set_visible(False)
        ax3.spines['bottom'].set_visible(False)
        ax3.spines['left'].set_visible(False)

    # Panel D: ARMS vs CEN ratio
    ax4 = fig.add_subplot(gs[1, 1])

    if all_data['marker_density'].sum() > 0:
        ratio_data = []
        for k in kmer_sizes:
            k_data = all_data[all_data['kmer_size'] == k]
            region_density = k_data.groupby(['parent', 'region_type'], group_keys=False).apply(
                lambda x: x['kmer_count'].sum() / x['region_size_mb'].sum() if x['region_size_mb'].sum() > 0 else 0,
                include_groups=False
            ).unstack(fill_value=0)

            if 'ARMS' in region_density.columns and 'CEN' in region_density.columns:
                avg_ratio = (region_density['ARMS'] / region_density['CEN']).mean()
                ratio_data.append(avg_ratio)
            else:
                ratio_data.append(0)

        bars = ax4.bar(x, ratio_data, color=colors['ARMS'], edgecolor='black', linewidth=1.5, alpha=0.8)

        for bar, val in zip(bars, ratio_data):
            ax4.text(bar.get_x() + bar.get_width()/2, val,
                    f'{val:.1f}x', ha='center', va='bottom',
                    fontsize=9, fontweight='bold')

        ax4.set_xticks(x)
        ax4.set_xticklabels([f'k={k}' for k in kmer_sizes])
        ax4.set_xlabel('K-mer Size', fontweight='bold')
        ax4.set_ylabel('ARMS:CEN Density Ratio', fontweight='bold')
        ax4.set_title('D. ARMS vs CEN Density Ratio', fontweight='bold', loc='left')
        ax4.spines['top'].set_visible(False)
        ax4.spines['right'].set_visible(False)
        ax4.grid(axis='y', alpha=0.3, linestyle='--')
    else:
        ax4.text(0.5, 0.5, 'Ratio data not available\n(provide BED files)',
                ha='center', va='center', transform=ax4.transAxes,
                fontsize=12, style='italic', color='gray')
        ax4.set_xticks([])
        ax4.set_yticks([])
        ax4.spines['top'].set_visible(False)
        ax4.spines['right'].set_visible(False)
        ax4.spines['bottom'].set_visible(False)
        ax4.spines['left'].set_visible(False)

    # Panel E: Chromosome variation (CV)
    ax5 = fig.add_subplot(gs[2, 0])
    cv_data = []

    for k in kmer_sizes:
        k_data = all_data[all_data['kmer_size'] == k]
        chr_totals = k_data.groupby(['parent', 'chromosome'])['kmer_count'].sum().unstack(fill_value=0)

        cvs = []
        for parent in chr_totals.index:
            counts = chr_totals.loc[parent].values
            cv = np.std(counts) / np.mean(counts) * 100 if np.mean(counts) > 0 else 0
            cvs.append(cv)

        cv_data.append(np.mean(cvs))

    bars = ax5.bar(x, cv_data, color='#9b59b6', edgecolor='black', linewidth=1.5, alpha=0.8)

    for bar, val in zip(bars, cv_data):
        ax5.text(bar.get_x() + bar.get_width()/2, val,
                f'{val:.1f}%', ha='center', va='bottom',
                fontsize=9, fontweight='bold')

    ax5.axhline(y=30, color='orange', linestyle='--', linewidth=2, label='30% threshold')

    ax5.set_xticks(x)
    ax5.set_xticklabels([f'k={k}' for k in kmer_sizes])
    ax5.set_xlabel('K-mer Size', fontweight='bold')
    ax5.set_ylabel('Coefficient of Variation (%)', fontweight='bold')
    ax5.set_title('E. Chromosome Variation vs K-mer Size', fontweight='bold', loc='left')
    ax5.legend(frameon=False, loc='upper right', fontsize=8)
    ax5.spines['top'].set_visible(False)
    ax5.spines['right'].set_visible(False)
    ax5.grid(axis='y', alpha=0.3, linestyle='--')

    # Panel F: Trend summary (line plot)
    ax6 = fig.add_subplot(gs[2, 1])

    # Normalize all metrics to 0-100 scale for comparison
    def normalize(data):
        data_arr = np.array(data)
        return (data_arr - data_arr.min()) / (data_arr.max() - data_arr.min()) * 100 if data_arr.max() > data_arr.min() else data_arr

    total_markers_norm = normalize([all_data[all_data['kmer_size'] == k]['kmer_count'].sum() for k in kmer_sizes])
    balance_norm = normalize([100 - b for b in balance_data])  # Invert so higher is better
    cv_norm = normalize([100 - c for c in cv_data])  # Invert so higher is better

    ax6.plot(x, total_markers_norm, 'o-', linewidth=2.5, markersize=8,
            label='Total Markers', color='#3498db', markeredgecolor='black', markeredgewidth=1.5)
    ax6.plot(x, balance_norm, 's-', linewidth=2.5, markersize=8,
            label='Parent Balance', color='#27ae60', markeredgecolor='black', markeredgewidth=1.5)
    ax6.plot(x, cv_norm, '^-', linewidth=2.5, markersize=8,
            label='Chr Consistency', color='#9b59b6', markeredgecolor='black', markeredgewidth=1.5)

    ax6.set_xticks(x)
    ax6.set_xticklabels([f'k={k}' for k in kmer_sizes])
    ax6.set_xlabel('K-mer Size', fontweight='bold')
    ax6.set_ylabel('Normalized Score (0-100)', fontweight='bold')
    ax6.set_title('F. Overall Quality Trends', fontweight='bold', loc='left')
    ax6.legend(frameon=False, loc='best')
    ax6.spines['top'].set_visible(False)
    ax6.spines['right'].set_visible(False)
    ax6.grid(axis='both', alpha=0.3, linestyle='--')
    ax6.set_ylim(-5, 105)

    plt.suptitle('K-mer Size Comparison for Cenhapmer Marker Database',
                fontweight='bold', fontsize=15, y=0.995)

    # Save figure
    formats_to_save = ['png', 'pdf', 'svg'] if format == 'all' else [format]

    for fmt in formats_to_save:
        output_file = f"{output_prefix}_kmer_comparison.{fmt}"
        dpi = 300 if fmt == 'png' else None
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight', format=fmt)
        print(f"  Saved: {output_file}")

    plt.close()

    # ========== FIGURE 2: Detailed region comparison ==========
    if all_data['marker_density'].sum() > 0:
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        # Top left: ARMS density by chromosome
        ax = axes[0, 0]
        for parent in all_data['parent'].unique():
            parent_data = all_data[(all_data['parent'] == parent) &
                                  (all_data['region_type'] == 'ARMS')]
            chr_density = parent_data.groupby(['kmer_size', 'chromosome'], group_keys=False).apply(
                lambda x: x['kmer_count'].sum() / x['region_size_mb'].sum() if x['region_size_mb'].sum() > 0 else 0,
                include_groups=False
            ).unstack(fill_value=0)

            for chr_num in sorted(chr_density.columns, key=int):
                ax.plot(chr_density.index, chr_density[chr_num],
                       'o-', linewidth=2, markersize=6,
                       label=f'{parent} Chr{chr_num}',
                       color=colors.get(parent, '#888'), alpha=0.3 + 0.1*int(chr_num))

        ax.set_xlabel('K-mer Size', fontweight='bold')
        ax.set_ylabel('Marker Density (markers/Mb)', fontweight='bold')
        ax.set_title('ARMS Marker Density by Chromosome', fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='both', alpha=0.3, linestyle='--')

        # Top right: CEN density by chromosome
        ax = axes[0, 1]
        for parent in all_data['parent'].unique():
            parent_data = all_data[(all_data['parent'] == parent) &
                                  (all_data['region_type'] == 'CEN')]
            chr_density = parent_data.groupby(['kmer_size', 'chromosome'], group_keys=False).apply(
                lambda x: x['kmer_count'].sum() / x['region_size_mb'].sum() if x['region_size_mb'].sum() > 0 else 0,
                include_groups=False
            ).unstack(fill_value=0)

            for chr_num in sorted(chr_density.columns, key=int):
                ax.plot(chr_density.index, chr_density[chr_num],
                       's-', linewidth=2, markersize=6,
                       label=f'{parent} Chr{chr_num}',
                       color=colors.get(parent, '#888'), alpha=0.3 + 0.1*int(chr_num))

        ax.set_xlabel('K-mer Size', fontweight='bold')
        ax.set_ylabel('Marker Density (markers/Mb)', fontweight='bold')
        ax.set_title('CEN Marker Density by Chromosome', fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='both', alpha=0.3, linestyle='--')

        # Bottom: Heatmap of density by k-mer size
        ax = axes[1, 0]
        heatmap_data = []
        for k in kmer_sizes:
            k_data = all_data[all_data['kmer_size'] == k]
            for parent in sorted(k_data['parent'].unique()):
                for region in ['ARMS', 'CEN']:
                    subset = k_data[(k_data['parent'] == parent) &
                                   (k_data['region_type'] == region)]
                    if len(subset) > 0:
                        density = subset['kmer_count'].sum() / subset['region_size_mb'].sum() if subset['region_size_mb'].sum() > 0 else 0
                        heatmap_data.append({
                            'K-mer Size': k,
                            'Group': f'{parent}\n{region}',
                            'Density': density
                        })

        heatmap_df = pd.DataFrame(heatmap_data)
        pivot = heatmap_df.pivot(index='Group', columns='K-mer Size', values='Density')

        sns.heatmap(pivot, annot=True, fmt='.0f', cmap='YlOrRd', ax=ax,
                   cbar_kws={'label': 'Marker Density (markers/Mb)'},
                   linewidths=1, linecolor='white',
                   annot_kws={'fontsize': 9, 'fontweight': 'bold'})

        ax.set_xlabel('K-mer Size', fontweight='bold')
        ax.set_ylabel('Parent / Region', fontweight='bold')
        ax.set_title('Marker Density Heatmap Across K-mer Sizes', fontweight='bold')

        # Bottom right: Summary statistics
        ax = axes[1, 1]
        ax.axis('off')

        summary_text = "K-MER SIZE SELECTION GUIDE\n" + "="*40 + "\n\n"

        # Find best k-mer size
        best_balance_idx = np.argmin(balance_data)
        best_cv_idx = np.argmin(cv_data)
        max_markers_idx = np.argmax([all_data[all_data['kmer_size'] == k]['kmer_count'].sum() for k in kmer_sizes])

        summary_text += f"Best parent balance: k={kmer_sizes[best_balance_idx]}\n"
        summary_text += f"  ({balance_data[best_balance_idx]:.1f}% imbalance)\n\n"

        summary_text += f"Most chromosome consistency: k={kmer_sizes[best_cv_idx]}\n"
        summary_text += f"  ({cv_data[best_cv_idx]:.1f}% CV)\n\n"

        summary_text += f"Maximum markers: k={kmer_sizes[max_markers_idx]}\n"
        summary_text += f"  ({all_data[all_data['kmer_size'] == kmer_sizes[max_markers_idx]]['kmer_count'].sum()/1e6:.1f}M markers)\n\n"

        summary_text += "\nRECOMMENDATIONS:\n" + "-"*40 + "\n\n"

        # Recommendation logic
        if all([b < 10 for b in balance_data]):
            summary_text += "✓ All k-mer sizes are well-balanced\n\n"
            summary_text += "For high-quality reads:\n"
            summary_text += f"  → Use k={kmer_sizes[-1]} (more markers)\n\n"
            summary_text += "For noisy reads:\n"
            summary_text += f"  → Use k={kmer_sizes[0]} (error tolerant)\n\n"
            summary_text += "Balanced choice:\n"
            summary_text += f"  → Use k={kmer_sizes[len(kmer_sizes)//2]}\n"
        else:
            good_k = [k for k, b in zip(kmer_sizes, balance_data) if b < 10]
            if good_k:
                summary_text += f"Recommended k-mer sizes: {', '.join(map(str, good_k))}\n"
                summary_text += f"  (< 10% parent imbalance)\n"
            else:
                summary_text += "⚠️ All k-mer sizes show imbalance\n"
                summary_text += f"  Best option: k={kmer_sizes[best_balance_idx]}\n"

        ax.text(0.05, 0.95, summary_text, transform=ax.transAxes,
               fontsize=9, verticalalignment='top', fontfamily='monospace',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

        plt.suptitle('Detailed K-mer Size Analysis by Genomic Region',
                    fontweight='bold', fontsize=14, y=0.995)

        plt.tight_layout()

        for fmt in formats_to_save:
            output_file = f"{output_prefix}_kmer_detailed.{fmt}"
            dpi = 300 if fmt == 'png' else None
            plt.savefig(output_file, dpi=dpi, bbox_inches='tight', format=fmt)
            print(f"  Saved: {output_file}")

        plt.close()


def save_comparison_summary(all_data, output_prefix):
    """Save text summary of k-mer size comparison."""

    kmer_sizes = sorted(all_data['kmer_size'].unique())

    summary_file = f"{output_prefix}_kmer_summary.txt"

    with open(summary_file, 'w') as f:
        f.write("K-MER SIZE COMPARISON SUMMARY\n")
        f.write("="*70 + "\n\n")

        for k in kmer_sizes:
            f.write(f"\nK-MER SIZE: {k}\n")
            f.write("-"*70 + "\n")

            k_data = all_data[all_data['kmer_size'] == k]

            # Parent totals
            parent_totals = k_data.groupby('parent')['kmer_count'].sum()
            f.write("\nParent Totals:\n")
            for parent, count in parent_totals.items():
                f.write(f"  {parent}: {count:,} markers\n")

            if len(parent_totals) == 2:
                ratio = parent_totals.max() / parent_totals.min()
                imbalance = (ratio - 1) * 100
                f.write(f"\n  Imbalance: {imbalance:.1f}% ({ratio:.2f}:1)\n")

                if imbalance < 10:
                    f.write("  Status: ✓ Well-balanced\n")
                elif imbalance < 20:
                    f.write("  Status: ⚡ Moderate imbalance\n")
                else:
                    f.write("  Status: ⚠️  High imbalance\n")

            # Region breakdown
            f.write("\nRegion Breakdown:\n")
            region_totals = k_data.groupby(['parent', 'region_type'])['kmer_count'].sum().unstack(fill_value=0)
            f.write(region_totals.to_string())
            f.write("\n")

            # Density if available
            if k_data['marker_density'].sum() > 0:
                f.write("\nMarker Density (markers/Mb):\n")
                parent_density = k_data.groupby('parent', group_keys=False).apply(
                    lambda x: x['kmer_count'].sum() / x['region_size_mb'].sum() if x['region_size_mb'].sum() > 0 else 0,
                    include_groups=False
                )
                for parent, density in parent_density.items():
                    f.write(f"  {parent}: {density:,.1f}\n")

        f.write("\n\n" + "="*70 + "\n")
        f.write("OVERALL RECOMMENDATIONS\n")
        f.write("="*70 + "\n\n")

        # Find best overall k-mer size
        balance_scores = []
        for k in kmer_sizes:
            k_data = all_data[all_data['kmer_size'] == k]
            parent_totals = k_data.groupby('parent')['kmer_count'].sum()
            if len(parent_totals) == 2:
                ratio = parent_totals.max() / parent_totals.min()
                imbalance = (ratio - 1) * 100
            else:
                imbalance = 0
            balance_scores.append(imbalance)

        best_idx = np.argmin(balance_scores)

        f.write(f"Best balanced k-mer size: k={kmer_sizes[best_idx]}\n")
        f.write(f"  (Lowest parent imbalance: {balance_scores[best_idx]:.1f}%)\n\n")

        f.write("General guidelines:\n")
        f.write("  - Higher k: More unique markers, less error-tolerant\n")
        f.write("  - Lower k: Fewer markers, more robust to sequencing errors\n")
        f.write("  - For long reads (ONT): Consider k=21 or k=25\n")
        f.write("  - For high-quality reads: Consider k=31 or higher\n")

    print(f"\nSaved comparison summary: {summary_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Compare marker balance across different k-mer sizes',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Compare all k-mer sizes found in cenhapmer directory
  python 17-compare_kmer_sizes.py \\
      --cenhapmer-dir /path/to/cenhapmers \\
      --output-prefix marker_analysis/kmer_comparison \\
      --col-bed index/Col-0.bed \\
      --ler-bed index/Ler-0.bed \\
      --format all

  # Compare specific k-mer sizes only
  python 17-compare_kmer_sizes.py \\
      --cenhapmer-dir /path/to/cenhapmers \\
      --kmer-sizes 21 31 41 \\
      --output-prefix marker_analysis/kmer_comparison \\
      --col-bed index/Col-0.bed \\
      --ler-bed index/Ler-0.bed

  # This will generate:
  #   - kmer_comparison_kmer_comparison.[png/pdf/svg] (6-panel overview)
  #   - kmer_comparison_kmer_detailed.[png/pdf/svg] (detailed analysis)
  #   - kmer_comparison_kmer_summary.txt (text summary)
  #   - kmer_comparison_all_data.tsv (raw data)
        """
    )

    parser.add_argument('--cenhapmer-dir', required=True,
                       help='Base directory containing k-mer subdirectories (k21, k31, etc.)')
    parser.add_argument('--kmer-sizes', nargs='+', type=int,
                       help='Specific k-mer sizes to compare (default: auto-detect all)')
    parser.add_argument('--output-prefix', required=True,
                       help='Output file prefix for results')
    parser.add_argument('--col-bed', required=False,
                       help='Col-0 BED file with genomic regions (for density calculation)')
    parser.add_argument('--ler-bed', required=False,
                       help='Ler-0 BED file with genomic regions (for density calculation)')
    parser.add_argument('--format', default='png', choices=['png', 'pdf', 'svg', 'all'],
                       help='Output format for plots (default: png)')

    args = parser.parse_args()

    # Create output directory
    output_dir = os.path.dirname(args.output_prefix)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # Auto-detect k-mer sizes if not specified
    if args.kmer_sizes:
        kmer_sizes = sorted(args.kmer_sizes)
    else:
        print("Auto-detecting k-mer sizes...")
        kmer_sizes = []
        for item in os.listdir(args.cenhapmer_dir):
            item_path = os.path.join(args.cenhapmer_dir, item)
            if os.path.isdir(item_path) and item.startswith('k'):
                try:
                    k = int(item[1:])
                    kmer_sizes.append(k)
                except ValueError:
                    continue
        kmer_sizes = sorted(kmer_sizes)

    if not kmer_sizes:
        print("Error: No k-mer sizes found or specified", file=sys.stderr)
        sys.exit(1)

    print(f"Comparing k-mer sizes: {', '.join(map(str, kmer_sizes))}")

    # Analyze each k-mer size
    all_data_frames = []

    for k in kmer_sizes:
        df = analyze_single_kmer_size(k, args.cenhapmer_dir, args.col_bed, args.ler_bed)
        if df is not None:
            all_data_frames.append(df)

    if not all_data_frames:
        print("Error: No data could be analyzed", file=sys.stderr)
        sys.exit(1)

    # Combine all data
    all_data = pd.concat(all_data_frames, ignore_index=True)

    # Save combined data
    data_file = f"{args.output_prefix}_all_data.tsv"
    all_data.to_csv(data_file, sep='\t', index=False)
    print(f"\nSaved combined data: {data_file}")

    # Generate comparison plots
    plot_kmer_size_comparison(all_data, args.output_prefix, format=args.format)

    # Save text summary
    save_comparison_summary(all_data, args.output_prefix)

    print("\n" + "="*70)
    print("K-MER SIZE COMPARISON COMPLETE")
    print("="*70)
    print(f"\nAll results saved with prefix: {args.output_prefix}")


if __name__ == '__main__':
    main()
