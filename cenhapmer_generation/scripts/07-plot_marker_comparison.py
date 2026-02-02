#!/usr/bin/env python3

"""
Script: 07-plot_marker_comparison.py
Description: Visualize comparison between RAW k-mers and UNIQUE hapmers
             Shows absolute counts and normalized counts (per Mbp)
             Helps understand information loss during uniqueness filtering

Usage: python3 07-plot_marker_comparison.py [options]

Options:
    --marker_dir DIR       Directory with marker density analysis (default: 04-marker_density)
    --output_dir DIR       Output directory for plots (default: 05-marker_plots)
    --parent1 NAME         Name of parent 1 (default: Col-0)
    --parent2 NAME         Name of parent 2 (default: Ler-0)
    --kmer_sizes SIZES     Comma-separated k-mer sizes (default: 21,31,41)
    --skip-cen             Skip centromere mode (chromosome-level hapmers)
    --help                 Show this help message

Requirements:
    - Python 3.6+
    - matplotlib
    - pandas
    - seaborn (optional)

Output:
    - Comparison plots (raw vs unique)
    - Retention rates (% of raw k-mers that are unique)
    - Per-chromosome breakdowns
    - Summary statistics
"""

import os
import sys
import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.filterwarnings('ignore')

# Try to import seaborn for better styling
try:
    import seaborn as sns
    HAS_SEABORN = True
    sns.set_style("whitegrid")
    sns.set_palette("husl")
except ImportError:
    HAS_SEABORN = False


def load_data(marker_dir, kmer_sizes):
    """Load raw and unique k-mer data for all k-mer sizes"""

    data = {}

    for k in kmer_sizes:
        k_dir = os.path.join(marker_dir, f"k{k}")

        raw_file = os.path.join(k_dir, "raw_kmers_summary.csv")
        unique_file = os.path.join(k_dir, "unique_hapmers_summary.csv")

        if os.path.exists(raw_file) and os.path.exists(unique_file):
            raw_df = pd.read_csv(raw_file)
            unique_df = pd.read_csv(unique_file)

            data[k] = {
                'raw': raw_df,
                'unique': unique_df
            }
        else:
            print(f"⚠️  Warning: Missing data for k={k}")

    return data


def calculate_retention_rates(data):
    """Calculate what % of raw k-mers are retained as unique hapmers"""

    retention_data = []

    for k, datasets in data.items():
        raw_df = datasets['raw']
        unique_df = datasets['unique']

        # Merge on common columns
        merged = pd.merge(
            raw_df[['parent', 'chromosome', 'unique_kmers', 'total_kmers']],
            unique_df[['parent', 'chromosome', 'unique_kmers', 'total_kmers']],
            on=['parent', 'chromosome'],
            suffixes=('_raw', '_unique')
        )

        # Calculate retention rates
        merged['unique_retention_pct'] = (merged['unique_kmers_unique'] / merged['unique_kmers_raw']) * 100
        merged['total_retention_pct'] = (merged['total_kmers_unique'] / merged['total_kmers_raw']) * 100
        merged['kmer_size'] = k

        retention_data.append(merged)

    return pd.concat(retention_data, ignore_index=True) if retention_data else pd.DataFrame()


def natural_sort_key(chr_name):
    """Sort chromosomes naturally: chr1, chr2, ..., chr9, chr10, chr11"""
    import re
    # Extract numeric part from chromosome name
    parts = re.split(r'(\d+)', str(chr_name))
    return [int(part) if part.isdigit() else part for part in parts]


def plot_raw_vs_unique_counts(data, output_dir, parent1, parent2):
    """Plot absolute counts: raw vs unique k-mers"""

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Raw K-mers vs Unique Hapmers (Absolute Counts)', fontsize=16, fontweight='bold')

    kmer_sizes = sorted(data.keys())

    for idx, parent in enumerate([parent1, parent2]):
        # Plot unique k-mers
        ax = axes[idx, 0]
        for k in kmer_sizes:
            raw_df = data[k]['raw']
            unique_df = data[k]['unique']

            parent_raw = raw_df[raw_df['parent'] == parent].copy()
            parent_raw = parent_raw.iloc[parent_raw['chromosome'].map(natural_sort_key).argsort()]

            parent_unique = unique_df[unique_df['parent'] == parent].copy()
            parent_unique = parent_unique.iloc[parent_unique['chromosome'].map(natural_sort_key).argsort()]

            if not parent_raw.empty:
                ax.plot(parent_raw['chromosome'], parent_raw['unique_kmers'],
                       marker='o', label=f'k={k} (raw)', linestyle='--', alpha=0.7)
            if not parent_unique.empty:
                ax.plot(parent_unique['chromosome'], parent_unique['unique_kmers'],
                       marker='s', label=f'k={k} (unique)', linewidth=2)

        ax.set_title(f'{parent} - Unique K-mers', fontsize=12, fontweight='bold')
        ax.set_xlabel('Chromosome', fontsize=10)
        ax.set_ylabel('Unique K-mer Count', fontsize=10)
        ax.legend(fontsize=8)
        ax.tick_params(axis='x', rotation=45)
        ax.grid(alpha=0.3)

        # Plot total k-mers
        ax = axes[idx, 1]
        for k in kmer_sizes:
            raw_df = data[k]['raw']
            unique_df = data[k]['unique']

            parent_raw = raw_df[raw_df['parent'] == parent].copy()
            parent_raw = parent_raw.iloc[parent_raw['chromosome'].map(natural_sort_key).argsort()]

            parent_unique = unique_df[unique_df['parent'] == parent].copy()
            parent_unique = parent_unique.iloc[parent_unique['chromosome'].map(natural_sort_key).argsort()]

            if not parent_raw.empty:
                ax.plot(parent_raw['chromosome'], parent_raw['total_kmers'],
                       marker='o', label=f'k={k} (raw)', linestyle='--', alpha=0.7)
            if not parent_unique.empty:
                ax.plot(parent_unique['chromosome'], parent_unique['total_kmers'],
                       marker='s', label=f'k={k} (unique)', linewidth=2)

        ax.set_title(f'{parent} - Total K-mers', fontsize=12, fontweight='bold')
        ax.set_xlabel('Chromosome', fontsize=10)
        ax.set_ylabel('Total K-mer Count', fontsize=10)
        ax.legend(fontsize=8)
        ax.tick_params(axis='x', rotation=45)
        ax.grid(alpha=0.3)

    plt.tight_layout()
    output_file = os.path.join(output_dir, 'raw_vs_unique_absolute_counts.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✅ Saved: {output_file}")
    plt.close()


def plot_normalized_counts(data, output_dir, parent1, parent2):
    """Plot normalized counts: per Mbp"""

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Raw K-mers vs Unique Hapmers (Normalized per Mbp)', fontsize=16, fontweight='bold')

    kmer_sizes = sorted(data.keys())

    for idx, parent in enumerate([parent1, parent2]):
        # Plot unique k-mers per Mbp
        ax = axes[idx, 0]
        for k in kmer_sizes:
            raw_df = data[k]['raw']
            unique_df = data[k]['unique']

            parent_raw = raw_df[raw_df['parent'] == parent].copy()
            parent_raw = parent_raw.iloc[parent_raw['chromosome'].map(natural_sort_key).argsort()]

            parent_unique = unique_df[unique_df['parent'] == parent].copy()
            parent_unique = parent_unique.iloc[parent_unique['chromosome'].map(natural_sort_key).argsort()]

            if not parent_raw.empty:
                ax.plot(parent_raw['chromosome'], parent_raw['unique_per_mbp'],
                       marker='o', label=f'k={k} (raw)', linestyle='--', alpha=0.7)
            if not parent_unique.empty:
                ax.plot(parent_unique['chromosome'], parent_unique['unique_per_mbp'],
                       marker='s', label=f'k={k} (unique)', linewidth=2)

        ax.set_title(f'{parent} - Unique K-mers per Mbp', fontsize=12, fontweight='bold')
        ax.set_xlabel('Chromosome', fontsize=10)
        ax.set_ylabel('Unique K-mers / Mbp', fontsize=10)
        ax.legend(fontsize=8)
        ax.tick_params(axis='x', rotation=45)
        ax.grid(alpha=0.3)

        # Plot total k-mers per Mbp
        ax = axes[idx, 1]
        for k in kmer_sizes:
            raw_df = data[k]['raw']
            unique_df = data[k]['unique']

            parent_raw = raw_df[raw_df['parent'] == parent].copy()
            parent_raw = parent_raw.iloc[parent_raw['chromosome'].map(natural_sort_key).argsort()]

            parent_unique = unique_df[unique_df['parent'] == parent].copy()
            parent_unique = parent_unique.iloc[parent_unique['chromosome'].map(natural_sort_key).argsort()]

            if not parent_raw.empty:
                ax.plot(parent_raw['chromosome'], parent_raw['total_per_mbp'],
                       marker='o', label=f'k={k} (raw)', linestyle='--', alpha=0.7)
            if not parent_unique.empty:
                ax.plot(parent_unique['chromosome'], parent_unique['total_per_mbp'],
                       marker='s', label=f'k={k} (unique)', linewidth=2)

        ax.set_title(f'{parent} - Total K-mers per Mbp', fontsize=12, fontweight='bold')
        ax.set_xlabel('Chromosome', fontsize=10)
        ax.set_ylabel('Total K-mers / Mbp', fontsize=10)
        ax.legend(fontsize=8)
        ax.tick_params(axis='x', rotation=45)
        ax.grid(alpha=0.3)

    plt.tight_layout()
    output_file = os.path.join(output_dir, 'raw_vs_unique_normalized.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✅ Saved: {output_file}")
    plt.close()


def plot_retention_rates(retention_df, output_dir, parent1, parent2):
    """Plot retention rates - what % of raw k-mers are kept as unique"""

    if retention_df.empty:
        print("  ⚠️  No retention data to plot")
        return

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    fig.suptitle('K-mer Retention Rates (% of Raw K-mers Retained as Unique)',
                 fontsize=16, fontweight='bold')

    kmer_sizes = sorted(retention_df['kmer_size'].unique())

    for idx, parent in enumerate([parent1, parent2]):
        ax = axes[idx]
        parent_df = retention_df[retention_df['parent'] == parent]

        for k in kmer_sizes:
            k_df = parent_df[parent_df['kmer_size'] == k].copy()
            k_df = k_df.iloc[k_df['chromosome'].map(natural_sort_key).argsort()]
            if not k_df.empty:
                ax.plot(k_df['chromosome'], k_df['unique_retention_pct'],
                       marker='o', label=f'k={k}', linewidth=2)

        ax.set_title(f'{parent} - Retention Rate', fontsize=12, fontweight='bold')
        ax.set_xlabel('Chromosome', fontsize=10)
        ax.set_ylabel('Retention Rate (%)', fontsize=10)
        ax.set_ylim([0, 100])
        ax.legend(fontsize=9)
        ax.tick_params(axis='x', rotation=45)
        ax.grid(alpha=0.3)

        # Add horizontal line at 50%
        ax.axhline(y=50, color='red', linestyle='--', alpha=0.5, label='50% threshold')

    plt.tight_layout()
    output_file = os.path.join(output_dir, 'retention_rates.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✅ Saved: {output_file}")
    plt.close()


def plot_kmer_size_comparison(data, output_dir):
    """Compare total markers across k-mer sizes"""

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle('K-mer Size Comparison', fontsize=16, fontweight='bold')

    # Collect totals for each k-mer size
    raw_totals = {}
    unique_totals = {}

    for k in sorted(data.keys()):
        raw_totals[k] = data[k]['raw']['unique_kmers'].sum()
        unique_totals[k] = data[k]['unique']['unique_kmers'].sum()

    # Bar plot
    x = np.arange(len(raw_totals))
    width = 0.35

    ax1.bar(x - width/2, list(raw_totals.values()), width,
            label='Raw k-mers', alpha=0.8, color='steelblue')
    ax1.bar(x + width/2, list(unique_totals.values()), width,
            label='Unique hapmers', alpha=0.8, color='coral')

    ax1.set_xlabel('K-mer Size', fontsize=12)
    ax1.set_ylabel('Total Unique K-mers', fontsize=12)
    ax1.set_title('Total Unique K-mers by Size', fontsize=12, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels([f'k={k}' for k in raw_totals.keys()])
    ax1.legend()
    ax1.grid(axis='y', alpha=0.3)

    # Retention rate
    retention_rates = {k: (unique_totals[k] / raw_totals[k] * 100)
                      for k in raw_totals.keys()}

    ax2.bar(range(len(retention_rates)), list(retention_rates.values()),
            color='mediumseagreen', alpha=0.8)
    ax2.set_xlabel('K-mer Size', fontsize=12)
    ax2.set_ylabel('Retention Rate (%)', fontsize=12)
    ax2.set_title('Overall Retention Rate by K-mer Size', fontsize=12, fontweight='bold')
    ax2.set_xticks(range(len(retention_rates)))
    ax2.set_xticklabels([f'k={k}' for k in retention_rates.keys()])
    ax2.set_ylim([0, 100])
    ax2.grid(axis='y', alpha=0.3)

    # Add value labels
    for i, (k, v) in enumerate(retention_rates.items()):
        ax2.text(i, v + 2, f'{v:.1f}%', ha='center', fontsize=10)

    plt.tight_layout()
    output_file = os.path.join(output_dir, 'kmer_size_comparison.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✅ Saved: {output_file}")
    plt.close()


def generate_summary_report(data, retention_df, output_dir, parent1, parent2):
    """Generate text summary report"""

    report_file = os.path.join(output_dir, 'marker_comparison_summary.txt')

    with open(report_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("  Raw K-mers vs Unique Hapmers - Comparison Report\n")
        f.write("="*80 + "\n\n")

        for k in sorted(data.keys()):
            f.write(f"\n{'='*80}\n")
            f.write(f"K-mer Size: {k}\n")
            f.write(f"{'='*80}\n\n")

            raw_df = data[k]['raw']
            unique_df = data[k]['unique']

            # Overall statistics
            f.write("Overall Statistics:\n")
            f.write("-" * 80 + "\n")

            raw_total = raw_df['unique_kmers'].sum()
            unique_total = unique_df['unique_kmers'].sum()
            retention = (unique_total / raw_total * 100) if raw_total > 0 else 0

            f.write(f"  Raw k-mers (unique):      {raw_total:,}\n")
            f.write(f"  Unique hapmers:           {unique_total:,}\n")
            f.write(f"  Retention rate:           {retention:.2f}%\n")
            f.write(f"  K-mers lost:              {raw_total - unique_total:,}\n\n")

            # Per-parent statistics
            f.write("Per-Parent Statistics:\n")
            f.write("-" * 80 + "\n")

            for parent in [parent1, parent2]:
                f.write(f"\n{parent}:\n")

                parent_raw = raw_df[raw_df['parent'] == parent]
                parent_unique = unique_df[unique_df['parent'] == parent]

                raw_sum = parent_raw['unique_kmers'].sum()
                unique_sum = parent_unique['unique_kmers'].sum()
                parent_retention = (unique_sum / raw_sum * 100) if raw_sum > 0 else 0

                f.write(f"  Raw k-mers:      {raw_sum:,}\n")
                f.write(f"  Unique hapmers:  {unique_sum:,}\n")
                f.write(f"  Retention rate:  {parent_retention:.2f}%\n")

                # Normalized values
                if not parent_raw.empty:
                    raw_per_mbp = parent_raw['unique_per_mbp'].mean()
                    unique_per_mbp = parent_unique['unique_per_mbp'].mean()
                    f.write(f"  Avg raw k-mers/Mbp:    {raw_per_mbp:.2f}\n")
                    f.write(f"  Avg unique hapmers/Mbp: {unique_per_mbp:.2f}\n")

        # Retention rate analysis
        if not retention_df.empty:
            f.write("\n" + "="*80 + "\n")
            f.write("Retention Rate Analysis\n")
            f.write("="*80 + "\n\n")

            for k in sorted(retention_df['kmer_size'].unique()):
                k_df = retention_df[retention_df['kmer_size'] == k]

                f.write(f"\nK-mer size {k}:\n")
                f.write(f"  Mean retention:   {k_df['unique_retention_pct'].mean():.2f}%\n")
                f.write(f"  Median retention: {k_df['unique_retention_pct'].median():.2f}%\n")
                f.write(f"  Min retention:    {k_df['unique_retention_pct'].min():.2f}%\n")
                f.write(f"  Max retention:    {k_df['unique_retention_pct'].max():.2f}%\n")

        # Recommendations
        f.write("\n" + "="*80 + "\n")
        f.write("Interpretation:\n")
        f.write("="*80 + "\n\n")
        f.write("Higher retention rates indicate:\n")
        f.write("  ✓ More k-mers are unique to each parent\n")
        f.write("  ✓ Better discriminatory power for haplotyping\n")
        f.write("  ✓ Less overlap between parental genomes\n\n")
        f.write("Lower retention rates indicate:\n")
        f.write("  ✗ Many k-mers are shared between parents\n")
        f.write("  ✗ Higher sequence similarity in those regions\n")
        f.write("  ✗ May need longer k-mers for better discrimination\n\n")
        f.write("Generally:\n")
        f.write("  - Higher k → More unique k-mers (higher specificity)\n")
        f.write("  - Lower k → Fewer unique k-mers (more sharing)\n")
        f.write("  - Retention rate of ~30-50% is typical for divergent accessions\n")

        f.write("\n" + "="*80 + "\n")

    print(f"  ✅ Saved: {report_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Visualize raw k-mers vs unique hapmers comparison",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument('--marker_dir', default='04-marker_density',
                       help='Directory with marker density analysis (default: 04-marker_density)')
    parser.add_argument('--output_dir', default='05-marker_plots',
                       help='Output directory for plots (default: 05-marker_plots)')
    parser.add_argument('--parent1', default='Col-0',
                       help='Name of parent 1 (default: Col-0)')
    parser.add_argument('--parent2', default='Ler-0',
                       help='Name of parent 2 (default: Ler-0)')
    parser.add_argument('--kmer_sizes', default='21,31,41',
                       help='Comma-separated k-mer sizes (default: 21,31,41)')
    parser.add_argument('--skip-cen', action='store_true',
                       help='Skip centromere mode (chromosome-level hapmers)')

    args = parser.parse_args()

    # Parse k-mer sizes
    kmer_sizes = [int(k) for k in args.kmer_sizes.split(',')]

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    print("\n" + "="*80)
    print("  Raw K-mers vs Unique Hapmers - Visualization")
    print("="*80)
    print(f"\nInput directory:  {args.marker_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Parents:          {args.parent1}, {args.parent2}")
    print(f"K-mer sizes:      {', '.join(map(str, kmer_sizes))}")
    print("="*80 + "\n")

    # Load data
    print("📊 Loading marker data...")
    data = load_data(args.marker_dir, kmer_sizes)

    if not data:
        print("\n❌ Error: No marker data found!")
        print("   Please run marker density analysis first (step 05)")
        return 1

    print(f"   Loaded data for {len(data)} k-mer sizes")
    print("")

    # Calculate retention rates
    print("📈 Calculating retention rates...")
    retention_df = calculate_retention_rates(data)
    print(f"   Calculated retention for {len(retention_df)} entries")
    print("")

    # Generate plots
    print("📊 Generating visualizations...")
    plot_raw_vs_unique_counts(data, args.output_dir, args.parent1, args.parent2)
    plot_normalized_counts(data, args.output_dir, args.parent1, args.parent2)
    plot_retention_rates(retention_df, args.output_dir, args.parent1, args.parent2)
    plot_kmer_size_comparison(data, args.output_dir)
    print("")

    # Generate summary report
    print("📝 Generating summary report...")
    generate_summary_report(data, retention_df, args.output_dir, args.parent1, args.parent2)
    print("")

    print("="*80)
    print("✅ Visualization complete!")
    print("="*80)
    print(f"\nOutput directory: {args.output_dir}")
    print("\nGenerated files:")
    print("  - raw_vs_unique_absolute_counts.png")
    print("  - raw_vs_unique_normalized.png")
    print("  - retention_rates.png")
    print("  - kmer_size_comparison.png")
    print("  - marker_comparison_summary.txt")
    print("")

    return 0


if __name__ == "__main__":
    sys.exit(main())
