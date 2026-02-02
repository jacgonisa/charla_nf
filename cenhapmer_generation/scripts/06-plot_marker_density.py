#!/usr/bin/env python3

"""
Script: 06-plot_marker_density.py
Description: Visualize marker (cenhapmer) density across chromosomes to help
             determine optimal k-mer size for CHARLA pipeline

Usage: python3 06-plot_marker_density.py [options]

Options:
    --marker_dir DIR       Directory with marker density analysis (default: 04-marker_density)
    --output_dir DIR       Output directory for plots (default: 05-marker_plots)
    --parent1 NAME         Name of parent 1 (default: Col-0)
    --parent2 NAME         Name of parent 2 (default: Ler-0)
    --kmer_sizes SIZES     Comma-separated k-mer sizes (default: 21,31,41)
    --num_chr N            Number of chromosomes (default: 5)
    --chr_prefix PREFIX    Chromosome prefix (default: Chr)
    --window_size N        Window size for binning in bp (default: 10000)
    --skip-cen             Skip centromere mode (chromosome-level hapmers)
    --help                 Show this help message

Requirements:
    - Python 3.6+
    - matplotlib
    - pandas
    - seaborn (optional, for better styling)

Output:
    - Marker density plots for each k-mer size
    - Comparison plots across k-mer sizes
    - Summary statistics
"""

import os
import sys
import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# Try to import seaborn for better styling
try:
    import seaborn as sns
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False


def read_histogram_file(hist_file):
    """Read KMC histogram file and return total marker count"""
    if not os.path.exists(hist_file):
        return 0

    try:
        total = 0
        with open(hist_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    count = int(parts[1])
                    total += count
        return total
    except Exception as e:
        print(f"  Warning: Error reading {hist_file}: {e}")
        return 0


def collect_marker_statistics(marker_dir, parent1, parent2, kmer_sizes, num_chr, chr_prefix, skip_cen=False):
    """Collect marker statistics from all cenhapmer databases"""

    results = []

    # Determine regions based on mode
    regions = ['CHR'] if skip_cen else ['CEN', 'ARMS']

    for k in kmer_sizes:
        k_dir = os.path.join(marker_dir, f"k{k}")

        if not os.path.exists(k_dir):
            print(f"Warning: Directory not found: {k_dir}")
            continue

        # Process each parent
        for parent in [parent1, parent2]:
            # Process each chromosome
            for chr_num in range(1, num_chr + 1):
                chr_name = f"{chr_prefix}{chr_num}"

                # Process regions (CEN/ARMS or CHR)
                for region in regions:
                    base_name = f"unique_{parent}_{region}_{chr_name}_k{k}"
                    hist_file = os.path.join(k_dir, f"{base_name}.counts.hist")

                    if os.path.exists(hist_file):
                        marker_count = read_histogram_file(hist_file)

                        results.append({
                            'kmer_size': k,
                            'parent': parent,
                            'chromosome': chr_name,
                            'region': region,
                            'marker_count': marker_count,
                            'database': base_name
                        })

    return pd.DataFrame(results)


def plot_marker_counts_per_ksize(df, output_dir, parent1, parent2):
    """Plot marker counts grouped by k-mer size"""

    if df.empty:
        print("Warning: No data to plot")
        return

    # Set style
    if HAS_SEABORN:
        sns.set_style("whitegrid")

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Plot for each parent
    for idx, parent in enumerate([parent1, parent2]):
        ax = axes[idx]
        parent_df = df[df['parent'] == parent]

        if parent_df.empty:
            continue

        # Group by kmer_size and chromosome
        pivot_data = parent_df.pivot_table(
            values='marker_count',
            index='chromosome',
            columns='kmer_size',
            aggfunc='sum',
            fill_value=0
        )

        # Plot grouped bar chart
        pivot_data.plot(kind='bar', ax=ax, width=0.8)
        ax.set_title(f'{parent} - Marker Counts by K-mer Size', fontsize=14, fontweight='bold')
        ax.set_xlabel('Chromosome', fontsize=12)
        ax.set_ylabel('Total Marker Count', fontsize=12)
        ax.legend(title='K-mer Size', fontsize=10)
        ax.tick_params(axis='x', rotation=45)
        ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    output_file = os.path.join(output_dir, 'marker_counts_per_chromosome.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✅ Saved: {output_file}")
    plt.close()


def plot_region_comparison(df, output_dir, parent1, parent2):
    """Plot CEN vs ARMS marker distribution"""

    if df.empty:
        return

    fig, axes = plt.subplots(len(df['kmer_size'].unique()), 2,
                            figsize=(14, 4 * len(df['kmer_size'].unique())))

    if len(df['kmer_size'].unique()) == 1:
        axes = axes.reshape(1, -1)

    for row_idx, k in enumerate(sorted(df['kmer_size'].unique())):
        k_df = df[df['kmer_size'] == k]

        for col_idx, parent in enumerate([parent1, parent2]):
            ax = axes[row_idx, col_idx]
            parent_df = k_df[k_df['parent'] == parent]

            if parent_df.empty:
                continue

            # Pivot for CEN vs ARMS
            pivot_data = parent_df.pivot_table(
                values='marker_count',
                index='chromosome',
                columns='region',
                aggfunc='sum',
                fill_value=0
            )

            pivot_data.plot(kind='bar', ax=ax, width=0.7, color=['#2E86AB', '#A23B72'])
            ax.set_title(f'{parent} - k={k}: CEN vs ARMS Markers', fontsize=12, fontweight='bold')
            ax.set_xlabel('Chromosome', fontsize=10)
            ax.set_ylabel('Marker Count', fontsize=10)
            ax.legend(title='Region', fontsize=9)
            ax.tick_params(axis='x', rotation=45)
            ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    output_file = os.path.join(output_dir, 'marker_counts_cen_vs_arms.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✅ Saved: {output_file}")
    plt.close()


def plot_ksize_comparison(df, output_dir):
    """Plot comparison across k-mer sizes"""

    if df.empty:
        return

    # Calculate total markers per k-mer size
    ksize_totals = df.groupby('kmer_size')['marker_count'].sum().sort_index()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Bar plot of totals
    ksize_totals.plot(kind='bar', ax=ax1, color='steelblue', width=0.6)
    ax1.set_title('Total Markers by K-mer Size', fontsize=14, fontweight='bold')
    ax1.set_xlabel('K-mer Size', fontsize=12)
    ax1.set_ylabel('Total Marker Count', fontsize=12)
    ax1.tick_params(axis='x', rotation=0)
    ax1.grid(axis='y', alpha=0.3)

    # Add value labels on bars
    for i, v in enumerate(ksize_totals):
        ax1.text(i, v, f'{int(v):,}', ha='center', va='bottom', fontsize=10)

    # Pie chart showing distribution
    colors = plt.cm.Set3(range(len(ksize_totals)))
    ax2.pie(ksize_totals, labels=[f'k={k}' for k in ksize_totals.index],
            autopct='%1.1f%%', colors=colors, startangle=90)
    ax2.set_title('Marker Distribution by K-mer Size', fontsize=14, fontweight='bold')

    plt.tight_layout()
    output_file = os.path.join(output_dir, 'kmer_size_comparison.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✅ Saved: {output_file}")
    plt.close()


def generate_summary_report(df, output_dir):
    """Generate a text summary report"""

    if df.empty:
        return

    report_file = os.path.join(output_dir, 'marker_density_summary.txt')

    with open(report_file, 'w') as f:
        f.write("="*70 + "\n")
        f.write("  CHARLA Cenhapmer Marker Density Analysis - Summary Report\n")
        f.write("="*70 + "\n\n")

        # Overall statistics
        f.write("Overall Statistics:\n")
        f.write("-" * 70 + "\n")
        for k in sorted(df['kmer_size'].unique()):
            k_df = df[df['kmer_size'] == k]
            total = k_df['marker_count'].sum()
            mean = k_df['marker_count'].mean()
            std = k_df['marker_count'].std()

            f.write(f"\nK-mer size: {k}\n")
            f.write(f"  Total markers:   {total:,}\n")
            f.write(f"  Mean per region: {mean:,.0f}\n")
            f.write(f"  Std deviation:   {std:,.0f}\n")

        # Per-parent statistics
        f.write("\n" + "="*70 + "\n")
        f.write("Per-Parent Statistics:\n")
        f.write("-" * 70 + "\n")

        for parent in df['parent'].unique():
            f.write(f"\n{parent}:\n")
            parent_df = df[df['parent'] == parent]

            for k in sorted(parent_df['kmer_size'].unique()):
                k_df = parent_df[parent_df['kmer_size'] == k]
                total = k_df['marker_count'].sum()
                f.write(f"  k={k}: {total:,} markers\n")

        # Region-specific statistics
        f.write("\n" + "="*70 + "\n")
        f.write("CEN vs ARMS Distribution:\n")
        f.write("-" * 70 + "\n")

        for k in sorted(df['kmer_size'].unique()):
            f.write(f"\nK-mer size: {k}\n")
            k_df = df[df['kmer_size'] == k]

            for region in ['CEN', 'ARMS']:
                region_total = k_df[k_df['region'] == region]['marker_count'].sum()
                percentage = (region_total / k_df['marker_count'].sum()) * 100
                f.write(f"  {region}: {region_total:,} ({percentage:.1f}%)\n")

        # Recommendations
        f.write("\n" + "="*70 + "\n")
        f.write("Recommendations:\n")
        f.write("-" * 70 + "\n")

        # Find k-mer size with highest total markers
        ksize_totals = df.groupby('kmer_size')['marker_count'].sum()
        best_k = ksize_totals.idxmax()

        f.write(f"\nHighest marker count: k={best_k}\n")
        f.write(f"  Total markers: {ksize_totals[best_k]:,}\n")
        f.write("\nK-mer Size Tradeoffs:\n")
        f.write("  - Higher k: MORE unique markers, higher specificity\n")
        f.write("              BUT less resilient to sequencing errors\n")
        f.write("  - Lower k:  FEWER unique markers, less specific\n")
        f.write("              BUT more tolerant of read errors (better for noisy reads)\n")
        f.write("\nKey Insight:\n")
        f.write("  A single sequencing error breaks a k=41 marker but may preserve\n")
        f.write("  multiple k=21 markers from the same region.\n")
        f.write("\nRecommendations:\n")
        f.write("  - k=21: Best for noisy long reads (ONT, PacBio CLR)\n")
        f.write("  - k=31: Good for high-quality reads (PacBio HiFi)\n")
        f.write("  - k=41: Only for very high-quality data\n")

        f.write("\n" + "="*70 + "\n")

    print(f"  ✅ Saved: {report_file}")

    # Also create a CSV summary
    summary_csv = os.path.join(output_dir, 'marker_density_summary.csv')
    summary_df = df.groupby(['kmer_size', 'parent', 'region']).agg({
        'marker_count': ['sum', 'mean', 'std', 'min', 'max']
    }).reset_index()
    summary_df.to_csv(summary_csv, index=False)
    print(f"  ✅ Saved: {summary_csv}")


def main():
    parser = argparse.ArgumentParser(
        description="Visualize marker density to determine optimal k-mer size",
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
    parser.add_argument('--num_chr', type=int, default=5,
                       help='Number of chromosomes (default: 5)')
    parser.add_argument('--chr_prefix', default='Chr',
                       help='Chromosome prefix (default: Chr)')
    parser.add_argument('--skip-cen', action='store_true',
                       help='Skip centromere mode (chromosome-level hapmers)')

    args = parser.parse_args()

    # Parse k-mer sizes
    kmer_sizes = [int(k) for k in args.kmer_sizes.split(',')]

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    print("\n" + "="*70)
    print("  Marker Density Visualization")
    print("="*70)
    print(f"\nInput directory:  {args.marker_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Parents:          {args.parent1}, {args.parent2}")
    print(f"K-mer sizes:      {', '.join(map(str, kmer_sizes))}")
    print(f"Skip centromeres: {args.skip_cen}")
    print("="*70 + "\n")

    # Collect marker statistics
    print("📊 Collecting marker statistics...")
    df = collect_marker_statistics(
        args.marker_dir,
        args.parent1,
        args.parent2,
        kmer_sizes,
        args.num_chr,
        args.chr_prefix,
        skip_cen=args.skip_cen
    )

    if df.empty:
        print("\n❌ Error: No marker data found!")
        print("   Please run marker density analysis first (step 05)")
        return 1

    print(f"   Collected {len(df)} records")
    print("")

    # Generate plots
    print("📈 Generating visualizations...")
    plot_marker_counts_per_ksize(df, args.output_dir, args.parent1, args.parent2)
    plot_region_comparison(df, args.output_dir, args.parent1, args.parent2)
    plot_ksize_comparison(df, args.output_dir)
    print("")

    # Generate summary report
    print("📝 Generating summary report...")
    generate_summary_report(df, args.output_dir)
    print("")

    print("="*70)
    print("✅ Marker density analysis complete!")
    print("="*70)
    print(f"\nOutput directory: {args.output_dir}")
    print("\nGenerated files:")
    print("  - marker_counts_per_chromosome.png")
    print("  - marker_counts_cen_vs_arms.png")
    print("  - kmer_size_comparison.png")
    print("  - marker_density_summary.txt")
    print("  - marker_density_summary.csv")
    print("")

    return 0


if __name__ == "__main__":
    sys.exit(main())
