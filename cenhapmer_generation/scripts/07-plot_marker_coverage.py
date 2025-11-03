#!/usr/bin/env python3
"""
Optimized script to generate genome-wide k-mer coverage plots from marker count CSV files.
Handles large datasets by intelligently downsampling for visualization.
Supports plotting multiple samples on the same plot for comparison.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import sys
from pathlib import Path


def get_chromosome_lengths_from_fasta(fasta_file):
    """
    Extract chromosome lengths from a FASTA file.
    
    Returns:
    - Dictionary mapping chromosome name to length
    """
    chrom_lengths = {}
    current_chrom = None
    current_length = 0
    
    try:
        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Save previous chromosome
                    if current_chrom:
                        chrom_lengths[current_chrom] = current_length
                    # Start new chromosome
                    current_chrom = line[1:].split()[0]  # Get first part of header
                    current_length = 0
                else:
                    current_length += len(line)
            
            # Save last chromosome
            if current_chrom:
                chrom_lengths[current_chrom] = current_length
                
        return chrom_lengths
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        sys.exit(1)


def load_marker_data(csv_files, sample_names=None, filter_samples=None):
    """
    Load marker count data from one or more CSV files.
    
    Parameters:
    - csv_files: List of CSV file paths
    - sample_names: Optional list of sample names to override CSV sample column
    - filter_samples: Optional list of sample names to filter from CSV (when samples are in the file)
    
    Returns:
    - DataFrame with all data, including 'sample_label' column for plotting
    """
    dfs = []
    
    for idx, csv_file in enumerate(csv_files):
        try:
            df = pd.read_csv(csv_file)
            required_cols = ['chrom', 'window_start', 'window_end', 'marker_count']
            
            if not all(col in df.columns for col in required_cols):
                print(f"Error: CSV must contain columns: {required_cols}")
                sys.exit(1)
            
            # Check if we should filter by samples from within the CSV
            if filter_samples and 'sample' in df.columns:
                # Filter to only specified samples
                df = df[df['sample'].isin(filter_samples)].copy()
                if len(df) == 0:
                    print(f"Warning: No data found for samples {filter_samples} in {csv_file}")
                    continue
                df['sample_label'] = df['sample']
                samples_found = df['sample'].unique()
                print(f"  Loaded {len(df):,} windows from {csv_file}")
                print(f"    Samples: {', '.join(samples_found)}")
            # Add sample label from command line (when each file is a different sample)
            elif sample_names and idx < len(sample_names):
                df['sample_label'] = sample_names[idx]
                print(f"  Loaded {len(df):,} windows from {csv_file} (sample: {sample_names[idx]})")
            # Use sample column from CSV if available
            elif 'sample' in df.columns:
                df['sample_label'] = df['sample']
                samples_found = df['sample'].unique()
                print(f"  Loaded {len(df):,} windows from {csv_file}")
                print(f"    Samples in file: {', '.join(samples_found)}")
            else:
                df['sample_label'] = f'Sample_{idx+1}'
                print(f"  Loaded {len(df):,} windows from {csv_file} (sample: Sample_{idx+1})")
            
            dfs.append(df)
            
        except Exception as e:
            print(f"Error loading CSV file {csv_file}: {e}")
            sys.exit(1)
    
    if not dfs:
        print("Error: No data loaded!")
        sys.exit(1)
    
    return pd.concat(dfs, ignore_index=True)


def downsample_chromosome_data(chrom_data, max_points=10000):
    """
    Downsample chromosome data for visualization while preserving peaks.
    Uses max pooling to preserve high values.
    """
    n_windows = len(chrom_data)
    
    if n_windows <= max_points:
        return chrom_data['window_start'].values, chrom_data['marker_count'].values
    
    bin_size = int(np.ceil(n_windows / max_points))
    x_coords = chrom_data['window_start'].values
    y_vals = chrom_data['marker_count'].values
    
    n_bins = int(np.ceil(n_windows / bin_size))
    x_binned = np.zeros(n_bins)
    y_binned = np.zeros(n_bins)
    
    for i in range(n_bins):
        start_idx = i * bin_size
        end_idx = min((i + 1) * bin_size, n_windows)
        x_binned[i] = x_coords[start_idx:end_idx].mean()
        y_binned[i] = y_vals[start_idx:end_idx].max()
    
    return x_binned, y_binned


def plot_chromosome_coverage(df, chrom_lengths, output_file, log_scale=False, max_points_per_chrom=10000, ncols=5):
    """
    Create per-chromosome coverage plots showing marker density.
    Each subplot shows that chromosome's k-mer markers ONLY within that chromosome's boundaries.

    For each chromosome:
    - Filter to windows within that chromosome's genomic coordinates
    - Plot those windows using their genome-wide coordinates
    - Set x-axis to show only that chromosome's range (0 to chrom_length)
    """
    # Get unique samples and chromosomes
    samples = sorted(df['sample_label'].unique())

    # Sort chromosomes numerically (chr1, chr2, ..., chr9, chr10)
    def chrom_sort_key(chrom_name):
        """Extract numeric part for proper sorting of chromosomes."""
        import re
        # Remove 'chr' or 'Chr' prefix
        chrom_clean = chrom_name.replace('chr', '').replace('Chr', '')
        # Try to convert to int, otherwise use string
        try:
            return (0, int(chrom_clean))  # Numeric chromosomes first
        except ValueError:
            return (1, chrom_clean)  # Non-numeric chromosomes after

    chroms = sorted(df['chrom'].unique(), key=chrom_sort_key)

    n_samples = len(samples)
    n_chroms = len(chroms)

    print(f"Plotting {n_samples} sample(s) across {n_chroms} chromosomes")
    print(f"Downsampling to max {max_points_per_chrom:,} points per chromosome")

    # Calculate cumulative chromosome positions in the genome
    chrom_positions = {}
    cumulative = 0

    # Use the same sorted order for cumulative positions
    for chrom in chroms:
        if chrom in chrom_lengths:
            chrom_positions[chrom] = {
                'start': cumulative,
                'end': cumulative + chrom_lengths[chrom],
                'length': chrom_lengths[chrom]
            }
            cumulative += chrom_lengths[chrom]

    print("\nChromosome positions in genome:")
    for chrom in chroms:
        if chrom in chrom_positions:
            pos = chrom_positions[chrom]
            print(f"  {chrom}: {pos['start']:,} - {pos['end']:,} bp (length: {pos['length']:,} bp)")

    # Define colors for multiple samples
    colors = ['#4A90E2', '#E24A4A', '#4AE290', '#E2904A', '#904AE2']

    # Calculate layout: multiple rows if needed
    nrows = int(np.ceil(n_chroms / ncols))
    fig_width = min(4 * ncols, 20)  # Cap at 20 inches wide
    fig_height = 4 * nrows

    print(f"Layout: {nrows} rows x {ncols} columns")

    # Create figure with subplots for each chromosome
    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_width, fig_height), sharey=True)

    # Flatten axes array for easier indexing
    if nrows == 1 and ncols == 1:
        axes = np.array([axes])
    elif nrows == 1 or ncols == 1:
        axes = axes.flatten()
    else:
        axes = axes.flatten()
    
    for idx, chrom in enumerate(chroms):
        ax = axes[idx]
        print(f"\n  Processing {chrom}:")

        # Get the chromosome boundaries
        if chrom not in chrom_positions:
            print(f"    Warning: {chrom} not found in genome file, skipping...")
            continue

        chrom_start = chrom_positions[chrom]['start']
        chrom_end = chrom_positions[chrom]['end']
        chrom_length = chrom_positions[chrom]['length']

        print(f"    Chromosome range: {chrom_start:,} - {chrom_end:,} bp (length: {chrom_length:,} bp)")

        # Plot each sample
        for sample_idx, sample in enumerate(samples):
            # Filter data for:
            # 1. This chromosome's k-mers (chrom column)
            # 2. Windows that fall within this chromosome's genomic boundaries
            mask = (
                (df['chrom'] == chrom) &
                (df['sample_label'] == sample) &
                (df['window_start'] >= chrom_start) &
                (df['window_start'] < chrom_end)
            )
            chrom_data = df[mask].sort_values('window_start').copy()

            if len(chrom_data) == 0:
                print(f"    {sample}: No windows found within chromosome boundaries")
                continue

            n_original = len(chrom_data)

            # Convert genome-wide coordinates to chromosome-relative coordinates
            chrom_data['chrom_relative_pos'] = chrom_data['window_start'] - chrom_start

            x_coords_raw = chrom_data['chrom_relative_pos'].values
            y_vals_raw = chrom_data['marker_count'].values

            # Downsample if needed
            if n_original <= max_points_per_chrom:
                x_coords = x_coords_raw
                y_vals = y_vals_raw
                n_plotted = n_original
            else:
                bin_size = int(np.ceil(n_original / max_points_per_chrom))
                n_bins = int(np.ceil(n_original / bin_size))
                x_coords = np.zeros(n_bins)
                y_vals = np.zeros(n_bins)

                for i in range(n_bins):
                    start_idx = i * bin_size
                    end_idx = min((i + 1) * bin_size, n_original)
                    x_coords[i] = x_coords_raw[start_idx:end_idx].mean()
                    y_vals[i] = y_vals_raw[start_idx:end_idx].max()

                n_plotted = n_bins

            print(f"    {sample}: {n_original:,} windows → {n_plotted:,} points")

            # Select color
            color = colors[sample_idx % len(colors)]

            # Plot as line only (no fill)
            linewidth = 0.5 if n_samples > 1 else 0.8

            ax.plot(x_coords, y_vals, color=color, linewidth=linewidth, alpha=0.9, label=sample)

        # Set x-axis limits to this chromosome's length (relative coordinates 0 to length)
        ax.set_xlim(0, chrom_length)

        # Styling
        ax.set_xlabel('Chromosome coordinate (bp)', fontsize=10)
        ax.set_title(f'{chrom}', fontsize=11, fontweight='bold')

        # Grid
        ax.grid(True, alpha=0.3, linewidth=0.5)
        ax.set_axisbelow(True)

        # Set log scale if requested
        if log_scale:
            ax.set_yscale('log')
            # Get all data for this chromosome to find min
            chrom_mask = (
                (df['chrom'] == chrom) &
                (df['window_start'] >= chrom_start) &
                (df['window_start'] < chrom_end)
            )
            chrom_all = df[chrom_mask]['marker_count'].values
            if len(chrom_all) > 0 and np.any(chrom_all > 0):
                y_min = chrom_all[chrom_all > 0].min()
                ax.set_ylim(bottom=max(0.1, y_min * 0.5))
            else:
                ax.set_ylim(bottom=0.1)

        # Format x-axis with scientific notation
        ax.ticklabel_format(style='scientific', axis='x', scilimits=(0,0))

        # Add legend only to first subplot if multiple samples
        if idx == 0 and n_samples > 1:
            ax.legend(loc='upper right', fontsize=8, framealpha=0.9)

    # Hide unused subplots
    for idx in range(n_chroms, len(axes)):
        axes[idx].set_visible(False)
    
    # Set y-label for the leftmost subplots in each row
    for row in range(nrows):
        if row * ncols < len(axes):
            if log_scale:
                axes[row * ncols].set_ylabel(f'Markers / 1kb (log scale)', fontsize=11)
            else:
                axes[row * ncols].set_ylabel(f'Markers / 1kb', fontsize=11)

    # Adjust layout
    plt.tight_layout()

    # Save figure in multiple formats
    print(f"\nSaving figure...")
    if isinstance(output_file, list):
        for out_file in output_file:
            plt.savefig(out_file, dpi=300, bbox_inches='tight')
            print(f"  ✓ Saved: {out_file}")
    else:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"  ✓ Saved: {output_file}")
    print(f"✓ Chromosome coverage plot saved successfully!")
    plt.close()


def plot_marker_distribution(df, chrom_lengths, k_value, output_file, percentile=50, max_bins=100):
    """
    Create histogram of marker count distribution with percentile indicator.
    Creates separate plots for each chromosome plus an overall distribution.
    
    Each chromosome plot shows the distribution of that chromosome's specific k-mers
    found within that chromosome's genomic region.
    """
    samples = sorted(df['sample_label'].unique())
    n_samples = len(samples)
    
    # Get chromosomes in proper order
    def chrom_sort_key(chrom_name):
        import re
        chrom_clean = chrom_name.replace('chr', '').replace('Chr', '')
        try:
            return (0, int(chrom_clean))
        except ValueError:
            return (1, chrom_clean)
    
    chroms = sorted(df['chrom'].unique(), key=chrom_sort_key)
    n_chroms = len(chroms)
    
    # Calculate chromosome boundaries
    chrom_positions = {}
    cumulative = 0
    for chrom in chroms:
        if chrom in chrom_lengths:
            chrom_positions[chrom] = {
                'start': cumulative,
                'end': cumulative + chrom_lengths[chrom]
            }
            cumulative += chrom_lengths[chrom]
    
    # Create figure: one column per sample, one row per chromosome + 1 for "All"
    fig, axes = plt.subplots(n_chroms + 1, n_samples, 
                            figsize=(6*n_samples, 3*(n_chroms+1)), 
                            squeeze=False)
    
    colors = ['#6BA3D4', '#D46B6B', '#6BD4A3', '#D4A36B', '#A36BD4']
    
    for sample_idx, sample in enumerate(samples):
        # Plot each chromosome
        for chrom_idx, chrom in enumerate(chroms):
            ax = axes[chrom_idx, sample_idx]
            
            if chrom not in chrom_positions:
                ax.text(0.5, 0.5, f'{chrom}\nNot in genome', 
                       ha='center', va='center', transform=ax.transAxes)
                ax.set_xticks([])
                ax.set_yticks([])
                continue
            
            chrom_start = chrom_positions[chrom]['start']
            chrom_end = chrom_positions[chrom]['end']
            
            # Filter for this chromosome's k-mers within its genomic region
            mask = (
                (df['chrom'] == chrom) & 
                (df['sample_label'] == sample) &
                (df['window_start'] >= chrom_start) &
                (df['window_start'] < chrom_end)
            )
            
            marker_counts = df[mask]['marker_count'].values
            
            if len(marker_counts) == 0:
                ax.text(0.5, 0.5, f'{chrom}\nNo data', 
                       ha='center', va='center', transform=ax.transAxes)
                ax.set_xticks([])
                ax.set_yticks([])
                continue
            
            # Calculate percentile
            percentile_val = np.percentile(marker_counts, percentile)
            
            # Create histogram
            n_bins = min(max_bins, max(10, int(np.sqrt(len(marker_counts)))))
            color = colors[sample_idx % len(colors)]
            
            ax.hist(marker_counts, bins=n_bins, color=color, 
                   edgecolor='#2E5C8A', alpha=0.8, linewidth=0.5)
            
            # Add percentile line
            ax.axvline(percentile_val, color='red', linestyle='--', linewidth=1.5,
                      label=f'{percentile}th: {int(percentile_val)}')
            
            # Styling
            ax.set_title(f'{sample} - {chrom}', fontsize=10, fontweight='bold')
            if sample_idx == 0:
                ax.set_ylabel('Frequency', fontsize=9)
            ax.legend(loc='upper right', fontsize=8)
            ax.grid(True, alpha=0.3, axis='y')
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{int(x):,}'))
        
        # Plot overall distribution (all chromosomes combined for this sample)
        ax = axes[n_chroms, sample_idx]
        
        # Get all marker counts for this sample (across all chromosomes in their regions)
        all_counts = []
        for chrom in chroms:
            if chrom not in chrom_positions:
                continue
            chrom_start = chrom_positions[chrom]['start']
            chrom_end = chrom_positions[chrom]['end']
            mask = (
                (df['chrom'] == chrom) & 
                (df['sample_label'] == sample) &
                (df['window_start'] >= chrom_start) &
                (df['window_start'] < chrom_end)
            )
            all_counts.extend(df[mask]['marker_count'].values)
        
        all_counts = np.array(all_counts)
        
        if len(all_counts) > 0:
            percentile_val = np.percentile(all_counts, percentile)
            n_bins = min(max_bins, max(10, int(np.sqrt(len(all_counts)))))
            color = colors[sample_idx % len(colors)]
            
            ax.hist(all_counts, bins=n_bins, color=color, 
                   edgecolor='#2E5C8A', alpha=0.8, linewidth=0.5)
            ax.axvline(percentile_val, color='red', linestyle='--', linewidth=2,
                      label=f'{percentile}th: {int(percentile_val)}')
            
            ax.set_title(f'{sample} - All Chromosomes', fontsize=11, fontweight='bold')
            ax.set_xlabel('Marker Count per 1kb Window', fontsize=10)
            if sample_idx == 0:
                ax.set_ylabel('Frequency', fontsize=9)
            ax.legend(loc='upper right', fontsize=9)
            ax.grid(True, alpha=0.3, axis='y')
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{int(x):,}'))
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved: {output_file}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Generate genome-wide k-mer coverage plots (optimized for large datasets)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # Single CSV file with multiple samples (e.g., B73 and Mo17 in same file)
  python plot_kmer_coverage_optimized_v2.py -i markers.csv -s B73 Mo17 -g B73.fasta -k 21 --log-scale
  
  # Multiple separate CSV files (one per sample)
  python plot_kmer_coverage_optimized_v2.py -i B73.csv Mo17.csv -s B73 Mo17 -g B73.fasta -k 21 --log-scale
  
  # Auto-detect all samples in CSV file
  python plot_kmer_coverage_optimized_v2.py -i markers.csv -g genome.fasta -k 21 --log-scale
        """
    )
    
    parser.add_argument('-i', '--input', nargs='+', required=True,
                        help='Input CSV file(s) with marker counts (space-separated for multiple)')
    parser.add_argument('-s', '--sample', nargs='+',
                        help='Sample name(s). Use to: (1) filter samples from CSV if they exist in "sample" column, OR (2) label separate input files')
    parser.add_argument('-g', '--genome', required=True,
                        help='Genome FASTA file to get chromosome lengths')
    parser.add_argument('-k', '--kmer', type=int, required=True,
                        help='K-mer size used')
    parser.add_argument('-w', '--window', type=int, default=1000,
                        help='Window size in bp (default: 1000)')
    parser.add_argument('-o', '--output', default='.',
                        help='Output directory (default: current directory)')
    parser.add_argument('-p', '--percentile', type=int, default=25,
                        help='Percentile to highlight in distribution plot (default: 25)')
    parser.add_argument('--log-scale', action='store_true',
                        help='Use logarithmic scale for y-axis in chromosome plots')
    parser.add_argument('--max-points', type=int, default=10000,
                        help='Maximum points to plot per chromosome (default: 10000)')
    parser.add_argument('--ncols', type=int, default=5,
                        help='Number of columns in the chromosome plot layout (default: 5)')
    parser.add_argument('--format', nargs='+', default=['png'],
                        choices=['png', 'pdf', 'svg'],
                        help='Output format(s) for plots (default: png). Can specify multiple formats')

    args = parser.parse_args()
    
    # Determine mode of operation
    if len(args.input) == 1 and args.sample:
        # Mode 1: Single CSV file with multiple samples specified - filter by sample names
        filter_samples = args.sample
        sample_names = None
        mode = "filter"
    elif len(args.input) > 1 and args.sample and len(args.sample) == len(args.input):
        # Mode 2: Multiple CSV files, each is a different sample
        filter_samples = None
        sample_names = args.sample
        mode = "multiple_files"
    elif not args.sample:
        # Mode 3: Use all samples from CSV file(s)
        filter_samples = None
        sample_names = None
        mode = "auto"
    else:
        print(f"Error: Number of sample names ({len(args.sample)}) must match number of input files ({len(args.input)})")
        print("\nUsage modes:")
        print("  1. Single CSV with samples inside: -i data.csv -s B73 Mo17")
        print("  2. Multiple CSV files: -i b73.csv mo17.csv -s B73 Mo17")
        print("  3. Auto-detect all samples: -i data.csv")
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load data
    print(f"Loading data from {len(args.input)} file(s)...")
    if mode == "filter":
        print(f"Filtering for samples: {', '.join(filter_samples)}")
    df = load_marker_data(args.input, sample_names, filter_samples)
    print(f"\nTotal: {len(df):,} windows across {df['chrom'].nunique()} chromosomes")
    
    # Load chromosome lengths from genome
    print(f"\nLoading chromosome lengths from {args.genome}...")
    chrom_lengths = get_chromosome_lengths_from_fasta(args.genome)
    print(f"Found {len(chrom_lengths)} chromosomes in genome file:")
    for chrom, length in sorted(chrom_lengths.items()):
        print(f"  {chrom}: {length:,} bp")
    
    # Generate output file names
    samples_in_data = sorted(df['sample_label'].unique())
    sample_str = '_vs_'.join(samples_in_data)
    base_name = f"{sample_str}_k{args.kmer}_{args.window}bp"

    # Create output files for each requested format
    if args.log_scale:
        chrom_plots = [output_dir / f"{base_name}_chromosomes_log.{fmt}" for fmt in args.format]
    else:
        chrom_plots = [output_dir / f"{base_name}_chromosomes.{fmt}" for fmt in args.format]
    dist_plots = [output_dir / f"{base_name}_distribution.{fmt}" for fmt in args.format]

    # Create plots
    print("\nGenerating chromosome coverage plot...")
    print(f"Output formats: {', '.join(args.format)}")
    plot_chromosome_coverage(df, chrom_lengths, chrom_plots, args.log_scale, args.max_points, args.ncols)

    print("\nGenerating marker distribution plot...")
    for dist_plot in dist_plots:
        plot_marker_distribution(df, chrom_lengths, args.kmer, dist_plot, args.percentile)
    
    print("\n✓ All plots generated successfully!")


if __name__ == '__main__':
    main()
