#!/usr/bin/env python3

"""
Script: 06-bin_marker_coverage.py
Description: Process per-base marker coverage into 1kb windows
             Converts comma-separated coverage values into binned marker counts

Usage: python3 06-bin_marker_coverage.py [options]

Options:
    --coverage_dir DIR     Directory with marker coverage files (default: 04-marker_coverage)
    --output_dir DIR       Output directory for binned data (default: 05-marker_bins)
    --window_size N        Window size in bp (default: 1000)
    --kmer_sizes SIZES     Comma-separated k-mer sizes (default: 21,31,41)
    --skip-cen             Skip centromere mode (chromosome-level hapmers)
    --help                 Show this help message

Requirements:
    - Python 3.6+
    - pandas

Input:
    - Coverage files from step 05 (*.counts)
    - Format: comma-separated values, one value per base position

Output:
    - ${output_dir}/k${K}_marker_counts_${WINDOW}bp.csv
    - Binned marker counts with columns:
      sample, region_type, chrom, window_start, window_end, marker_count, window_len
"""

import os
import sys
import argparse
import pandas as pd
from pathlib import Path


def process_counts_file(filepath, window_size=1000, skip_cen=False):
    """
    Process a coverage file and bin into windows

    Args:
        filepath: Path to .counts file
        window_size: Size of bins in bp
        skip_cen: Whether in skip-cen mode (affects filename parsing)

    Returns:
        DataFrame with binned marker counts
    """
    basename = os.path.basename(filepath)

    # Parse filename: unique_{PARENT}_{REGION}_{CHR}_k{K}.counts
    parts = basename.replace(".counts", "").split("_")

    try:
        # Expected format: unique_PARENT_REGION_CHR_kNN
        sample = parts[1]           # B73, Mo17, Col-0, Ler-0, etc.
        region_type = parts[2]      # CHR, CEN, or ARMS
        chrom = parts[3]            # chr1, Chr1, etc.
        k_value = parts[4]          # k21, k31, etc.
    except IndexError:
        print(f"  ⚠️  Warning: Could not parse filename: {basename}")
        print(f"     Expected format: unique_PARENT_REGION_CHR_kNN.counts")
        return None

    # Read file and collect all non-header lines
    print(f"  Processing: {basename}")

    with open(filepath, 'r') as f:
        lines = f.readlines()
        values_raw = []

        for line in lines:
            line = line.strip()
            # Skip FASTA headers and empty lines
            if line.startswith(">") or not line:
                continue
            # Split by comma and extend
            values_raw.extend(line.split(","))

    # Convert to binary presence/absence (1 if marker present, 0 if not)
    values = [int(x) > 0 for x in values_raw if x != '']

    if len(values) == 0:
        print(f"  ⚠️  Warning: No coverage values found in {basename}")
        return None

    print(f"     Found {len(values)} base positions")
    print(f"     Total markers: {sum(values)}")

    # Bin into windows and count markers
    bins = []
    for i in range(0, len(values), window_size):
        window = values[i:i+window_size]
        marker_count = sum(window)
        bins.append({
            "sample": sample,
            "region_type": region_type,
            "chrom": chrom,
            "window_start": i,
            "window_end": min(i + window_size, len(values)),
            "marker_count": marker_count,
            "window_len": len(window)
        })

    print(f"     Created {len(bins)} windows of {window_size}bp")

    return pd.DataFrame(bins)


def main():
    parser = argparse.ArgumentParser(
        description="Bin per-base marker coverage into windows",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument('--coverage_dir', default='04-marker_coverage',
                       help='Directory with coverage files (default: 04-marker_coverage)')
    parser.add_argument('--output_dir', default='05-marker_bins',
                       help='Output directory for binned data (default: 05-marker_bins)')
    parser.add_argument('--window_size', type=int, default=1000,
                       help='Window size in bp (default: 1000)')
    parser.add_argument('--kmer_sizes', default='21,31,41',
                       help='Comma-separated k-mer sizes (default: 21,31,41)')
    parser.add_argument('--skip-cen', action='store_true',
                       help='Skip centromere mode (chromosome-level hapmers)')

    args = parser.parse_args()

    # Parse k-mer sizes
    kmer_sizes = [int(k) for k in args.kmer_sizes.split(',')]

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    print("\n" + "="*70)
    print("  Marker Coverage Binning")
    print("="*70)
    print(f"\nInput directory:  {args.coverage_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Window size:      {args.window_size} bp")
    print(f"K-mer sizes:      {', '.join(map(str, kmer_sizes))}")
    print(f"Skip centromeres: {args.skip_cen}")
    print("="*70 + "\n")

    # Process each k-mer size
    for k in kmer_sizes:
        k_dir = os.path.join(args.coverage_dir, f"k{k}")

        if not os.path.exists(k_dir):
            print(f"⚠️  Warning: Directory not found: {k_dir}")
            continue

        print(f"📊 Processing k={k}...")

        # Find all .counts files
        counts_files = list(Path(k_dir).glob("*.counts"))

        if len(counts_files) == 0:
            print(f"  ⚠️  Warning: No .counts files found in {k_dir}")
            continue

        print(f"  Found {len(counts_files)} coverage files")

        # Process all files
        all_data = []
        for filepath in counts_files:
            df = process_counts_file(
                str(filepath),
                window_size=args.window_size,
                skip_cen=args.skip_cen
            )
            if df is not None:
                all_data.append(df)

        if len(all_data) == 0:
            print(f"  ⚠️  Warning: No data processed for k={k}")
            continue

        # Combine all results
        final_df = pd.concat(all_data, ignore_index=True)

        # Save to CSV
        output_file = os.path.join(
            args.output_dir,
            f"k{k}_marker_counts_{args.window_size}bp.csv"
        )
        final_df.to_csv(output_file, index=False)

        print(f"  ✅ Saved: {output_file}")
        print(f"     Total rows: {len(final_df)}")
        print(f"     Total markers: {final_df['marker_count'].sum()}")
        print("")

    print("="*70)
    print("✅ Binning complete!")
    print("="*70)
    print(f"\nOutput directory: {args.output_dir}")
    print("\nGenerated files:")
    for k in kmer_sizes:
        csv_file = f"k{k}_marker_counts_{args.window_size}bp.csv"
        if os.path.exists(os.path.join(args.output_dir, csv_file)):
            print(f"  - {csv_file}")
    print("\nNext step: Visualize binned marker density")
    print(f"  python3 scripts/07-plot_marker_coverage.py --bins_dir {args.output_dir}")
    print("")

    return 0


if __name__ == "__main__":
    sys.exit(main())
