#!/usr/bin/env python3
"""
Normalize k-mer counts to correct for marker database imbalances.

This script applies normalization to k-mer count matrices to account for:
1. Unequal marker counts between parents
2. Variable marker density across regions
3. Chromosome-specific imbalances

Normalization methods:
- Probability scoring: Convert counts to probabilities based on marker availability
- Weighted scoring: Weight counts by inverse marker frequency
- Min-max scaling: Scale to 0-1 range per region

Input: Raw k-mer count matrix from process_cenhapmer_counts
Output: Normalized k-mer count matrix
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
from collections import defaultdict


def load_marker_stats(marker_stats_file):
    """
    Load marker statistics from analysis output.

    Expected format (TSV):
    database    parent  region_code  region_type  chromosome  kmer_count
    Col-0_CA1   Col-0   CA          ARMS         1           12345
    ...
    """
    if not os.path.exists(marker_stats_file):
        print(f"Warning: Marker stats file not found: {marker_stats_file}", file=sys.stderr)
        return None

    df = pd.read_csv(marker_stats_file, sep='\t')
    return df


def calculate_normalization_weights(marker_stats_df, method='probability'):
    """
    Calculate normalization weights based on marker statistics.

    Parameters
    ----------
    marker_stats_df : pd.DataFrame
        Marker statistics with columns: database, parent, region_code, kmer_count
    method : str
        Normalization method: 'probability', 'weighted', or 'minmax'

    Returns
    -------
    dict
        Dictionary mapping region_code to normalization weight
    """
    weights = {}

    if method == 'probability':
        # Convert to probability: weight = 1 / marker_count
        # So count/weight gives count/total_markers = probability
        for _, row in marker_stats_df.iterrows():
            region_code = row['region_code']
            kmer_count = row['kmer_count']
            if kmer_count > 0:
                weights[region_code] = kmer_count  # Store total, will divide counts by this

    elif method == 'weighted':
        # Weight by inverse frequency relative to median
        median_markers = marker_stats_df['kmer_count'].median()
        for _, row in marker_stats_df.iterrows():
            region_code = row['region_code']
            kmer_count = row['kmer_count']
            if kmer_count > 0:
                # Regions with more markers get lower weight
                weights[region_code] = median_markers / kmer_count

    elif method == 'global_mean':
        # Normalize by global mean across all regions
        mean_markers = marker_stats_df['kmer_count'].mean()
        for _, row in marker_stats_df.iterrows():
            region_code = row['region_code']
            kmer_count = row['kmer_count']
            if kmer_count > 0:
                weights[region_code] = mean_markers / kmer_count

    else:
        raise ValueError(f"Unknown normalization method: {method}")

    return weights


def normalize_kmer_matrix(input_tsv, output_tsv, marker_stats_file, method='probability', scale_factor=1000):
    """
    Apply normalization to k-mer count matrix.

    Parameters
    ----------
    input_tsv : str
        Input k-mer count matrix (read_id x region columns)
    output_tsv : str
        Output normalized matrix
    marker_stats_file : str
        Marker statistics file from 16-analyze_marker_balance.py
    method : str
        Normalization method
    scale_factor : float
        Scaling factor for probability method (to avoid tiny numbers)
    """

    print(f"Loading k-mer count matrix: {input_tsv}")
    df = pd.read_csv(input_tsv, sep='\t')

    print(f"Matrix shape: {df.shape[0]} reads x {df.shape[1]} columns")

    # Load marker statistics
    if marker_stats_file and os.path.exists(marker_stats_file):
        print(f"Loading marker statistics: {marker_stats_file}")
        marker_stats_df = load_marker_stats(marker_stats_file)

        if marker_stats_df is not None:
            print(f"Calculating normalization weights (method: {method})...")
            weights = calculate_normalization_weights(marker_stats_df, method=method)

            # Apply normalization to each region column
            normalized_count = 0
            for col in df.columns:
                if col == 'read_id' or col == 'read':
                    continue  # Skip read ID column

                # Extract region code from column name
                # Column format might be: CA1, CC1, LA1, LC1, etc.
                region_code = col[:2] if len(col) >= 2 else None

                if region_code and region_code in weights:
                    if method == 'probability':
                        # Convert to probability: count / total_markers * scale_factor
                        total_markers = weights[region_code]
                        df[col] = (df[col] / total_markers) * scale_factor
                    elif method in ['weighted', 'global_mean']:
                        # Multiply by weight
                        df[col] = df[col] * weights[region_code]

                    normalized_count += 1

            print(f"Normalized {normalized_count} columns")
        else:
            print("Warning: Could not load marker statistics, skipping normalization", file=sys.stderr)
    else:
        print(f"Warning: Marker stats file not found, outputting original counts", file=sys.stderr)

    # Save normalized matrix
    print(f"Saving normalized matrix: {output_tsv}")
    df.to_csv(output_tsv, sep='\t', index=False)

    # Print summary statistics
    print("\nNormalization Summary:")
    print(f"  Method: {method}")
    if method == 'probability':
        print(f"  Scale factor: {scale_factor}")
    print(f"  Normalized columns: {normalized_count}")

    # Show before/after statistics for a few columns
    print("\nExample column statistics (after normalization):")
    numeric_cols = df.select_dtypes(include=[np.number]).columns[:5]
    for col in numeric_cols:
        mean_val = df[col].mean()
        median_val = df[col].median()
        max_val = df[col].max()
        print(f"  {col}: mean={mean_val:.2f}, median={median_val:.2f}, max={max_val:.2f}")


def main():
    parser = argparse.ArgumentParser(
        description='Normalize k-mer counts to correct for marker imbalances',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Normalization Methods:
  probability    - Convert to probability scores (count / total_markers)
                   Best for comparing across very imbalanced regions
  weighted       - Weight by inverse frequency relative to median
                   Good for moderate imbalances (10-30%)
  global_mean    - Normalize by global mean across all regions
                   Simple and robust for most cases

Examples:
  # Probability normalization (for high imbalance)
  python 17-normalize_kmer_counts.py \\
      --input kmer_matrix.tsv \\
      --output kmer_matrix_normalized.tsv \\
      --marker-stats marker_analysis/k31_balance_marker_counts.tsv \\
      --method probability \\
      --scale-factor 10000

  # Weighted normalization (for moderate imbalance)
  python 17-normalize_kmer_counts.py \\
      --input kmer_matrix.tsv \\
      --output kmer_matrix_normalized.tsv \\
      --marker-stats marker_analysis/k31_balance_marker_counts.tsv \\
      --method weighted

  # No normalization (just copy)
  python 17-normalize_kmer_counts.py \\
      --input kmer_matrix.tsv \\
      --output kmer_matrix.tsv
        """
    )

    parser.add_argument('--input', required=True,
                       help='Input k-mer count matrix (TSV)')
    parser.add_argument('--output', required=True,
                       help='Output normalized matrix (TSV)')
    parser.add_argument('--marker-stats', required=False,
                       help='Marker statistics file from 16-analyze_marker_balance.py')
    parser.add_argument('--method', choices=['probability', 'weighted', 'global_mean'],
                       default='weighted',
                       help='Normalization method (default: weighted)')
    parser.add_argument('--scale-factor', type=float, default=1000,
                       help='Scale factor for probability method (default: 1000)')

    args = parser.parse_args()

    # Normalize
    normalize_kmer_matrix(
        args.input,
        args.output,
        args.marker_stats,
        method=args.method,
        scale_factor=args.scale_factor
    )

    print("\n✓ Normalization complete!")


if __name__ == '__main__':
    main()
