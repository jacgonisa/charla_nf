#!/usr/bin/env python3
"""
Read Segment Visualization from BED Files

Visualizes read structures showing:
- Assigned haplotype segments (Col, Ler)
- Unassigned regions (gaps between segments)
- Crossover positions
- Read length distribution

Input: Segmentation BED file from CHARLA pipeline
Output: Read structure plots with clear unassigned regions

Author: Claude Code
Date: 2025-11-28
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import seaborn as sns
import argparse
import sys
from collections import defaultdict

# Publication settings
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.size'] = 10
plt.rcParams['figure.dpi'] = 300
sns.set_style("white")

# Color scheme for haplotypes
COLORS = {
    # Col haplotypes (blue shades)
    'CA1': '#0066CC', 'CA2': '#0077DD', 'CA3': '#0088EE', 'CA4': '#0099FF', 'CA5': '#00AAFF',
    # Ler haplotypes (orange shades)
    'LA1': '#FF6600', 'LA2': '#FF7711', 'LA3': '#FF8822', 'LA4': '#FF9933', 'LA5': '#FFAA44',
    # Complex (gray shades)
    'LC1': '#666666', 'LC2': '#777777', 'LC3': '#888888', 'LC4': '#999999', 'LC5': '#AAAAAA',
    # Unassigned (light gray)
    'UNASSIGNED': '#E0E0E0'
}


def load_segments_bed(bed_file):
    """
    Load segmentation BED file.

    Format:
    read_id    start    end    segment_type
    """
    df = pd.read_csv(bed_file, sep='\t', header=None,
                     names=['read_id', 'start', 'end', 'segment_type'])

    print(f"Loaded {len(df)} segments from {len(df['read_id'].unique())} reads")

    return df


def get_read_length(segments_df, read_id):
    """Get total read length from rightmost end coordinate."""
    read_segments = segments_df[segments_df['read_id'] == read_id]
    return read_segments['end'].max()


def identify_unassigned_regions(segments_df, read_id):
    """
    Find unassigned regions (gaps between segments) for a read.

    Returns list of (start, end) tuples for unassigned regions.
    """
    read_segments = segments_df[segments_df['read_id'] == read_id].sort_values('start')
    read_length = get_read_length(segments_df, read_id)

    unassigned = []

    # Check if there's a gap at the start
    if read_segments.iloc[0]['start'] > 0:
        unassigned.append((0, read_segments.iloc[0]['start']))

    # Check gaps between consecutive segments
    for i in range(len(read_segments) - 1):
        gap_start = read_segments.iloc[i]['end']
        gap_end = read_segments.iloc[i + 1]['start']

        if gap_end > gap_start:
            unassigned.append((gap_start, gap_end))

    # Check if there's a gap at the end
    if read_segments.iloc[-1]['end'] < read_length:
        unassigned.append((read_segments.iloc[-1]['end'], read_length))

    return unassigned


def classify_read_type(segments_df, read_id):
    """
    Classify read as SCO, NCO, or Complex based on segment patterns.

    SCO: Single crossover (exactly 2 different parental haplotypes)
    NCO: Non-crossover or gene conversion (small segments, complex patterns)
    Complex: Multiple crossovers or unclear patterns
    """
    read_segments = segments_df[segments_df['read_id'] == read_id].sort_values('start')

    # Get unique haplotype types (ignore chromosome number)
    haplotypes = set()
    for seg_type in read_segments['segment_type']:
        if seg_type.startswith('C'):
            haplotypes.add('Col')
        elif seg_type.startswith('L'):
            haplotypes.add('Ler')

    n_segments = len(read_segments)
    n_crossovers = n_segments - 1  # Transitions between segments

    # Classification logic
    if len(haplotypes) == 1:
        return 'PARENTAL'  # No crossover
    elif len(haplotypes) == 2 and n_segments == 2:
        return 'SCO'  # Single crossover
    elif len(haplotypes) == 2 and n_segments > 2:
        # Check if it's a small NCO embedded in larger segments
        segment_lengths = read_segments['end'] - read_segments['start']
        if segment_lengths.min() < 500:  # Small segment suggests NCO
            return 'NCO'
        else:
            return 'COMPLEX'
    else:
        return 'COMPLEX'


def plot_read_structures(segments_df, output_prefix, max_reads=200,
                         read_type_filter=None, confidence_scores=None):
    """
    Plot read structures showing segments and unassigned regions.

    Args:
        segments_df: DataFrame with segment information
        output_prefix: Output file prefix
        max_reads: Maximum number of reads to plot
        read_type_filter: Filter for specific read type (SCO, NCO, etc.)
        confidence_scores: Optional TSV with confidence scores for filtering
    """

    # Load confidence scores if provided
    high_conf_reads = None
    if confidence_scores:
        try:
            conf_df = pd.read_csv(confidence_scores, sep='\t')
            if 'confidence_score' in conf_df.columns:
                high_conf_reads = set(conf_df[conf_df['confidence_score'] >= 0.8]['read_id'])
                print(f"Filtering for {len(high_conf_reads)} high-confidence reads")
        except Exception as e:
            print(f"Warning: Could not load confidence scores: {e}")

    # Classify all reads
    read_classifications = {}
    for read_id in segments_df['read_id'].unique():
        read_classifications[read_id] = classify_read_type(segments_df, read_id)

    # Filter by read type if specified
    if read_type_filter:
        filtered_reads = [rid for rid, rtype in read_classifications.items()
                         if rtype == read_type_filter]
        print(f"Filtered to {len(filtered_reads)} {read_type_filter} reads")
    else:
        filtered_reads = list(read_classifications.keys())

    # Filter by confidence if provided
    if high_conf_reads:
        filtered_reads = [rid for rid in filtered_reads if rid in high_conf_reads]
        print(f"After confidence filter: {len(filtered_reads)} reads")

    # Sort by read length for better visualization
    filtered_reads = sorted(filtered_reads,
                           key=lambda x: get_read_length(segments_df, x),
                           reverse=True)

    # Limit number of reads
    reads_to_plot = filtered_reads[:max_reads]
    print(f"Plotting {len(reads_to_plot)} reads")

    # Create figure
    fig, ax = plt.subplots(figsize=(14, min(20, len(reads_to_plot) * 0.15)))

    y_position = 0
    yticks = []
    yticklabels = []

    for read_id in reads_to_plot:
        read_segments = segments_df[segments_df['read_id'] == read_id].sort_values('start')
        read_length = get_read_length(segments_df, read_id)
        read_type = read_classifications[read_id]

        # Plot each segment
        for _, segment in read_segments.iterrows():
            start = segment['start']
            end = segment['end']
            seg_type = segment['segment_type']

            color = COLORS.get(seg_type, '#CCCCCC')

            rect = patches.Rectangle((start, y_position), end - start, 0.8,
                                    linewidth=0.5, edgecolor='black',
                                    facecolor=color)
            ax.add_patch(rect)

        # Plot unassigned regions
        unassigned_regions = identify_unassigned_regions(segments_df, read_id)
        for unassigned_start, unassigned_end in unassigned_regions:
            rect = patches.Rectangle((unassigned_start, y_position),
                                    unassigned_end - unassigned_start, 0.8,
                                    linewidth=0.5, edgecolor='red', linestyle='--',
                                    facecolor=COLORS['UNASSIGNED'], alpha=0.5)
            ax.add_patch(rect)

        # Add read label
        yticks.append(y_position + 0.4)
        label = f"{read_id[:8]}... ({read_length/1000:.1f}kb, {read_type})"
        yticklabels.append(label)

        y_position += 1

    # Format plot
    ax.set_xlim(0, max(get_read_length(segments_df, rid) for rid in reads_to_plot))
    ax.set_ylim(-0.5, y_position)
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels, fontsize=6)
    ax.set_xlabel('Position (bp)', fontsize=12)
    ax.set_title(f'Read Structures - Segment Assignments\n' +
                f'(n={len(reads_to_plot)} reads, ' +
                (f'type={read_type_filter}, ' if read_type_filter else '') +
                'red dashes = unassigned regions)',
                fontsize=14, fontweight='bold')

    # Add legend
    legend_elements = [
        patches.Patch(facecolor='#0088EE', edgecolor='black', label='Col haplotypes'),
        patches.Patch(facecolor='#FF8822', edgecolor='black', label='Ler haplotypes'),
        patches.Patch(facecolor='#888888', edgecolor='black', label='Complex'),
        patches.Patch(facecolor=COLORS['UNASSIGNED'], edgecolor='red',
                     linestyle='--', alpha=0.5, label='Unassigned regions')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='x', alpha=0.3, linestyle=':')

    plt.tight_layout()

    # Save
    output_file = f"{output_prefix}_read_structures.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()


def plot_segment_statistics(segments_df, output_prefix):
    """Plot summary statistics about segments."""

    # Classify all reads
    read_types = []
    for read_id in segments_df['read_id'].unique():
        read_types.append({
            'read_id': read_id,
            'type': classify_read_type(segments_df, read_id),
            'n_segments': len(segments_df[segments_df['read_id'] == read_id]),
            'read_length': get_read_length(segments_df, read_id)
        })

    types_df = pd.DataFrame(read_types)

    # Create summary figure
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Panel 1: Read type distribution
    ax = axes[0, 0]
    type_counts = types_df['type'].value_counts()
    colors_pie = {'SCO': '#2ecc71', 'NCO': '#e74c3c',
                  'PARENTAL': '#3498db', 'COMPLEX': '#9b59b6'}
    ax.pie(type_counts.values, labels=type_counts.index, autopct='%1.1f%%',
           colors=[colors_pie.get(t, '#95a5a6') for t in type_counts.index],
           startangle=90)
    ax.set_title('Read Type Distribution', fontweight='bold')

    # Panel 2: Segments per read
    ax = axes[0, 1]
    for read_type in ['SCO', 'NCO', 'PARENTAL', 'COMPLEX']:
        data = types_df[types_df['type'] == read_type]['n_segments']
        if len(data) > 0:
            ax.hist(data, alpha=0.6, label=read_type, bins=range(1, 11),
                   color=colors_pie.get(read_type, '#95a5a6'))
    ax.set_xlabel('Number of segments per read')
    ax.set_ylabel('Frequency')
    ax.set_title('Segments per Read by Type', fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)

    # Panel 3: Read length distribution
    ax = axes[1, 0]
    for read_type in ['SCO', 'NCO', 'PARENTAL', 'COMPLEX']:
        data = types_df[types_df['type'] == read_type]['read_length'] / 1000
        if len(data) > 0:
            ax.hist(data, alpha=0.6, label=read_type, bins=30,
                   color=colors_pie.get(read_type, '#95a5a6'))
    ax.set_xlabel('Read length (kb)')
    ax.set_ylabel('Frequency')
    ax.set_title('Read Length Distribution by Type', fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)

    # Panel 4: Segment length distribution
    ax = axes[1, 1]
    segment_lengths = (segments_df['end'] - segments_df['start']) / 1000

    # Separate by Col/Ler/Complex
    col_segs = segments_df[segments_df['segment_type'].str.startswith('C')]
    ler_segs = segments_df[segments_df['segment_type'].str.startswith('L')]

    col_lengths = (col_segs['end'] - col_segs['start']) / 1000
    ler_lengths = (ler_segs['end'] - ler_segs['start']) / 1000

    ax.hist([col_lengths, ler_lengths], bins=50, alpha=0.6,
           label=['Col segments', 'Ler segments'],
           color=['#0088EE', '#FF8822'], stacked=False)
    ax.set_xlabel('Segment length (kb)')
    ax.set_ylabel('Frequency')
    ax.set_title('Segment Length Distribution', fontweight='bold')
    ax.set_yscale('log')
    ax.legend()
    ax.grid(alpha=0.3)

    plt.tight_layout()

    output_file = f"{output_prefix}_statistics.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()

    # Write summary table
    summary_file = f"{output_prefix}_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("Read Segment Analysis Summary\n")
        f.write("=" * 60 + "\n\n")

        f.write(f"Total reads: {len(types_df)}\n")
        f.write(f"Total segments: {len(segments_df)}\n\n")

        f.write("Read type breakdown:\n")
        for read_type in ['SCO', 'NCO', 'PARENTAL', 'COMPLEX']:
            count = len(types_df[types_df['type'] == read_type])
            pct = 100 * count / len(types_df)
            f.write(f"  {read_type}: {count} ({pct:.1f}%)\n")

        f.write(f"\nAverage segments per read: {types_df['n_segments'].mean():.2f}\n")
        f.write(f"Average read length: {types_df['read_length'].mean()/1000:.1f} kb\n")
        f.write(f"Average segment length: {segment_lengths.mean():.1f} kb\n")

    print(f"Saved: {summary_file}")


def main():
    parser = argparse.ArgumentParser(description='Visualize read segments from BED file')
    parser.add_argument('--bed', required=True, help='Segmentation BED file')
    parser.add_argument('--output-prefix', required=True, help='Output file prefix')
    parser.add_argument('--max-reads', type=int, default=200,
                       help='Maximum reads to plot (default: 200)')
    parser.add_argument('--read-type', choices=['SCO', 'NCO', 'PARENTAL', 'COMPLEX'],
                       help='Filter for specific read type')
    parser.add_argument('--confidence-scores', help='Optional confidence scores TSV for filtering')

    args = parser.parse_args()

    print("=" * 60)
    print("Read Segment Visualization")
    print("=" * 60)

    # Load data
    segments_df = load_segments_bed(args.bed)

    # Plot structures
    plot_read_structures(segments_df, args.output_prefix,
                        max_reads=args.max_reads,
                        read_type_filter=args.read_type,
                        confidence_scores=args.confidence_scores)

    # Plot statistics
    plot_segment_statistics(segments_df, args.output_prefix)

    print("=" * 60)
    print("Complete!")
    print("=" * 60)


if __name__ == '__main__':
    main()
