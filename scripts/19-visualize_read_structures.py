#!/usr/bin/env python3
"""
Read Structure Visualization and Recombination Analysis

Creates visual representations of crossover reads showing:
- Parental haplotype segments (Col=blue, Ler=orange, undetermined=grey)
- Read structure sorted by length
- Parental proportions
- MAPQ quality per segment
- Recombination rate (events per Mb)

Author: Claude Code
Date: 2025-11-27
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import seaborn as sns
import argparse
import sys
import os
import gzip
from pathlib import Path

# Publication settings
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.size'] = 9
plt.rcParams['figure.dpi'] = 300
sns.set_style("white")

# Colors
COL_COLOR = '#1f77b4'  # Blue
LER_COLOR = '#ff7f0e'  # Orange
UNDETERMINED_COLOR = '#808080'  # Grey


def load_paf_files(paf_dir):
    """
    Load all PAF files from directory with full alignment information.
    """
    paf_records = []

    if not os.path.exists(paf_dir):
        print(f"[WARNING] PAF directory not found: {paf_dir}", file=sys.stderr)
        return pd.DataFrame()

    # Find all PAF files
    paf_files = list(Path(paf_dir).rglob('*.paf'))

    print(f"[INFO] Found {len(paf_files)} PAF files", file=sys.stderr)

    for paf_file in paf_files:
        try:
            open_fn = gzip.open if str(paf_file).endswith('.gz') else open
            with open_fn(paf_file, 'rt') as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue

                    parts = line.strip().split('\t')
                    if len(parts) < 12:
                        continue

                    # Extract base read ID and segment info from query name
                    query_name = parts[0]

                    # Parse segment info: read_id_segment_parent
                    # Example: 7ffd019e-57a9-4faf-bb55-a5a0be8ccc70_1_CA5
                    parts_split = query_name.rsplit('_', 2)
                    if len(parts_split) >= 3:
                        base_read_id = parts_split[0]
                        segment_idx = parts_split[1]
                        segment_label = parts_split[2]  # e.g., CA5, LA5
                    else:
                        base_read_id = query_name
                        segment_idx = '0'
                        segment_label = 'unknown'

                    # Determine parent from target or segment label
                    target = parts[5]
                    if 'Col' in target:
                        parent = 'Col'
                    elif 'Ler' in target:
                        parent = 'Ler'
                    elif 'C' in segment_label[0]:
                        parent = 'Col'
                    elif 'L' in segment_label[0]:
                        parent = 'Ler'
                    else:
                        parent = 'undetermined'

                    paf_record = {
                        'query_name': query_name,
                        'base_read_id': base_read_id,
                        'segment_idx': int(segment_idx) if segment_idx.isdigit() else 0,
                        'segment_label': segment_label,
                        'parent': parent,
                        'qlen': int(parts[1]),
                        'qstart': int(parts[2]),
                        'qend': int(parts[3]),
                        'strand': parts[4],
                        'target': parts[5],
                        'tlen': int(parts[6]),
                        'tstart': int(parts[7]),
                        'tend': int(parts[8]),
                        'matches': int(parts[9]),
                        'block_len': int(parts[10]),
                        'mapq': int(parts[11])
                    }

                    # Parse optional tags
                    for tag in parts[12:]:
                        if tag.startswith('AS:i:'):
                            paf_record['AS'] = int(tag[5:])
                        elif tag.startswith('NM:i:'):
                            paf_record['NM'] = int(tag[5:])
                        elif tag.startswith('de:f:'):
                            paf_record['divergence'] = float(tag[5:])

                    # Calculate identity
                    paf_record['identity'] = paf_record['matches'] / paf_record['block_len'] if paf_record['block_len'] > 0 else 0
                    paf_record['length'] = paf_record['qend'] - paf_record['qstart']

                    paf_records.append(paf_record)

        except Exception as e:
            print(f"[WARNING] Error loading {paf_file}: {e}", file=sys.stderr)

    if not paf_records:
        print(f"[WARNING] No PAF records loaded from {paf_dir}", file=sys.stderr)
        return pd.DataFrame()

    df = pd.DataFrame(paf_records)
    print(f"[INFO] Loaded {len(df)} PAF alignments for {df['base_read_id'].nunique()} reads", file=sys.stderr)

    return df


def calculate_read_structures(paf_df):
    """
    Calculate fictional positions within each read for visualization.
    Similar to the R code you provided.
    """
    # Sort by read and segment
    paf_sorted = paf_df.sort_values(['base_read_id', 'segment_idx']).copy()

    # Calculate cumulative positions
    paf_sorted['frag_start'] = paf_sorted.groupby('base_read_id')['length'].apply(
        lambda x: x.shift(1, fill_value=0).cumsum()
    ).values

    paf_sorted['frag_end'] = paf_sorted['frag_start'] + paf_sorted['length']

    # Calculate total read length
    read_lengths = paf_sorted.groupby('base_read_id')['length'].sum().reset_index()
    read_lengths.columns = ['base_read_id', 'read_length']

    # Merge
    paf_sorted = paf_sorted.merge(read_lengths, on='base_read_id')

    # Sort reads by length (longest first)
    read_order = read_lengths.sort_values('read_length', ascending=False)['base_read_id'].tolist()

    return paf_sorted, read_order


def calculate_parental_proportions(paf_df):
    """
    Calculate proportion of each read from Col vs Ler.
    """
    props = []

    for read_id, group in paf_df.groupby('base_read_id'):
        total_length = group['length'].sum()

        col_length = group[group['parent'] == 'Col']['length'].sum()
        ler_length = group[group['parent'] == 'Ler']['length'].sum()
        und_length = group[group['parent'] == 'undetermined']['length'].sum()

        props.append({
            'base_read_id': read_id,
            'total_length': total_length,
            'col_length': col_length,
            'ler_length': ler_length,
            'undetermined_length': und_length,
            'col_prop': col_length / total_length if total_length > 0 else 0,
            'ler_prop': ler_length / total_length if total_length > 0 else 0,
            'undetermined_prop': und_length / total_length if total_length > 0 else 0,
            'n_segments': len(group),
            'n_switches': (group['parent'] != group['parent'].shift()).sum() - 1,  # -1 for first
            'mean_mapq': group['mapq'].mean(),
            'mean_identity': group['identity'].mean()
        })

    return pd.DataFrame(props)


def calculate_recombination_rate(props_df):
    """
    Calculate recombination events per Mb sequenced.
    """
    total_mb = props_df['total_length'].sum() / 1_000_000
    total_switches = props_df['n_switches'].sum()

    rate = total_switches / total_mb if total_mb > 0 else 0

    return rate, total_switches, total_mb


def plot_read_structures(paf_df, read_order, output_file, max_reads=200):
    """
    Create the read structure visualization like your R ggplot code.
    Reads are sorted by length, colored by parent.
    """
    # Limit to top N reads if too many
    if len(read_order) > max_reads:
        print(f"[INFO] Limiting visualization to top {max_reads} longest reads", file=sys.stderr)
        read_order = read_order[:max_reads]

    # Filter to selected reads
    plot_df = paf_df[paf_df['base_read_id'].isin(read_order)].copy()

    # Create figure
    fig, ax = plt.subplots(figsize=(16, max(8, len(read_order) * 0.12)))

    # Assign y positions (reverse order so longest at top)
    read_to_y = {read_id: i for i, read_id in enumerate(read_order)}
    plot_df['y_pos'] = plot_df['base_read_id'].map(read_to_y)

    # Plot each segment as a rectangle
    for _, row in plot_df.iterrows():
        y = row['y_pos']

        # Color by parent
        if row['parent'] == 'Col':
            color = COL_COLOR
        elif row['parent'] == 'Ler':
            color = LER_COLOR
        else:
            color = UNDETERMINED_COLOR

        # Rectangle
        rect = mpatches.Rectangle(
            (row['frag_start'], y - 0.4),
            row['length'],
            0.8,
            facecolor=color,
            edgecolor='black',
            linewidth=0.5,
            alpha=0.8
        )
        ax.add_patch(rect)

    # Set limits
    ax.set_xlim(0, plot_df.groupby('base_read_id')['frag_end'].max().max() * 1.02)
    ax.set_ylim(-1, len(read_order))

    # Labels
    ax.set_xlabel('Fictional read position (bp)', fontweight='bold', fontsize=12)
    ax.set_ylabel('Read ID (sorted by length, longest at top)', fontweight='bold', fontsize=12)
    ax.set_title(f'Structure of Crossover Reads (n={len(read_order)})', fontweight='bold', fontsize=14)

    # Remove y-axis ticks
    ax.set_yticks([])

    # Legend
    col_patch = mpatches.Patch(color=COL_COLOR, label='Col-0', alpha=0.8)
    ler_patch = mpatches.Patch(color=LER_COLOR, label='Ler-0', alpha=0.8)
    und_patch = mpatches.Patch(color=UNDETERMINED_COLOR, label='Undetermined', alpha=0.8)
    ax.legend(handles=[col_patch, ler_patch, und_patch], loc='upper right', fontsize=10)

    ax.grid(True, alpha=0.2, axis='x')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"[SUCCESS] Saved read structure plot: {output_file}", file=sys.stderr)


def plot_parental_proportions(props_df, output_file):
    """
    Plot parental proportion analysis.
    """
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('Parental Haplotype Analysis', fontsize=16, fontweight='bold')

    # Panel 1: Stacked bar of proportions
    ax = axes[0, 0]
    props_sorted = props_df.sort_values('total_length', ascending=False).head(100)

    x = np.arange(len(props_sorted))
    ax.bar(x, props_sorted['col_prop'], color=COL_COLOR, alpha=0.7, label='Col-0')
    ax.bar(x, props_sorted['ler_prop'], bottom=props_sorted['col_prop'],
           color=LER_COLOR, alpha=0.7, label='Ler-0')
    ax.bar(x, props_sorted['undetermined_prop'],
           bottom=props_sorted['col_prop'] + props_sorted['ler_prop'],
           color=UNDETERMINED_COLOR, alpha=0.7, label='Undetermined')

    ax.set_xlabel('Read (sorted by length)', fontweight='bold')
    ax.set_ylabel('Proportion', fontweight='bold')
    ax.set_title(f'Parental Proportions per Read (top 100)', fontweight='bold')
    ax.set_ylim(0, 1)
    ax.legend()

    # Panel 2: Col vs Ler proportion scatter
    ax = axes[0, 1]
    ax.scatter(props_df['col_prop'], props_df['ler_prop'], alpha=0.3, s=20,
              c=props_df['n_switches'], cmap='viridis', edgecolor='black', linewidth=0.3)
    ax.plot([0, 1], [1, 0], 'r--', alpha=0.5, label='Col+Ler=100%')
    ax.set_xlabel('Col-0 Proportion', fontweight='bold')
    ax.set_ylabel('Ler-0 Proportion', fontweight='bold')
    ax.set_title('Col vs Ler Balance', fontweight='bold')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend()
    cbar = plt.colorbar(ax.collections[0], ax=ax)
    cbar.set_label('N Switches', fontweight='bold')

    # Panel 3: Number of segments
    ax = axes[0, 2]
    seg_counts = props_df['n_segments'].value_counts().sort_index()
    ax.bar(seg_counts.index, seg_counts.values, color='steelblue', alpha=0.7,
          edgecolor='black', linewidth=1)
    ax.set_xlabel('Number of Segments', fontweight='bold')
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_title('Segment Count Distribution', fontweight='bold')

    # Panel 4: Number of switches (recombination events)
    ax = axes[1, 0]
    switch_counts = props_df['n_switches'].value_counts().sort_index()
    ax.bar(switch_counts.index, switch_counts.values, color='#2ecc71', alpha=0.7,
          edgecolor='black', linewidth=1)
    ax.set_xlabel('Number of Parental Switches', fontweight='bold')
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_title('Recombination Events per Read', fontweight='bold')

    # Panel 5: MAPQ distribution
    ax = axes[1, 1]
    ax.hist(props_df['mean_mapq'], bins=50, color='purple', alpha=0.7,
           edgecolor='black', linewidth=0.5)
    ax.axvline(props_df['mean_mapq'].median(), color='red', linestyle='--',
              linewidth=2, label=f'Median: {props_df["mean_mapq"].median():.1f}')
    ax.set_xlabel('Mean MAPQ per Read', fontweight='bold')
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_title('Mapping Quality Distribution', fontweight='bold')
    ax.legend()

    # Panel 6: Identity distribution
    ax = axes[1, 2]
    ax.hist(props_df['mean_identity'], bins=50, color='orange', alpha=0.7,
           edgecolor='black', linewidth=0.5)
    ax.axvline(props_df['mean_identity'].median(), color='red', linestyle='--',
              linewidth=2, label=f'Median: {props_df["mean_identity"].median():.3f}')
    ax.set_xlabel('Mean Identity per Read', fontweight='bold')
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_title('Alignment Identity Distribution', fontweight='bold')
    ax.legend()

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"[SUCCESS] Saved parental proportions plot: {output_file}", file=sys.stderr)


def create_summary_report(props_df, recomb_rate, total_switches, total_mb, output_file):
    """
    Create text summary of the analysis.
    """
    with open(output_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("READ STRUCTURE AND RECOMBINATION ANALYSIS SUMMARY\n")
        f.write("="*80 + "\n\n")

        f.write(f"Total Reads Analyzed: {len(props_df):,}\n")
        f.write(f"Total Sequence Length: {total_mb:.2f} Mb\n\n")

        f.write("PARENTAL PROPORTIONS:\n")
        f.write(f"  Mean Col-0 proportion: {props_df['col_prop'].mean():.3f}\n")
        f.write(f"  Mean Ler-0 proportion: {props_df['ler_prop'].mean():.3f}\n")
        f.write(f"  Mean undetermined: {props_df['undetermined_prop'].mean():.3f}\n\n")

        f.write("RECOMBINATION EVENTS:\n")
        f.write(f"  Total switches: {total_switches:,}\n")
        f.write(f"  Recombination rate: {recomb_rate:.2f} events/Mb\n")
        f.write(f"  Mean switches per read: {props_df['n_switches'].mean():.2f}\n")
        f.write(f"  Median switches per read: {props_df['n_switches'].median():.0f}\n\n")

        f.write("SEGMENT STATISTICS:\n")
        f.write(f"  Mean segments per read: {props_df['n_segments'].mean():.2f}\n")
        f.write(f"  Median segments per read: {props_df['n_segments'].median():.0f}\n\n")

        f.write("QUALITY METRICS:\n")
        f.write(f"  Mean MAPQ: {props_df['mean_mapq'].mean():.1f}\n")
        f.write(f"  Median MAPQ: {props_df['mean_mapq'].median():.1f}\n")
        f.write(f"  Mean Identity: {props_df['mean_identity'].mean():.4f}\n")
        f.write(f"  Median Identity: {props_df['mean_identity'].median():.4f}\n\n")

        f.write("="*80 + "\n")

    print(f"[SUCCESS] Saved summary report: {output_file}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='Visualize read structures and analyze recombination'
    )

    parser.add_argument('--paf-dir', required=True, help='Directory containing PAF files')
    parser.add_argument('--output-prefix', required=True, help='Output prefix for plots')
    parser.add_argument('--max-reads', type=int, default=200, help='Max reads to plot (default: 200)')
    parser.add_argument('--crossover-type', default='all', choices=['all', 'SCO', 'NCO'],
                       help='Filter by crossover type (requires confidence scores)')
    parser.add_argument('--confidence-scores', help='Optional: confidence scores TSV to filter by type')

    args = parser.parse_args()

    print("[INFO] Loading PAF files...", file=sys.stderr)
    paf_df = load_paf_files(args.paf_dir)

    if len(paf_df) == 0:
        print("[ERROR] No PAF data loaded!", file=sys.stderr)
        sys.exit(1)

    # Filter by crossover type if requested
    if args.crossover_type != 'all' and args.confidence_scores:
        print(f"[INFO] Filtering for {args.crossover_type} reads...", file=sys.stderr)
        conf_df = pd.read_csv(args.confidence_scores, sep='\t')
        if args.crossover_type == 'SCO':
            selected_reads = conf_df[conf_df['crossover_type'] == 'SCO']['read_id'].tolist()
        else:
            selected_reads = conf_df[conf_df['crossover_type'] != 'SCO']['read_id'].tolist()

        paf_df = paf_df[paf_df['base_read_id'].isin(selected_reads)]
        print(f"[INFO] After filtering: {paf_df['base_read_id'].nunique()} reads", file=sys.stderr)

    print("[INFO] Calculating read structures...", file=sys.stderr)
    paf_structured, read_order = calculate_read_structures(paf_df)

    print("[INFO] Calculating parental proportions...", file=sys.stderr)
    props_df = calculate_parental_proportions(paf_structured)

    print("[INFO] Calculating recombination rate...", file=sys.stderr)
    recomb_rate, total_switches, total_mb = calculate_recombination_rate(props_df)

    print(f"\n[INFO] Recombination rate: {recomb_rate:.2f} events/Mb", file=sys.stderr)
    print(f"[INFO] Total switches: {total_switches:,}", file=sys.stderr)
    print(f"[INFO] Total sequence: {total_mb:.2f} Mb\n", file=sys.stderr)

    print("[INFO] Creating read structure visualization...", file=sys.stderr)
    plot_read_structures(paf_structured, read_order,
                        f'{args.output_prefix}_read_structures.png',
                        args.max_reads)

    print("[INFO] Creating parental proportion plots...", file=sys.stderr)
    plot_parental_proportions(props_df, f'{args.output_prefix}_parental_analysis.png')

    print("[INFO] Creating summary report...", file=sys.stderr)
    create_summary_report(props_df, recomb_rate, total_switches, total_mb,
                         f'{args.output_prefix}_summary.txt')

    # Save detailed table
    props_df.to_csv(f'{args.output_prefix}_parental_proportions.tsv', sep='\t', index=False, float_format='%.4f')
    print(f"[SUCCESS] Saved parental proportions table", file=sys.stderr)

    print("\n" + "="*80, file=sys.stderr)
    print("READ STRUCTURE VISUALIZATION COMPLETE", file=sys.stderr)
    print("="*80, file=sys.stderr)


if __name__ == '__main__':
    main()
