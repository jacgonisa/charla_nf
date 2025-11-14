#!/usr/bin/env python3
"""
Analyze MAFFT alignment output for satellite blocks.
Updates the satellite_analysis.tsv with conservation statistics.
Creates visualization of MSA showing conservation patterns.
"""

import os
import sys
import pandas as pd
from Bio import AlignIO
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns


def analyze_alignment_quality(alignment_file):
    """
    Analyze quality and conservation of MSA alignment.
    """
    if not os.path.exists(alignment_file) or os.path.getsize(alignment_file) == 0:
        print(f"[WARNING] Alignment file not found or empty: {alignment_file}")
        return None

    try:
        # Try reading as FASTA format (MAFFT default output)
        alignment = AlignIO.read(alignment_file, "fasta")

        # Calculate alignment statistics
        num_sequences = len(alignment)
        alignment_length = alignment.get_alignment_length()

        # Calculate conservation score (simple identity-based)
        conservation_scores = []
        for i in range(alignment_length):
            column = alignment[:, i]
            most_common = Counter(column).most_common(1)[0]
            conservation = most_common[1] / num_sequences
            conservation_scores.append(conservation)

        avg_conservation = np.mean(conservation_scores)
        highly_conserved_positions = sum(1 for x in conservation_scores if x >= 0.8)

        return {
            'num_sequences': num_sequences,
            'alignment_length': alignment_length,
            'avg_conservation': round(avg_conservation, 4),
            'highly_conserved_positions': highly_conserved_positions,
            'conservation_percentage': round(highly_conserved_positions / alignment_length, 4) if alignment_length > 0 else 0
        }
    except Exception as e:
        print(f"[ERROR] Could not analyze alignment: {e}")
        return None


def create_msa_visualization(alignment_file, output_dir):
    """
    Create comprehensive MSA visualization showing:
    1. Alignment heatmap
    2. Conservation plot
    3. Sequence identity matrix
    """
    if not os.path.exists(alignment_file) or os.path.getsize(alignment_file) == 0:
        print(f"[WARNING] Cannot create visualization - alignment file not found or empty")
        return

    try:
        alignment = AlignIO.read(alignment_file, "fasta")
        num_sequences = len(alignment)
        alignment_length = alignment.get_alignment_length()

        print(f"\n[INFO] Creating MSA visualization...")
        print(f"  Sequences: {num_sequences}, Length: {alignment_length} bp")

        # Convert alignment to numerical matrix for visualization
        # A=0, C=1, G=2, T=3, gap=-1, N=-2
        base_to_num = {'A': 0, 'C': 1, 'G': 2, 'T': 3, '-': -1, 'N': -2, 'a': 0, 'c': 1, 'g': 2, 't': 3}
        alignment_matrix = np.zeros((num_sequences, alignment_length))

        for i, record in enumerate(alignment):
            for j, base in enumerate(str(record.seq)):
                alignment_matrix[i, j] = base_to_num.get(base.upper(), -2)

        # Calculate conservation scores
        conservation_scores = []
        for i in range(alignment_length):
            column = alignment[:, i]
            most_common = Counter(column).most_common(1)[0]
            conservation = most_common[1] / num_sequences
            conservation_scores.append(conservation)

        # Create figure with multiple subplots
        fig = plt.figure(figsize=(20, 12))
        gs = fig.add_gridspec(3, 2, height_ratios=[1, 3, 2], hspace=0.3, wspace=0.3)

        # Subplot 1: Conservation plot (top, spanning both columns)
        ax1 = fig.add_subplot(gs[0, :])
        positions = np.arange(alignment_length)
        ax1.fill_between(positions, conservation_scores, alpha=0.7, color='green', edgecolor='darkgreen', linewidth=1)
        ax1.axhline(y=0.8, color='red', linestyle='--', linewidth=2, label='80% conservation threshold')
        ax1.set_xlabel('Alignment Position (bp)', fontweight='bold', fontsize=12)
        ax1.set_ylabel('Conservation', fontweight='bold', fontsize=12)
        ax1.set_title('Sequence Conservation Across Alignment', fontweight='bold', fontsize=14)
        ax1.set_ylim(0, 1)
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # Add statistics text
        avg_cons = np.mean(conservation_scores)
        highly_conserved = sum(1 for x in conservation_scores if x >= 0.8)
        ax1.text(0.02, 0.95,
                f'Avg conservation: {avg_cons:.2%}\n'
                f'Highly conserved (≥80%): {highly_conserved} / {alignment_length} ({highly_conserved/alignment_length:.1%})',
                transform=ax1.transAxes, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
                fontsize=10, fontweight='bold')

        # Subplot 2: Alignment heatmap (middle left)
        ax2 = fig.add_subplot(gs[1, 0])

        # Create custom colormap for bases
        colors_list = ['white', 'green', 'blue', 'orange', 'red', 'gray']  # gap, A, C, G, T, N
        n_bins = 6
        cmap = LinearSegmentedColormap.from_list('bases',
                                                 [(i/(n_bins-1), colors_list[i]) for i in range(n_bins)],
                                                 N=n_bins)

        # Subsample if alignment is too long for visualization
        max_display_width = 500
        if alignment_length > max_display_width:
            step = alignment_length // max_display_width
            display_matrix = alignment_matrix[:, ::step]
            display_positions = positions[::step]
            title_suffix = f' (subsampled, showing every {step}th position)'
        else:
            display_matrix = alignment_matrix
            display_positions = positions
            title_suffix = ''

        im = ax2.imshow(display_matrix, cmap=cmap, aspect='auto', interpolation='nearest', vmin=-1, vmax=3)
        ax2.set_xlabel('Alignment Position (bp)', fontweight='bold')
        ax2.set_ylabel('Sequence ID', fontweight='bold')
        ax2.set_title(f'Multiple Sequence Alignment Heatmap{title_suffix}', fontweight='bold')

        # Add colorbar with labels
        cbar = plt.colorbar(im, ax=ax2, ticks=[-1, 0, 1, 2, 3])
        cbar.set_label('Base', fontweight='bold')
        cbar.ax.set_yticklabels(['Gap', 'A', 'C', 'G', 'T'])

        # Subplot 3: Base composition plot (middle right)
        ax3 = fig.add_subplot(gs[1, 1])

        # Calculate base composition per position
        bases = ['A', 'C', 'G', 'T', 'Gap']
        base_counts = {base: [] for base in bases}

        for i in range(alignment_length):
            column = str(alignment[:, i])
            total = len(column)
            base_counts['A'].append(column.count('A') + column.count('a'))
            base_counts['C'].append(column.count('C') + column.count('c'))
            base_counts['G'].append(column.count('G') + column.count('g'))
            base_counts['T'].append(column.count('T') + column.count('t'))
            base_counts['Gap'].append(column.count('-'))

        # Convert to percentages
        base_percentages = {base: np.array(counts) / num_sequences * 100
                           for base, counts in base_counts.items()}

        # Create stacked area plot
        ax3.stackplot(positions,
                     base_percentages['A'], base_percentages['C'],
                     base_percentages['G'], base_percentages['T'],
                     base_percentages['Gap'],
                     labels=['A', 'C', 'G', 'T', 'Gap'],
                     colors=['green', 'blue', 'orange', 'red', 'lightgray'],
                     alpha=0.8)

        ax3.set_xlabel('Alignment Position (bp)', fontweight='bold')
        ax3.set_ylabel('Percentage (%)', fontweight='bold')
        ax3.set_title('Base Composition Across Alignment', fontweight='bold')
        ax3.set_ylim(0, 100)
        ax3.legend(loc='upper right', framealpha=0.9)
        ax3.grid(True, alpha=0.3)

        # Subplot 4: Pairwise identity matrix (bottom left)
        ax4 = fig.add_subplot(gs[2, 0])

        # Calculate pairwise identity
        identity_matrix = np.zeros((num_sequences, num_sequences))
        for i in range(num_sequences):
            for j in range(num_sequences):
                if i == j:
                    identity_matrix[i, j] = 100.0
                else:
                    seq_i = str(alignment[i].seq)
                    seq_j = str(alignment[j].seq)
                    matches = sum(a == b and a != '-' for a, b in zip(seq_i, seq_j))
                    valid_positions = sum(a != '-' and b != '-' for a, b in zip(seq_i, seq_j))
                    if valid_positions > 0:
                        identity_matrix[i, j] = (matches / valid_positions) * 100

        sns.heatmap(identity_matrix, annot=True if num_sequences <= 10 else False,
                   fmt='.1f', cmap='RdYlGn', center=80, vmin=50, vmax=100,
                   cbar_kws={'label': 'Identity (%)'}, ax=ax4, square=True)
        ax4.set_xlabel('Sequence ID', fontweight='bold')
        ax4.set_ylabel('Sequence ID', fontweight='bold')
        ax4.set_title('Pairwise Sequence Identity Matrix', fontweight='bold')

        # Subplot 5: Gap distribution (bottom right)
        ax5 = fig.add_subplot(gs[2, 1])

        gap_counts = [str(record.seq).count('-') for record in alignment]
        gap_percentages = [(count / alignment_length) * 100 for count in gap_counts]

        bars = ax5.bar(range(num_sequences), gap_percentages, color='lightcoral', edgecolor='black')
        ax5.set_xlabel('Sequence ID', fontweight='bold')
        ax5.set_ylabel('Gap Percentage (%)', fontweight='bold')
        ax5.set_title('Gap Distribution Per Sequence', fontweight='bold')
        ax5.grid(True, alpha=0.3, axis='y')

        # Add mean line
        mean_gap = np.mean(gap_percentages)
        ax5.axhline(y=mean_gap, color='red', linestyle='--', linewidth=2,
                   label=f'Mean: {mean_gap:.1f}%')
        ax5.legend()

        # Add summary text
        ax5.text(0.02, 0.98,
                f'Total sequences: {num_sequences}\n'
                f'Alignment length: {alignment_length} bp\n'
                f'Mean gaps: {mean_gap:.1f}%',
                transform=ax5.transAxes, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8),
                fontsize=10)

        # Overall title
        fig.suptitle('Multiple Sequence Alignment Analysis - 178bp Satellite Blocks',
                    fontsize=16, fontweight='bold', y=0.995)

        # Save figure
        output_file = os.path.join(output_dir, 'msa_visualization.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

        print(f"[SUCCESS] Saved MSA visualization: {output_file}")

    except Exception as e:
        print(f"[ERROR] Could not create MSA visualization: {e}")
        import traceback
        traceback.print_exc()


def main():
    if len(sys.argv) != 3:
        print("Usage: python3 13.4-analyze_mafft_alignment.py <alignment_file> <satellite_analysis.tsv>")
        sys.exit(1)

    alignment_file = sys.argv[1]
    analysis_file = sys.argv[2]

    # Get output directory from alignment file path
    output_dir = os.path.dirname(alignment_file)

    print("\n" + "="*80)
    print("=== ANALYZING MAFFT ALIGNMENT ===")
    print("="*80)
    print(f"[CONFIG] Alignment file: {alignment_file}")
    print(f"[CONFIG] Analysis file: {analysis_file}")
    print(f"[CONFIG] Output directory: {output_dir}")
    print("="*80)

    # Analyze alignment
    print(f"\n[INFO] Analyzing alignment...")
    alignment_stats = analyze_alignment_quality(alignment_file)

    if alignment_stats:
        print(f"\n[SUCCESS] Alignment Statistics:")
        print(f"  Sequences: {alignment_stats['num_sequences']}")
        print(f"  Length: {alignment_stats['alignment_length']} bp")
        print(f"  Avg conservation: {alignment_stats['avg_conservation']:.2%}")
        print(f"  Highly conserved positions: {alignment_stats['highly_conserved_positions']} ({alignment_stats['conservation_percentage']:.2%})")

        # Create MSA visualization
        create_msa_visualization(alignment_file, output_dir)

        # Update analysis file with conservation statistics
        if os.path.exists(analysis_file):
            print(f"\n[INFO] Updating {analysis_file} with alignment statistics...")
            df = pd.read_csv(analysis_file, sep='\t')

            # Add alignment quality columns
            df['alignment_included'] = True
            df['avg_block_conservation'] = alignment_stats['avg_conservation']
            df['alignment_num_sequences'] = alignment_stats['num_sequences']
            df['alignment_length'] = alignment_stats['alignment_length']

            # Save updated file
            df.to_csv(analysis_file, sep='\t', index=False)
            print(f"[SUCCESS] Updated {analysis_file}")
        else:
            print(f"[WARNING] Analysis file not found: {analysis_file}")
    else:
        print(f"[WARNING] Could not analyze alignment")

    print("\n" + "="*80)
    print("=== ALIGNMENT ANALYSIS COMPLETE ===")
    print("="*80)


if __name__ == "__main__":
    main()
