#!/usr/bin/env python3
"""
SCO vs NCO Comparison Analysis

Compares Single Crossovers (SCO) vs Non-Crossovers/Complex (NCO) reads
to understand quality differences and features.

Key metrics:
- Crossover width (resolution)
- Segment lengths
- Mapping consensus
- Confidence scores
- Genomic distribution

Author: Claude Code
Date: 2025-11-27
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import sys

# Try to import scipy, but make it optional
try:
    from scipy import stats
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    print("WARNING: scipy not available - statistical tests will be skipped")

# Publication settings
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.size'] = 10
plt.rcParams['figure.dpi'] = 300
sns.set_style("white")

# Color scheme
SCO_COLOR = '#2ecc71'  # Green
NCO_COLOR = '#e74c3c'  # Red


def load_data(confidence_file, raw_metrics_file, cowidth_file):
    """Load all data files."""

    # Confidence scores (TSV)
    conf_df = pd.read_csv(confidence_file, sep='\t')

    # Raw metrics (CSV - from ARMS summary)
    # Detect separator automatically
    if raw_metrics_file.endswith('.csv'):
        raw_df = pd.read_csv(raw_metrics_file, sep=',')
    else:
        raw_df = pd.read_csv(raw_metrics_file, sep='\t')

    # Cowidths (TSV)
    cowidth_df = pd.read_csv(cowidth_file, sep='\t')

    # Merge confidence with raw metrics
    merged = pd.merge(conf_df, raw_df, on='read_id', how='left', suffixes=('', '_raw'))

    # Calculate min crossover width per read from cowidth data
    # Cowidth file has columns: read, segment_pair, width
    if 'read' in cowidth_df.columns and 'width' in cowidth_df.columns:
        min_cowidths = cowidth_df.groupby('read')['width'].min().reset_index()
        min_cowidths.columns = ['read_id', 'min_co_width']

        # Merge with main dataframe
        merged = pd.merge(merged, min_cowidths, on='read_id', how='left')

    return merged, cowidth_df


def plot_sco_nco_overview(df, output_prefix):
    """
    Create comprehensive 8-panel comparison of SCO vs NCO.
    """
    fig, axes = plt.subplots(3, 3, figsize=(18, 16))
    fig.suptitle('SCO vs NCO Comprehensive Comparison', fontsize=16, fontweight='bold', y=0.995)

    # Separate SCO and NCO
    sco = df[df['crossover_type'] == 'SCO']
    nco = df[df['crossover_type'] != 'SCO']

    # Panel 1: Count comparison
    ax = axes[0, 0]
    counts = df['crossover_type'].value_counts()
    colors = [SCO_COLOR if 'SCO' in x else NCO_COLOR for x in counts.index]
    ax.bar(range(len(counts)), counts.values, color=colors, alpha=0.7, edgecolor='black', linewidth=2)
    ax.set_xticks(range(len(counts)))
    ax.set_xticklabels(counts.index, rotation=0)
    ax.set_ylabel('Count', fontweight='bold', fontsize=12)
    ax.set_title(f'Read Counts\n(Total: {len(df):,})', fontweight='bold')
    for i, v in enumerate(counts.values):
        ax.text(i, v + max(counts.values)*0.02, f'{v:,}\n({v/len(df)*100:.1f}%)',
               ha='center', va='bottom', fontweight='bold')

    # Panel 2: Confidence score comparison
    ax = axes[0, 1]
    data_to_plot = [sco['confidence_score'].dropna(), nco['confidence_score'].dropna()]
    bp = ax.boxplot(data_to_plot, tick_labels=['SCO', 'NCO/Complex'],
                    patch_artist=True, widths=0.6)
    for patch, color in zip(bp['boxes'], [SCO_COLOR, NCO_COLOR]):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    ax.set_ylabel('Confidence Score', fontweight='bold', fontsize=12)
    ax.set_title('Confidence Score Distribution', fontweight='bold')
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3, axis='y')

    # Add stats
    if HAS_SCIPY and len(data_to_plot[0]) > 0 and len(data_to_plot[1]) > 0:
        t_stat, p_val = stats.ttest_ind(data_to_plot[0], data_to_plot[1])
        ax.text(0.5, 0.95, f'p = {p_val:.2e}' if p_val < 0.001 else f'p = {p_val:.3f}',
               transform=ax.transAxes, ha='center', va='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Panel 3: Crossover width (only for SCO)
    ax = axes[0, 2]
    if 'min_co_width' in sco.columns and sco['min_co_width'].notna().any():
        widths = sco['min_co_width'].dropna()
        ax.hist(widths, bins=50, color=SCO_COLOR, alpha=0.7, edgecolor='black', linewidth=0.5)
        ax.axvline(widths.median(), color='red', linestyle='--', linewidth=2,
                  label=f'Median: {widths.median():.0f} bp')
        ax.axvline(widths.quantile(0.25), color='orange', linestyle=':', linewidth=2,
                  label=f'Q25: {widths.quantile(0.25):.0f} bp')
        ax.axvline(widths.quantile(0.75), color='green', linestyle=':', linewidth=2,
                  label=f'Q75: {widths.quantile(0.75):.0f} bp')
        ax.set_xlabel('Crossover Width (bp)', fontweight='bold', fontsize=12)
        ax.set_ylabel('Count', fontweight='bold', fontsize=12)
        ax.set_title(f'SCO Resolution\n(n={len(widths):,})', fontweight='bold')
        ax.legend(fontsize=9)
        ax.set_xlim(0, min(5000, widths.quantile(0.99)))
    else:
        ax.text(0.5, 0.5, 'No crossover width data', ha='center', va='center',
               transform=ax.transAxes, fontsize=12)
        ax.set_title('SCO Resolution', fontweight='bold')

    # Panel 4: Min segment length comparison
    ax = axes[1, 0]
    if 'min_segment_length' in df.columns:
        data_to_plot = [sco['min_segment_length'].dropna(), nco['min_segment_length'].dropna()]
        bp = ax.boxplot(data_to_plot, tick_labels=['SCO', 'NCO/Complex'],
                        patch_artist=True, widths=0.6)
        for patch, color in zip(bp['boxes'], [SCO_COLOR, NCO_COLOR]):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        ax.set_ylabel('Min Segment Length (bp)', fontweight='bold', fontsize=12)
        ax.set_title('Minimum Segment Length\n(smaller = worse quality)', fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')

        # Add stats
        if HAS_SCIPY and len(data_to_plot[0]) > 0 and len(data_to_plot[1]) > 0:
            t_stat, p_val = stats.ttest_ind(data_to_plot[0], data_to_plot[1])
            ax.text(0.5, 0.95, f'p = {p_val:.2e}' if p_val < 0.001 else f'p = {p_val:.3f}',
                   transform=ax.transAxes, ha='center', va='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Panel 5: Mean segment length comparison
    ax = axes[1, 1]
    if 'mean_segment_length' in df.columns:
        data_to_plot = [sco['mean_segment_length'].dropna(), nco['mean_segment_length'].dropna()]
        bp = ax.boxplot(data_to_plot, tick_labels=['SCO', 'NCO/Complex'],
                        patch_artist=True, widths=0.6)
        for patch, color in zip(bp['boxes'], [SCO_COLOR, NCO_COLOR]):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        ax.set_ylabel('Mean Segment Length (bp)', fontweight='bold', fontsize=12)
        ax.set_title('Average Segment Length\n(larger = better coverage)', fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')

        if len(data_to_plot[0]) > 0 and len(data_to_plot[1]) > 0:
            t_stat, p_val = stats.ttest_ind(data_to_plot[0], data_to_plot[1])
            ax.text(0.5, 0.95, f'p = {p_val:.2e}' if p_val < 0.001 else f'p = {p_val:.3f}',
                   transform=ax.transAxes, ha='center', va='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Panel 6: Number of classification segments (from BED)
    ax = axes[1, 2]
    if 'n_classification_segments' in df.columns:
        sco_counts = sco['n_classification_segments'].value_counts().sort_index()
        nco_counts = nco['n_classification_segments'].value_counts().sort_index()

        # Only show up to 10 segments for clarity
        x = np.arange(1, min(11, max(sco_counts.index.max(), nco_counts.index.max()) + 2))

        sco_vals = [sco_counts.get(i, 0) for i in x]
        nco_vals = [nco_counts.get(i, 0) for i in x]

        width = 0.35
        ax.bar(x - width/2, sco_vals, width, label='SCO', color=SCO_COLOR, alpha=0.7, edgecolor='black')
        ax.bar(x + width/2, nco_vals, width, label='NCO/Complex', color=NCO_COLOR, alpha=0.7, edgecolor='black')

        ax.set_xlabel('Number of Parental Segments', fontweight='bold', fontsize=12)
        ax.set_ylabel('Count', fontweight='bold', fontsize=12)
        ax.set_title('Classification Segment Distribution\n(SCO=2, NCO=3, DCO=4, Complex=5+)', fontweight='bold')
        ax.legend()
        ax.set_xticks(x)

        # Add expected classification markers
        ax.axvline(2.5, color='red', linestyle='--', alpha=0.3, linewidth=1)
        ax.text(2, max(max(sco_vals), max(nco_vals)) * 0.9, 'SCO', ha='center', fontsize=8, color='red')
        ax.text(3, max(max(sco_vals), max(nco_vals)) * 0.9, 'NCO', ha='center', fontsize=8, color='red')
        ax.text(4, max(max(sco_vals), max(nco_vals)) * 0.9, 'DCO', ha='center', fontsize=8, color='red')

    # Panel 7: Mapping consensus comparison
    ax = axes[2, 0]
    if 'n_modes_best' in df.columns:
        data_to_plot = [sco['n_modes_best'].dropna(), nco['n_modes_best'].dropna()]
        bp = ax.boxplot(data_to_plot, tick_labels=['SCO', 'NCO/Complex'],
                        patch_artist=True, widths=0.6)
        for patch, color in zip(bp['boxes'], [SCO_COLOR, NCO_COLOR]):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        ax.set_ylabel('N Modes Calling "Best"', fontweight='bold', fontsize=12)
        ax.set_title('Mapping Consensus\n(5 = all modes agree)', fontweight='bold')
        ax.set_ylim(0, 5.5)
        ax.grid(True, alpha=0.3, axis='y')

        if len(data_to_plot[0]) > 0 and len(data_to_plot[1]) > 0:
            t_stat, p_val = stats.ttest_ind(data_to_plot[0], data_to_plot[1])
            ax.text(0.5, 0.95, f'p = {p_val:.2e}' if p_val < 0.001 else f'p = {p_val:.3f}',
                   transform=ax.transAxes, ha='center', va='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Panel 8: Region type (ARMS vs CEN)
    ax = axes[2, 1]
    if 'region_type' in df.columns:
        region_counts = pd.crosstab(df['crossover_type'], df['region_type'])
        region_counts.plot(kind='bar', stacked=True, ax=ax, color=['#3498db', '#9b59b6'],
                          alpha=0.7, edgecolor='black', linewidth=1)
        ax.set_xlabel('Crossover Type', fontweight='bold', fontsize=12)
        ax.set_ylabel('Count', fontweight='bold', fontsize=12)
        ax.set_title('Genomic Distribution\n(ARMS vs CEN)', fontweight='bold')
        ax.legend(title='Region', fontsize=9)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=0)

    # Panel 9: Summary statistics table
    ax = axes[2, 2]
    ax.axis('off')

    # Calculate summary stats
    sco_class_segs = sco['n_classification_segments'].mean() if 'n_classification_segments' in sco.columns else 'N/A'
    nco_class_segs = nco['n_classification_segments'].mean() if 'n_classification_segments' in nco.columns else 'N/A'
    sco_bed_segs = sco['n_bed_segments'].mean() if 'n_bed_segments' in sco.columns else 'N/A'
    nco_bed_segs = nco['n_bed_segments'].mean() if 'n_bed_segments' in nco.columns else 'N/A'

    summary_text = f"""
    Summary Statistics

    SCO Reads: {len(sco):,} ({len(sco)/len(df)*100:.1f}%)
    - Mean confidence: {sco['confidence_score'].mean():.3f}
    - Median CO width: {sco['min_co_width'].median():.0f} bp
    - Parental segments: {sco_class_segs if isinstance(sco_class_segs, str) else f'{sco_class_segs:.1f}'}
    - BED segments: {sco_bed_segs if isinstance(sco_bed_segs, str) else f'{sco_bed_segs:.1f}'}

    NCO/Complex: {len(nco):,} ({len(nco)/len(df)*100:.1f}%)
    - Mean confidence: {nco['confidence_score'].mean():.3f}
    - Parental segments: {nco_class_segs if isinstance(nco_class_segs, str) else f'{nco_class_segs:.1f}'}
    - BED segments: {nco_bed_segs if isinstance(nco_bed_segs, str) else f'{nco_bed_segs:.1f}'}

    High Confidence (≥0.8):
    - SCO: {(sco['confidence_score'] >= 0.8).sum():,} ({(sco['confidence_score'] >= 0.8).sum()/len(sco)*100:.1f}%)
    - NCO: {(nco['confidence_score'] >= 0.8).sum():,} ({(nco['confidence_score'] >= 0.8).sum()/len(nco)*100:.1f}%)
    """

    ax.text(0.1, 0.9, summary_text, transform=ax.transAxes,
           fontsize=10, verticalalignment='top', fontfamily='monospace',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"[SUCCESS] Saved SCO vs NCO comparison plot")


def plot_resolution_analysis(df, cowidth_df, output_prefix):
    """
    Detailed analysis of crossover resolution (width).
    """
    # Merge cowidth data
    sco = df[df['crossover_type'] == 'SCO'].copy()

    if len(sco) == 0:
        print("[WARNING] No SCO reads found for resolution analysis")
        return

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle('SCO Resolution Analysis (Crossover Width)', fontsize=16, fontweight='bold')

    # Panel 1: Crossover width distribution
    ax = axes[0, 0]
    widths = sco['min_co_width'].dropna()
    # Filter out zero or negative widths for log scale
    widths = widths[widths > 0]

    if len(widths) > 10:  # Need at least 10 points for meaningful histogram
        # Histogram with log bins
        bins = np.logspace(np.log10(widths.min()), np.log10(widths.quantile(0.99)), 50)
        ax.hist(widths, bins=bins, color=SCO_COLOR, alpha=0.7, edgecolor='black', linewidth=0.5)
        ax.set_xscale('log')
        ax.axvline(widths.median(), color='red', linestyle='--', linewidth=2,
                  label=f'Median: {widths.median():.0f} bp')
        ax.set_xlabel('Crossover Width (bp, log scale)', fontweight='bold')
        ax.set_ylabel('Count', fontweight='bold')
        ax.set_title(f'Width Distribution (n={len(widths):,})', fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
    else:
        ax.text(0.5, 0.5, f'Insufficient data\n(n={len(widths)})',
                ha='center', va='center', transform=ax.transAxes)

    # Panel 2: Width vs confidence
    ax = axes[0, 1]
    if 'min_co_width' in sco.columns and 'confidence_score' in sco.columns:
        valid = sco[sco['min_co_width'].notna() & sco['confidence_score'].notna()]
        # Filter out zero/negative widths for log scale
        valid = valid[valid['min_co_width'] > 0]

        if len(valid) > 10:
            ax.scatter(valid['min_co_width'], valid['confidence_score'],
                      alpha=0.3, s=20, color=SCO_COLOR, edgecolor='black', linewidth=0.5)

            # Add trend line (with error handling)
            try:
                log_widths = np.log10(valid['min_co_width'])
                # Check for inf/nan values
                valid_for_fit = ~(np.isinf(log_widths) | np.isnan(log_widths))
                if valid_for_fit.sum() > 10:
                    z = np.polyfit(log_widths[valid_for_fit],
                                  valid['confidence_score'][valid_for_fit], 1)
                    p = np.poly1d(z)
                    x_trend = np.logspace(np.log10(valid['min_co_width'].min()),
                                         np.log10(valid['min_co_width'].max()), 100)
                    ax.plot(x_trend, p(np.log10(x_trend)), "r--", alpha=0.8, linewidth=2,
                           label=f'Trend: y={z[0]:.3f}*log(x)+{z[1]:.3f}')
            except (np.linalg.LinAlgError, ValueError, RuntimeWarning) as e:
                print(f"[WARNING] Could not fit trend line: {e}", file=sys.stderr)

            ax.set_xscale('log')
            ax.set_xlabel('Crossover Width (bp, log scale)', fontweight='bold')
            ax.set_ylabel('Confidence Score', fontweight='bold')
            ax.set_title('Resolution vs Confidence', fontweight='bold')
            ax.set_ylim(0, 1)
            ax.legend()
            ax.grid(True, alpha=0.3)
        else:
            ax.text(0.5, 0.5, f'Insufficient data\n(n={len(valid)})',
                    ha='center', va='center', transform=ax.transAxes)

            # Correlation
            corr = valid[['min_co_width', 'confidence_score']].corr().iloc[0, 1]
            ax.text(0.05, 0.95, f'Correlation: {corr:.3f}',
                   transform=ax.transAxes, ha='left', va='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Panel 3: Width by region
    ax = axes[1, 0]
    if 'region_type' in sco.columns:
        arms_widths = sco[sco['region_type'] == 'ARMS']['min_co_width'].dropna()
        cen_widths = sco[sco['region_type'] == 'CEN']['min_co_width'].dropna()

        bp = ax.boxplot([arms_widths, cen_widths], tick_labels=['ARMS', 'CEN'],
                        patch_artist=True, widths=0.6)
        for patch in bp['boxes']:
            patch.set_facecolor(SCO_COLOR)
            patch.set_alpha(0.7)

        ax.set_ylabel('Crossover Width (bp)', fontweight='bold')
        ax.set_title('Resolution by Genomic Region', fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')

        if HAS_SCIPY and len(arms_widths) > 0 and len(cen_widths) > 0:
            t_stat, p_val = stats.ttest_ind(arms_widths, cen_widths)
            ax.text(0.5, 0.95, f'p = {p_val:.2e}' if p_val < 0.001 else f'p = {p_val:.3f}',
                   transform=ax.transAxes, ha='center', va='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Panel 4: Cumulative distribution
    ax = axes[1, 1]
    if len(widths) > 0:
        sorted_widths = np.sort(widths)
        cumulative = np.arange(1, len(sorted_widths) + 1) / len(sorted_widths) * 100

        ax.plot(sorted_widths, cumulative, color=SCO_COLOR, linewidth=2)
        ax.axvline(widths.median(), color='red', linestyle='--', linewidth=2,
                  label=f'Median: {widths.median():.0f} bp')
        ax.axhline(50, color='gray', linestyle=':', alpha=0.5)

        # Mark percentiles
        for p in [25, 75, 90]:
            val = widths.quantile(p/100)
            ax.axvline(val, color='orange', linestyle=':', alpha=0.5)
            ax.text(val, 5, f'P{p}', rotation=90, va='bottom', ha='right')

        ax.set_xlabel('Crossover Width (bp)', fontweight='bold')
        ax.set_ylabel('Cumulative Percentage (%)', fontweight='bold')
        ax.set_title('Cumulative Distribution', fontweight='bold')
        ax.set_xlim(0, min(5000, widths.quantile(0.99)))
        ax.set_ylim(0, 100)
        ax.legend()
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_resolution_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"[SUCCESS] Saved resolution analysis plot")


def create_summary_table(df, output_file):
    """Create detailed summary statistics table."""

    summary_stats = []

    for co_type in ['SCO', 'NCO/Complex']:
        subset = df[df['crossover_type'] == co_type] if co_type == 'SCO' else df[df['crossover_type'] != 'SCO']

        stats_dict = {
            'Crossover_Type': co_type,
            'Count': len(subset),
            'Percentage': len(subset) / len(df) * 100,
            'Mean_Confidence': subset['confidence_score'].mean(),
            'Median_Confidence': subset['confidence_score'].median(),
            'High_Confidence_Count': (subset['confidence_score'] >= 0.8).sum(),
            'High_Confidence_Pct': (subset['confidence_score'] >= 0.8).sum() / len(subset) * 100 if len(subset) > 0 else 0,
        }

        if 'min_co_width' in subset.columns:
            stats_dict['Median_CO_Width'] = subset['min_co_width'].median()
            stats_dict['Mean_CO_Width'] = subset['min_co_width'].mean()

        if 'min_segment_length' in subset.columns:
            stats_dict['Mean_Min_Seg_Length'] = subset['min_segment_length'].mean()
            stats_dict['Mean_Avg_Seg_Length'] = subset['mean_segment_length'].mean()

        if 'n_segments' in subset.columns:
            stats_dict['Mean_N_Segments'] = subset['n_segments'].mean()

        if 'n_modes_best' in subset.columns:
            stats_dict['Mean_Consensus'] = subset['n_modes_best'].mean()

        summary_stats.append(stats_dict)

    summary_df = pd.DataFrame(summary_stats)
    summary_df.to_csv(output_file, sep='\t', index=False, float_format='%.4f')

    print(f"[SUCCESS] Saved summary table to {output_file}")

    return summary_df


def extract_segment_info_from_bed(bed_file):
    """
    Extract segment information from BED file.

    Returns dict with:
    - n_bed_segments: number of segments per read (from BED)
    - n_classification_segments: collapsed parental segments (for classification)
    - segment_pattern: parental pattern (e.g., "CL", "CLC")
    """
    import re

    # Read BED file
    bed_data = []
    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                bed_data.append({
                    'read_id': parts[0],
                    'start': int(parts[1]),
                    'end': int(parts[2]),
                    'segment_type': parts[3]  # e.g., "LA4", "CA3", "LC2"
                })

    bed_df = pd.DataFrame(bed_data)

    # Group by read and extract info
    segment_info = {}
    for read_id, group in bed_df.groupby('read_id'):
        # Sort by start position
        group = group.sort_values('start')

        # Count total BED segments
        n_bed_segments = len(group)

        # Extract parental types (C or L) from segment types
        parental_types = []
        for seg_type in group['segment_type']:
            # Find first C or L in segment type (e.g., "LA4" → "L", "CA3" → "C")
            for char in seg_type:
                if char in ['C', 'L']:
                    parental_types.append(char)
                    break

        # Create parental pattern
        parental_pattern = ''.join(parental_types)

        # Collapse consecutive same parental type (CCC → C, LLL → L)
        collapsed = re.sub(r'C+', 'C', parental_pattern)
        collapsed = re.sub(r'L+', 'L', collapsed)

        # Count classification segments (length of collapsed pattern)
        n_classification_segments = len(collapsed)

        segment_info[read_id] = {
            'n_bed_segments': n_bed_segments,
            'n_classification_segments': n_classification_segments,
            'segment_pattern': collapsed
        }

    return pd.DataFrame.from_dict(segment_info, orient='index').reset_index().rename(columns={'index': 'read_id'})


def main():
    parser = argparse.ArgumentParser(
        description='Analyze SCO vs NCO reads comparison'
    )

    parser.add_argument('--confidence', required=True, help='Confidence scores TSV')
    parser.add_argument('--raw-metrics', required=True, help='Raw metrics TSV')
    parser.add_argument('--cowidth', required=True, help='Cowidth file (combined ARMS+CEN)')
    parser.add_argument('--segments-bed', required=True, help='Segments BED file from SEGMENT_READS')
    parser.add_argument('--output-prefix', required=True, help='Output prefix for plots')

    args = parser.parse_args()

    print("[INFO] Loading data...")
    df, cowidth_df = load_data(args.confidence, args.raw_metrics, args.cowidth)

    # Extract segment information from BED file
    print("[INFO] Extracting segment information from BED file...")
    segment_info_df = extract_segment_info_from_bed(args.segments_bed)

    # Merge segment info into main dataframe
    df = df.merge(segment_info_df, on='read_id', how='left')
    print(f"[INFO] Added BED segment info for {segment_info_df['read_id'].nunique():,} reads")

    print(f"[INFO] Total reads: {len(df):,}")
    print(f"[INFO] SCO reads: {(df['crossover_type'] == 'SCO').sum():,}")
    print(f"[INFO] NCO/Complex reads: {(df['crossover_type'] != 'SCO').sum():,}")

    print("[INFO] Creating SCO vs NCO comparison plots...")
    plot_sco_nco_overview(df, args.output_prefix)

    print("[INFO] Creating resolution analysis plots...")
    plot_resolution_analysis(df, cowidth_df, args.output_prefix)

    print("[INFO] Creating summary table...")
    summary_df = create_summary_table(df, f'{args.output_prefix}_summary.tsv')

    print("\n" + "="*80)
    print("SCO vs NCO ANALYSIS COMPLETE")
    print("="*80)
    print(summary_df.to_string(index=False))
    print("="*80)


if __name__ == '__main__':
    main()
