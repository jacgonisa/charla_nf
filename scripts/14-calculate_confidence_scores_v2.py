#!/usr/bin/env python3
"""
Improved Confidence Scoring with Diagnostics (v2)

This version:
1. Creates diagnostic plots of RAW metric distributions
2. Properly loads FULL PAF files from mapping directories
3. Fixes normalization based on actual data distributions
4. Separates SCO vs NCO confidence scores
5. Provides detailed diagnostic output

Author: Claude Code
Date: 2025-11-27
"""

import pandas as pd
import numpy as np
import argparse
import sys
import os
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

# Set publication-quality plotting
plt.rcParams['pdf.fonttype'] = 42  # TrueType fonts
plt.rcParams['font.size'] = 10
sns.set_style("white")


def load_full_paf_data(mapping_base_dir, read_ids):
    """
    Load FULL PAF files from all mapping modes in the mapping directories.
    This gives us complete alignment information including AS, NM, MAPQ.

    Returns: dict of {read_id: [paf_records]}
    """
    paf_data = defaultdict(list)

    if not os.path.exists(mapping_base_dir):
        return paf_data

    # Find all best_alignments.paf files in summary_* directories
    paf_files = list(Path(mapping_base_dir).rglob('summary_*/best_alignments.paf'))

    print(f"[INFO] Found {len(paf_files)} PAF files in {mapping_base_dir}", file=sys.stderr)

    for paf_file in paf_files:
        try:
            with open(paf_file, 'r') as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue

                    parts = line.strip().split('\t')
                    if len(parts) < 12:
                        continue

                    segment_id = parts[0]
                    # Extract base read ID (remove segment suffix like _1_2000_2500)
                    base_read_id = '_'.join(segment_id.split('_')[:-3]) if segment_id.count('_') >= 3 else segment_id

                    if base_read_id not in read_ids:
                        continue

                    paf_record = {
                        'segment_id': segment_id,
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
                        'mapq': int(parts[11]),
                        'mapping_mode': paf_file.parent.name  # e.g., summary_mm_sr
                    }

                    # Parse optional tags (AS, NM, etc.)
                    for tag in parts[12:]:
                        if tag.startswith('AS:i:'):
                            paf_record['AS'] = int(tag[5:])
                        elif tag.startswith('NM:i:'):
                            paf_record['NM'] = int(tag[5:])
                        elif tag.startswith('ms:i:'):
                            paf_record['ms'] = int(tag[5:])
                        elif tag.startswith('nn:i:'):
                            paf_record['nn'] = int(tag[5:])

                    paf_data[base_read_id].append(paf_record)

        except Exception as e:
            print(f"[WARNING] Error loading {paf_file}: {e}", file=sys.stderr)

    print(f"[INFO] Loaded alignment data for {len(paf_data)} reads", file=sys.stderr)
    return dict(paf_data)


def load_bed_data(bed_path, read_ids):
    """Load BED file with segment information."""
    bed_data = defaultdict(list)

    if not os.path.exists(bed_path):
        return bed_data

    try:
        with open(bed_path, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue

                parts = line.strip().split('\t')
                if len(parts) < 4:
                    continue

                read_id = parts[0]
                if read_id not in read_ids:
                    continue

                bed_data[read_id].append({
                    'start': int(parts[1]),
                    'end': int(parts[2]),
                    'region': parts[3],
                    'length': int(parts[2]) - int(parts[1])
                })

    except Exception as e:
        print(f"[WARNING] Error loading BED {bed_path}: {e}", file=sys.stderr)

    return dict(bed_data)


def load_cowidth_data(cowidth_path, read_ids):
    """Load crossover width data."""
    cowidth_data = {}

    if not os.path.exists(cowidth_path):
        return cowidth_data

    try:
        df = pd.read_csv(cowidth_path, sep='\t')
        for read_id in read_ids:
            if read_id in df['read'].values:
                records = df[df['read'] == read_id]
                cowidth_data[read_id] = records['width'].tolist()

    except Exception as e:
        print(f"[WARNING] Error loading cowidth {cowidth_path}: {e}", file=sys.stderr)

    return cowidth_data


def calculate_raw_metrics(read_id, summary_row, bed_data, paf_data, cowidth_data):
    """
    Calculate RAW metrics (before normalization).

    Returns: dict with all raw metric values
    """
    metrics = {
        'read_id': read_id,
        # Raw mapping consensus metrics
        'n_modes_best': 0,
        'n_modes_total': 5,
        'has_problems': False,
    }

    # 1. Mapping consensus (count "best" across modes)
    modes = ['summary_lrhq', 'summary_lrhqae', 'summary_mm_ont', 'summary_mm_sr', 'summary_wm']
    for mode in modes:
        if mode in summary_row:
            if summary_row[mode] in ['best', 'best;other']:
                metrics['n_modes_best'] += 1
            elif summary_row[mode] in ['chimeric', 'inversion', 'inconclusive_tie', 'wrong_ecotype']:
                metrics['has_problems'] = True

    # 2. Segment quality metrics (from BED)
    if read_id in bed_data:
        segments = bed_data[read_id]
        lengths = [s['length'] for s in segments]
        metrics['n_segments'] = len(segments)
        metrics['min_segment_length'] = min(lengths)
        metrics['mean_segment_length'] = np.mean(lengths)
        metrics['total_segment_length'] = sum(lengths)
    else:
        metrics['n_segments'] = 0
        metrics['min_segment_length'] = 0
        metrics['mean_segment_length'] = 0
        metrics['total_segment_length'] = 0

    # 3. Mapping quality metrics (from PAF)
    if read_id in paf_data:
        paf_records = paf_data[read_id]

        # MAPQ values
        mapqs = [r['mapq'] for r in paf_records]
        metrics['n_alignments'] = len(mapqs)
        metrics['mean_mapq'] = np.mean(mapqs)
        metrics['min_mapq'] = min(mapqs)
        metrics['max_mapq'] = max(mapqs)

        # Alignment identity (matches / block_len)
        identities = [r['matches'] / r['block_len'] if r['block_len'] > 0 else 0 for r in paf_records]
        metrics['mean_identity'] = np.mean(identities)
        metrics['min_identity'] = min(identities)

        # AS-NM scores (if available)
        as_nm_scores = []
        for r in paf_records:
            if 'AS' in r and 'NM' in r:
                aligned_len = r['qend'] - r['qstart']
                if aligned_len > 0:
                    as_nm_scores.append((r['AS'] - r['NM']) / aligned_len)

        if as_nm_scores:
            metrics['mean_as_nm_norm'] = np.mean(as_nm_scores)
            metrics['min_as_nm_norm'] = min(as_nm_scores)
        else:
            metrics['mean_as_nm_norm'] = None
            metrics['min_as_nm_norm'] = None

        # 4. Parent discrimination (Col vs Ler alignment quality difference)
        col_records = [r for r in paf_records if 'Col' in r['target']]
        ler_records = [r for r in paf_records if 'Ler' in r['target']]

        if col_records and ler_records:
            # Get best alignment score from each parent
            def get_aln_score(rec):
                if 'AS' in rec and 'NM' in rec:
                    return rec['AS'] - rec['NM']
                else:
                    return rec['matches']  # Fallback to raw matches

            best_col_score = max(get_aln_score(r) for r in col_records)
            best_ler_score = max(get_aln_score(r) for r in ler_records)

            # Discrimination is the relative difference
            if max(best_col_score, best_ler_score) > 0:
                metrics['parent_discrimination'] = abs(best_col_score - best_ler_score) / max(best_col_score, best_ler_score)
            else:
                metrics['parent_discrimination'] = 0.0

            metrics['best_col_score'] = best_col_score
            metrics['best_ler_score'] = best_ler_score
        else:
            metrics['parent_discrimination'] = None
            metrics['best_col_score'] = None
            metrics['best_ler_score'] = None

    else:
        # No PAF data available
        metrics['n_alignments'] = 0
        metrics['mean_mapq'] = None
        metrics['min_mapq'] = None
        metrics['max_mapq'] = None
        metrics['mean_identity'] = None
        metrics['min_identity'] = None
        metrics['mean_as_nm_norm'] = None
        metrics['min_as_nm_norm'] = None
        metrics['parent_discrimination'] = None
        metrics['best_col_score'] = None
        metrics['best_ler_score'] = None

    # 5. Crossover resolution (minimum gap width)
    if read_id in cowidth_data:
        widths = cowidth_data[read_id]
        metrics['n_crossovers'] = len(widths)
        metrics['min_co_width'] = min(widths)
        metrics['mean_co_width'] = np.mean(widths)
    else:
        metrics['n_crossovers'] = 0
        metrics['min_co_width'] = None
        metrics['mean_co_width'] = None

    return metrics


def normalize_metrics_adaptive(raw_metrics_df):
    """
    Normalize metrics based on their ACTUAL distributions in the data.

    Uses percentile-based normalization to handle skewed distributions.
    """
    normalized = raw_metrics_df.copy()

    # 1. Mapping consensus: simple ratio (already 0-1)
    normalized['mapping_consensus_score'] = raw_metrics_df['n_modes_best'] / 5.0
    # Penalize if problems detected
    normalized.loc[raw_metrics_df['has_problems'], 'mapping_consensus_score'] *= 0.7

    # 2. Segment quality: composite score
    # Min segment length (percentile-based)
    if raw_metrics_df['min_segment_length'].notna().any():
        p25 = raw_metrics_df['min_segment_length'].quantile(0.25)
        p75 = raw_metrics_df['min_segment_length'].quantile(0.75)
        normalized['min_seg_len_score'] = np.clip((raw_metrics_df['min_segment_length'] - p25) / (p75 - p25), 0, 1)
    else:
        normalized['min_seg_len_score'] = 0.5

    # Mean segment length
    if raw_metrics_df['mean_segment_length'].notna().any():
        p25 = raw_metrics_df['mean_segment_length'].quantile(0.25)
        p75 = raw_metrics_df['mean_segment_length'].quantile(0.75)
        normalized['mean_seg_len_score'] = np.clip((raw_metrics_df['mean_segment_length'] - p25) / (p75 - p25), 0, 1)
    else:
        normalized['mean_seg_len_score'] = 0.5

    # Number of segments (2 is ideal for simple crossover)
    normalized['n_seg_score'] = raw_metrics_df['n_segments'].apply(
        lambda x: 1.0 if x == 2 else 0.7 if x == 3 else 0.3
    )

    # Composite segment quality
    normalized['segment_quality_score'] = (
        0.4 * normalized['min_seg_len_score'] +
        0.3 * normalized['mean_seg_len_score'] +
        0.3 * normalized['n_seg_score']
    )

    # 3. Mapping quality: percentile-based normalization of MAPQ
    if raw_metrics_df['mean_mapq'].notna().any():
        p10 = raw_metrics_df['mean_mapq'].quantile(0.10)
        p90 = raw_metrics_df['mean_mapq'].quantile(0.90)
        normalized['mapq_score'] = np.clip((raw_metrics_df['mean_mapq'] - p10) / (p90 - p10), 0, 1)
    else:
        normalized['mapq_score'] = 0.5

    # Identity score
    if raw_metrics_df['mean_identity'].notna().any():
        # Identity should be high (>0.95 is good)
        normalized['identity_score'] = np.clip((raw_metrics_df['mean_identity'] - 0.90) / 0.08, 0, 1)
    else:
        normalized['identity_score'] = 0.5

    # Composite mapping quality
    normalized['mapping_quality_score'] = 0.6 * normalized['mapq_score'] + 0.4 * normalized['identity_score']

    # 4. Parent discrimination: percentile-based
    if raw_metrics_df['parent_discrimination'].notna().any():
        p25 = raw_metrics_df['parent_discrimination'].quantile(0.25)
        p75 = raw_metrics_df['parent_discrimination'].quantile(0.75)
        normalized['parent_discrimination_score'] = np.clip(
            (raw_metrics_df['parent_discrimination'] - p25) / (p75 - p25), 0, 1
        )
    else:
        normalized['parent_discrimination_score'] = 0.5

    # 5. Crossover resolution: smaller is better
    if raw_metrics_df['min_co_width'].notna().any():
        p90 = raw_metrics_df['min_co_width'].quantile(0.90)
        p10 = raw_metrics_df['min_co_width'].quantile(0.10)
        # Invert: smaller widths = higher scores
        normalized['crossover_resolution_score'] = np.clip((p90 - raw_metrics_df['min_co_width']) / (p90 - p10), 0, 1)
    else:
        normalized['crossover_resolution_score'] = 0.5

    # Fill NaN values with 0.5 (neutral)
    score_columns = [
        'mapping_consensus_score', 'segment_quality_score', 'mapping_quality_score',
        'parent_discrimination_score', 'crossover_resolution_score'
    ]
    for col in score_columns:
        normalized[col] = normalized[col].fillna(0.5)

    return normalized


def calculate_confidence_scores(normalized_df):
    """
    Calculate final confidence scores using weighted combination.
    """
    weights = {
        'mapping_consensus': 0.30,
        'segment_quality': 0.25,
        'mapping_quality': 0.20,
        'parent_discrimination': 0.15,
        'crossover_resolution': 0.10
    }

    confidence = (
        weights['mapping_consensus'] * normalized_df['mapping_consensus_score'] +
        weights['segment_quality'] * normalized_df['segment_quality_score'] +
        weights['mapping_quality'] * normalized_df['mapping_quality_score'] +
        weights['parent_discrimination'] * normalized_df['parent_discrimination_score'] +
        weights['crossover_resolution'] * normalized_df['crossover_resolution_score']
    )

    return confidence


def create_diagnostic_plots(raw_metrics_df, normalized_df, output_dir):
    """
    Create comprehensive diagnostic plots showing RAW and normalized distributions.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Figure 1: Raw metric distributions
    fig, axes = plt.subplots(3, 3, figsize=(18, 16))
    fig.suptitle('Raw Metric Distributions (Before Normalization)', fontsize=16, fontweight='bold')

    # Plot 1: Mapping consensus
    ax = axes[0, 0]
    if raw_metrics_df['n_modes_best'].notna().any():
        ax.hist(raw_metrics_df['n_modes_best'], bins=6, alpha=0.7, color='steelblue', edgecolor='black')
        ax.axvline(raw_metrics_df['n_modes_best'].mean(), color='red', linestyle='--', linewidth=2,
                  label=f'Mean: {raw_metrics_df["n_modes_best"].mean():.2f}')
    ax.set_xlabel('Number of modes calling "best"', fontweight='bold')
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_title('Mapping Consensus (0-5 modes)', fontweight='bold')
    ax.legend()

    # Plot 2: Min segment length
    ax = axes[0, 1]
    if raw_metrics_df['min_segment_length'].notna().any():
        data = raw_metrics_df['min_segment_length'].dropna()
        ax.hist(data, bins=50, alpha=0.7, color='steelblue', edgecolor='black')
        ax.axvline(data.mean(), color='red', linestyle='--', linewidth=2,
                  label=f'Mean: {data.mean():.0f} bp')
        ax.axvline(data.quantile(0.25), color='orange', linestyle=':', linewidth=2,
                  label=f'Q25: {data.quantile(0.25):.0f} bp')
        ax.axvline(data.quantile(0.75), color='green', linestyle=':', linewidth=2,
                  label=f'Q75: {data.quantile(0.75):.0f} bp')
    ax.set_xlabel('Minimum segment length (bp)', fontweight='bold')
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_title('Segment Quality: Min Length', fontweight='bold')
    ax.legend(fontsize=8)

    # Plot 3: Mean segment length
    ax = axes[0, 2]
    if raw_metrics_df['mean_segment_length'].notna().any():
        data = raw_metrics_df['mean_segment_length'].dropna()
        ax.hist(data, bins=50, alpha=0.7, color='steelblue', edgecolor='black')
        ax.axvline(data.mean(), color='red', linestyle='--', linewidth=2,
                  label=f'Mean: {data.mean():.0f} bp')
    ax.set_xlabel('Mean segment length (bp)', fontweight='bold')
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_title('Segment Quality: Mean Length', fontweight='bold')
    ax.legend()

    # Plot 4: Number of segments
    ax = axes[1, 0]
    if raw_metrics_df['n_segments'].notna().any():
        counts = raw_metrics_df['n_segments'].value_counts().sort_index()
        ax.bar(counts.index, counts.values, alpha=0.7, color='steelblue', edgecolor='black')
    ax.set_xlabel('Number of segments', fontweight='bold')
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_title('Segment Count Distribution', fontweight='bold')

    # Plot 5: MAPQ distribution
    ax = axes[1, 1]
    if raw_metrics_df['mean_mapq'].notna().any():
        data = raw_metrics_df['mean_mapq'].dropna()
        ax.hist(data, bins=50, alpha=0.7, color='steelblue', edgecolor='black')
        ax.axvline(data.mean(), color='red', linestyle='--', linewidth=2,
                  label=f'Mean: {data.mean():.1f}')
        ax.axvline(data.quantile(0.10), color='orange', linestyle=':', linewidth=2,
                  label=f'P10: {data.quantile(0.10):.1f}')
        ax.axvline(data.quantile(0.90), color='green', linestyle=':', linewidth=2,
                  label=f'P90: {data.quantile(0.90):.1f}')
    ax.set_xlabel('Mean MAPQ', fontweight='bold')
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_title('Mapping Quality: MAPQ', fontweight='bold')
    ax.legend(fontsize=8)

    # Plot 6: Alignment identity
    ax = axes[1, 2]
    if raw_metrics_df['mean_identity'].notna().any():
        data = raw_metrics_df['mean_identity'].dropna()
        ax.hist(data, bins=50, alpha=0.7, color='steelblue', edgecolor='black')
        ax.axvline(data.mean(), color='red', linestyle='--', linewidth=2,
                  label=f'Mean: {data.mean():.3f}')
    ax.set_xlabel('Mean alignment identity', fontweight='bold')
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_title('Mapping Quality: Identity', fontweight='bold')
    ax.legend()

    # Plot 7: Parent discrimination
    ax = axes[2, 0]
    if raw_metrics_df['parent_discrimination'].notna().any():
        data = raw_metrics_df['parent_discrimination'].dropna()
        ax.hist(data, bins=50, alpha=0.7, color='steelblue', edgecolor='black')
        ax.axvline(data.mean(), color='red', linestyle='--', linewidth=2,
                  label=f'Mean: {data.mean():.3f}')
        ax.axvline(data.quantile(0.25), color='orange', linestyle=':', linewidth=2,
                  label=f'Q25: {data.quantile(0.25):.3f}')
        ax.axvline(data.quantile(0.75), color='green', linestyle=':', linewidth=2,
                  label=f'Q75: {data.quantile(0.75):.3f}')
    ax.set_xlabel('Parent discrimination ratio', fontweight='bold')
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_title('Parent Discrimination (Col vs Ler)', fontweight='bold')
    ax.legend(fontsize=8)

    # Plot 8: Crossover width
    ax = axes[2, 1]
    if raw_metrics_df['min_co_width'].notna().any():
        data = raw_metrics_df['min_co_width'].dropna()
        ax.hist(data, bins=50, alpha=0.7, color='steelblue', edgecolor='black')
        ax.axvline(data.mean(), color='red', linestyle='--', linewidth=2,
                  label=f'Mean: {data.mean():.0f} bp')
        ax.axvline(data.quantile(0.10), color='orange', linestyle=':', linewidth=2,
                  label=f'P10: {data.quantile(0.10):.0f} bp')
        ax.axvline(data.quantile(0.90), color='green', linestyle=':', linewidth=2,
                  label=f'P90: {data.quantile(0.90):.0f} bp')
    ax.set_xlabel('Minimum crossover width (bp)', fontweight='bold')
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_title('Crossover Resolution', fontweight='bold')
    ax.legend(fontsize=8)

    # Plot 9: Number of crossovers
    ax = axes[2, 2]
    if raw_metrics_df['n_crossovers'].notna().any():
        counts = raw_metrics_df['n_crossovers'].value_counts().sort_index()
        ax.bar(counts.index, counts.values, alpha=0.7, color='steelblue', edgecolor='black')
    ax.set_xlabel('Number of crossovers', fontweight='bold')
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_title('Crossover Count per Read', fontweight='bold')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'diagnostic_raw_metrics.png'), dpi=300, bbox_inches='tight')
    plt.close()

    print(f"[SUCCESS] Saved raw metrics diagnostic plot", file=sys.stderr)

    # Figure 2: Normalized score distributions
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('Normalized Score Distributions (After Adaptive Normalization)', fontsize=16, fontweight='bold')

    score_cols = [
        ('mapping_consensus_score', 'Mapping Consensus'),
        ('segment_quality_score', 'Segment Quality'),
        ('mapping_quality_score', 'Mapping Quality'),
        ('parent_discrimination_score', 'Parent Discrimination'),
        ('crossover_resolution_score', 'Crossover Resolution')
    ]

    for idx, (col, title) in enumerate(score_cols):
        row = idx // 3
        col_idx = idx % 3
        ax = axes[row, col_idx]

        if col in normalized_df.columns:
            data = normalized_df[col].dropna()
            ax.hist(data, bins=50, alpha=0.7, color='steelblue', edgecolor='black')
            ax.axvline(data.mean(), color='red', linestyle='--', linewidth=2,
                      label=f'Mean: {data.mean():.3f}')
            ax.axvline(data.median(), color='orange', linestyle=':', linewidth=2,
                      label=f'Median: {data.median():.3f}')

        ax.set_xlabel('Score (0-1)', fontweight='bold')
        ax.set_ylabel('Count', fontweight='bold')
        ax.set_title(title, fontweight='bold')
        ax.set_xlim(0, 1)
        ax.legend()

    # Final confidence score
    ax = axes[1, 2]
    if 'confidence_score' in normalized_df.columns:
        data = normalized_df['confidence_score'].dropna()
        ax.hist(data, bins=50, alpha=0.7, color='green', edgecolor='black')
        ax.axvline(data.mean(), color='red', linestyle='--', linewidth=2,
                  label=f'Mean: {data.mean():.3f}')
        ax.axvline(data.median(), color='orange', linestyle=':', linewidth=2,
                  label=f'Median: {data.median():.3f}')

    ax.set_xlabel('Score (0-1)', fontweight='bold')
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_title('Final Confidence Score', fontweight='bold')
    ax.set_xlim(0, 1)
    ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'diagnostic_normalized_scores.png'), dpi=300, bbox_inches='tight')
    plt.close()

    print(f"[SUCCESS] Saved normalized scores diagnostic plot", file=sys.stderr)


def identify_sco_nco(raw_metrics_df):
    """
    Identify reads as SCO (single crossover) or NCO (non-crossover/complex).

    Simple heuristic:
    - SCO: n_crossovers == 1 AND n_segments == 2
    - NCO/Complex: anything else
    """
    sco_mask = (raw_metrics_df['n_crossovers'] == 1) & (raw_metrics_df['n_segments'] == 2)

    co_type = ['NCO/Complex'] * len(raw_metrics_df)
    for idx in sco_mask[sco_mask].index:
        co_type[idx] = 'SCO'

    return co_type


def main():
    parser = argparse.ArgumentParser(
        description='Improved confidence scoring with diagnostics (v2)',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--summary-arms', required=True, help='Summary table CSV for ARMS')
    parser.add_argument('--summary-cen', required=True, help='Summary table CSV for CEN')
    parser.add_argument('--bed-arms', required=True, help='BED file for ARMS')
    parser.add_argument('--bed-cen', required=True, help='BED file for CEN')
    parser.add_argument('--mapping-base-arms', required=True, help='Mapping base directory for ARMS (contains summary_*/ subdirs)')
    parser.add_argument('--mapping-base-cen', required=True, help='Mapping base directory for CEN (contains summary_*/ subdirs)')
    parser.add_argument('--cowidth-arms', required=True, help='Cowidth file for ARMS')
    parser.add_argument('--cowidth-cen', required=True, help='Cowidth file for CEN')
    parser.add_argument('--output', required=True, help='Output TSV file')
    parser.add_argument('--diagnostic-dir', required=True, help='Directory for diagnostic plots')

    args = parser.parse_args()

    # Create output directory
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    os.makedirs(args.diagnostic_dir, exist_ok=True)

    all_raw_metrics = []

    # Process ARMS
    print(f"[INFO] Processing ARMS reads...", file=sys.stderr)
    if os.path.exists(args.summary_arms):
        summary = pd.read_csv(args.summary_arms)
        read_ids = set(summary['read_id'].values)

        # Load all data
        bed_data = load_bed_data(args.bed_arms, read_ids)
        paf_data = load_full_paf_data(args.mapping_base_arms, read_ids)
        cowidth_data = load_cowidth_data(args.cowidth_arms, read_ids)

        # Calculate raw metrics for each read
        for _, row in summary.iterrows():
            read_id = row['read_id']
            metrics = calculate_raw_metrics(read_id, row, bed_data, paf_data, cowidth_data)
            metrics['region_type'] = 'ARMS'
            all_raw_metrics.append(metrics)

        print(f"[INFO] Processed {len(all_raw_metrics)} ARMS reads", file=sys.stderr)

    # Process CEN
    print(f"[INFO] Processing CEN reads...", file=sys.stderr)
    if os.path.exists(args.summary_cen):
        summary = pd.read_csv(args.summary_cen)
        read_ids = set(summary['read_id'].values)

        # Load all data
        bed_data = load_bed_data(args.bed_cen, read_ids)
        paf_data = load_full_paf_data(args.mapping_base_cen, read_ids)
        cowidth_data = load_cowidth_data(args.cowidth_cen, read_ids)

        arms_count = len(all_raw_metrics)

        # Calculate raw metrics for each read
        for _, row in summary.iterrows():
            read_id = row['read_id']
            metrics = calculate_raw_metrics(read_id, row, bed_data, paf_data, cowidth_data)
            metrics['region_type'] = 'CEN'
            all_raw_metrics.append(metrics)

        print(f"[INFO] Processed {len(all_raw_metrics) - arms_count} CEN reads", file=sys.stderr)

    if not all_raw_metrics:
        print("[ERROR] No reads processed!", file=sys.stderr)
        sys.exit(1)

    # Create DataFrame of raw metrics
    raw_metrics_df = pd.DataFrame(all_raw_metrics)

    # Identify SCO vs NCO
    raw_metrics_df['crossover_type'] = identify_sco_nco(raw_metrics_df)

    # Save raw metrics for inspection
    raw_metrics_file = args.output.replace('.tsv', '_raw_metrics.tsv')
    raw_metrics_df.to_csv(raw_metrics_file, sep='\t', index=False, float_format='%.4f')
    print(f"[SUCCESS] Saved raw metrics to {raw_metrics_file}", file=sys.stderr)

    # Normalize metrics
    print(f"[INFO] Normalizing metrics...", file=sys.stderr)
    normalized_df = normalize_metrics_adaptive(raw_metrics_df)

    # Calculate confidence scores
    normalized_df['confidence_score'] = calculate_confidence_scores(normalized_df)

    # Create diagnostic plots
    print(f"[INFO] Creating diagnostic plots...", file=sys.stderr)
    create_diagnostic_plots(raw_metrics_df, normalized_df, args.diagnostic_dir)

    # Prepare output
    output_cols = [
        'read_id', 'confidence_score', 'region_type', 'crossover_type',
        'mapping_consensus_score', 'segment_quality_score', 'mapping_quality_score',
        'parent_discrimination_score', 'crossover_resolution_score'
    ]

    output_df = normalized_df[output_cols].copy()
    output_df = output_df.sort_values('confidence_score', ascending=False)

    # Save
    output_df.to_csv(args.output, sep='\t', index=False, float_format='%.4f')

    # Print summary
    print("\n" + "="*80, file=sys.stderr)
    print("CONFIDENCE SCORING SUMMARY (v2 - Improved)", file=sys.stderr)
    print("="*80, file=sys.stderr)
    print(f"Total reads scored: {len(output_df)}", file=sys.stderr)
    print(f"  - ARMS: {(output_df['region_type'] == 'ARMS').sum()}", file=sys.stderr)
    print(f"  - CEN: {(output_df['region_type'] == 'CEN').sum()}", file=sys.stderr)
    print(f"\nCrossover types:", file=sys.stderr)
    print(f"  - SCO (single crossover): {(output_df['crossover_type'] == 'SCO').sum()}", file=sys.stderr)
    print(f"  - NCO/Complex: {(output_df['crossover_type'] == 'NCO/Complex').sum()}", file=sys.stderr)
    print(f"\nConfidence score distribution:", file=sys.stderr)
    print(f"  High (≥0.8): {(output_df['confidence_score'] >= 0.8).sum()} ({(output_df['confidence_score'] >= 0.8).sum() / len(output_df) * 100:.1f}%)", file=sys.stderr)
    print(f"  Medium (0.5-0.8): {((output_df['confidence_score'] >= 0.5) & (output_df['confidence_score'] < 0.8)).sum()} ({((output_df['confidence_score'] >= 0.5) & (output_df['confidence_score'] < 0.8)).sum() / len(output_df) * 100:.1f}%)", file=sys.stderr)
    print(f"  Low (<0.5): {(output_df['confidence_score'] < 0.5).sum()} ({(output_df['confidence_score'] < 0.5).sum() / len(output_df) * 100:.1f}%)", file=sys.stderr)

    print(f"\nScore statistics:", file=sys.stderr)
    print(f"  Mean: {output_df['confidence_score'].mean():.4f}", file=sys.stderr)
    print(f"  Median: {output_df['confidence_score'].median():.4f}", file=sys.stderr)
    print(f"  Std: {output_df['confidence_score'].std():.4f}", file=sys.stderr)
    print(f"  Min: {output_df['confidence_score'].min():.4f}", file=sys.stderr)
    print(f"  Max: {output_df['confidence_score'].max():.4f}", file=sys.stderr)

    print(f"\nComponent score means:", file=sys.stderr)
    print(f"  Mapping consensus: {output_df['mapping_consensus_score'].mean():.4f}", file=sys.stderr)
    print(f"  Segment quality: {output_df['segment_quality_score'].mean():.4f}", file=sys.stderr)
    print(f"  Mapping quality: {output_df['mapping_quality_score'].mean():.4f}", file=sys.stderr)
    print(f"  Parent discrimination: {output_df['parent_discrimination_score'].mean():.4f}", file=sys.stderr)
    print(f"  Crossover resolution: {output_df['crossover_resolution_score'].mean():.4f}", file=sys.stderr)

    # SCO vs NCO breakdown
    print(f"\nSCO vs NCO confidence breakdown:", file=sys.stderr)
    for co_type in ['SCO', 'NCO/Complex']:
        subset = output_df[output_df['crossover_type'] == co_type]
        if len(subset) > 0:
            print(f"  {co_type}:", file=sys.stderr)
            print(f"    Count: {len(subset)}", file=sys.stderr)
            print(f"    Mean confidence: {subset['confidence_score'].mean():.4f}", file=sys.stderr)
            print(f"    High confidence: {(subset['confidence_score'] >= 0.8).sum()} ({(subset['confidence_score'] >= 0.8).sum() / len(subset) * 100:.1f}%)", file=sys.stderr)

    print("="*80, file=sys.stderr)
    print(f"\nOutput files:", file=sys.stderr)
    print(f"  - Confidence scores: {args.output}", file=sys.stderr)
    print(f"  - Raw metrics: {raw_metrics_file}", file=sys.stderr)
    print(f"  - Diagnostic plots: {args.diagnostic_dir}/", file=sys.stderr)
    print("="*80, file=sys.stderr)


if __name__ == '__main__':
    main()
