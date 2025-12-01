#!/usr/bin/env python3
"""
Confidence Scoring V3 - Redesigned from Scratch

NEW METRICS (no more penalizing small segments!):
1. MAPQ Score (30%) - Mapping quality from PAF files
2. Identity Score (25%) - Alignment identity
3. Mapping Consensus (20%) - Agreement across algorithms
4. Parental Balance (15%) - Col/Ler ratio (for SCOs)
5. Coverage Completeness (10%) - Fraction that's NOT undetermined

Author: Claude Code
Date: 2025-11-27
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import sys
import os
from pathlib import Path

# Publication settings
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.size'] = 10
sns.set_style("white")


def load_parental_proportions(parental_file):
    """
    Load parental proportions table from read structure analysis.
    Contains: col_prop, ler_prop, undetermined_prop, mean_mapq, mean_identity, n_switches
    """
    if not os.path.exists(parental_file):
        print(f"[WARNING] Parental proportions file not found: {parental_file}", file=sys.stderr)
        return pd.DataFrame()

    df = pd.read_csv(parental_file, sep='\t')
    print(f"[INFO] Loaded parental data for {len(df)} reads", file=sys.stderr)
    return df


def load_summary_table(summary_file):
    """Load summary table with mapping consensus."""
    if not os.path.exists(summary_file):
        print(f"[WARNING] Summary table not found: {summary_file}", file=sys.stderr)
        return pd.DataFrame()

    df = pd.read_csv(summary_file)
    print(f"[INFO] Loaded summary table for {len(df)} reads", file=sys.stderr)
    return df


def load_bed_data(bed_file):
    """Load BED file for segment information."""
    if not os.path.exists(bed_file):
        return pd.DataFrame()

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
                    'region': parts[3],
                    'length': int(parts[2]) - int(parts[1])
                })

    return pd.DataFrame(bed_data)


def calculate_mapq_score(mean_mapq):
    """
    Score based on MAPQ (Phred-scaled mapping quality).

    MAPQ scale:
    - 60: 1 in 1,000,000 error (perfect!)
    - 40: 1 in 10,000 error (excellent)
    - 20: 1 in 100 error (good)
    - 10: 1 in 10 error (poor)
    - 0: unmapped or very poor

    Returns: 0-1 score
    """
    if pd.isna(mean_mapq):
        return 0.5  # Neutral if no data

    # Sigmoid-like scoring
    if mean_mapq >= 50:
        return 1.0
    elif mean_mapq >= 30:
        return 0.8 + 0.2 * (mean_mapq - 30) / 20
    elif mean_mapq >= 20:
        return 0.6 + 0.2 * (mean_mapq - 20) / 10
    elif mean_mapq >= 10:
        return 0.3 + 0.3 * (mean_mapq - 10) / 10
    else:
        return 0.3 * (mean_mapq / 10)


def calculate_identity_score(mean_identity):
    """
    Score based on alignment identity.

    Identity scale:
    - >0.995: Excellent
    - 0.99-0.995: Very good
    - 0.98-0.99: Good
    - 0.95-0.98: OK
    - <0.95: Poor

    Returns: 0-1 score
    """
    if pd.isna(mean_identity):
        return 0.5  # Neutral if no data

    if mean_identity >= 0.995:
        return 1.0
    elif mean_identity >= 0.99:
        return 0.8 + 0.2 * (mean_identity - 0.99) / 0.005
    elif mean_identity >= 0.98:
        return 0.6 + 0.2 * (mean_identity - 0.98) / 0.01
    elif mean_identity >= 0.95:
        return 0.3 + 0.3 * (mean_identity - 0.95) / 0.03
    else:
        return 0.3 * (mean_identity / 0.95)


def calculate_consensus_score(summary_row):
    """
    Score based on mapping consensus (agreement across algorithms).
    Same as before - this one works well!
    """
    modes = ['summary_lrhq', 'summary_lrhqae', 'summary_mm_ont', 'summary_mm_sr', 'summary_wm']

    best_count = 0
    for mode in modes:
        if mode in summary_row and summary_row[mode] in ['best', 'best;other']:
            best_count += 1

    consensus_score = best_count / 5.0

    # Penalize problematic flags
    problematic_flags = ['chimeric', 'inversion', 'inconclusive_tie', 'wrong_ecotype']
    for mode in modes:
        if mode in summary_row and summary_row[mode] in problematic_flags:
            consensus_score *= 0.7
            break

    return consensus_score


def calculate_parental_balance_score(col_prop, ler_prop, undetermined_prop, crossover_type):
    """
    Score based on parental balance.

    For SCOs:
    - Ideal: ~50/50 Col/Ler
    - Good: 40/60 to 60/40
    - Poor: Very skewed (>80/20)

    For NCOs:
    - Expected to be skewed (small tract in one parent background)
    - Less stringent scoring

    Returns: 0-1 score
    """
    if pd.isna(col_prop) or pd.isna(ler_prop):
        return 0.5

    # Calculate balance (how close to 50/50)
    balance = 1.0 - abs(col_prop - ler_prop)  # 1.0 = perfect 50/50, 0.0 = 100/0

    if crossover_type == 'SCO':
        # SCOs should be balanced
        if balance >= 0.8:  # 45/55 to 55/45
            return 1.0
        elif balance >= 0.6:  # 40/60 to 60/40
            return 0.8 + 0.2 * (balance - 0.6) / 0.2
        elif balance >= 0.4:  # 30/70 to 70/30
            return 0.5 + 0.3 * (balance - 0.4) / 0.2
        else:
            return 0.5 * (balance / 0.4)
    else:
        # NCOs can be skewed (gene conversion)
        # Less stringent - just check it's not ALL one parent
        if balance >= 0.2:  # At least some of each parent
            return 0.8 + 0.2 * (balance / 0.8)
        else:
            return 0.5 + 0.3 * (balance / 0.2)


def calculate_coverage_completeness_score(undetermined_prop):
    """
    Score based on how much of the read is determined (NOT grey).

    - Low undetermined (<5%) = excellent k-mer coverage
    - Medium undetermined (5-20%) = good
    - High undetermined (>20%) = repetitive/ambiguous regions

    Returns: 0-1 score
    """
    if pd.isna(undetermined_prop):
        return 0.5

    determined_prop = 1.0 - undetermined_prop

    if determined_prop >= 0.95:
        return 1.0
    elif determined_prop >= 0.8:
        return 0.7 + 0.3 * (determined_prop - 0.8) / 0.15
    elif determined_prop >= 0.5:
        return 0.4 + 0.3 * (determined_prop - 0.5) / 0.3
    else:
        return 0.4 * (determined_prop / 0.5)


def identify_crossover_type_mapping(n_segments, n_switches):
    """
    Classify crossover type based on MAPPING structure (for comparison).

    This is less reliable than k-mer-based classification.

    SCO: 2 segments, 1 switch (classic single crossover)
    NCO/Complex: anything else
    """
    if n_segments == 2 and n_switches == 1:
        return 'SCO'
    else:
        return 'NCO/Complex'


def classify_from_bed_pattern(bed_file, read_id):
    """
    Classify read type from BED file haplotype pattern.

    This gives us classification AFTER mapping validation.
    Reads in BED file have passed mapping quality checks.
    """
    # This will be loaded once and cached
    # For now, return Unknown - will be filled from merged data
    return 'Unknown'


def classify_read_type_from_segments(segments_str):
    """
    Classify read based on segment pattern string.

    Input: comma-separated segments like "325LA4,319LA3"

    Logic:
    - Extract parental type (C or L) from each segment
    - Collapse consecutive same parent
    - Classify pattern:
      * C or L → Parental
      * CL or LC → SCO
      * CLC or LCL → NCO
      * CLCL or LCLC → DCO
      * More switches → Complex
    """
    if not segments_str or pd.isna(segments_str):
        return 'Unknown'

    # Extract parental types from segments
    # Segment format: <number><parental><region><chr> e.g., "325LA4" = 325 bp, Ler, ARMS, Chr4
    # Need to find the first letter (C or L) in the segment
    segments = segments_str.split(',')
    parental_types = []
    for seg in segments:
        seg = seg.strip()
        # Find first letter that is C or L
        for char in seg:
            if char in ['C', 'L']:
                parental_types.append(char)
                break
        else:
            parental_types.append('X')

    parental_pattern = ''.join(parental_types)

    # Collapse consecutive same parental type
    import re
    collapsed = re.sub(r'C+', 'C', parental_pattern)
    collapsed = re.sub(r'L+', 'L', collapsed)

    # Classify based on number of segments (switches between parents)
    # User clarification: 3 segments = NCO, >3 segments = Complex
    num_segments = len(collapsed)

    if num_segments == 1:
        # Single parent (C or L)
        return 'Parental'
    elif num_segments == 2:
        # One switch (CL or LC)
        return 'SCO'
    elif num_segments == 3:
        # Two switches (CLC or LCL) - this is NCO (gene conversion)
        return 'NCO'
    elif num_segments == 4:
        # Three switches (CLCL or LCLC) - double crossover
        return 'DCO'
    else:
        # More than 4 segments = Complex
        return 'Complex'


def normalize_hybrid_type(hybrid_type):
    """
    Normalize hybrid_type from curated profiles to standard categories.

    From curated profiles (BEFORE mapping):
    - SCO → SCO
    - NCO → NCO
    - GC-CO → NCO (gene conversion)
    - None → Parental
    - Unknown → Unknown
    """
    if hybrid_type == 'SCO':
        return 'SCO'
    elif hybrid_type in ['NCO', 'GC-CO']:
        return 'NCO'
    elif hybrid_type == 'None':
        return 'Parental'
    else:
        return 'Unknown'


def calculate_confidence_scores(merged_df):
    """
    Calculate confidence scores for all reads using new metrics.
    """
    results = []

    for _, row in merged_df.iterrows():
        # Classification AFTER mapping (from ultra_curated_profile segments)
        # This is what passed mapping validation!
        ultra_profile = row.get('ultra_curated_profile', '')
        crossover_type_after = classify_read_type_from_segments(ultra_profile)

        # Classification BEFORE mapping (from hybrid_type in curated profiles)
        crossover_type_before = normalize_hybrid_type(row.get('hybrid_type', 'Unknown'))

        # Mapping-based classification (from PAF alignment patterns)
        n_segments = row.get('n_segments', 2)
        n_switches = row.get('n_switches', 1)
        crossover_type_mapping = identify_crossover_type_mapping(n_segments, n_switches)

        # Check concordance between AFTER mapping and PAF-based
        concordant = (crossover_type_after == crossover_type_mapping) or \
                    (crossover_type_after in ['NCO', 'DCO', 'Complex'] and crossover_type_mapping == 'NCO/Complex')

        # Use AFTER MAPPING classification for scoring (the validated one!)
        crossover_type = crossover_type_after

        # Calculate component scores
        mapq_score = calculate_mapq_score(row.get('mean_mapq'))
        identity_score = calculate_identity_score(row.get('mean_identity'))
        consensus_score = calculate_consensus_score(row)
        balance_score = calculate_parental_balance_score(
            row.get('col_prop'), row.get('ler_prop'),
            row.get('undetermined_prop'), crossover_type
        )
        coverage_score = calculate_coverage_completeness_score(row.get('undetermined_prop'))

        # Weighted combination
        weights = {
            'mapq': 0.30,
            'identity': 0.25,
            'consensus': 0.20,
            'balance': 0.15,
            'coverage': 0.10
        }

        confidence = (
            weights['mapq'] * mapq_score +
            weights['identity'] * identity_score +
            weights['consensus'] * consensus_score +
            weights['balance'] * balance_score +
            weights['coverage'] * coverage_score
        )

        results.append({
            'read_id': row['read_id'],
            'confidence_score': confidence,
            'mapq_score': mapq_score,
            'identity_score': identity_score,
            'consensus_score': consensus_score,
            'balance_score': balance_score,
            'coverage_score': coverage_score,
            'crossover_type': crossover_type,  # AFTER mapping (validated) - USED FOR ANALYSIS
            'crossover_type_before': crossover_type_before,  # BEFORE mapping (k-mer only)
            'crossover_type_mapping': crossover_type_mapping,  # PAF-based (for comparison)
            'classification_concordant': concordant,  # Does AFTER match PAF?
            'region_type': row.get('region_type', 'unknown'),
            # Raw values for reference
            'mean_mapq': row.get('mean_mapq'),
            'mean_identity': row.get('mean_identity'),
            'col_prop': row.get('col_prop'),
            'ler_prop': row.get('ler_prop'),
            'undetermined_prop': row.get('undetermined_prop'),
            'n_segments': n_segments,
            'n_switches': n_switches
        })

    return pd.DataFrame(results)


def create_diagnostic_plots(results_df, output_dir):
    """
    Create diagnostic plots showing score distributions.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Figure 1: Component score distributions
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('Confidence Score Components (V3 - Redesigned)', fontsize=16, fontweight='bold')

    components = [
        ('mapq_score', 'MAPQ Score (30%)', 'Based on mapping quality'),
        ('identity_score', 'Identity Score (25%)', 'Based on alignment identity'),
        ('consensus_score', 'Consensus Score (20%)', 'Agreement across mappers'),
        ('balance_score', 'Balance Score (15%)', 'Parental balance (SCO: 50/50)'),
        ('coverage_score', 'Coverage Score (10%)', 'Fraction determined (not grey)')
    ]

    for idx, (col, title, desc) in enumerate(components):
        ax = axes[idx // 3, idx % 3]

        data = results_df[col].dropna()
        ax.hist(data, bins=50, alpha=0.7, color='steelblue', edgecolor='black', linewidth=0.5)
        ax.axvline(data.mean(), color='red', linestyle='--', linewidth=2,
                  label=f'Mean: {data.mean():.3f}')
        ax.axvline(data.median(), color='orange', linestyle=':', linewidth=2,
                  label=f'Median: {data.median():.3f}')

        ax.set_xlabel('Score', fontweight='bold')
        ax.set_ylabel('Count', fontweight='bold')
        ax.set_title(f'{title}\n{desc}', fontweight='bold', fontsize=10)
        ax.set_xlim(0, 1)
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')

    # Final confidence score
    ax = axes[1, 2]
    data = results_df['confidence_score'].dropna()
    ax.hist(data, bins=50, alpha=0.7, color='green', edgecolor='black', linewidth=0.5)
    ax.axvline(data.mean(), color='red', linestyle='--', linewidth=2,
              label=f'Mean: {data.mean():.3f}')
    ax.axvline(data.median(), color='orange', linestyle=':', linewidth=2,
              label=f'Median: {data.median():.3f}')

    ax.set_xlabel('Score', fontweight='bold')
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_title('Final Confidence Score\n(Weighted combination)', fontweight='bold', fontsize=10)
    ax.set_xlim(0, 1)
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'confidence_v3_component_scores.png'), dpi=300, bbox_inches='tight')
    plt.close()

    print(f"[SUCCESS] Saved component scores plot", file=sys.stderr)

    # Figure 2: Raw metric distributions
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('Raw Metric Distributions (V3)', fontsize=16, fontweight='bold')

    # MAPQ
    ax = axes[0, 0]
    data = results_df['mean_mapq'].dropna()
    ax.hist(data, bins=50, alpha=0.7, color='purple', edgecolor='black', linewidth=0.5)
    ax.axvline(data.median(), color='red', linestyle='--', linewidth=2,
              label=f'Median: {data.median():.1f}')
    ax.set_xlabel('Mean MAPQ', fontweight='bold')
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_title('MAPQ Distribution\n(60 = 1 in 1M error)', fontweight='bold')
    ax.legend()

    # Identity
    ax = axes[0, 1]
    data = results_df['mean_identity'].dropna()
    ax.hist(data, bins=50, alpha=0.7, color='orange', edgecolor='black', linewidth=0.5)
    ax.axvline(data.median(), color='red', linestyle='--', linewidth=2,
              label=f'Median: {data.median():.4f}')
    ax.set_xlabel('Mean Identity', fontweight='bold')
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_title('Alignment Identity\n(1.0 = 100% match)', fontweight='bold')
    ax.legend()

    # Parental balance (Col vs Ler)
    ax = axes[0, 2]
    sco = results_df[results_df['crossover_type'] == 'SCO']
    if len(sco) > 0:
        ax.scatter(sco['col_prop'], sco['ler_prop'], alpha=0.3, s=20,
                  c=sco['confidence_score'], cmap='viridis', edgecolor='black', linewidth=0.3)
        ax.plot([0, 1], [1, 0], 'r--', alpha=0.5, label='Col+Ler=100%')
        ax.plot([0.5], [0.5], 'ro', markersize=10, label='Ideal (50/50)')
        ax.set_xlabel('Col Proportion', fontweight='bold')
        ax.set_ylabel('Ler Proportion', fontweight='bold')
        ax.set_title('Parental Balance (SCO only)\nIdeal = 50/50', fontweight='bold')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.legend()
        plt.colorbar(ax.collections[0], ax=ax, label='Confidence')

    # Undetermined proportion
    ax = axes[1, 0]
    data = results_df['undetermined_prop'].dropna()
    ax.hist(data, bins=50, alpha=0.7, color='grey', edgecolor='black', linewidth=0.5)
    ax.axvline(data.median(), color='red', linestyle='--', linewidth=2,
              label=f'Median: {data.median():.3f}')
    ax.set_xlabel('Undetermined Proportion', fontweight='bold')
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_title('Grey Regions\n(Lower = better k-mer coverage)', fontweight='bold')
    ax.legend()

    # Number of switches
    ax = axes[1, 1]
    switch_counts = results_df['n_switches'].value_counts().sort_index()
    ax.bar(switch_counts.index, switch_counts.values, alpha=0.7, color='green',
          edgecolor='black', linewidth=1)
    ax.set_xlabel('Number of Switches', fontweight='bold')
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_title('Recombination Events\n(1 = SCO, >1 = complex)', fontweight='bold')

    # Confidence by crossover type
    ax = axes[1, 2]
    sco_conf = results_df[results_df['crossover_type'] == 'SCO']['confidence_score']
    nco_conf = results_df[results_df['crossover_type'] != 'SCO']['confidence_score']

    bp = ax.boxplot([sco_conf.dropna(), nco_conf.dropna()],
                    labels=['SCO', 'NCO/Complex'],
                    patch_artist=True, widths=0.6)
    for patch, color in zip(bp['boxes'], ['#2ecc71', '#e74c3c']):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax.set_ylabel('Confidence Score', fontweight='bold')
    ax.set_title('Confidence by Type\n(V3 scoring)', fontweight='bold')
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'confidence_v3_raw_metrics.png'), dpi=300, bbox_inches='tight')
    plt.close()

    print(f"[SUCCESS] Saved raw metrics plot", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='Calculate confidence scores V3 (redesigned metrics)'
    )

    parser.add_argument('--parental-proportions', required=True,
                       help='Parental proportions TSV from read structure analysis')
    parser.add_argument('--curated-profiles', required=True,
                       help='Curated hybrid profiles TSV (for TRUE classification)')
    parser.add_argument('--threshold', type=int, required=True,
                       help='Threshold value to filter curated profiles')
    parser.add_argument('--summary-arms', required=True, help='ARMS summary table')
    parser.add_argument('--summary-cen', required=True, help='CEN summary table')
    parser.add_argument('--bed-arms', help='ARMS BED file (optional)')
    parser.add_argument('--bed-cen', help='CEN BED file (optional)')
    parser.add_argument('--output', required=True, help='Output confidence scores TSV')
    parser.add_argument('--diagnostic-dir', required=True, help='Diagnostic plots directory')

    args = parser.parse_args()

    print("[INFO] Loading data...", file=sys.stderr)

    # Load curated profiles for classification comparison
    print(f"[INFO] Loading curated profiles (threshold={args.threshold})...", file=sys.stderr)
    curated_df = pd.read_csv(args.curated_profiles, sep='\t')

    # Filter for correct threshold and select relevant columns
    curated_df = curated_df[curated_df['threshold'] == args.threshold][
        ['readname', 'ultra_curated_profile', 'hybrid_type', 'ectopic']
    ]
    curated_df.columns = ['read_id', 'ultra_curated_profile', 'hybrid_type', 'ectopic']

    print(f"[INFO] Loaded {len(curated_df)} curated profiles", file=sys.stderr)

    # Count classification BEFORE mapping
    before_counts = curated_df['hybrid_type'].value_counts()
    print(f"\n[INFO] Classification BEFORE mapping (k-mer only):", file=sys.stderr)
    for htype, count in before_counts.items():
        print(f"  {htype}: {count}", file=sys.stderr)

    # Load parental proportions (has MAPQ, identity, col/ler props)
    parental_df = load_parental_proportions(args.parental_proportions)

    if len(parental_df) == 0:
        print("[ERROR] No parental data loaded!", file=sys.stderr)
        sys.exit(1)

    # Load summary tables (has mapping consensus)
    arms_summary = load_summary_table(args.summary_arms)
    cen_summary = load_summary_table(args.summary_cen)

    # Combine summaries
    arms_summary['region_type'] = 'ARMS'
    cen_summary['region_type'] = 'CEN'
    summary_df = pd.concat([arms_summary, cen_summary], ignore_index=True)

    # Merge parental + summary
    merged_df = parental_df.merge(summary_df, left_on='base_read_id', right_on='read_id', how='left')

    # Merge with curated profiles to get TRUE classification
    merged_df = merged_df.merge(curated_df, left_on='base_read_id', right_on='read_id', how='left', suffixes=('', '_curated'))

    # If read_id_y exists (from merge), use it; otherwise use base_read_id
    if 'read_id_y' in merged_df.columns:
        merged_df['read_id'] = merged_df['read_id_y'].fillna(merged_df['base_read_id'])
    elif 'read_id' not in merged_df.columns:
        merged_df['read_id'] = merged_df['base_read_id']

    print(f"[INFO] Merged data for {len(merged_df)} reads", file=sys.stderr)

    # Calculate confidence scores
    print("[INFO] Calculating confidence scores V3...", file=sys.stderr)
    results_df = calculate_confidence_scores(merged_df)

    # Create diagnostic plots
    print("[INFO] Creating diagnostic plots...", file=sys.stderr)
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    create_diagnostic_plots(results_df, args.diagnostic_dir)

    # Save results
    results_df.to_csv(args.output, sep='\t', index=False, float_format='%.4f')
    print(f"[SUCCESS] Saved confidence scores to {args.output}", file=sys.stderr)

    # Print summary
    print("\n" + "="*80, file=sys.stderr)
    print("CONFIDENCE SCORING V3 SUMMARY", file=sys.stderr)
    print("="*80, file=sys.stderr)
    print(f"Total reads: {len(results_df)}", file=sys.stderr)

    # Classification AFTER mapping (FINAL, used for analysis)
    print(f"\nClassification AFTER mapping (validated by mappers):", file=sys.stderr)
    after_counts = results_df['crossover_type'].value_counts()
    for ctype, count in after_counts.items():
        pct = 100 * count / len(results_df)
        print(f"  {ctype}: {count} ({pct:.1f}%)", file=sys.stderr)

    # Comparison with BEFORE mapping
    print(f"\nClassification concordance (BEFORE vs AFTER mapping):", file=sys.stderr)
    concordant_classification = (results_df['crossover_type'] == results_df['crossover_type_before']).sum()
    print(f"  Matching: {concordant_classification} ({100*concordant_classification/len(results_df):.1f}%)", file=sys.stderr)
    print(f"  Changed: {len(results_df)-concordant_classification} ({100*(len(results_df)-concordant_classification)/len(results_df):.1f}%)", file=sys.stderr)
    print(f"\nConfidence distribution:", file=sys.stderr)
    print(f"  High (≥0.8): {(results_df['confidence_score'] >= 0.8).sum()} ({(results_df['confidence_score'] >= 0.8).sum()/len(results_df)*100:.1f}%)", file=sys.stderr)
    print(f"  Medium (0.5-0.8): {((results_df['confidence_score'] >= 0.5) & (results_df['confidence_score'] < 0.8)).sum()} ({((results_df['confidence_score'] >= 0.5) & (results_df['confidence_score'] < 0.8)).sum()/len(results_df)*100:.1f}%)", file=sys.stderr)
    print(f"  Low (<0.5): {(results_df['confidence_score'] < 0.5).sum()} ({(results_df['confidence_score'] < 0.5).sum()/len(results_df)*100:.1f}%)", file=sys.stderr)
    print(f"\nMean scores:", file=sys.stderr)
    print(f"  Overall confidence: {results_df['confidence_score'].mean():.3f}", file=sys.stderr)
    print(f"  MAPQ score: {results_df['mapq_score'].mean():.3f}", file=sys.stderr)
    print(f"  Identity score: {results_df['identity_score'].mean():.3f}", file=sys.stderr)
    print(f"  Consensus score: {results_df['consensus_score'].mean():.3f}", file=sys.stderr)
    print(f"  Balance score: {results_df['balance_score'].mean():.3f}", file=sys.stderr)
    print(f"  Coverage score: {results_df['coverage_score'].mean():.3f}", file=sys.stderr)
    print("="*80, file=sys.stderr)


if __name__ == '__main__':
    main()
