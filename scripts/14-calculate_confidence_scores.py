#!/usr/bin/env python3
"""
Calculate confidence scores for detected crossovers.

This script assigns a confidence score (0.0-1.0) to each crossover based on:
1. Mapping consensus (30%): Agreement across multiple mapping algorithms
2. Segment quality (25%): Segment length and k-mer support
3. Mapping quality (20%): MAPQ and alignment scores
4. Parent discrimination (15%): Clarity of Col vs Ler assignment
5. Crossover resolution (10%): Precision of breakpoint localization

Output: TSV file with confidence scores and component scores per read.
"""

import pandas as pd
import numpy as np
import argparse
import sys
import os
from pathlib import Path


def mapping_consensus_score(row):
    """
    Score based on agreement across 5 mapping modes.

    Returns: Score 0.0-1.0
    """
    modes = ['summary_lrhq', 'summary_lrhqae', 'summary_mm_ont', 'summary_mm_sr', 'summary_wm']

    # Count how many modes say "best" or "best;other"
    best_count = 0
    for mode in modes:
        if mode in row and row[mode] in ['best', 'best;other']:
            best_count += 1

    # Base consensus score
    consensus_score = best_count / 5.0

    # Penalize if any mode flagged as problematic
    problematic_flags = ['chimeric', 'inversion', 'inconclusive_tie', 'wrong_ecotype']
    has_problems = False
    for mode in modes:
        if mode in row and row[mode] in problematic_flags:
            has_problems = True
            break

    if has_problems:
        consensus_score *= 0.7  # Penalize if ANY mode found issues

    return consensus_score


def segment_quality_score(bed_records):
    """
    Score based on segment lengths and k-mer support.

    Returns: Score 0.0-1.0
    """
    if not bed_records or len(bed_records) == 0:
        return 0.0

    # Calculate segment lengths
    segment_lengths = [rec['end'] - rec['start'] for rec in bed_records]

    if not segment_lengths:
        return 0.0

    min_length = min(segment_lengths)
    mean_length = np.mean(segment_lengths)
    num_segments = len(segment_lengths)

    # A. Minimum segment length score
    if min_length >= 200:
        min_length_score = 1.0
    elif min_length >= 100:
        min_length_score = 0.7
    elif min_length >= 50:
        min_length_score = 0.4
    else:
        min_length_score = 0.1

    # B. Average segment length score
    if mean_length >= 2000:
        mean_length_score = 1.0
    elif mean_length >= 1000:
        mean_length_score = 0.8
    elif mean_length >= 500:
        mean_length_score = 0.6
    else:
        mean_length_score = 0.3

    # C. Number of segments (2 is ideal for simple crossover)
    if num_segments == 2:
        segment_count_score = 1.0
    elif num_segments == 3:
        segment_count_score = 0.7
    else:
        segment_count_score = 0.3

    # Weighted average
    segment_quality = (0.4 * min_length_score +
                      0.3 * mean_length_score +
                      0.3 * segment_count_score)

    return segment_quality


def mapping_quality_score(paf_records):
    """
    Score based on MAPQ and alignment scores (AS-NM).

    Returns: Score 0.0-1.0
    """
    if not paf_records or len(paf_records) == 0:
        return 0.0

    mapqs = [rec['mapq'] for rec in paf_records if 'mapq' in rec]

    if not mapqs:
        return 0.5  # Neutral if no MAPQ data

    mean_mapq = np.mean(mapqs)

    # MAPQ score
    if mean_mapq >= 50:
        mapq_score = 1.0
    elif mean_mapq >= 30:
        mapq_score = 0.8
    elif mean_mapq >= 10:
        mapq_score = 0.5
    else:
        mapq_score = 0.2

    # Alignment score quality (AS-NM)
    as_nm_scores = []
    lengths = []
    for rec in paf_records:
        if 'AS' in rec and 'NM' in rec and 'qend' in rec and 'qstart' in rec:
            as_nm_scores.append(rec['AS'] - rec['NM'])
            lengths.append(rec['qend'] - rec['qstart'])

    if as_nm_scores and lengths:
        mean_as_nm = np.mean(as_nm_scores)
        mean_length = np.mean(lengths)
        normalized_as_nm = mean_as_nm / mean_length if mean_length > 0 else 0

        if normalized_as_nm >= 0.9:
            alignment_score = 1.0
        elif normalized_as_nm >= 0.7:
            alignment_score = 0.8
        elif normalized_as_nm >= 0.5:
            alignment_score = 0.6
        else:
            alignment_score = 0.3
    else:
        alignment_score = 0.5  # Neutral if no AS-NM data

    # Weighted average
    mapping_quality = 0.6 * mapq_score + 0.4 * alignment_score

    return mapping_quality


def parent_discrimination_score(paf_records):
    """
    Score based on how clearly segments discriminate between parents.

    Returns: Score 0.0-1.0
    """
    if not paf_records or len(paf_records) == 0:
        return 0.5  # Neutral if no data

    # Group by segment (read_id contains segment info)
    segment_groups = {}
    for rec in paf_records:
        read_id = rec.get('read_id', '')
        # Simple grouping: assume segment info is in read_id
        segment_groups.setdefault(read_id, []).append(rec)

    segment_scores = []

    for segment_id, records in segment_groups.items():
        # Separate Col and Ler alignments
        col_alignments = [r for r in records if 'target' in r and 'Col' in r['target']]
        ler_alignments = [r for r in records if 'target' in r and 'Ler' in r['target']]

        if not col_alignments or not ler_alignments:
            segment_scores.append(0.5)  # Neutral if only one parent
            continue

        # Get best from each parent
        def get_score(rec):
            if 'AS' in rec and 'NM' in rec:
                return rec['AS'] - rec['NM']
            return 0

        best_col = max(col_alignments, key=get_score)
        best_ler = max(ler_alignments, key=get_score)

        col_score = get_score(best_col)
        ler_score = get_score(best_ler)

        # Calculate discrimination ratio
        if col_score > ler_score and col_score > 0:
            discrimination_ratio = (col_score - ler_score) / col_score
        elif ler_score > col_score and ler_score > 0:
            discrimination_ratio = (ler_score - col_score) / ler_score
        else:
            discrimination_ratio = 0.0

        # Score discrimination
        if discrimination_ratio >= 0.20:
            disc_score = 1.0
        elif discrimination_ratio >= 0.10:
            disc_score = 0.7
        elif discrimination_ratio >= 0.05:
            disc_score = 0.4
        else:
            disc_score = 0.1

        segment_scores.append(disc_score)

    return np.mean(segment_scores) if segment_scores else 0.5


def crossover_resolution_score(cowidth_records):
    """
    Score based on crossover width (gap between segments).

    Returns: Score 0.0-1.0
    """
    if not cowidth_records or len(cowidth_records) == 0:
        return 0.5  # Neutral if no data

    widths = [rec['width'] for rec in cowidth_records if 'width' in rec]

    if not widths:
        return 0.5

    min_width = min(widths)

    # Score based on minimum gap
    if min_width <= 100:
        return 1.0
    elif min_width <= 500:
        return 0.8
    elif min_width <= 1000:
        return 0.6
    elif min_width <= 5000:
        return 0.4
    else:
        return 0.2


def load_bed_file(bed_path, read_ids):
    """Load BED file and organize by read ID."""
    bed_data = {}

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
                if read_id in read_ids:
                    bed_data.setdefault(read_id, []).append({
                        'read_id': read_id,
                        'start': int(parts[1]),
                        'end': int(parts[2]),
                        'region': parts[3]
                    })
    except Exception as e:
        print(f"Warning: Error loading BED file {bed_path}: {e}", file=sys.stderr)

    return bed_data


def load_paf_file(paf_path, read_ids):
    """Load PAF file and organize by read ID (base read, not segment)."""
    paf_data = {}

    if not os.path.exists(paf_path):
        return paf_data

    try:
        with open(paf_path, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 12:
                    continue

                read_id = parts[0]
                # Extract base read ID (remove segment info)
                base_read_id = '_'.join(read_id.split('_')[:-2]) if '_' in read_id else read_id

                if base_read_id in read_ids:
                    paf_record = {
                        'read_id': read_id,
                        'qlen': int(parts[1]),
                        'qstart': int(parts[2]),
                        'qend': int(parts[3]),
                        'strand': parts[4],
                        'target': parts[5],
                        'tlen': int(parts[6]),
                        'tstart': int(parts[7]),
                        'tend': int(parts[8]),
                        'mapq': int(parts[11])
                    }

                    # Parse optional tags
                    for tag in parts[12:]:
                        if tag.startswith('AS:i:'):
                            paf_record['AS'] = int(tag[5:])
                        elif tag.startswith('NM:i:'):
                            paf_record['NM'] = int(tag[5:])

                    paf_data.setdefault(base_read_id, []).append(paf_record)
    except Exception as e:
        print(f"Warning: Error loading PAF file {paf_path}: {e}", file=sys.stderr)

    return paf_data


def load_cowidth_file(cowidth_path, read_ids):
    """Load crossover width file and organize by read ID."""
    cowidth_data = {}

    if not os.path.exists(cowidth_path):
        return cowidth_data

    try:
        df = pd.read_csv(cowidth_path, sep='\t')
        for read_id in read_ids:
            if read_id in df['read'].values:
                records = df[df['read'] == read_id].to_dict('records')
                cowidth_data[read_id] = records
    except Exception as e:
        print(f"Warning: Error loading cowidth file {cowidth_path}: {e}", file=sys.stderr)

    return cowidth_data


def calculate_confidence_score(read_id, summary_row, bed_data, paf_data, cowidth_data):
    """
    Calculate overall confidence score for a crossover call.

    Returns: dict with confidence score and component scores
    """
    # Component scores
    mapping_consensus = mapping_consensus_score(summary_row)

    segment_quality = 0.5
    if read_id in bed_data:
        segment_quality = segment_quality_score(bed_data[read_id])

    mapping_quality = 0.5
    if read_id in paf_data:
        mapping_quality = mapping_quality_score(paf_data[read_id])

    parent_discrimination = 0.5
    if read_id in paf_data:
        parent_discrimination = parent_discrimination_score(paf_data[read_id])

    crossover_resolution = 0.5
    if read_id in cowidth_data:
        crossover_resolution = crossover_resolution_score(cowidth_data[read_id])

    # Weights
    weights = {
        'mapping_consensus': 0.30,
        'segment_quality': 0.25,
        'mapping_quality': 0.20,
        'parent_discrimination': 0.15,
        'crossover_resolution': 0.10
    }

    # Weighted sum
    confidence = (
        weights['mapping_consensus'] * mapping_consensus +
        weights['segment_quality'] * segment_quality +
        weights['mapping_quality'] * mapping_quality +
        weights['parent_discrimination'] * parent_discrimination +
        weights['crossover_resolution'] * crossover_resolution
    )

    return {
        'read_id': read_id,
        'confidence_score': confidence,
        'mapping_consensus': mapping_consensus,
        'segment_quality': segment_quality,
        'mapping_quality': mapping_quality,
        'parent_discrimination': parent_discrimination,
        'crossover_resolution': crossover_resolution
    }


def main():
    parser = argparse.ArgumentParser(
        description='Calculate confidence scores for detected crossovers',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
    python calculate_confidence_scores.py \\
        --summary-arms summary_table_ARMS.csv \\
        --summary-cen summary_table_CEN.csv \\
        --bed-arms CANDIDATE_reads_nonectopic_ARMS.bed \\
        --bed-cen CANDIDATE_reads_nonectopic_CEN.bed \\
        --paf-dir mapping_results/ \\
        --cowidth-arms CANDIDATE_reads_nonectopic_ARMS.bed_of_best.cowidths \\
        --cowidth-cen CANDIDATE_reads_nonectopic_CEN.bed_of_best.cowidths \\
        --output confidence_scores.tsv
        """
    )

    parser.add_argument('--summary-arms', help='Summary table CSV for ARMS reads')
    parser.add_argument('--summary-cen', help='Summary table CSV for CEN reads')
    parser.add_argument('--bed-arms', help='BED file for ARMS reads')
    parser.add_argument('--bed-cen', help='BED file for CEN reads')
    parser.add_argument('--paf-dir', help='Directory containing PAF files')
    parser.add_argument('--cowidth-arms', help='Cowidth file for ARMS reads')
    parser.add_argument('--cowidth-cen', help='Cowidth file for CEN reads')
    parser.add_argument('--output', required=True, help='Output TSV file')

    args = parser.parse_args()

    results = []

    # Process ARMS reads
    if args.summary_arms and os.path.exists(args.summary_arms):
        print(f"Processing ARMS reads from {args.summary_arms}...", file=sys.stderr)
        summary = pd.read_csv(args.summary_arms)
        read_ids = set(summary['read_id'].values)

        # Load supporting data
        bed_data = load_bed_file(args.bed_arms, read_ids) if args.bed_arms else {}

        paf_data = {}
        if args.paf_dir and os.path.isdir(args.paf_dir):
            # Load all PAF files in directory
            for paf_file in Path(args.paf_dir).glob('*.paf'):
                paf_data.update(load_paf_file(paf_file, read_ids))

        cowidth_data = load_cowidth_file(args.cowidth_arms, read_ids) if args.cowidth_arms else {}

        # Calculate scores
        for _, row in summary.iterrows():
            read_id = row['read_id']
            score_data = calculate_confidence_score(read_id, row, bed_data, paf_data, cowidth_data)
            score_data['region_type'] = 'ARMS'
            results.append(score_data)

        print(f"  Processed {len(results)} ARMS reads", file=sys.stderr)

    # Process CEN reads
    if args.summary_cen and os.path.exists(args.summary_cen):
        print(f"Processing CEN reads from {args.summary_cen}...", file=sys.stderr)
        summary = pd.read_csv(args.summary_cen)
        read_ids = set(summary['read_id'].values)

        # Load supporting data
        bed_data = load_bed_file(args.bed_cen, read_ids) if args.bed_cen else {}

        paf_data = {}
        if args.paf_dir and os.path.isdir(args.paf_dir):
            for paf_file in Path(args.paf_dir).glob('*.paf'):
                paf_data.update(load_paf_file(paf_file, read_ids))

        cowidth_data = load_cowidth_file(args.cowidth_cen, read_ids) if args.cowidth_cen else {}

        # Calculate scores
        for _, row in summary.iterrows():
            read_id = row['read_id']
            score_data = calculate_confidence_score(read_id, row, bed_data, paf_data, cowidth_data)
            score_data['region_type'] = 'CEN'
            results.append(score_data)

        print(f"  Processed {len(results) - len([r for r in results if r['region_type'] == 'ARMS'])} CEN reads", file=sys.stderr)

    if not results:
        print("Warning: No reads processed. Check input files.", file=sys.stderr)
        # Create empty output file
        with open(args.output, 'w') as f:
            f.write("read_id\tconfidence_score\tmapping_consensus\tsegment_quality\tmapping_quality\tparent_discrimination\tcrossover_resolution\tregion_type\n")
        return

    # Create DataFrame and save
    df = pd.DataFrame(results)

    # Reorder columns
    cols = ['read_id', 'confidence_score', 'mapping_consensus', 'segment_quality',
            'mapping_quality', 'parent_discrimination', 'crossover_resolution', 'region_type']
    df = df[cols]

    # Sort by confidence score (descending)
    df = df.sort_values('confidence_score', ascending=False)

    df.to_csv(args.output, sep='\t', index=False, float_format='%.4f')

    # Print summary statistics
    print("\n" + "="*60, file=sys.stderr)
    print("CONFIDENCE SCORING SUMMARY", file=sys.stderr)
    print("="*60, file=sys.stderr)
    print(f"Total reads scored: {len(df)}", file=sys.stderr)
    print(f"  - ARMS: {(df['region_type'] == 'ARMS').sum()}", file=sys.stderr)
    print(f"  - CEN: {(df['region_type'] == 'CEN').sum()}", file=sys.stderr)
    print(f"\nConfidence score distribution:", file=sys.stderr)
    print(f"  High confidence (≥0.8): {(df['confidence_score'] >= 0.8).sum()} ({(df['confidence_score'] >= 0.8).sum() / len(df) * 100:.1f}%)", file=sys.stderr)
    print(f"  Medium confidence (0.5-0.8): {((df['confidence_score'] >= 0.5) & (df['confidence_score'] < 0.8)).sum()} ({((df['confidence_score'] >= 0.5) & (df['confidence_score'] < 0.8)).sum() / len(df) * 100:.1f}%)", file=sys.stderr)
    print(f"  Low confidence (<0.5): {(df['confidence_score'] < 0.5).sum()} ({(df['confidence_score'] < 0.5).sum() / len(df) * 100:.1f}%)", file=sys.stderr)
    print(f"\nScore statistics:", file=sys.stderr)
    print(f"  Mean: {df['confidence_score'].mean():.4f}", file=sys.stderr)
    print(f"  Median: {df['confidence_score'].median():.4f}", file=sys.stderr)
    print(f"  Std: {df['confidence_score'].std():.4f}", file=sys.stderr)
    print(f"  Min: {df['confidence_score'].min():.4f}", file=sys.stderr)
    print(f"  Max: {df['confidence_score'].max():.4f}", file=sys.stderr)
    print(f"\nComponent score means:", file=sys.stderr)
    print(f"  Mapping consensus: {df['mapping_consensus'].mean():.4f}", file=sys.stderr)
    print(f"  Segment quality: {df['segment_quality'].mean():.4f}", file=sys.stderr)
    print(f"  Mapping quality: {df['mapping_quality'].mean():.4f}", file=sys.stderr)
    print(f"  Parent discrimination: {df['parent_discrimination'].mean():.4f}", file=sys.stderr)
    print(f"  Crossover resolution: {df['crossover_resolution'].mean():.4f}", file=sys.stderr)
    print("="*60, file=sys.stderr)
    print(f"\nOutput written to: {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()
