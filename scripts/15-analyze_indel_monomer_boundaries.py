#!/usr/bin/env python3
"""
Analyze centromeric indel boundaries relative to the 178bp satellite monomer.

This script addresses the key biological question: Do indels represent whole-monomer
deletions (in-frame) or random breakpoints?

Key analyses:
1. Align indels to doubled monomer reference (handles circularity)
2. Create metaprofile of breakpoint positions across the monomer
3. Test if breakpoints cluster at specific positions (biological signal)
4. Visualize genome-wide/centromere-wide patterns
5. Dimensional reduction to find "active" monomer variants

Author: Claude Code
"""

import os
import sys
import subprocess
import tempfile
from pathlib import Path
import pandas as pd
import numpy as np
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import substitution_matrices
from collections import Counter, defaultdict
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist

# Try importing sklearn for dimensionality reduction
try:
    from sklearn.decomposition import PCA
    from sklearn.manifold import TSNE
    from sklearn.cluster import KMeans
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    print("[WARNING] sklearn not available - dimensionality reduction will be skipped")

# Try importing umap
try:
    import umap
    HAS_UMAP = True
except ImportError:
    HAS_UMAP = False


# Arabidopsis CEN178 consensus monomer
CEN178_CONSENSUS = (
    "AGTATAAGAACTTAAACCGCAACCGATCTTA"
    "AAAGCCTAAGTAGTGTTTCCTTGTTAGAAGA"
    "CACAAAGCCAAAGACTCATATGGACTTTGGC"
    "TACACCATGAAAGCTTTGAGAAGCAAGAAGA"
    "AGGTTGGTTAGTGTTTTGGAGTCGAATATGA"
    "CTTGATGTCATGTGTATGATTG"
)

#CEN178_CONSENSUS_alternative = (
#    "TACACCATGAAAGCTTTGAGAAGCAAGAAGA"
#    "AGGTTGGTTAGTGTTTTGGAGTCGAATATGA"
#    "CTTGATGTCATGTGTATGATTG"
#   "AGTATAAGAACTTAAACCGCAACCGATCTTA"
#    "AAAGCCTAAGTAGTGTTTCCTTGTTAGAAGA"
#    "CACAAAGCCAAAGACTCATATGGACTTTGGC"
#    "TACACCATGAAAGCTTTGAGAAGCAAGAAGA"
#    "AGGTTGGTTAGTGTTTTGGAGTCGAATATGA"
#    "CTTGATGTCATGTGTATGATTG"
#)


MONOMER_LENGTH = 178


def create_doubled_monomer():
    """Create doubled monomer for circular alignment."""
    return CEN178_CONSENSUS + CEN178_CONSENSUS


def get_reverse_complement(seq):
    """Get reverse complement of a sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq.upper()))


def align_to_monomer_pairwise(sequence, monomer_ref):
    """
    Align a sequence to the doubled monomer using local alignment.
    Returns alignment details including start/end positions in monomer coordinates.
    """
    # Use Smith-Waterman local alignment
    # Parameters tuned for satellite DNA (allow gaps for indels)
    alignments = pairwise2.align.localms(
        monomer_ref,      # Reference (doubled monomer)
        sequence.upper(), # Query (indel sequence)
        2,                # Match score
        -1,               # Mismatch penalty
        -2,               # Gap open penalty
        -0.5,             # Gap extend penalty
        one_alignment_only=True
    )

    if not alignments:
        return None

    best = alignments[0]
    ref_aligned, query_aligned, score, ref_start, ref_end = best

    # Find actual alignment boundaries (excluding leading/trailing gaps)
    query_start = 0
    query_end = len(sequence)

    # Calculate coverage
    aligned_bases = sum(1 for i, c in enumerate(query_aligned) if c != '-')
    coverage = aligned_bases / len(sequence) if len(sequence) > 0 else 0

    # Calculate identity
    matches = sum(1 for r, q in zip(ref_aligned, query_aligned) if r == q and r != '-')
    identity = matches / aligned_bases if aligned_bases > 0 else 0

    # Normalize position to single monomer (0-177)
    monomer_start = ref_start % MONOMER_LENGTH
    monomer_end = (ref_start + len(sequence)) % MONOMER_LENGTH

    return {
        'ref_start': ref_start,
        'ref_end': ref_start + len(sequence),
        'monomer_start': monomer_start,
        'monomer_end': monomer_end,
        'score': score,
        'identity': identity,
        'coverage': coverage,
        'strand': '+'
    }


def align_with_minimap2(indel_sequences_file, output_dir):
    """
    Align indel sequences to doubled monomer using minimap2.
    Returns alignment results as a dictionary.
    """
    doubled_monomer = create_doubled_monomer()

    # Create temporary reference file
    ref_file = os.path.join(output_dir, "doubled_monomer_ref.fa")
    with open(ref_file, 'w') as f:
        f.write(">doubled_CEN178\n")
        f.write(doubled_monomer + "\n")
        # Also add reverse complement
        f.write(">doubled_CEN178_rc\n")
        f.write(get_reverse_complement(doubled_monomer) + "\n")

    # Run minimap2
    paf_file = os.path.join(output_dir, "indels_vs_monomer.paf")
    cmd = [
        'minimap2',
        '-c',           # Generate CIGAR
        '-x', 'map-ont', # ONT preset
        '-N', '1',      # Only output best alignment
        '--secondary=no',
        ref_file,
        indel_sequences_file
    ]

    print(f"[INFO] Running minimap2 alignment...")
    try:
        with open(paf_file, 'w') as outfile:
            result = subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True)
        print(f"[SUCCESS] Alignment complete: {paf_file}")
        return parse_paf_file(paf_file)
    except FileNotFoundError:
        print("[WARNING] minimap2 not found, using pairwise alignment fallback")
        return None


def parse_paf_file(paf_file):
    """Parse PAF alignment file and extract boundary information."""
    alignments = {}

    with open(paf_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue

            query_name = fields[0]
            query_len = int(fields[1])
            query_start = int(fields[2])
            query_end = int(fields[3])
            strand = fields[4]
            target_name = fields[5]
            target_len = int(fields[6])
            target_start = int(fields[7])
            target_end = int(fields[8])

            # Calculate alignment quality
            matches = int(fields[9])
            block_len = int(fields[10])
            mapq = int(fields[11])

            identity = matches / block_len if block_len > 0 else 0
            coverage = (query_end - query_start) / query_len if query_len > 0 else 0

            # Normalize to single monomer coordinates (0-177)
            if 'rc' in target_name:
                # Reverse complement alignment
                monomer_start = (MONOMER_LENGTH * 2 - target_end) % MONOMER_LENGTH
                monomer_end = (MONOMER_LENGTH * 2 - target_start) % MONOMER_LENGTH
                actual_strand = '-'
            else:
                monomer_start = target_start % MONOMER_LENGTH
                monomer_end = target_end % MONOMER_LENGTH
                actual_strand = strand

            alignments[query_name] = {
                'ref_start': target_start,
                'ref_end': target_end,
                'monomer_start': monomer_start,
                'monomer_end': monomer_end,
                'identity': identity,
                'coverage': coverage,
                'mapq': mapq,
                'strand': actual_strand,
                'query_len': query_len
            }

    return alignments


def analyze_boundary_distribution(alignments, output_dir):
    """
    Analyze the distribution of indel boundaries across the monomer.
    Tests if breakpoints cluster at specific positions (suggesting biological mechanism).
    """
    start_positions = []
    end_positions = []

    for indel_id, aln in alignments.items():
        if aln['coverage'] >= 0.5 and aln['identity'] >= 0.6:  # Quality filter
            start_positions.append(aln['monomer_start'])
            end_positions.append(aln['monomer_end'])

    if not start_positions:
        print("[WARNING] No quality alignments found for boundary analysis")
        return None

    # Create metaprofile
    start_counts = Counter(start_positions)
    end_counts = Counter(end_positions)
    all_boundaries = start_positions + end_positions
    boundary_counts = Counter(all_boundaries)

    # Statistical test: Are boundaries uniformly distributed?
    # Use chi-square test against uniform distribution
    observed = [boundary_counts.get(i, 0) for i in range(MONOMER_LENGTH)]
    expected = [len(all_boundaries) / MONOMER_LENGTH] * MONOMER_LENGTH

    # Only test if we have enough data
    if len(all_boundaries) >= MONOMER_LENGTH:
        chi2, p_value = stats.chisquare(observed, expected)
        uniform_test = {'chi2': chi2, 'p_value': p_value}
    else:
        uniform_test = {'chi2': None, 'p_value': None}

    # Find hotspots (positions with significantly more boundaries)
    mean_count = np.mean(observed)
    std_count = np.std(observed)
    threshold = mean_count + 2 * std_count  # 2 SD above mean

    hotspots = [(pos, count) for pos, count in enumerate(observed) if count > threshold]
    hotspots.sort(key=lambda x: x[1], reverse=True)

    return {
        'start_positions': start_positions,
        'end_positions': end_positions,
        'boundary_counts': boundary_counts,
        'uniform_test': uniform_test,
        'hotspots': hotspots[:10],  # Top 10 hotspots
        'total_boundaries': len(all_boundaries)
    }


def create_metaprofile_plot(boundary_analysis, output_dir):
    """
    Create metaprofile visualization showing where boundaries occur in the monomer.
    """
    fig, axes = plt.subplots(3, 1, figsize=(16, 14))

    # Full metaprofile
    ax = axes[0]
    positions = list(range(MONOMER_LENGTH))
    counts = [boundary_analysis['boundary_counts'].get(p, 0) for p in positions]

    # Smooth the profile for visualization
    window = 5
    smoothed = np.convolve(counts + counts[:window], np.ones(window)/window, mode='valid')[:MONOMER_LENGTH]

    ax.bar(positions, counts, alpha=0.5, color='steelblue', label='Raw counts')
    ax.plot(positions, smoothed, color='red', linewidth=2, label=f'Smoothed (window={window})')

    # Mark hotspots
    for pos, count in boundary_analysis['hotspots'][:5]:
        ax.axvline(pos, color='green', linestyle='--', alpha=0.7, linewidth=1)
        ax.text(pos, max(counts) * 1.05, f'{pos}', rotation=90,
                ha='center', va='bottom', fontsize=8, color='green')

    ax.set_xlabel('Position in Monomer (bp)', fontweight='bold', fontsize=12)
    ax.set_ylabel('Number of Boundaries', fontweight='bold', fontsize=12)
    ax.set_title('Metaprofile: Indel Boundaries Across the 178bp Monomer',
                 fontweight='bold', fontsize=14)
    ax.legend()
    ax.set_xlim(0, MONOMER_LENGTH)

    # Add statistical annotation
    if boundary_analysis['uniform_test']['p_value'] is not None:
        p_val = boundary_analysis['uniform_test']['p_value']
        significance = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else "ns"
        ax.text(0.02, 0.95, f"Chi-square test for uniformity:\np = {p_val:.2e} ({significance})",
               transform=ax.transAxes, fontsize=10, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Start vs End positions
    ax = axes[1]
    start_counts = [boundary_analysis['start_positions'].count(p) for p in positions]
    end_counts = [boundary_analysis['end_positions'].count(p) for p in positions]

    ax.bar(positions, start_counts, alpha=0.6, color='blue', label='Start boundaries')
    ax.bar(positions, end_counts, alpha=0.6, color='orange', label='End boundaries', bottom=start_counts)

    ax.set_xlabel('Position in Monomer (bp)', fontweight='bold', fontsize=12)
    ax.set_ylabel('Number of Boundaries', fontweight='bold', fontsize=12)
    ax.set_title('Start vs End Boundary Positions', fontweight='bold', fontsize=14)
    ax.legend()
    ax.set_xlim(0, MONOMER_LENGTH)

    # Circular representation
    ax = axes[2]

    # Create polar plot for circular visualization
    theta = np.linspace(0, 2 * np.pi, MONOMER_LENGTH, endpoint=False)
    radii = np.array(counts)

    # Normalize for visualization
    radii_norm = radii / max(radii) if max(radii) > 0 else radii

    # Plot as bar chart wrapped around
    width = 2 * np.pi / MONOMER_LENGTH
    bars = ax.bar(theta, radii_norm, width=width, alpha=0.7, color='steelblue')

    # Color hotspots
    for pos, _ in boundary_analysis['hotspots'][:5]:
        if pos < len(bars):
            bars[pos].set_color('red')

    ax.set_xlabel('Position in Monomer (polar representation)', fontweight='bold', fontsize=12)
    ax.set_ylabel('Normalized Boundary Count', fontweight='bold', fontsize=12)
    ax.set_title('Circular Metaprofile of Indel Boundaries', fontweight='bold', fontsize=14)

    # Add monomer sequence annotation
    ax.text(0.5, -0.15, f'Total boundaries: {boundary_analysis["total_boundaries"]} | '
            f'Top hotspot: position {boundary_analysis["hotspots"][0][0] if boundary_analysis["hotspots"] else "N/A"}',
            transform=ax.transAxes, ha='center', fontsize=10)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'indel_boundary_metaprofile.png'), dpi=300, bbox_inches='tight')
    plt.close()

    print(f"[SUCCESS] Saved metaprofile plot")


def analyze_inframe_alignment(catalog_df, alignments, output_dir):
    """
    Analyze whether indels represent in-frame (whole monomer) deletions.

    An indel is considered "in-frame" if:
    1. Its size is a multiple of 178bp
    2. Its boundaries align with monomer boundaries (start ≈ end position)
    """
    results = []

    for _, row in catalog_df.iterrows():
        indel_id = row['indel_id']
        indel_size = row['size']

        # Check if size is multiple of 178
        is_size_multiple = (indel_size % MONOMER_LENGTH == 0)
        closest_multiple = round(indel_size / MONOMER_LENGTH)
        size_deviation = abs(indel_size - closest_multiple * MONOMER_LENGTH)

        # Check alignment boundaries
        if indel_id in alignments:
            aln = alignments[indel_id]
            boundary_distance = abs(aln['monomer_end'] - aln['monomer_start'])
            # Account for wrapping around
            boundary_distance = min(boundary_distance, MONOMER_LENGTH - boundary_distance)

            # In-frame if boundaries are close (within 10bp tolerance)
            is_boundary_aligned = boundary_distance <= 10

            results.append({
                'indel_id': indel_id,
                'size': indel_size,
                'is_size_multiple': is_size_multiple,
                'closest_multiple': closest_multiple,
                'size_deviation': size_deviation,
                'monomer_start': aln['monomer_start'],
                'monomer_end': aln['monomer_end'],
                'boundary_distance': boundary_distance,
                'is_boundary_aligned': is_boundary_aligned,
                'is_inframe': is_size_multiple and is_boundary_aligned,
                'alignment_identity': aln['identity'],
                'alignment_coverage': aln['coverage'],
                'strand': aln['strand']
            })
        else:
            results.append({
                'indel_id': indel_id,
                'size': indel_size,
                'is_size_multiple': is_size_multiple,
                'closest_multiple': closest_multiple,
                'size_deviation': size_deviation,
                'monomer_start': None,
                'monomer_end': None,
                'boundary_distance': None,
                'is_boundary_aligned': None,
                'is_inframe': None,
                'alignment_identity': None,
                'alignment_coverage': None,
                'strand': None
            })

    return pd.DataFrame(results)


def create_inframe_analysis_plot(inframe_df, output_dir):
    """
    Visualize in-frame analysis results.
    """
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))

    # Filter to aligned indels
    aligned_df = inframe_df.dropna(subset=['is_inframe'])

    if len(aligned_df) == 0:
        print("[WARNING] No aligned indels for plotting")
        plt.close()
        return

    # Plot 1: Size deviation from 178bp multiples
    ax = axes[0, 0]
    deviations = aligned_df['size_deviation'].values
    ax.hist(deviations, bins=30, alpha=0.7, color='steelblue', edgecolor='black')
    ax.axvline(0, color='red', linestyle='--', linewidth=2, label='Exact multiple')
    ax.set_xlabel('Deviation from 178bp Multiple (bp)', fontweight='bold')
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_title('Size Deviation from Monomer Multiples', fontweight='bold')
    ax.legend()

    # Plot 2: Boundary distance distribution
    ax = axes[0, 1]
    boundary_dist = aligned_df['boundary_distance'].dropna().values
    ax.hist(boundary_dist, bins=30, alpha=0.7, color='purple', edgecolor='black')
    ax.axvline(10, color='red', linestyle='--', linewidth=2, label='Alignment threshold (10bp)')
    ax.set_xlabel('Start-End Boundary Distance (bp)', fontweight='bold')
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_title('Indel Boundary Alignment', fontweight='bold')
    ax.legend()

    # Plot 3: In-frame classification
    ax = axes[0, 2]
    inframe_counts = aligned_df['is_inframe'].value_counts()
    colors = ['#2ecc71' if x else '#e74c3c' for x in [True, False]]
    labels = ['In-frame\n(biological)', 'Out-of-frame\n(random)']
    values = [inframe_counts.get(True, 0), inframe_counts.get(False, 0)]
    ax.bar(labels, values, color=colors, alpha=0.7, edgecolor='black')
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_title('In-Frame vs Out-of-Frame Indels', fontweight='bold')

    # Add percentages
    total = sum(values)
    for i, v in enumerate(values):
        pct = v / total * 100 if total > 0 else 0
        ax.text(i, v + 1, f'{pct:.1f}%', ha='center', fontweight='bold')

    # Plot 4: Size vs boundary alignment scatter
    ax = axes[1, 0]
    scatter = ax.scatter(
        aligned_df['size_deviation'],
        aligned_df['boundary_distance'],
        c=aligned_df['is_inframe'].map({True: 'green', False: 'red'}),
        alpha=0.6,
        s=50
    )
    ax.axhline(10, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('Size Deviation from Multiple (bp)', fontweight='bold')
    ax.set_ylabel('Boundary Distance (bp)', fontweight='bold')
    ax.set_title('Size Deviation vs Boundary Alignment', fontweight='bold')

    # Add quadrant labels
    ax.text(0.02, 0.98, 'In-frame\n(green)', transform=ax.transAxes,
            va='top', fontsize=10, color='green')

    # Plot 5: Start vs End positions
    ax = axes[1, 1]
    ax.scatter(
        aligned_df['monomer_start'],
        aligned_df['monomer_end'],
        c=aligned_df['is_inframe'].map({True: 'green', False: 'red', None: 'gray'}),
        alpha=0.6,
        s=50
    )
    # Add diagonal (in-frame would be on diagonal)
    ax.plot([0, MONOMER_LENGTH], [0, MONOMER_LENGTH], 'k--', alpha=0.3, label='Perfect alignment')
    ax.set_xlabel('Monomer Start Position', fontweight='bold')
    ax.set_ylabel('Monomer End Position', fontweight='bold')
    ax.set_title('Start vs End Positions in Monomer', fontweight='bold')
    ax.set_xlim(0, MONOMER_LENGTH)
    ax.set_ylim(0, MONOMER_LENGTH)
    ax.legend()

    # Plot 6: Number of monomers per indel
    ax = axes[1, 2]
    monomer_counts = aligned_df['closest_multiple'].value_counts().sort_index()
    bars = ax.bar(monomer_counts.index, monomer_counts.values, alpha=0.7, color='steelblue', edgecolor='black')

    # Color in-frame vs out-of-frame
    for mult in monomer_counts.index:
        subset = aligned_df[aligned_df['closest_multiple'] == mult]
        inframe_pct = subset['is_inframe'].sum() / len(subset) if len(subset) > 0 else 0
        if inframe_pct > 0.5:
            bars[list(monomer_counts.index).index(mult)].set_color('#2ecc71')

    ax.set_xlabel('Number of Monomers', fontweight='bold')
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_title('Indel Size by Monomer Units', fontweight='bold')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'inframe_analysis.png'), dpi=300, bbox_inches='tight')
    plt.close()

    print(f"[SUCCESS] Saved in-frame analysis plot")


def extract_monomer_variants(indel_sequences_file, alignments, output_dir):
    """
    Extract individual monomer units from indels and analyze their variation.
    Uses alignments to properly phase the extraction.
    """
    monomer_variants = []

    for record in SeqIO.parse(indel_sequences_file, 'fasta'):
        indel_id = record.id
        seq = str(record.seq).upper()

        if indel_id not in alignments:
            continue

        aln = alignments[indel_id]
        if aln['coverage'] < 0.5 or aln['identity'] < 0.6:
            continue

        # Use alignment to phase the sequence
        start_offset = aln['monomer_start']

        # Adjust sequence to start at monomer boundary
        if start_offset > 0:
            # Trim to nearest boundary
            seq = seq[MONOMER_LENGTH - start_offset:]

        # Extract complete monomers
        num_monomers = len(seq) // MONOMER_LENGTH
        for i in range(num_monomers):
            monomer_seq = seq[i * MONOMER_LENGTH:(i + 1) * MONOMER_LENGTH]
            if len(monomer_seq) == MONOMER_LENGTH:
                monomer_variants.append({
                    'source_indel': indel_id,
                    'monomer_index': i,
                    'sequence': monomer_seq,
                    'strand': aln['strand']
                })

    return monomer_variants


def create_sequence_matrix(sequences, reference=None):
    """
    Create a numerical matrix from sequences for dimensionality reduction.
    Each sequence is encoded as presence/absence of k-mers or base composition.
    """
    if not sequences:
        return None

    # Simple encoding: base composition at each position
    # This works well for aligned sequences of same length
    max_len = max(len(s) for s in sequences)

    matrix = []
    for seq in sequences:
        # Pad sequence if needed
        padded = seq.ljust(max_len, 'N')
        # Encode: A=0, C=1, G=2, T=3, N=4
        encoding = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
        row = [encoding.get(base, 4) for base in padded]
        matrix.append(row)

    return np.array(matrix)


def run_dimensionality_reduction(monomer_variants, output_dir):
    """
    Perform dimensionality reduction on monomer variants to find clusters.
    """
    if not HAS_SKLEARN:
        print("[WARNING] sklearn not available, skipping dimensionality reduction")
        return None

    if len(monomer_variants) < 10:
        print(f"[WARNING] Not enough monomer variants for clustering ({len(monomer_variants)})")
        return None

    sequences = [v['sequence'] for v in monomer_variants]

    # Create sequence matrix
    seq_matrix = create_sequence_matrix(sequences)
    if seq_matrix is None:
        return None

    print(f"[INFO] Running dimensionality reduction on {len(sequences)} monomer variants...")

    # PCA
    pca = PCA(n_components=min(10, len(sequences) - 1))
    pca_coords = pca.fit_transform(seq_matrix)

    # K-means clustering
    n_clusters = min(5, len(sequences) // 3)
    if n_clusters >= 2:
        kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
        clusters = kmeans.fit_predict(seq_matrix)
    else:
        clusters = np.zeros(len(sequences))

    # UMAP if available
    if HAS_UMAP and len(sequences) >= 15:
        reducer = umap.UMAP(n_components=2, random_state=42)
        umap_coords = reducer.fit_transform(seq_matrix)
    else:
        umap_coords = pca_coords[:, :2]

    return {
        'pca_coords': pca_coords,
        'umap_coords': umap_coords,
        'clusters': clusters,
        'pca_variance': pca.explained_variance_ratio_
    }


def create_dimensionality_reduction_plot(monomer_variants, dim_red_results, output_dir):
    """
    Visualize monomer variant clustering.
    """
    if dim_red_results is None:
        return

    fig, axes = plt.subplots(2, 2, figsize=(14, 14))

    # Plot 1: PCA
    ax = axes[0, 0]
    scatter = ax.scatter(
        dim_red_results['pca_coords'][:, 0],
        dim_red_results['pca_coords'][:, 1],
        c=dim_red_results['clusters'],
        cmap='tab10',
        alpha=0.6,
        s=50
    )
    ax.set_xlabel(f'PC1 ({dim_red_results["pca_variance"][0]*100:.1f}%)', fontweight='bold')
    ax.set_ylabel(f'PC2 ({dim_red_results["pca_variance"][1]*100:.1f}%)', fontweight='bold')
    ax.set_title('PCA of Monomer Variants', fontweight='bold')
    plt.colorbar(scatter, ax=ax, label='Cluster')

    # Plot 2: UMAP/t-SNE
    ax = axes[0, 1]
    scatter = ax.scatter(
        dim_red_results['umap_coords'][:, 0],
        dim_red_results['umap_coords'][:, 1],
        c=dim_red_results['clusters'],
        cmap='tab10',
        alpha=0.6,
        s=50
    )
    ax.set_xlabel('UMAP 1', fontweight='bold')
    ax.set_ylabel('UMAP 2', fontweight='bold')
    ax.set_title('UMAP of Monomer Variants', fontweight='bold')
    plt.colorbar(scatter, ax=ax, label='Cluster')

    # Plot 3: Cluster sizes
    ax = axes[1, 0]
    cluster_counts = Counter(dim_red_results['clusters'])
    clusters_sorted = sorted(cluster_counts.keys())
    counts = [cluster_counts[c] for c in clusters_sorted]
    bars = ax.bar(clusters_sorted, counts, alpha=0.7, color='steelblue', edgecolor='black')
    ax.set_xlabel('Cluster', fontweight='bold')
    ax.set_ylabel('Number of Monomers', fontweight='bold')
    ax.set_title('Cluster Sizes', fontweight='bold')

    # Plot 4: PCA variance explained
    ax = axes[1, 1]
    n_components = len(dim_red_results['pca_variance'])
    ax.bar(range(1, n_components + 1), dim_red_results['pca_variance'] * 100,
           alpha=0.7, color='purple', edgecolor='black')
    ax.set_xlabel('Principal Component', fontweight='bold')
    ax.set_ylabel('Variance Explained (%)', fontweight='bold')
    ax.set_title('PCA Scree Plot', fontweight='bold')

    # Add cumulative line
    cumsum = np.cumsum(dim_red_results['pca_variance']) * 100
    ax2 = ax.twinx()
    ax2.plot(range(1, n_components + 1), cumsum, 'r-o', linewidth=2, markersize=4)
    ax2.set_ylabel('Cumulative Variance (%)', fontweight='bold', color='red')
    ax2.tick_params(axis='y', labelcolor='red')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'monomer_variant_clustering.png'), dpi=300, bbox_inches='tight')
    plt.close()

    print(f"[SUCCESS] Saved monomer variant clustering plot")


def create_alignment_heatmap(indel_sequences_file, alignments, output_dir):
    """
    Create a heatmap showing how each indel aligns across the monomer positions.
    Each row is an indel, each column is a position in the monomer (0-177).
    This visualization clearly shows if boundaries cluster at specific positions.
    """
    # Build alignment matrix
    indel_ids = []
    indel_sizes = []
    alignment_matrix = []

    for record in SeqIO.parse(indel_sequences_file, 'fasta'):
        indel_id = record.id
        if indel_id not in alignments:
            continue

        aln = alignments[indel_id]
        if aln['coverage'] < 0.3:  # Skip poor alignments
            continue

        # Create row: mark aligned positions
        row = np.zeros(MONOMER_LENGTH)
        start = aln['monomer_start']
        end = aln['monomer_end']

        # Handle wrap-around
        if start <= end:
            row[start:end] = 1
        else:
            row[start:] = 1
            row[:end] = 1

        alignment_matrix.append(row)
        indel_ids.append(indel_id)
        indel_sizes.append(len(record.seq))

    if len(alignment_matrix) < 5:
        print("[WARNING] Not enough alignments for heatmap visualization")
        return

    alignment_matrix = np.array(alignment_matrix)

    # Sort by start position for better visualization
    start_positions = [alignments[id]['monomer_start'] for id in indel_ids]
    sort_idx = np.argsort(start_positions)
    alignment_matrix = alignment_matrix[sort_idx]
    indel_ids = [indel_ids[i] for i in sort_idx]
    indel_sizes = [indel_sizes[i] for i in sort_idx]

    # Create figure
    fig = plt.figure(figsize=(16, max(8, len(indel_ids) * 0.1)))
    gs = fig.add_gridspec(1, 2, width_ratios=[20, 1], wspace=0.02)

    # Main heatmap
    ax_main = fig.add_subplot(gs[0])

    # Custom colormap: white for no alignment, blue for aligned
    cmap = plt.cm.Blues
    im = ax_main.imshow(alignment_matrix, aspect='auto', cmap=cmap, interpolation='nearest')

    ax_main.set_xlabel('Position in Monomer (bp)', fontweight='bold', fontsize=12)
    ax_main.set_ylabel('Indels (sorted by start position)', fontweight='bold', fontsize=12)
    ax_main.set_title('Indel Alignment Across the 178bp Monomer\n(Each row = one indel, colored region = aligned)',
                     fontweight='bold', fontsize=14)

    # Add vertical lines at notable positions (every 30bp)
    for pos in range(0, MONOMER_LENGTH, 30):
        ax_main.axvline(pos, color='red', linestyle=':', alpha=0.3, linewidth=0.5)

    # Size annotation on right
    ax_size = fig.add_subplot(gs[1])
    ax_size.barh(range(len(indel_sizes)), indel_sizes, color='steelblue', alpha=0.7)
    ax_size.set_xlabel('Size\n(bp)', fontweight='bold', fontsize=10)
    ax_size.set_yticks([])
    ax_size.set_ylim(ax_main.get_ylim())

    # Add boundary density plot on top
    boundary_density = alignment_matrix.sum(axis=0)

    # Create inset for density
    ax_density = ax_main.inset_axes([0, 1.02, 1, 0.15])
    ax_density.fill_between(range(MONOMER_LENGTH), boundary_density, alpha=0.7, color='steelblue')
    ax_density.set_xlim(0, MONOMER_LENGTH)
    ax_density.set_ylabel('Coverage', fontsize=8)
    ax_density.set_xticks([])
    ax_density.spines['bottom'].set_visible(False)

    plt.savefig(os.path.join(output_dir, 'indel_alignment_heatmap.png'), dpi=300, bbox_inches='tight')
    plt.close()

    print(f"[SUCCESS] Saved alignment heatmap")


def create_genome_wide_plot(catalog_df, output_dir, centromere_sizes=None):
    """
    Create genome-wide visualization of indels across chromosomes.
    Each centromere gets its own subplot with proper coordinate scaling.

    Properly handles haplotype naming: CC1 = Col Centromere 1, LC5 = Ler Centromere 5

    Args:
        centromere_sizes: dict of {(ecotype, chr_num): size} from BED file
    """
    # Use provided centromere sizes or fall back to defaults
    if not centromere_sizes:
        print("[WARNING] No BED file provided - using default centromere sizes for filtering")
        # Default sizes (approximate) - these should come from BED file
        centromere_sizes = {
            ('C', 1): 2_500_000, ('L', 1): 2_500_000,
            ('C', 2): 2_300_000, ('L', 2): 2_300_000,
            ('C', 3): 2_300_000, ('L', 3): 2_300_000,
            ('C', 4): 2_900_000, ('L', 4): 2_900_000,
            ('C', 5): 2_900_000, ('L', 5): 2_900_000,
        }

    # Parse haplotype into ecotype and chromosome
    # Format: [C/L][C/A][1-5] e.g., CC1 = Col Centromere 1, LC5 = Ler Centromere 5
    catalog_df = catalog_df.copy()

    # Extract ecotype (first char: C=Col, L=Ler)
    catalog_df['ecotype'] = catalog_df['haplotype'].str[0]
    # Extract region type (second char: C=Centromere, A=Arms)
    catalog_df['region'] = catalog_df['haplotype'].str[1]
    # Extract chromosome number
    catalog_df['chromosome'] = catalog_df['haplotype'].str.extract(r'(\d+)').astype(float)

    # Filter out outliers that mapped outside expected centromere range
    # These are likely reads that mapped to chromosome arms instead
    original_count = len(catalog_df)
    filtered_rows = []
    for idx, row in catalog_df.iterrows():
        ecotype = row['ecotype']
        chrom = int(row['chromosome']) if pd.notna(row['chromosome']) else 0
        # Look up max size using (ecotype, chromosome) key
        max_size = centromere_sizes.get((ecotype, chrom), 3_000_000)
        if row['ref_pos'] <= max_size:
            filtered_rows.append(idx)
        else:
            print(f"[WARNING] Filtering outlier: {row['haplotype']} ref_pos={row['ref_pos']} > {max_size} (likely mapped to arms)")

    catalog_df = catalog_df.loc[filtered_rows]
    filtered_count = original_count - len(catalog_df)
    if filtered_count > 0:
        print(f"[INFO] Filtered {filtered_count} outliers with ref_pos exceeding expected centromere size")

    # Create a unique identifier for each centromere (ecotype + chromosome)
    catalog_df['centromere_id'] = catalog_df['ecotype'] + '_Chr' + catalog_df['chromosome'].astype(int).astype(str)

    # Get unique centromeres
    centromeres = sorted(catalog_df['centromere_id'].dropna().unique())
    n_cen = len(centromeres)

    if n_cen == 0:
        print("[WARNING] No centromere data found for genome-wide plot")
        return

    # Use different color schemes for Col vs Ler
    col_colors = plt.cm.Blues(np.linspace(0.4, 0.9, 5))
    ler_colors = plt.cm.Oranges(np.linspace(0.4, 0.9, 5))

    # Get centromere boundaries for each unique centromere
    centromere_info = {}
    for cen_id in centromeres:
        cen_data = catalog_df[catalog_df['centromere_id'] == cen_id]
        if len(cen_data) > 0:
            min_pos = cen_data['ref_pos'].min()
            max_pos = cen_data['ref_pos'].max()
            cen_size = max_pos - min_pos
            ecotype = cen_data['ecotype'].iloc[0]
            chrom = int(cen_data['chromosome'].iloc[0])
            centromere_info[cen_id] = {
                'start': min_pos,
                'end': max_pos,
                'size': cen_size,
                'n_indels': len(cen_data),
                'ecotype': ecotype,
                'chromosome': chrom,
                'color': col_colors[chrom-1] if ecotype == 'C' else ler_colors[chrom-1]
            }

    # Create figure with subplots for each centromere
    # Each row: position scatter + density histogram
    fig = plt.figure(figsize=(18, 4 * n_cen))

    for i, cen_id in enumerate(centromeres):
        cen_data = catalog_df[catalog_df['centromere_id'] == cen_id]

        if len(cen_data) == 0:
            continue

        cen_info = centromere_info[cen_id]
        ecotype_name = 'Col' if cen_info['ecotype'] == 'C' else 'Ler'

        # Left: Position scatter with size encoding
        ax_scatter = fig.add_subplot(n_cen, 2, i * 2 + 1)

        # Scatter plot: x=position, y=indel size, color by type
        scatter = ax_scatter.scatter(
            cen_data['ref_pos'],
            cen_data['size'],
            c=cen_data['type'].map({'Insertion': 'blue', 'Deletion': 'red'}).fillna('gray'),
            alpha=0.6,
            s=50,
            edgecolors='black',
            linewidth=0.5
        )

        # Set x-axis to actual centromere coordinates
        padding = max(cen_info['size'] * 0.05, 10000)  # At least 10kb padding
        ax_scatter.set_xlim(cen_info['start'] - padding,
                           cen_info['end'] + padding)

        # Add 178bp reference lines
        max_size = cen_data['size'].max() if not cen_data.empty else 178
        for mult in range(1, int(max_size / 178) + 2):
            ax_scatter.axhline(mult * 178, color='green', linestyle='--', alpha=0.3, linewidth=1)
            if mult <= 5:
                ax_scatter.text(cen_info['start'], mult * 178 + 10, f'{mult}×178',
                               fontsize=7, va='bottom', color='green', alpha=0.7)

        # Format x-axis in Mb
        ax_scatter.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e6:.2f}'))
        ax_scatter.set_xlabel('Reference Position (Mb)', fontweight='bold')
        ax_scatter.set_ylabel('Indel Size (bp)', fontweight='bold')

        # Title with centromere size info
        cen_size_kb = cen_info['size'] / 1000
        ax_scatter.set_title(f'{ecotype_name} Centromere {cen_info["chromosome"]} - Indel Positions\n'
                            f'(n={cen_info["n_indels"]}, span={cen_size_kb:.1f} kb, '
                            f'{cen_info["start"]/1e6:.2f}-{cen_info["end"]/1e6:.2f} Mb)',
                            fontweight='bold', fontsize=11, color=cen_info['color'])

        # Add legend for insertion/deletion
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=8, label='Insertion'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=8, label='Deletion')
        ]
        ax_scatter.legend(handles=legend_elements, loc='upper right', fontsize=8)

        # Right: Coverage/density histogram along centromere
        ax_hist = fig.add_subplot(n_cen, 2, i * 2 + 2)

        # Create histogram of indel positions with proper binning
        positions = cen_data['ref_pos'].values
        if len(positions) > 0:
            # Adaptive binning based on centromere size
            n_bins = min(50, max(10, int(cen_info['size'] / 5000))) if cen_info['size'] > 0 else 10

            # Separate insertions and deletions
            ins_pos = cen_data[cen_data['type'] == 'Insertion']['ref_pos'].values
            del_pos = cen_data[cen_data['type'] == 'Deletion']['ref_pos'].values

            # Use same coordinate range as scatter
            bins = np.linspace(cen_info['start'], cen_info['end'], n_bins + 1)

            if len(ins_pos) > 0:
                ax_hist.hist(ins_pos, bins=bins, alpha=0.6, color='blue', label='Insertions', edgecolor='black')
            if len(del_pos) > 0:
                ax_hist.hist(del_pos, bins=bins, alpha=0.6, color='red', label='Deletions', edgecolor='black')

            # Set same x-axis limits
            ax_hist.set_xlim(cen_info['start'] - padding,
                            cen_info['end'] + padding)

            # Format x-axis in Mb
            ax_hist.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e6:.2f}'))
            ax_hist.set_xlabel('Reference Position (Mb)', fontweight='bold')
            ax_hist.set_ylabel('Number of Indels', fontweight='bold')
            ax_hist.set_title(f'{ecotype_name} Centromere {cen_info["chromosome"]} - Indel Density',
                             fontweight='bold', fontsize=11, color=cen_info['color'])
            ax_hist.legend(loc='upper right', fontsize=8)

            # Add summary statistics
            mean_pos = np.mean(positions)
            ax_hist.axvline(mean_pos, color='black', linestyle='-', linewidth=2, alpha=0.5)

            # Calculate indel density (indels per kb)
            density = len(positions) / (cen_info['size'] / 1000) if cen_info['size'] > 0 else 0
            ax_hist.text(0.02, 0.95, f'Mean: {mean_pos/1e6:.3f} Mb\n'
                                    f'Density: {density:.2f} indels/kb',
                        transform=ax_hist.transAxes, fontsize=9, va='top',
                        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'centromere_indel_coverage.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # Create a normalized comparison plot (relative position within each centromere)
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # Create labels and colors for each centromere
    labels = []
    bar_colors = []
    for cen_id in centromeres:
        info = centromere_info[cen_id]
        ecotype = 'Col' if info['ecotype'] == 'C' else 'Ler'
        labels.append(f'{ecotype}\nChr{info["chromosome"]}')
        bar_colors.append(info['color'])

    # Plot 1: Indel count by centromere
    ax = axes[0, 0]
    counts = [centromere_info[c]['n_indels'] for c in centromeres]
    bars = ax.bar(range(len(centromeres)), counts, color=bar_colors, alpha=0.7, edgecolor='black')
    ax.set_xticks(range(len(centromeres)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_xlabel('Centromere', fontweight='bold')
    ax.set_ylabel('Number of Indels', fontweight='bold')
    ax.set_title('Total Indels per Centromere', fontweight='bold')

    for bar, count in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
               str(count), ha='center', va='bottom', fontsize=9)

    # Plot 2: Centromere sizes (span where indels were found)
    ax = axes[0, 1]
    sizes_kb = [centromere_info[c]['size'] / 1000 for c in centromeres]
    bars = ax.bar(range(len(centromeres)), sizes_kb, color=bar_colors, alpha=0.7, edgecolor='black')
    ax.set_xticks(range(len(centromeres)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_xlabel('Centromere', fontweight='bold')
    ax.set_ylabel('Indel Span (kb)', fontweight='bold')
    ax.set_title('Span of Indel Positions (from data)', fontweight='bold')

    for bar, size in zip(bars, sizes_kb):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
               f'{size:.1f}', ha='center', va='bottom', fontsize=8)

    # Plot 3: Indel density per centromere
    ax = axes[1, 0]
    densities = [centromere_info[c]['n_indels'] / (centromere_info[c]['size'] / 1000)
                 if centromere_info[c]['size'] > 0 else 0 for c in centromeres]
    bars = ax.bar(range(len(centromeres)), densities, color=bar_colors, alpha=0.7, edgecolor='black')
    ax.set_xticks(range(len(centromeres)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_xlabel('Centromere', fontweight='bold')
    ax.set_ylabel('Indel Density (indels/kb)', fontweight='bold')
    ax.set_title('Indel Density per Centromere', fontweight='bold')

    for bar, dens in zip(bars, densities):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
               f'{dens:.2f}', ha='center', va='bottom', fontsize=8)

    # Plot 4: Normalized position distribution (0-100% of each centromere)
    ax = axes[1, 1]

    for i, cen_id in enumerate(centromeres):
        cen_data = catalog_df[catalog_df['centromere_id'] == cen_id]
        if len(cen_data) == 0:
            continue

        cen_info = centromere_info[cen_id]
        ecotype = 'Col' if cen_info['ecotype'] == 'C' else 'Ler'

        # Normalize positions to 0-100%
        if cen_info['size'] > 0:
            norm_pos = (cen_data['ref_pos'] - cen_info['start']) / cen_info['size'] * 100
        else:
            norm_pos = cen_data['ref_pos'] * 0

        # Create density estimate
        if len(norm_pos) > 1:
            from scipy.stats import gaussian_kde
            try:
                kde = gaussian_kde(norm_pos)
                x_range = np.linspace(0, 100, 200)
                density = kde(x_range)
                ax.plot(x_range, density, color=cen_info['color'], linewidth=2,
                       label=f'{ecotype} Chr{cen_info["chromosome"]} (n={len(norm_pos)})')
                ax.fill_between(x_range, density, alpha=0.2, color=cen_info['color'])
            except:
                # Fall back to histogram if KDE fails
                ax.hist(norm_pos, bins=20, alpha=0.3, color=cen_info['color'], density=True,
                       label=f'{ecotype} Chr{cen_info["chromosome"]} (n={len(norm_pos)})')

    ax.set_xlabel('Relative Position in Centromere (%)', fontweight='bold')
    ax.set_ylabel('Density', fontweight='bold')
    ax.set_title('Normalized Indel Distribution\n(All centromeres aligned 0-100%)', fontweight='bold')
    ax.set_xlim(0, 100)
    ax.legend(loc='upper right', fontsize=7)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'centromere_comparison.png'), dpi=300, bbox_inches='tight')
    plt.close()

    print(f"[SUCCESS] Saved centromere coverage plots")


def parse_centromere_sizes_from_bed(bed_file):
    """
    Parse centromere sizes from BED file.
    Returns dict: {(ecotype, chromosome): max_size}

    Expected format: Chr1_CEN_Col, Chr2_CEN_Ler, etc.
    """
    centromere_sizes = {}

    if not bed_file or not os.path.exists(bed_file) or 'NO_FILE' in bed_file:
        print(f"[WARNING] BED file not provided or not found: {bed_file}")
        return centromere_sizes

    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue

            name = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            size = end - start

            # Parse: Chr1_CEN_Col, Chr2_CEN_Ler, etc.
            parts = name.split('_')
            if len(parts) >= 3 and parts[1] == 'CEN':
                chr_num = int(parts[0].replace('Chr', ''))
                ecotype = 'C' if 'Col' in parts[2] else 'L' if 'Ler' in parts[2] else None

                if ecotype:
                    centromere_sizes[(ecotype, chr_num)] = size
                    print(f"[INFO] Parsed centromere size: {ecotype} Chr{chr_num} = {size/1e6:.2f} Mb")

    return centromere_sizes


def main():
    if len(sys.argv) < 5:
        print("Usage: python3 14-analyze_indel_monomer_boundaries.py <indel_sequences.fa> <catalog.tsv> <output_dir> <threads> [bed_file]")
        print("       bed_file is optional but recommended for filtering outliers")
        sys.exit(1)

    indel_sequences_file = sys.argv[1]
    catalog_file = sys.argv[2]
    output_dir = sys.argv[3]
    threads = int(sys.argv[4])
    bed_file = sys.argv[5] if len(sys.argv) > 5 else None

    os.makedirs(output_dir, exist_ok=True)

    print("=" * 80)
    print("CENTROMERIC INDEL BOUNDARY ANALYSIS")
    print("=" * 80)
    print(f"Indel sequences: {indel_sequences_file}")
    print(f"Catalog: {catalog_file}")
    print(f"Output: {output_dir}")
    print(f"BED file: {bed_file if bed_file else 'Not provided'}")
    print(f"Monomer reference: CEN178 ({MONOMER_LENGTH}bp)")
    print("=" * 80)

    # Parse centromere sizes from BED file
    centromere_sizes = parse_centromere_sizes_from_bed(bed_file) if bed_file else None

    # Load catalog
    catalog_df = pd.read_csv(catalog_file, sep='\t')
    print(f"\n[INFO] Loaded {len(catalog_df)} indels from catalog")

    # Check if sequences file exists and has content
    if not os.path.exists(indel_sequences_file) or os.path.getsize(indel_sequences_file) == 0:
        print("[ERROR] No indel sequences found")
        sys.exit(1)

    # Count sequences
    n_seqs = sum(1 for _ in SeqIO.parse(indel_sequences_file, 'fasta'))
    print(f"[INFO] Found {n_seqs} indel sequences")

    # Step 1: Align indels to doubled monomer
    print("\n" + "=" * 80)
    print("STEP 1: Aligning indels to doubled monomer reference")
    print("=" * 80)

    # Try minimap2 first, fall back to pairwise alignment
    alignments = align_with_minimap2(indel_sequences_file, output_dir)

    if alignments is None:
        # Fallback to pairwise alignment
        print("[INFO] Using pairwise alignment fallback...")
        doubled_monomer = create_doubled_monomer()
        alignments = {}

        for record in SeqIO.parse(indel_sequences_file, 'fasta'):
            aln_result = align_to_monomer_pairwise(str(record.seq), doubled_monomer)
            if aln_result:
                alignments[record.id] = aln_result

        print(f"[SUCCESS] Aligned {len(alignments)} indels")

    # Step 2: Analyze boundary distribution
    print("\n" + "=" * 80)
    print("STEP 2: Analyzing boundary distribution (metaprofile)")
    print("=" * 80)

    boundary_analysis = analyze_boundary_distribution(alignments, output_dir)

    if boundary_analysis:
        create_metaprofile_plot(boundary_analysis, output_dir)

        # Create alignment heatmap
        create_alignment_heatmap(indel_sequences_file, alignments, output_dir)

        # Print hotspot summary
        print(f"\n[RESULTS] Total boundaries analyzed: {boundary_analysis['total_boundaries']}")
        if boundary_analysis['uniform_test']['p_value'] is not None:
            p_val = boundary_analysis['uniform_test']['p_value']
            if p_val < 0.05:
                print(f"[RESULTS] Boundaries are NOT uniformly distributed (p = {p_val:.2e})")
                print(f"[RESULTS] This suggests biological clustering of breakpoints!")
            else:
                print(f"[RESULTS] Boundaries appear uniformly distributed (p = {p_val:.2e})")

        print(f"\n[RESULTS] Top 5 hotspot positions in monomer:")
        for pos, count in boundary_analysis['hotspots'][:5]:
            print(f"  Position {pos}: {count} boundaries")

    # Step 3: Analyze in-frame alignment
    print("\n" + "=" * 80)
    print("STEP 3: Analyzing in-frame (whole-monomer) indels")
    print("=" * 80)

    inframe_df = analyze_inframe_alignment(catalog_df, alignments, output_dir)

    # Save results
    inframe_file = os.path.join(output_dir, 'inframe_analysis_results.tsv')
    inframe_df.to_csv(inframe_file, sep='\t', index=False)
    print(f"[SUCCESS] Saved in-frame analysis to {inframe_file}")

    # Create visualization
    create_inframe_analysis_plot(inframe_df, output_dir)

    # Print summary
    aligned_df = inframe_df.dropna(subset=['is_inframe'])
    if len(aligned_df) > 0:
        inframe_count = aligned_df['is_inframe'].sum()
        total = len(aligned_df)
        print(f"\n[RESULTS] In-frame indels: {inframe_count}/{total} ({inframe_count/total*100:.1f}%)")

        # Size-based vs boundary-based analysis
        size_multiple = aligned_df['is_size_multiple'].sum()
        boundary_aligned = aligned_df['is_boundary_aligned'].sum()
        print(f"[RESULTS] Exact 178bp size multiples: {size_multiple}/{total} ({size_multiple/total*100:.1f}%)")
        print(f"[RESULTS] Boundary-aligned: {boundary_aligned}/{total} ({boundary_aligned/total*100:.1f}%)")

    # Step 4: Extract and cluster monomer variants
    print("\n" + "=" * 80)
    print("STEP 4: Extracting and clustering monomer variants")
    print("=" * 80)

    monomer_variants = extract_monomer_variants(indel_sequences_file, alignments, output_dir)
    print(f"[INFO] Extracted {len(monomer_variants)} monomer variants from indels")

    if len(monomer_variants) >= 10:
        # Run dimensionality reduction
        dim_red_results = run_dimensionality_reduction(monomer_variants, output_dir)

        if dim_red_results:
            create_dimensionality_reduction_plot(monomer_variants, dim_red_results, output_dir)

            # Save monomer variants
            variants_df = pd.DataFrame(monomer_variants)
            variants_df['cluster'] = dim_red_results['clusters']
            variants_file = os.path.join(output_dir, 'monomer_variants.tsv')
            variants_df.to_csv(variants_file, sep='\t', index=False)

            # Find consensus for each cluster
            print(f"\n[RESULTS] Found {len(set(dim_red_results['clusters']))} monomer variant clusters")
    else:
        print(f"[WARNING] Not enough monomer variants for clustering (need at least 10)")

    # Step 5: Genome-wide visualization
    print("\n" + "=" * 80)
    print("STEP 5: Creating genome-wide visualization")
    print("=" * 80)

    create_genome_wide_plot(catalog_df, output_dir, centromere_sizes)

    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
    print(f"\nOutput files in {output_dir}:")
    print("  - indel_boundary_metaprofile.png: Metaprofile of boundary positions")
    print("  - inframe_analysis.png: In-frame vs out-of-frame analysis")
    print("  - inframe_analysis_results.tsv: Detailed results table")
    print("  - genome_wide_indels.png: Genome-wide distribution")
    if len(monomer_variants) >= 10:
        print("  - monomer_variant_clustering.png: Monomer variant clusters")
        print("  - monomer_variants.tsv: Extracted monomer sequences")
    print("=" * 80)


if __name__ == "__main__":
    main()
