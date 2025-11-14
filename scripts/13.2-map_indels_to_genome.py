#!/usr/bin/env python3
"""
Map indel sequences back to reference genomes to understand their origin.

This script:
1. Takes indel sequences in FASTA format and catalog with metadata
2. Maps each indel ONLY to its corresponding parental genome:
   - CC* haplotypes (Col centromeres) → map to Col-0 only
   - LC* haplotypes (Ler centromeres) → map to Ler-0 only
3. Compares indel mapping location with original read mapping location
4. Analyzes if indels map near where the read was originally mapping
5. Generates summary statistics of mapping results
"""

import os
import sys
import subprocess
from pathlib import Path
import pandas as pd
from Bio import SeqIO
import pysam
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import LineCollection
import numpy as np


def find_reference_genomes(reference_dir):
    """
    Find all FASTA reference genome files in the reference directory.
    Returns a dictionary of genome_name: genome_path
    """
    genomes = {}

    for fasta_file in Path(reference_dir).glob("*.fa*"):
        if fasta_file.is_file() and not fasta_file.name.endswith('.fai'):
            genome_name = fasta_file.stem
            # Remove common suffixes to get clean genome names
            for suffix in ['_renamed', '.renamed', '.TAIR10', '_genome']:
                genome_name = genome_name.replace(suffix, '')
            genomes[genome_name] = str(fasta_file)

    return genomes


def get_parent_genome(haplotype_code):
    """
    Determine which parent genome based on haplotype code.
    CC* = Col-0, LC* = Ler-0
    Returns 'Col-0' or 'Ler-0'
    """
    if haplotype_code.startswith('C'):
        return 'Col-0'
    elif haplotype_code.startswith('L'):
        return 'Ler-0'
    else:
        return None


def get_original_read_mapping(bam_file, read_id):
    """
    Get the original mapping location of a read from the BAM file.
    Returns dict with chromosome, start, end, or None if not found.
    """
    if not os.path.exists(bam_file):
        return None

    try:
        bamfile = pysam.AlignmentFile(bam_file, "rb")
    except Exception as e:
        print(f"  [WARNING] Could not open BAM file {bam_file}: {e}")
        return None

    for read in bamfile.fetch(until_eof=True):
        if read.query_name == read_id:
            if not read.is_unmapped:
                result = {
                    'chromosome': read.reference_name,
                    'start': read.reference_start,
                    'end': read.reference_end,
                    'mapping_quality': read.mapping_quality
                }
                bamfile.close()
                return result

    bamfile.close()
    return None


def run_minimap2(query_fasta, reference_fasta, output_paf, threads=8):
    """
    Run minimap2 to map query sequences to reference genome.
    Uses parameters suitable for mapping short indel sequences.
    """
    cmd = [
        'minimap2',
        '-x', 'map-ont',  # Oxford Nanopore preset
        '-c',             # Generate CIGAR in PAF
        '-t', str(threads),
        '--secondary=no', # No secondary alignments
        '-N', '10',       # Keep up to 10 best alignments
        reference_fasta,
        query_fasta
    ]

    print(f"Running: {' '.join(cmd)}")

    try:
        with open(output_paf, 'w') as outfile:
            result = subprocess.run(
                cmd,
                stdout=outfile,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
        print(f"Mapping complete: {output_paf}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running minimap2: {e}")
        print(f"stderr: {e.stderr}")
        return False


def parse_paf_file(paf_file):
    """
    Parse PAF file and extract mapping information.
    Returns a list of mapping records.
    """
    mappings = []

    if not os.path.exists(paf_file) or os.path.getsize(paf_file) == 0:
        return mappings

    with open(paf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue

            mapping = {
                'query_name': fields[0],
                'query_length': int(fields[1]),
                'query_start': int(fields[2]),
                'query_end': int(fields[3]),
                'strand': fields[4],
                'target_name': fields[5],
                'target_length': int(fields[6]),
                'target_start': int(fields[7]),
                'target_end': int(fields[8]),
                'num_matches': int(fields[9]),
                'alignment_length': int(fields[10]),
                'mapping_quality': int(fields[11])
            }

            # Calculate additional metrics
            mapping['query_coverage'] = (mapping['query_end'] - mapping['query_start']) / mapping['query_length']
            mapping['identity'] = mapping['num_matches'] / mapping['alignment_length'] if mapping['alignment_length'] > 0 else 0

            mappings.append(mapping)

    return mappings


def parse_bed_files_for_chromosome_structure(parent1_bed, parent2_bed):
    """
    Parse BED files to reconstruct full chromosome structure.
    Returns dictionaries mapping fragment names to (chr_num, offset, length).

    Chromosome order: ARMS_1 -> CEN -> ARMS_2
    """
    chromosome_structure = {}
    chr_lengths = {}
    centromere_positions = {}

    for bed_file, parent in [(parent1_bed, 'Col'), (parent2_bed, 'Ler')]:
        if not bed_file or not os.path.exists(bed_file):
            print(f"[WARNING] BED file not found: {bed_file}")
            continue

        with open(bed_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 3:
                    continue

                fragment_name = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                length = end - start

                # Parse fragment name: Chr1_ARMS_Col_1, Chr1_CEN_Col, etc.
                parts = fragment_name.split('_')
                if len(parts) < 3:
                    continue

                chr_name = parts[0]  # Chr1, Chr2, etc.
                region_type = parts[1]  # ARMS or CEN

                if chr_name not in chr_lengths:
                    chr_lengths[chr_name] = {'ARMS_1': 0, 'CEN': 0, 'ARMS_2': 0}

                if region_type == 'ARMS':
                    arm_num = parts[3] if len(parts) > 3 else '1'  # 1 or 2
                    chr_lengths[chr_name][f'ARMS_{arm_num}'] = length
                elif region_type == 'CEN':
                    chr_lengths[chr_name]['CEN'] = length

    # Calculate offsets and full chromosome lengths
    for chr_name in sorted(chr_lengths.keys()):
        arms1_len = chr_lengths[chr_name]['ARMS_1']
        cen_len = chr_lengths[chr_name]['CEN']
        arms2_len = chr_lengths[chr_name]['ARMS_2']

        # Offsets
        arms1_offset = 0
        cen_offset = arms1_len
        arms2_offset = arms1_len + cen_len
        total_length = arms1_len + cen_len + arms2_len

        # Store chromosome info
        for parent in ['Col', 'Ler']:
            chromosome_structure[f'{chr_name}_ARMS_{parent}_1'] = (chr_name, arms1_offset, arms1_len)
            chromosome_structure[f'{chr_name}_CEN_{parent}'] = (chr_name, cen_offset, cen_len)
            chromosome_structure[f'{chr_name}_ARMS_{parent}_2'] = (chr_name, arms2_offset, arms2_len)

        # Store centromere positions for plotting
        centromere_positions[chr_name] = (cen_offset, cen_offset + cen_len)

        # Store total chromosome length
        chromosome_structure[f'{chr_name}_total_length'] = total_length

    return chromosome_structure, centromere_positions


def convert_fragment_to_full_chromosome_coords(fragment_name, fragment_pos, chromosome_structure):
    """
    Convert position on a fragment to full chromosome coordinates.

    Args:
        fragment_name: e.g., "Chr1_CEN_Col"
        fragment_pos: position within the fragment
        chromosome_structure: dict from parse_bed_files_for_chromosome_structure

    Returns:
        (chr_name, full_chr_position) or (None, None) if conversion fails
    """
    if fragment_name not in chromosome_structure:
        return None, None

    chr_name, offset, length = chromosome_structure[fragment_name]
    full_position = offset + fragment_pos

    return chr_name, full_position


def calculate_distance_to_original(indel_mapping_chr, indel_mapping_start, indel_mapping_end,
                                   original_chr, original_start, original_end):
    """
    Calculate distance between indel mapping and original read mapping.
    Returns distance in bp, or -1 if different chromosomes.
    """
    if indel_mapping_chr != original_chr:
        return -1  # Different chromosomes

    # Calculate distance between regions
    # If they overlap, distance is 0
    if (indel_mapping_start <= original_end and indel_mapping_end >= original_start):
        return 0

    # Otherwise, find minimum distance
    if indel_mapping_end < original_start:
        return original_start - indel_mapping_end
    else:
        return indel_mapping_start - original_end


def analyze_mapping_patterns(mappings_by_indel, catalog_df, original_mappings_by_indel, nonhybrid_dir):
    """
    Analyze mapping patterns for each indel.
    Compare with original read mapping locations.
    """
    analysis = []

    for indel_id, mappings in mappings_by_indel.items():
        # Get indel metadata from catalog
        indel_info = catalog_df[catalog_df['indel_id'] == indel_id]
        if indel_info.empty:
            continue

        indel_row = indel_info.iloc[0]
        haplotype = indel_row['haplotype']
        read_id = indel_row['read_id']
        source_bam = indel_row['source_bam']

        # Get original read mapping location
        original_mapping = original_mappings_by_indel.get(indel_id)

        if not mappings:
            # No mappings found
            analysis.append({
                'indel_id': indel_id,
                'haplotype': haplotype,
                'mapping_mode': source_bam.split('_')[1] if '_' in source_bam else 'unknown',
                'read_id': read_id,
                'mapped': False,
                'mapping_chr': 'N/A',
                'mapping_start': 0,
                'mapping_end': 0,
                'mapping_quality': 0,
                'identity': 0,
                'query_coverage': 0,
                'original_read_chr': original_mapping['chromosome'] if original_mapping else 'N/A',
                'original_read_start': original_mapping['start'] if original_mapping else 0,
                'original_read_end': original_mapping['end'] if original_mapping else 0,
                'distance_to_original': -999,
                'same_chromosome': False,
                'near_original': False
            })
            continue

        # Sort by mapping quality and identity
        mappings.sort(key=lambda x: (x['mapping_quality'], x['identity']), reverse=True)
        best_mapping = mappings[0]

        # Calculate distance to original mapping
        distance_to_original = -999
        same_chromosome = False
        near_original = False

        if original_mapping:
            distance_to_original = calculate_distance_to_original(
                best_mapping['target_name'],
                best_mapping['target_start'],
                best_mapping['target_end'],
                original_mapping['chromosome'],
                original_mapping['start'],
                original_mapping['end']
            )
            same_chromosome = (distance_to_original >= 0)
            near_original = (same_chromosome and distance_to_original < 100000)  # Within 100kb

        # Extract mapping mode from source_bam (e.g., "CC1_mm_ont.bam" -> "mm_ont")
        mapping_mode = source_bam.replace('.bam', '').split('_', 1)[1] if '_' in source_bam else 'unknown'

        analysis.append({
            'indel_id': indel_id,
            'haplotype': haplotype,
            'mapping_mode': mapping_mode,
            'read_id': read_id,
            'mapped': True,
            'mapping_chr': best_mapping['target_name'],
            'mapping_start': best_mapping['target_start'],
            'mapping_end': best_mapping['target_end'],
            'mapping_quality': best_mapping['mapping_quality'],
            'identity': round(best_mapping['identity'], 4),
            'query_coverage': round(best_mapping['query_coverage'], 4),
            'original_read_chr': original_mapping['chromosome'] if original_mapping else 'N/A',
            'original_read_start': original_mapping['start'] if original_mapping else 0,
            'original_read_end': original_mapping['end'] if original_mapping else 0,
            'distance_to_original': distance_to_original,
            'same_chromosome': same_chromosome,
            'near_original': near_original
        })

    return analysis


def create_genomic_plot(analysis_df, catalog_df, output_dir):
    """
    Create genomic plot showing distribution of indels across centromeres.
    Shows both original read location and where indel sequences map.
    """
    print(f"\n[INFO] Creating genomic distribution plots...")

    # Merge analysis with catalog to get size info
    plot_df = analysis_df.merge(catalog_df[['indel_id', 'size', 'type', 'multiple_of_178']],
                                  on='indel_id', how='left')

    # Filter to only mapped indels
    mapped_df = plot_df[plot_df['mapped'] == True].copy()

    if len(mapped_df) == 0:
        print("[WARNING] No mapped indels to plot")
        return

    # Extract chromosome number from chromosome name (e.g., "Chr1" -> 1)
    def extract_chr_num(chr_name):
        if pd.isna(chr_name) or chr_name == 'N/A':
            return None
        import re
        match = re.search(r'(\d+)', str(chr_name))
        return int(match.group(1)) if match else None

    mapped_df['orig_chr_num'] = mapped_df['original_read_chr'].apply(extract_chr_num)
    mapped_df['indel_chr_num'] = mapped_df['mapping_chr'].apply(extract_chr_num)

    # Remove rows where we couldn't extract chromosome number
    mapped_df = mapped_df.dropna(subset=['orig_chr_num', 'indel_chr_num'])

    if len(mapped_df) == 0:
        print("[WARNING] No valid chromosome data to plot")
        return

    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Genomic Distribution of Large Indels in Centromeres',
                 fontsize=16, fontweight='bold')

    # Plot 1: Original read positions across chromosomes
    ax1 = axes[0, 0]
    for mode in mapped_df['mapping_mode'].unique():
        mode_df = mapped_df[mapped_df['mapping_mode'] == mode]
        ax1.scatter(mode_df['original_read_start'], mode_df['orig_chr_num'],
                   alpha=0.6, s=30, label=mode)
    ax1.set_xlabel('Genomic Position (bp)', fontweight='bold')
    ax1.set_ylabel('Chromosome', fontweight='bold')
    ax1.set_title('Original Read Mapping Positions', fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Plot 2: Indel mapping positions across chromosomes
    ax2 = axes[0, 1]
    colors = ['red' if x else 'blue' for x in mapped_df['multiple_of_178']]
    sizes = mapped_df['size'] / 10  # Scale size for visibility
    scatter = ax2.scatter(mapped_df['mapping_start'], mapped_df['indel_chr_num'],
                         c=colors, s=sizes, alpha=0.6, edgecolors='black', linewidth=0.5)
    ax2.set_xlabel('Genomic Position (bp)', fontweight='bold')
    ax2.set_ylabel('Chromosome', fontweight='bold')
    ax2.set_title('Indel Sequence Mapping Positions', fontweight='bold')

    # Add legend for colors
    red_patch = mpatches.Patch(color='red', label='178bp multiple')
    blue_patch = mpatches.Patch(color='blue', label='Not 178bp multiple')
    ax2.legend(handles=[red_patch, blue_patch])
    ax2.grid(True, alpha=0.3)

    # Plot 3: Comparison - same chromosome vs different chromosome
    ax3 = axes[1, 0]
    same_chr_df = mapped_df[mapped_df['same_chromosome'] == True]
    diff_chr_df = mapped_df[mapped_df['same_chromosome'] == False]

    # Plot lines connecting original to indel position for same chromosome
    for _, row in same_chr_df.iterrows():
        ax3.plot([row['original_read_start'], row['mapping_start']],
                [row['orig_chr_num'], row['indel_chr_num']],
                'g-', alpha=0.3, linewidth=1)

    ax3.scatter(same_chr_df['original_read_start'], same_chr_df['orig_chr_num'],
               c='blue', s=50, alpha=0.7, label='Original position', marker='o')
    ax3.scatter(same_chr_df['mapping_start'], same_chr_df['indel_chr_num'],
               c='red', s=50, alpha=0.7, label='Indel maps to', marker='^')

    ax3.set_xlabel('Genomic Position (bp)', fontweight='bold')
    ax3.set_ylabel('Chromosome', fontweight='bold')
    ax3.set_title('Same Chromosome Mappings\n(lines connect read→indel)', fontweight='bold')
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    # Plot 4: Distance distribution for same chromosome mappings
    ax4 = axes[1, 1]
    same_chr_distances = mapped_df[mapped_df['same_chromosome'] == True]['distance_to_original']
    same_chr_distances = same_chr_distances[same_chr_distances >= 0]  # Remove -1 values

    if len(same_chr_distances) > 0:
        bins = np.logspace(0, np.log10(max(same_chr_distances.max(), 1)), 30)
        ax4.hist(same_chr_distances, bins=bins, alpha=0.7, color='green', edgecolor='black')
        ax4.set_xscale('log')
        ax4.axvline(100000, color='red', linestyle='--', linewidth=2,
                   label='100kb threshold')
        ax4.set_xlabel('Distance to Original Position (bp, log scale)', fontweight='bold')
        ax4.set_ylabel('Count', fontweight='bold')
        ax4.set_title('Distance Distribution\n(for same chromosome)', fontweight='bold')
        ax4.legend()
        ax4.grid(True, alpha=0.3, which='both')

        # Add statistics
        median_dist = same_chr_distances.median()
        near_count = sum(same_chr_distances < 100000)
        ax4.text(0.05, 0.95,
                f'Median: {median_dist:.0f} bp\n<100kb: {near_count}/{len(same_chr_distances)}',
                transform=ax4.transAxes, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    else:
        ax4.text(0.5, 0.5, 'No same-chromosome mappings',
                ha='center', va='center', transform=ax4.transAxes)

    plt.tight_layout()
    plot_file = os.path.join(output_dir, 'genomic_distribution_plot.png')
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"[SUCCESS] Saved genomic distribution plot: {plot_file}")

    # Create per-chromosome detail plots
    create_per_chromosome_plots(mapped_df, output_dir)


def create_chromosome_ideogram(analysis_df, catalog_df, chromosome_structure, centromere_positions, output_dir):
    """
    Create chromosome-scale ideogram showing full reconstructed chromosomes with centromeres
    and indel positions mapped along them.

    Uses BED file data to reconstruct full chromosomes from ARMS_1 + CEN + ARMS_2 fragments.
    """
    print(f"\n[INFO] Creating chromosome ideogram visualization...")

    # Extract chromosome lengths from structure
    chr_lengths = {}
    chromosomes = []
    for key in chromosome_structure.keys():
        if key.endswith('_total_length'):
            chr_name = key.replace('_total_length', '')
            chr_lengths[chr_name] = chromosome_structure[key]
            chromosomes.append(chr_name)

    chromosomes = sorted(chromosomes)

    if not chromosomes:
        print("[WARNING] No chromosome structure found from BED files")
        return

    print(f"[INFO] Reconstructed chromosomes: {', '.join(chromosomes)}")
    for chr_name in chromosomes:
        cen_start, cen_end = centromere_positions.get(chr_name, (0, 0))
        print(f"  {chr_name}: {chr_lengths[chr_name]/1e6:.2f} Mb (centromere: {cen_start/1e6:.2f}-{cen_end/1e6:.2f} Mb)")

    # Merge analysis with catalog
    plot_df = analysis_df.merge(catalog_df[['indel_id', 'size', 'type', 'multiple_of_178']],
                                  on='indel_id', how='left')

    # Filter to mapped indels only
    mapped_df = plot_df[plot_df['mapped'] == True].copy()

    if len(mapped_df) == 0:
        print("[WARNING] No mapped indels for chromosome ideogram")
        return

    # Convert fragment coordinates to full chromosome coordinates
    print(f"[INFO] Converting fragment coordinates to full chromosome coordinates...")
    converted_coords = []
    for _, row in mapped_df.iterrows():
        # Convert indel mapping position
        indel_chr, indel_pos = convert_fragment_to_full_chromosome_coords(
            row['mapping_chr'], row['mapping_start'], chromosome_structure
        )

        # Convert original read mapping position
        orig_chr, orig_pos = convert_fragment_to_full_chromosome_coords(
            row['original_read_chr'], row['original_read_start'], chromosome_structure
        )

        if indel_chr and orig_chr:
            converted_coords.append({
                **row.to_dict(),
                'full_chr_indel': indel_chr,
                'full_pos_indel': indel_pos,
                'full_chr_original': orig_chr,
                'full_pos_original': orig_pos
            })

    if not converted_coords:
        print("[WARNING] No coordinates could be converted")
        return

    converted_df = pd.DataFrame(converted_coords)
    print(f"[SUCCESS] Converted {len(converted_df)} indel coordinates")

    # Create figure with chromosome ideograms
    n_chr = len(chromosomes)
    fig, axes = plt.subplots(n_chr, 1, figsize=(20, 4*n_chr), sharex=False)
    if n_chr == 1:
        axes = [axes]

    fig.suptitle('Chromosome Ideogram - Large Indel Distribution\n(Reconstructed from ARMS + CEN fragments)',
                 fontsize=18, fontweight='bold')

    for idx, chr_name in enumerate(chromosomes):
        ax = axes[idx]

        chr_length = chr_lengths[chr_name]
        cen_start, cen_end = centromere_positions.get(chr_name, (0, 0))

        # Draw chromosome as a horizontal bar
        chromosome_height = 0.5
        ax.add_patch(mpatches.Rectangle((0, 0.25), chr_length, chromosome_height,
                                       facecolor='lightgray', edgecolor='black', linewidth=2))

        # Draw centromere region (darker color)
        if cen_start > 0 and cen_end > 0:
            ax.add_patch(mpatches.Rectangle((cen_start, 0.15), cen_end - cen_start, chromosome_height + 0.2,
                                           facecolor='darkgray', edgecolor='red', linewidth=2, alpha=0.7))

        # Get indels for this chromosome
        chr_indels = converted_df[converted_df['full_chr_indel'] == chr_name]
        chr_original = converted_df[converted_df['full_chr_original'] == chr_name]

        # Plot original read positions (below chromosome)
        if len(chr_original) > 0:
            for _, row in chr_original.iterrows():
                pos = row['full_pos_original']
                ax.plot([pos, pos], [0, 0.2], 'b-', linewidth=1.5, alpha=0.6)
                ax.scatter([pos], [0.1], c='blue', s=80, alpha=0.8, marker='v',
                          edgecolors='black', linewidths=0.5)

        # Plot indel mapping positions (above chromosome)
        if len(chr_indels) > 0:
            for _, row in chr_indels.iterrows():
                pos = row['full_pos_indel']
                color = 'red' if row['multiple_of_178'] else 'orange'
                size = min(row['size'] / 5, 150)  # Scale size, cap at 150

                ax.plot([pos, pos], [0.8, 1.1], color=color, linewidth=1.5, alpha=0.6)
                ax.scatter([pos], [1.0], c=color, s=size, alpha=0.8, marker='^',
                          edgecolors='black', linewidths=0.5)

                # If same chromosome as original, draw connecting line
                if row['full_chr_indel'] == row['full_chr_original']:
                    orig_pos = row['full_pos_original']
                    ax.plot([orig_pos, pos], [0.1, 1.0], 'g-', linewidth=0.8, alpha=0.3)

        # Formatting
        ax.set_xlim(0, chr_length)
        ax.set_ylim(-0.1, 1.3)
        ax.set_ylabel(f'{chr_name}\n({chr_length/1e6:.1f} Mb)', fontweight='bold', fontsize=11)
        ax.set_yticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)

        # Add scale bar (in Mb)
        ax.set_xlabel('Position (Mb)', fontweight='bold') if idx == n_chr - 1 else ax.set_xlabel('')
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e6:.1f}'))

        # Add statistics text
        n_original = len(chr_original)
        n_mapped = len(chr_indels)
        n_same_chr = sum(chr_indels['full_chr_indel'] == chr_indels['full_chr_original'])
        n_178_multiple = sum(chr_indels['multiple_of_178'])

        stats_text = f'Reads: {n_original} | Indels: {n_mapped} | Same chr: {n_same_chr} | 178bp mult: {n_178_multiple}'
        ax.text(0.02, 0.95, stats_text,
               transform=ax.transAxes, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
               fontsize=9, fontweight='bold')

    # Add legend
    legend_elements = [
        mpatches.Patch(facecolor='darkgray', edgecolor='red', label='Centromere'),
        plt.Line2D([0], [0], marker='v', color='w', markerfacecolor='blue', markersize=10,
                  label='Original read position'),
        plt.Line2D([0], [0], marker='^', color='w', markerfacecolor='red', markersize=10,
                  label='Indel mapping (178bp multiple)'),
        plt.Line2D([0], [0], marker='^', color='w', markerfacecolor='orange', markersize=10,
                  label='Indel mapping (not 178bp multiple)'),
        plt.Line2D([0], [0], color='green', linewidth=2, alpha=0.5,
                  label='Connection (same chromosome)')
    ]
    fig.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(0.98, 0.98),
              fontsize=11, framealpha=0.9)

    plt.tight_layout(rect=[0, 0, 1, 0.98])
    plot_file = os.path.join(output_dir, 'chromosome_ideogram.png')
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"[SUCCESS] Saved chromosome ideogram: {plot_file}")


def create_per_chromosome_plots(mapped_df, output_dir):
    """
    Create detailed plots for each chromosome showing indel positions.
    """
    print(f"\n[INFO] Creating per-chromosome detail plots...")

    chromosomes = sorted(mapped_df['orig_chr_num'].unique())

    if len(chromosomes) == 0:
        return

    n_chrs = len(chromosomes)
    n_cols = min(3, n_chrs)
    n_rows = (n_chrs + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 4*n_rows))
    if n_rows == 1 and n_cols == 1:
        axes = np.array([[axes]])
    elif n_rows == 1 or n_cols == 1:
        axes = axes.reshape(n_rows, n_cols)

    fig.suptitle('Per-Chromosome Indel Distribution', fontsize=16, fontweight='bold')

    for idx, chr_num in enumerate(chromosomes):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes[row, col]

        chr_df = mapped_df[mapped_df['orig_chr_num'] == chr_num]

        # Plot original positions
        ax.scatter(chr_df['original_read_start'],
                  [1]*len(chr_df),
                  c='blue', s=100, alpha=0.7,
                  label='Original read', marker='o')

        # Plot indel mapping positions
        same_chr = chr_df[chr_df['same_chromosome'] == True]
        diff_chr = chr_df[chr_df['same_chromosome'] == False]

        ax.scatter(same_chr['mapping_start'],
                  [2]*len(same_chr),
                  c='green', s=100, alpha=0.7,
                  label='Indel (same chr)', marker='^')

        ax.scatter(diff_chr['mapping_start'],
                  [2]*len(diff_chr),
                  c='red', s=100, alpha=0.7,
                  label='Indel (diff chr)', marker='x')

        # Draw lines connecting same chromosome mappings
        for _, row_data in same_chr.iterrows():
            ax.plot([row_data['original_read_start'], row_data['mapping_start']],
                   [1, 2], 'g-', alpha=0.3, linewidth=1)

        ax.set_xlabel('Genomic Position (bp)', fontweight='bold')
        ax.set_yticks([1, 2])
        ax.set_yticklabels(['Original\nRead', 'Indel\nMaps to'])
        ax.set_title(f'Chromosome {int(chr_num)}', fontweight='bold')
        ax.legend(fontsize=8, loc='upper right')
        ax.grid(True, alpha=0.3, axis='x')

        # Add statistics
        total = len(chr_df)
        same = len(same_chr)
        near = sum(same_chr['near_original'] == True)
        ax.text(0.02, 0.98,
               f'Total: {total}\nSame chr: {same}\nNear (<100kb): {near}',
               transform=ax.transAxes, verticalalignment='top',
               fontsize=8,
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Hide unused subplots
    for idx in range(len(chromosomes), n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        axes[row, col].axis('off')

    plt.tight_layout()
    plot_file = os.path.join(output_dir, 'per_chromosome_distribution.png')
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"[SUCCESS] Saved per-chromosome plot: {plot_file}")


def main():
    if len(sys.argv) != 9:
        print("Usage: python3 13.2-map_indels_to_genome.py <indel_sequences.fa> <catalog.tsv> <nonhybrid_dir> <reference_genomes_dir> <output_dir> <threads> <parent1_bed> <parent2_bed>")
        sys.exit(1)

    indel_fasta = sys.argv[1]
    catalog_file = sys.argv[2]
    nonhybrid_dir = sys.argv[3]
    reference_dir = sys.argv[4]
    output_dir = sys.argv[5]
    threads = int(sys.argv[6])
    parent1_bed = sys.argv[7]
    parent2_bed = sys.argv[8]

    print("\n" + "="*80)
    print("=== MAPPING INDELS TO PARENTAL GENOMES ===")
    print("="*80)
    print(f"[CONFIG] Indel sequences: {indel_fasta}")
    print(f"[CONFIG] Catalog file: {catalog_file}")
    print(f"[CONFIG] Nonhybrid directory: {nonhybrid_dir}")
    print(f"[CONFIG] Reference directory: {reference_dir}")
    print(f"[CONFIG] Output directory: {output_dir}")
    print(f"[CONFIG] Threads: {threads}")
    print(f"[CONFIG] Parent1 BED: {parent1_bed}")
    print(f"[CONFIG] Parent2 BED: {parent2_bed}")
    print("="*80)

    # Parse BED files to get chromosome structure
    print(f"\n[INFO] Parsing BED files to reconstruct chromosome structure...")
    chromosome_structure, centromere_positions = parse_bed_files_for_chromosome_structure(parent1_bed, parent2_bed)
    print(f"[SUCCESS] Parsed chromosome structure from BED files")

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Check if indel sequences file exists and has content
    if not os.path.exists(indel_fasta) or os.path.getsize(indel_fasta) == 0:
        print("[WARNING] No indel sequences to map. Exiting.")
        return

    # Load catalog
    print(f"\n[INFO] Loading indel catalog...")
    catalog_df = pd.read_csv(catalog_file, sep='\t')
    print(f"[SUCCESS] Loaded {len(catalog_df)} indels from catalog")

    # Find reference genomes
    genomes = find_reference_genomes(reference_dir)
    print(f"\n[INFO] Found {len(genomes)} reference genomes:")
    for name, path in genomes.items():
        print(f"  - {name}: {path}")

    if not genomes:
        print("[ERROR] No reference genomes found!")
        sys.exit(1)

    # Group indels by parent genome and mapping mode
    print(f"\n[INFO] Grouping indels by parent genome and mapping mode...")
    indels_by_parent_mode = {}

    for _, row in catalog_df.iterrows():
        indel_id = row['indel_id']
        haplotype = row['haplotype']
        source_bam = row['source_bam']

        # Determine parent genome
        parent = get_parent_genome(haplotype)
        if not parent:
            print(f"  [WARNING] Could not determine parent for {haplotype}")
            continue

        # Extract mapping mode (e.g., "CC1_mm_ont.bam" -> "mm_ont")
        mapping_mode = source_bam.replace('.bam', '').split('_', 1)[1] if '_' in source_bam else 'unknown'

        # Create key for grouping
        key = (parent, mapping_mode)
        if key not in indels_by_parent_mode:
            indels_by_parent_mode[key] = []
        indels_by_parent_mode[key].append(row)

    print(f"[INFO] Grouped indels into {len(indels_by_parent_mode)} parent-mode combinations:")
    for (parent, mode), indels in indels_by_parent_mode.items():
        print(f"  - {parent} + {mode}: {len(indels)} indels")

    # Get original read mappings for all indels
    print(f"\n[INFO] Retrieving original read mapping locations...")
    mappings_dir = os.path.join(nonhybrid_dir, "mappings")
    original_mappings_by_indel = {}

    for _, row in catalog_df.iterrows():
        indel_id = row['indel_id']
        read_id = row['read_id']
        source_bam = row['source_bam']

        bam_path = os.path.join(mappings_dir, source_bam)
        original_mapping = get_original_read_mapping(bam_path, read_id)

        if original_mapping:
            original_mappings_by_indel[indel_id] = original_mapping

    print(f"[SUCCESS] Retrieved original mappings for {len(original_mappings_by_indel)} indels")

    # Map each group to its corresponding parent genome
    mappings_by_indel = {}

    for (parent, mode), indels_list in indels_by_parent_mode.items():
        if parent not in genomes:
            print(f"\n[WARNING] Reference genome {parent} not found, skipping...")
            continue

        genome_path = genomes[parent]

        print(f"\n{'='*80}")
        print(f"[PROGRESS] Mapping {len(indels_list)} indels from {parent} ({mode}) to {parent} genome")
        print(f"{'='*80}")

        # Create temporary FASTA with only these indels
        temp_fasta = os.path.join(output_dir, f"temp_{parent}_{mode}.fa")
        temp_records = []
        for indel_row in indels_list:
            indel_id = indel_row['indel_id']
            # Find sequence in original fasta
            for record in SeqIO.parse(indel_fasta, 'fasta'):
                if record.id == indel_id:
                    temp_records.append(record)
                    break

        if not temp_records:
            print(f"  [WARNING] No sequences found for this group")
            continue

        SeqIO.write(temp_records, temp_fasta, 'fasta')
        print(f"  [INFO] Created temporary FASTA with {len(temp_records)} sequences")

        # Run minimap2
        output_paf = os.path.join(output_dir, f"indels_{parent}_{mode}.paf")
        success = run_minimap2(temp_fasta, genome_path, output_paf, threads)

        # Clean up temp file
        os.remove(temp_fasta)

        if success:
            # Parse PAF file
            mappings = parse_paf_file(output_paf)
            print(f"  [SUCCESS] Found {len(mappings)} mappings")

            # Group mappings by indel ID
            for mapping in mappings:
                indel_id = mapping['query_name']
                if indel_id not in mappings_by_indel:
                    mappings_by_indel[indel_id] = []
                mappings_by_indel[indel_id].append(mapping)

    # Analyze mapping patterns
    print(f"\n{'='*80}")
    print("[PROGRESS] Analyzing Mapping Patterns")
    print(f"{'='*80}")
    analysis = analyze_mapping_patterns(mappings_by_indel, catalog_df, original_mappings_by_indel, nonhybrid_dir)

    # Save analysis results
    analysis_file = os.path.join(output_dir, "indel_mapping_analysis.tsv")
    df = pd.DataFrame(analysis)
    df.to_csv(analysis_file, sep='\t', index=False)
    print(f"[SUCCESS] Saved mapping analysis to: {analysis_file}")

    # Create genomic plots
    if len(analysis) > 0:
        try:
            # Create chromosome ideogram (main visualization with proper scaling)
            create_chromosome_ideogram(df, catalog_df, chromosome_structure, centromere_positions, output_dir)
        except Exception as e:
            print(f"[WARNING] Could not create chromosome ideogram: {e}")
            import traceback
            traceback.print_exc()

    # Print summary
    print(f"\n{'='*80}")
    print("=== MAPPING SUMMARY ===")
    print(f"{'='*80}")
    print(f"Total indels analyzed: {len(analysis)}")

    if len(analysis) == 0:
        print("[WARNING] No indels were analyzed. Check that reference genomes were found correctly.")
    else:
        mapped = sum(1 for x in analysis if x['mapped'])
        print(f"Mapped indels: {mapped} ({mapped/len(analysis)*100:.1f}%)")
        print(f"Unmapped indels: {len(analysis) - mapped}")

        if mapped > 0:
            same_chr = sum(1 for x in analysis if x['same_chromosome'])
            near_original = sum(1 for x in analysis if x['near_original'])

            print(f"\nIndels mapping to same chromosome as original read: {same_chr} ({same_chr/mapped*100:.1f}%)")
            print(f"Indels mapping near original location (<100kb): {near_original} ({near_original/mapped*100:.1f}%)")

            # Distribution by mapping mode
            print(f"\nMapping statistics by mode:")
            for mode in df['mapping_mode'].unique():
                mode_df = df[df['mapping_mode'] == mode]
                mode_mapped = sum(mode_df['mapped'])
                mode_near = sum(mode_df['near_original'])
                print(f"  {mode}: {mode_mapped}/{len(mode_df)} mapped ({mode_mapped/len(mode_df)*100:.1f}%), {mode_near} near original")

    print(f"\n{'='*80}")
    print("=== MAPPING COMPLETE ===")
    print(f"{'='*80}")


if __name__ == "__main__":
    main()
