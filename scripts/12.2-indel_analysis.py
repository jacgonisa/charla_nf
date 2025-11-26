#!/usr/bin/env python3
"""
Indel Analysis Script for Non-Hybrid Reads

This script analyzes indels >50bp from BAM alignment files and generates
comprehensive visualizations with standardized axes for cross-condition comparison.

Features:
- Filters for large indels (>50bp)
- Standardized x and y axes across all subplots
- Full-range histograms showing complete size distribution
- Zoomed histograms (0-1000bp) for detailed view
- 178bp satellite repeat reference lines
- Publication-ready figures with clean styling

Visualization approach inspired by sv-satellite-analyzer:
- Consistent color schemes (#3498db for Col, #e74c3c for Ler)
- Clean white background with subtle grid
- Hatching patterns to distinguish insertions from deletions
"""
import os
import sys
import pysam
import csv
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from pathlib import Path
import numpy as np


def parse_haplotype(haplotype_code):
    """
    Parse haplotype code to extract components.
    E.g., 'CA1' -> {'background': 'Col', 'region': 'ARMS', 'chromosome': '1'}
    """
    if len(haplotype_code) != 3:
        return None
    
    background = 'Col' if haplotype_code[0] == 'C' else 'Ler' if haplotype_code[0] == 'L' else 'Unknown'
    region = 'ARMS' if haplotype_code[1] == 'A' else 'CEN' if haplotype_code[1] == 'C' else 'Unknown'
    chromosome = haplotype_code[2]
    
    return {
        'background': background,
        'region': region,
        'chromosome': chromosome
    }


def analyze_bam_file(bam_path, output_prefix):
    """
    Analyze a BAM file for indels and generate statistics.
    Returns dictionary with stats and indel details for combined plotting.
    """
    if not os.path.exists(bam_path):
        print(f"Warning: BAM file not found: {bam_path}")
        return None
    
    try:
        bamfile = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        print(f"Error opening BAM file {bam_path}: {e}")
        return None
    
    # Collect data for analysis
    indel_details = []
    total_reads = 0
    mapped_reads = 0
    reads_with_indels = set()
    
    # Parse all reads
    for read in bamfile.fetch(until_eof=True):
        total_reads += 1
        
        if read.is_unmapped:
            continue
        mapped_reads += 1
        
        if read.cigartuples is None:
            continue
            
        query_pos = 0
        ref_pos = read.reference_start
        
        for op, length in read.cigartuples:
            if op == 0:  # Match or mismatch (M)
                query_pos += length
                ref_pos += length
            elif op == 1 and length > 50:  # Insertion > 50bp
                reads_with_indels.add(read.query_name)
                indel_details.append({
                    "Read_ID": read.query_name,
                    "Type": "Insertion",
                    "Size": length,
                    "Read_Pos": query_pos,
                    "Ref_Pos": ref_pos
                })
                query_pos += length
            elif op == 2 and length > 50:  # Deletion > 50bp
                reads_with_indels.add(read.query_name)
                indel_details.append({
                    "Read_ID": read.query_name,
                    "Type": "Deletion",
                    "Size": length,
                    "Read_Pos": query_pos,
                    "Ref_Pos": ref_pos
                })
                ref_pos += length
            elif op in [3, 4, 5]:
                if op in [4, 5]:
                    query_pos += length
                if op == 3:
                    ref_pos += length
    
    bamfile.close()
    
    # Generate detailed TSV output
    tsv_output = f"{output_prefix}_indel_positions.tsv"
    with open(tsv_output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["Read_ID", "Type", "Size", "Read_Pos", "Ref_Pos"], delimiter="\t")
        writer.writeheader()
        for row in indel_details:
            writer.writerow(row)
    
    # Calculate stats
    insertion_lengths = [x["Size"] for x in indel_details if x["Type"] == "Insertion"]
    deletion_lengths = [x["Size"] for x in indel_details if x["Type"] == "Deletion"]
    num_indels = len(indel_details)
    avg_indel_size = np.mean([x["Size"] for x in indel_details]) if num_indels > 0 else 0
    indel_density = num_indels / mapped_reads if mapped_reads > 0 else 0
    reads_with_indels_count = len(reads_with_indels)
    
    print(f"Analysis complete for {output_prefix}:")
    print(f"  Total reads: {total_reads}")
    print(f"  Mapped reads: {mapped_reads}")
    print(f"  Indels (>50bp): {num_indels}")
    print(f"  Reads with indels: {reads_with_indels_count}")
    print(f"  Insertions: {len(insertion_lengths)}")
    print(f"  Deletions: {len(deletion_lengths)}")
    print(f"  Average indel size: {avg_indel_size:.2f} bp")
    print(f"  Indel density: {indel_density:.4f} indels/read")
    print(f"  Detailed results: {tsv_output}")
    
    return {
        'total_reads': total_reads,
        'mapped_reads': mapped_reads,
        'indels': num_indels,
        'insertions': len(insertion_lengths),
        'deletions': len(deletion_lengths),
        'avg_indel_size': round(avg_indel_size, 2),
        'indel_density': round(indel_density, 4),
        'reads_with_indels': reads_with_indels_count,
        'insertion_lengths': insertion_lengths,
        'deletion_lengths': deletion_lengths
    }


def create_combined_plots_by_mapping_mode(results_data, output_dir):
    """
    Create comprehensive plots for indel size distribution (>50bp).

    Generates:
    - Version 1: One plot per haplotype (separate Col and Ler) - FULL RANGE
    - Version 2: Combined plot with 4 distinct visual styles - FULL RANGE
    - Version 3: Zoomed plots (0-1000bp) for detailed view

    All plots use standardized axes for comparison across conditions.
    """
    # Set seaborn style for cleaner plots
    sns.set_style("white")

    # Group results by mapping mode
    mapping_modes = {}
    for item in results_data:
        mode = item['mapping_mode']
        if mode not in mapping_modes:
            mapping_modes[mode] = []
        mapping_modes[mode].append(item)

    # ========================================================================
    # CALCULATE GLOBAL AXIS LIMITS FOR STANDARDIZATION
    # ========================================================================
    global_max_size = 0
    global_max_count_full = 0
    global_max_count_zoom = 0

    for mode, mode_data in mapping_modes.items():
        for item in mode_data:
            all_indels = item['insertion_lengths'] + item['deletion_lengths']
            if all_indels:
                global_max_size = max(global_max_size, max(all_indels))

                # Calculate max counts for full range
                hist_full, _ = np.histogram(all_indels, bins=50, range=(0, max(all_indels)))
                global_max_count_full = max(global_max_count_full, max(hist_full))

                # Calculate max counts for zoomed range (0-1000bp)
                indels_zoom = [x for x in all_indels if x <= 1000]
                if indels_zoom:
                    hist_zoom, _ = np.histogram(indels_zoom, bins=50, range=(0, 1000))
                    global_max_count_zoom = max(global_max_count_zoom, max(hist_zoom))

    # Add 10% padding to y-axis
    global_max_count_full = int(global_max_count_full * 1.1) + 1
    global_max_count_zoom = int(global_max_count_zoom * 1.1) + 1

    print(f"\n[INFO] Global axis limits calculated:")
    print(f"  Max indel size: {global_max_size} bp")
    print(f"  Max count (full range): {global_max_count_full}")
    print(f"  Max count (0-1000bp zoom): {global_max_count_zoom}")

    # ========================================================================
    # VERSION 1: Separate plots per haplotype (background) - FULL RANGE
    # ========================================================================
    for mode, mode_data in mapping_modes.items():
        backgrounds = {}
        for item in mode_data:
            bg = item['background']
            if bg not in backgrounds:
                backgrounds[bg] = []
            backgrounds[bg].append(item)

        for background, bg_data in backgrounds.items():
            fig, axes = plt.subplots(2, 5, figsize=(20, 8))
            fig.suptitle(f'Indel Size Distribution (>50bp) - {mode} - {background}',
                        fontsize=16, fontweight='bold')

            for item in bg_data:
                region = item['region']
                chrom = item['chromosome']
                row = 0 if region == 'ARMS' else 1
                col = int(chrom) - 1

                if col < 0 or col >= 5:
                    continue

                ax = axes[row, col]
                insertions = item['insertion_lengths']
                deletions = item['deletion_lengths']
                all_indels = insertions + deletions

                if all_indels:
                    max_size = max(all_indels)
                    bins = np.linspace(0, max_size, 50)

                    ax.hist(insertions, bins=bins, color='#3498db', alpha=0.6,
                           label=f'Ins (n={len(insertions)})', edgecolor='darkblue', linewidth=0.5)
                    ax.hist(deletions, bins=bins, color='#e74c3c', alpha=0.6,
                           label=f'Del (n={len(deletions)})', edgecolor='darkred', linewidth=0.5)

                    # Add 178bp multiple reference lines
                    for i in range(1, int(max_size / 178) + 2):
                        ref_line = i * 178
                        if ref_line <= max_size + 100:
                            ax.axvline(ref_line, color='gray', linestyle='--',
                                     alpha=0.3, linewidth=0.8)

                    ax.legend(fontsize=8, loc='upper right', framealpha=0.9)
                else:
                    ax.text(0.5, 0.5, 'No indels >50bp', ha='center', va='center',
                           transform=ax.transAxes, fontsize=10, style='italic', color='gray')

                # STANDARDIZED AXES
                ax.set_xlim(0, global_max_size * 1.02)
                ax.set_ylim(0, global_max_count_full)

                ax.set_title(f'Chr{chrom} - {region}', fontsize=10, fontweight='bold')
                if col == 0:
                    ax.set_ylabel('Count', fontsize=9)
                if row == 1:
                    ax.set_xlabel('Indel size (bp)', fontsize=9)
                ax.tick_params(labelsize=8)
                ax.grid(True, alpha=0.2, axis='y')

            plt.tight_layout()
            plot_filename = os.path.join(output_dir, f'indel_histogram_{mode}_{background}_full.png')
            plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"Generated plot: {plot_filename}")

    # ========================================================================
    # VERSION 2: Combined plots (all backgrounds) - FULL RANGE
    # ========================================================================
    for mode, mode_data in mapping_modes.items():
        fig, axes = plt.subplots(2, 5, figsize=(20, 8))
        fig.suptitle(f'Indel Size Distribution (>50bp, Combined) - {mode}',
                    fontsize=16, fontweight='bold')

        aggregated_data = {}
        for item in mode_data:
            region = item['region']
            chrom = item['chromosome']
            background = item['background']
            key = (region, chrom)

            if key not in aggregated_data:
                aggregated_data[key] = {
                    'Col': {'insertions': [], 'deletions': []},
                    'Ler': {'insertions': [], 'deletions': []}
                }

            aggregated_data[key][background]['insertions'].extend(item['insertion_lengths'])
            aggregated_data[key][background]['deletions'].extend(item['deletion_lengths'])

        for (region, chrom), data in aggregated_data.items():
            row = 0 if region == 'ARMS' else 1
            col = int(chrom) - 1

            if col < 0 or col >= 5:
                continue

            ax = axes[row, col]
            all_indels = (data['Col']['insertions'] + data['Col']['deletions'] +
                         data['Ler']['insertions'] + data['Ler']['deletions'])

            if all_indels:
                max_size = max(all_indels)
                bins = np.linspace(0, max_size, 50)

                # Plot with distinct colors (sv-satellite-analyzer style)
                if data['Col']['insertions']:
                    ax.hist(data['Col']['insertions'], bins=bins,
                           color='#3498db', alpha=0.5, edgecolor='darkblue', linewidth=0.5,
                           label=f"Col-Ins (n={len(data['Col']['insertions'])})")

                if data['Col']['deletions']:
                    ax.hist(data['Col']['deletions'], bins=bins,
                           color='#85c1e9', alpha=0.5, edgecolor='darkblue', linewidth=0.5,
                           hatch='///', label=f"Col-Del (n={len(data['Col']['deletions'])})")

                if data['Ler']['insertions']:
                    ax.hist(data['Ler']['insertions'], bins=bins,
                           color='#e74c3c', alpha=0.5, edgecolor='darkred', linewidth=0.5,
                           label=f"Ler-Ins (n={len(data['Ler']['insertions'])})")

                if data['Ler']['deletions']:
                    ax.hist(data['Ler']['deletions'], bins=bins,
                           color='#f1948a', alpha=0.5, edgecolor='darkred', linewidth=0.5,
                           hatch='///', label=f"Ler-Del (n={len(data['Ler']['deletions'])})")

                # Add 178bp multiple reference lines
                for i in range(1, int(max_size / 178) + 2):
                    ref_line = i * 178
                    if ref_line <= max_size + 100:
                        ax.axvline(ref_line, color='gray', linestyle='--',
                                 alpha=0.3, linewidth=0.8)

                ax.legend(fontsize=6, loc='upper right', framealpha=0.9)
            else:
                ax.text(0.5, 0.5, 'No indels >50bp', ha='center', va='center',
                       transform=ax.transAxes, fontsize=10, style='italic', color='gray')

            # STANDARDIZED AXES
            ax.set_xlim(0, global_max_size * 1.02)
            ax.set_ylim(0, global_max_count_full)

            ax.set_title(f'Chr{chrom} - {region}', fontsize=10, fontweight='bold')
            if col == 0:
                ax.set_ylabel('Count', fontsize=9)
            if row == 1:
                ax.set_xlabel('Indel size (bp)', fontsize=9)
            ax.tick_params(labelsize=8)
            ax.grid(True, alpha=0.2, axis='y')

        plt.tight_layout()
        plot_filename = os.path.join(output_dir, f'indel_histogram_{mode}_combined_full.png')
        plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Generated combined plot: {plot_filename}")

    # ========================================================================
    # VERSION 3: ZOOMED plots (0-1000bp) for detailed view
    # ========================================================================
    for mode, mode_data in mapping_modes.items():
        fig, axes = plt.subplots(2, 5, figsize=(20, 8))
        fig.suptitle(f'Indel Size Distribution (>50bp, 0-1000bp zoom) - {mode}',
                    fontsize=16, fontweight='bold')

        aggregated_data = {}
        for item in mode_data:
            region = item['region']
            chrom = item['chromosome']
            background = item['background']
            key = (region, chrom)

            if key not in aggregated_data:
                aggregated_data[key] = {
                    'Col': {'insertions': [], 'deletions': []},
                    'Ler': {'insertions': [], 'deletions': []}
                }

            # Filter for 0-1000bp range
            aggregated_data[key][background]['insertions'].extend(
                [x for x in item['insertion_lengths'] if x <= 1000])
            aggregated_data[key][background]['deletions'].extend(
                [x for x in item['deletion_lengths'] if x <= 1000])

        for (region, chrom), data in aggregated_data.items():
            row = 0 if region == 'ARMS' else 1
            col = int(chrom) - 1

            if col < 0 or col >= 5:
                continue

            ax = axes[row, col]
            all_indels = (data['Col']['insertions'] + data['Col']['deletions'] +
                         data['Ler']['insertions'] + data['Ler']['deletions'])

            if all_indels:
                bins = np.linspace(0, 1000, 50)

                if data['Col']['insertions']:
                    ax.hist(data['Col']['insertions'], bins=bins,
                           color='#3498db', alpha=0.5, edgecolor='darkblue', linewidth=0.5,
                           label=f"Col-Ins (n={len(data['Col']['insertions'])})")

                if data['Col']['deletions']:
                    ax.hist(data['Col']['deletions'], bins=bins,
                           color='#85c1e9', alpha=0.5, edgecolor='darkblue', linewidth=0.5,
                           hatch='///', label=f"Col-Del (n={len(data['Col']['deletions'])})")

                if data['Ler']['insertions']:
                    ax.hist(data['Ler']['insertions'], bins=bins,
                           color='#e74c3c', alpha=0.5, edgecolor='darkred', linewidth=0.5,
                           label=f"Ler-Ins (n={len(data['Ler']['insertions'])})")

                if data['Ler']['deletions']:
                    ax.hist(data['Ler']['deletions'], bins=bins,
                           color='#f1948a', alpha=0.5, edgecolor='darkred', linewidth=0.5,
                           hatch='///', label=f"Ler-Del (n={len(data['Ler']['deletions'])})")

                # Add 178bp, 356bp, 534bp, 712bp, 890bp reference lines
                for i in range(1, 6):
                    ref_line = i * 178
                    if ref_line <= 1000:
                        ax.axvline(ref_line, color='gray', linestyle='--',
                                 alpha=0.4, linewidth=0.8)
                        ax.text(ref_line, ax.get_ylim()[1] * 0.95, f'{i}×178',
                               fontsize=7, ha='center', va='top', color='gray')

                ax.legend(fontsize=6, loc='upper right', framealpha=0.9)
            else:
                ax.text(0.5, 0.5, 'No indels 50-1000bp', ha='center', va='center',
                       transform=ax.transAxes, fontsize=10, style='italic', color='gray')

            # STANDARDIZED AXES (zoomed)
            ax.set_xlim(0, 1000)
            ax.set_ylim(0, global_max_count_zoom)

            ax.set_title(f'Chr{chrom} - {region}', fontsize=10, fontweight='bold')
            if col == 0:
                ax.set_ylabel('Count', fontsize=9)
            if row == 1:
                ax.set_xlabel('Indel size (bp)', fontsize=9)
            ax.tick_params(labelsize=8)
            ax.grid(True, alpha=0.2, axis='y')

        plt.tight_layout()
        plot_filename = os.path.join(output_dir, f'indel_histogram_{mode}_combined_zoom.png')
        plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Generated zoomed plot (0-1000bp): {plot_filename}")


def create_comparative_plots(results_df, output_dir):
    """
    Create comprehensive comparative plots for indel analysis (>50bp).
    Uses clean styling inspired by sv-satellite-analyzer.
    """
    sns.set_style("white")
    sns.set_palette(["#3498db", "#e74c3c"])
    
    fig = plt.figure(figsize=(18, 10))
    
    # Indel density by region
    ax1 = plt.subplot(2, 4, 1)
    region_data = results_df.groupby('region')['indel_density'].agg(['mean', 'std']).reset_index()
    ax1.bar(region_data['region'], region_data['mean'], yerr=region_data['std'], capsize=5, alpha=0.7)
    ax1.set_title('Average Indel Density: ARMS vs CEN', fontweight='bold')
    ax1.set_ylabel('Indels per read (>50bp)')
    ax1.grid(True, alpha=0.3, axis='y')

    # Average indel size
    ax2 = plt.subplot(2, 4, 2)
    region_size_data = results_df.groupby('region')['avg_indel_size'].agg(['mean', 'std']).reset_index()
    ax2.bar(region_size_data['region'], region_size_data['mean'], yerr=region_size_data['std'], capsize=5, alpha=0.7)
    ax2.set_title('Average Indel Size: ARMS vs CEN', fontweight='bold')
    ax2.set_ylabel('Average indel size (bp)')
    ax2.grid(True, alpha=0.3, axis='y')

    # Background comparison
    ax3 = plt.subplot(2, 4, 3)
    background_data = results_df.groupby('background')['indel_density'].agg(['mean', 'std']).reset_index()
    ax3.bar(background_data['background'], background_data['mean'], yerr=background_data['std'], capsize=5, alpha=0.7)
    ax3.set_title('Average Indel Density: Col vs Ler', fontweight='bold')
    ax3.set_ylabel('Indels per read (>50bp)')
    ax3.grid(True, alpha=0.3, axis='y')

    # Chromosome comparison
    ax4 = plt.subplot(2, 4, 4)
    chr_data = results_df.groupby('chromosome')['indel_density'].agg(['mean', 'std']).reset_index()
    ax4.bar(chr_data['chromosome'], chr_data['mean'], yerr=chr_data['std'], capsize=5, alpha=0.7)
    ax4.set_title('Average Indel Density by Chromosome', fontweight='bold')
    ax4.set_ylabel('Indels per read (>50bp)')
    ax4.set_xlabel('Chromosome')
    ax4.grid(True, alpha=0.3, axis='y')
    
    # Heatmap of indel density by region and chromosome
    ax5 = plt.subplot(2, 4, 5)
    pivot_data = results_df.pivot_table(values='indel_density', index='region', columns='chromosome', aggfunc='mean')
    sns.heatmap(pivot_data, annot=True, fmt='.4f', cmap='YlOrRd', ax=ax5, cbar_kws={'label': 'Indels per read'})
    ax5.set_title('Indel Density Heatmap\n(Region vs Chromosome)', fontweight='bold')
    
    # Insertion vs Deletion proportions
    ax6 = plt.subplot(2, 4, 6)
    total_ins = results_df['insertions'].sum()
    total_del = results_df['deletions'].sum()
    ax6.pie([total_ins, total_del], labels=['Insertions', 'Deletions'], autopct='%1.1f%%',
           startangle=90, colors=['#3498db', '#e74c3c'])
    ax6.set_title('Overall Insertion vs Deletion (>50bp)', fontweight='bold')
    
    # Indel density by mapping mode
    ax7 = plt.subplot(2, 4, 7)
    if 'mapping_mode' in results_df.columns:
        mode_data = results_df.groupby('mapping_mode')['indel_density'].agg(['mean', 'std']).reset_index()
        ax7.bar(mode_data['mapping_mode'], mode_data['mean'], yerr=mode_data['std'], capsize=5, alpha=0.7)
        ax7.set_title('Indel Density by Mapping Mode', fontweight='bold')
        ax7.set_ylabel('Indels per read (>50bp)')
        ax7.set_xlabel('Mapping Mode')
        ax7.grid(True, alpha=0.3, axis='y')
        plt.setp(ax7.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    # Total indels by haplotype
    ax8 = plt.subplot(2, 4, 8)
    hap_indels = results_df.groupby('haplotype')['indels'].sum().sort_values(ascending=False).head(10)
    ax8.barh(range(len(hap_indels)), hap_indels.values, alpha=0.7)
    ax8.set_yticks(range(len(hap_indels)))
    ax8.set_yticklabels(hap_indels.index)
    ax8.set_xlabel('Total Indels')
    ax8.set_title('Top 10 Haplotypes by Indel Count', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'comparative_indel_analysis.png'), dpi=300, bbox_inches='tight')
    plt.close()

    print("Generated comparative analysis plots.")


def main():
    if len(sys.argv) != 3:
        print("Usage: python3 12.2-indel_analysis.py <sequences_dir> <output_results.tsv>")
        sys.exit(1)
    
    sequences_dir = sys.argv[1]
    output_results = sys.argv[2]
    
    mappings_dir = os.path.join(sequences_dir, "mappings")
    
    if not os.path.exists(mappings_dir):
        print(f"Error: Mappings directory not found: {mappings_dir}")
        with open(output_results, 'w') as f:
            f.write("haplotype\tmapping_mode\ttotal_reads\tmapped_reads\tindels\tinsertions\tdeletions\tavg_indel_size\tindel_density\treads_with_indels\tbackground\tregion\tchromosome\n")
        sys.exit(0)
    
    print("=== Starting comprehensive indel analysis ===")
    results = []
    results_with_lengths = []  # Keep indel length data for plotting
    
    for bam_file in Path(mappings_dir).glob("*.bam"):
        if bam_file.is_file():
            basename = bam_file.stem
            parts = basename.split('_', 1)
            haplotype = parts[0]
            mapping_mode = parts[1] if len(parts) > 1 else "unknown"
            hap_info = parse_haplotype(haplotype)
            output_prefix = os.path.join(mappings_dir, basename)
            
            print(f"\nProcessing {bam_file.name} (haplotype: {haplotype}, mode: {mapping_mode})")
            stats = analyze_bam_file(str(bam_file), output_prefix)
            
            if stats and hap_info:
                # Store for TSV output (without length arrays)
                result_entry = {**hap_info, 'haplotype': haplotype, 'mapping_mode': mapping_mode}
                for key in ['total_reads', 'mapped_reads', 'indels', 'insertions', 'deletions', 
                           'avg_indel_size', 'indel_density', 'reads_with_indels']:
                    result_entry[key] = stats[key]
                results.append(result_entry)
                
                # Store with length data for plotting
                results_with_lengths.append({**hap_info, **stats, 'haplotype': haplotype, 'mapping_mode': mapping_mode})
    
    # Save TSV results
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_results, sep='\t', index=False)
    
    # Generate plots
    if results_with_lengths:
        create_combined_plots_by_mapping_mode(results_with_lengths, mappings_dir)
        
        if not results_df.empty:
            create_comparative_plots(results_df, mappings_dir)
    
    print(f"\n=== Comprehensive indel analysis complete ===")
    print(f"Summary results written to: {output_results}")


if __name__ == "__main__":
    main()
