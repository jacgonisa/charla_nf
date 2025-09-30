#!/usr/bin/env python3
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
    Analyze a BAM file for indels and generate both statistics and plots.
    Returns dictionary with stats for the summary table.
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
    indel_lengths = []
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
            
        query_pos = 0  # position in read
        ref_pos = read.reference_start  # position on reference
        
        for op, length in read.cigartuples:
            if op == 0:  # Match or mismatch (M)
                query_pos += length
                ref_pos += length
            elif op == 1 and length > 10:  # Insertion > 10bp
                indel_lengths.append(length)
                reads_with_indels.add(read.query_name)
                indel_details.append({
                    "Read_ID": read.query_name,
                    "Type": "Insertion",
                    "Size": length,
                    "Read_Pos": query_pos,
                    "Ref_Pos": ref_pos
                })
                query_pos += length
            elif op == 2 and length > 10:  # Deletion > 10bp
                indel_lengths.append(length)
                reads_with_indels.add(read.query_name)
                indel_details.append({
                    "Read_ID": read.query_name,
                    "Type": "Deletion",
                    "Size": length,
                    "Read_Pos": query_pos,
                    "Ref_Pos": ref_pos
                })
                ref_pos += length
            elif op in [3, 4, 5]:  # Skip (N), soft clip (S), hard clip (H)
                if op in [4, 5]:  # clipping affects read position
                    query_pos += length
                if op == 3:  # reference skip
                    ref_pos += length
    
    bamfile.close()
    
    # Generate detailed TSV output
    tsv_output = f"{output_prefix}_indel_positions.tsv"
    with open(tsv_output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["Read_ID", "Type", "Size", "Read_Pos", "Ref_Pos"], delimiter="\t")
        writer.writeheader()
        for row in indel_details:
            writer.writerow(row)
    
    # Generate histogram plot if we have indels
    if indel_lengths:
        plt.figure(figsize=(12, 8))
        plt.hist(indel_lengths, bins=50, color="skyblue", edgecolor="black", alpha=0.7)
        plt.title(f"Histogram of Indel Sizes (>10 bp) - {output_prefix}")
        plt.xlabel("Indel size (bp)")
        plt.ylabel("Count")
        
        # Add reference lines at 178, 356, 534... (centromeric repeat units)
        max_length = max(indel_lengths) if indel_lengths else 0
        reference_lines_added = False
        for i in range(1, max_length // 178 + 2):
            ref_line = i * 178
            if ref_line <= max_length + 100:  # Add a bit of margin
                plt.axvline(ref_line, color="red", linestyle="--", alpha=0.7, 
                           label=f"178bp multiples" if not reference_lines_added else "")
                reference_lines_added = True
        
        if reference_lines_added:
            plt.legend()
        
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        plot_output = f"{output_prefix}_indel_histogram.png"
        plt.savefig(plot_output, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Generated plot: {plot_output}")
    
    # Calculate statistics
    num_indels = len(indel_lengths)
    avg_indel_size = sum(indel_lengths) / num_indels if num_indels > 0 else 0
    indel_density = num_indels / mapped_reads if mapped_reads > 0 else 0
    reads_with_indels_count = len(reads_with_indels)
    
    print(f"Analysis complete for {output_prefix}:")
    print(f"  Total reads: {total_reads}")
    print(f"  Mapped reads: {mapped_reads}")
    print(f"  Indels (>10bp): {num_indels}")
    print(f"  Reads with indels: {reads_with_indels_count}")
    print(f"  Average indel size: {avg_indel_size:.2f} bp")
    print(f"  Indel density: {indel_density:.4f} indels/read")
    print(f"  Detailed results: {tsv_output}")
    
    return {
        'total_reads': total_reads,
        'mapped_reads': mapped_reads,
        'indels': num_indels,
        'avg_indel_size': round(avg_indel_size, 2),
        'indel_density': round(indel_density, 4),
        'reads_with_indels': reads_with_indels_count
    }

def create_comparative_plots(results_df, output_dir):
    """
    Create comprehensive comparative plots for indel analysis.
    """
    plt.style.use('default')
    
    # Set up the plotting style
    sns.set_palette("husl")
    
    # 1. ARMS vs CEN comparison
    plt.figure(figsize=(15, 10))
    
    # Subplot 1: Indel density by region
    plt.subplot(2, 3, 1)
    region_data = results_df.groupby('region')['indel_density'].agg(['mean', 'std']).reset_index()
    plt.bar(region_data['region'], region_data['mean'], yerr=region_data['std'], capsize=5, alpha=0.7)
    plt.title('Average Indel Density: ARMS vs CEN')
    plt.ylabel('Indels per read')
    plt.xticks(rotation=45)
    
    # Subplot 2: Average indel size by region
    plt.subplot(2, 3, 2)
    region_size_data = results_df.groupby('region')['avg_indel_size'].agg(['mean', 'std']).reset_index()
    plt.bar(region_size_data['region'], region_size_data['mean'], yerr=region_size_data['std'], capsize=5, alpha=0.7)
    plt.title('Average Indel Size: ARMS vs CEN')
    plt.ylabel('Average indel size (bp)')
    plt.xticks(rotation=45)
    
    # Subplot 3: Col vs Ler comparison
    plt.subplot(2, 3, 3)
    background_data = results_df.groupby('background')['indel_density'].agg(['mean', 'std']).reset_index()
    plt.bar(background_data['background'], background_data['mean'], yerr=background_data['std'], capsize=5, alpha=0.7)
    plt.title('Average Indel Density: Col vs Ler')
    plt.ylabel('Indels per read')
    plt.xticks(rotation=45)
    
    # Subplot 4: Chromosome comparison
    plt.subplot(2, 3, 4)
    chr_data = results_df.groupby('chromosome')['indel_density'].agg(['mean', 'std']).reset_index()
    plt.bar(chr_data['chromosome'], chr_data['mean'], yerr=chr_data['std'], capsize=5, alpha=0.7)
    plt.title('Average Indel Density by Chromosome')
    plt.ylabel('Indels per read')
    plt.xlabel('Chromosome')
    
    # Subplot 5: Heatmap of indel density
    plt.subplot(2, 3, 5)
    pivot_data = results_df.pivot_table(values='indel_density', 
                                       index='region', 
                                       columns='chromosome', 
                                       aggfunc='mean')
    sns.heatmap(pivot_data, annot=True, fmt='.4f', cmap='YlOrRd')
    plt.title('Indel Density Heatmap\n(Region vs Chromosome)')
    
    # Subplot 6: Box plot of indel densities
    plt.subplot(2, 3, 6)
    region_chr = results_df['region'] + '_Chr' + results_df['chromosome']
    plt.boxplot([results_df[region_chr == cat]['indel_density'].values 
                for cat in sorted(region_chr.unique())], 
                labels=sorted(region_chr.unique()))
    plt.title('Indel Density Distribution')
    plt.ylabel('Indels per read')
    plt.xticks(rotation=90)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'comparative_indel_analysis.png'), 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Detailed statistical comparison
    comparison_stats = []
    
    # ARMS vs CEN
    arms_data = results_df[results_df['region'] == 'ARMS']['indel_density']
    cen_data = results_df[results_df['region'] == 'CEN']['indel_density']
    
    comparison_stats.append({
        'comparison': 'ARMS_vs_CEN',
        'group1': 'ARMS',
        'group2': 'CEN',
        'group1_mean': arms_data.mean(),
        'group1_std': arms_data.std(),
        'group1_count': len(arms_data),
        'group2_mean': cen_data.mean(),
        'group2_std': cen_data.std(),
        'group2_count': len(cen_data),
        'fold_change': cen_data.mean() / arms_data.mean() if arms_data.mean() > 0 else float('inf')
    })
    
    # Col vs Ler
    col_data = results_df[results_df['background'] == 'Col']['indel_density']
    ler_data = results_df[results_df['background'] == 'Ler']['indel_density']
    
    comparison_stats.append({
        'comparison': 'Col_vs_Ler',
        'group1': 'Col',
        'group2': 'Ler', 
        'group1_mean': col_data.mean(),
        'group1_std': col_data.std(),
        'group1_count': len(col_data),
        'group2_mean': ler_data.mean(),
        'group2_std': ler_data.std(),
        'group2_count': len(ler_data),
        'fold_change': ler_data.mean() / col_data.mean() if col_data.mean() > 0 else float('inf')
    })
    
    # Save comparison statistics
    comparison_df = pd.DataFrame(comparison_stats)
    comparison_df.to_csv(os.path.join(output_dir, 'statistical_comparisons.tsv'), 
                        sep='\t', index=False)
    
    print("Generated comparative analysis plots and statistics:")
    print(f"  - Comparative plots: {os.path.join(output_dir, 'comparative_indel_analysis.png')}")
    print(f"  - Statistical comparisons: {os.path.join(output_dir, 'statistical_comparisons.tsv')}")

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 12.2-indel_analysis.py <sequences_dir> <output_results.tsv>")
        sys.exit(1)
    
    sequences_dir = sys.argv[1]
    output_results = sys.argv[2]
    
    mappings_dir = os.path.join(sequences_dir, "mappings")
    
    if not os.path.exists(mappings_dir):
        print(f"Error: Mappings directory not found: {mappings_dir}")
        # Create empty results file
        with open(output_results, 'w') as f:
            f.write("haplotype\tmapping_mode\ttotal_reads\tmapped_reads\tindels\tavg_indel_size\tindel_density\treads_with_indels\tbackground\tregion\tchromosome\n")
        sys.exit(0)
    
    print("=== Starting comprehensive indel analysis ===")
    
    # Initialize results
    results = []
    
    # Process each BAM file
    for bam_file in Path(mappings_dir).glob("*.bam"):
        if bam_file.is_file():
            # Extract haplotype and mapping mode from filename
            basename = bam_file.stem  # removes .bam
            parts = basename.split('_', 1)  # Split on first underscore
            if len(parts) >= 2:
                haplotype = parts[0]
                mapping_mode = parts[1]
            else:
                haplotype = basename
                mapping_mode = "unknown"
            
            # Parse haplotype components
            hap_info = parse_haplotype(haplotype)
            
            output_prefix = os.path.join(mappings_dir, basename)
            
            print(f"\nProcessing {bam_file.name} (haplotype: {haplotype}, mode: {mapping_mode})")
            
            # Analyze this BAM file
            stats = analyze_bam_file(str(bam_file), output_prefix)
            
            if stats and hap_info:
                results.append({
                    'haplotype': haplotype,
                    'mapping_mode': mapping_mode,
                    'background': hap_info['background'],
                    'region': hap_info['region'],
                    'chromosome': hap_info['chromosome'],
                    **stats
                })
            else:
                # Add empty result for failed analysis
                if hap_info:
                    results.append({
                        'haplotype': haplotype,
                        'mapping_mode': mapping_mode,
                        'background': hap_info['background'],
                        'region': hap_info['region'],
                        'chromosome': hap_info['chromosome'],
                        'total_reads': 0,
                        'mapped_reads': 0,
                        'indels': 0,
                        'avg_indel_size': 0,
                        'indel_density': 0,
                        'reads_with_indels': 0
                    })
    
    # Convert to DataFrame for analysis
    results_df = pd.DataFrame(results)
    
    # Write detailed summary results
    with open(output_results, 'w') as f:
        f.write("haplotype\tmapping_mode\ttotal_reads\tmapped_reads\tindels\tavg_indel_size\tindel_density\treads_with_indels\tbackground\tregion\tchromosome\n")
        for result in results:
            f.write(f"{result['haplotype']}\t{result['mapping_mode']}\t{result['total_reads']}\t"
                   f"{result['mapped_reads']}\t{result['indels']}\t{result['avg_indel_size']}\t"
                   f"{result['indel_density']}\t{result['reads_with_indels']}\t{result['background']}\t"
                   f"{result['region']}\t{result['chromosome']}\n")
    
    # Create comparative analysis plots
    if not results_df.empty:
        create_comparative_plots(results_df, mappings_dir)
    
    print(f"\n=== Comprehensive indel analysis complete ===")
    print(f"Summary results written to: {output_results}")
    print(f"Individual plots and detailed results in: {mappings_dir}")
    print(f"Comparative analysis plots generated in: {mappings_dir}")

if __name__ == "__main__":
    main()
