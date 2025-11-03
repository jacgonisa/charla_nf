process PLOT_CROSSOVER_MAP {
    tag "crossover_plot"
    publishDir "${params.outdir}/crossover_plots", mode: 'symlink'

    input:
    val parent1_bed                 // Parent 1 BED file (or "NA" for centromere-unaware mode)
    val parent2_bed                 // Parent 2 BED file (or "NA" for centromere-unaware mode)
    path parent1_paf                // Parent 1 PAF file(s) - ARMS+CEN or CHR
    path parent1_cen_paf            // Parent 1 CEN PAF file (or empty for centromere-unaware)
    path parent2_paf                // Parent 2 PAF file(s) - ARMS+CEN or CHR
    path parent2_cen_paf            // Parent 2 CEN PAF file (or empty for centromere-unaware)
    val threshold                   // Threshold value for naming
    val sample_name                 // Sample name for output files
    
    output:
    path "*.svg", emit: crossover_plots
    path "*.txt", emit: coverage_data
   // path "*.log", emit: log_files
    path "summary_mm_sr_ARMSandCEN_againstCol.paf", emit: col_merged_paf
    path "summary_mm_sr_ARMSandCEN_againstLer.paf", emit: ler_merged_paf
    path "summary_mm_sr_ARMSandCEN_againstCol.filtered.paf", emit: col_filtered_paf
    path "summary_mm_sr_ARMSandCEN_againstLer.filtered.paf", emit: ler_filtered_paf
    path "high_coverage_regions_*againstCol*.bed", emit: col_high_coverage_bed
    path "high_coverage_regions_*againstLer*.bed", emit: ler_high_coverage_bed
    path "filtering_stats.txt", emit: filtering_stats

script:
def cluster_tolerance = params.cluster_tolerance ?: 50
def max_cluster_count = params.max_cluster_count ?: 3

"""

# First, concatenate PAF files (handle missing CEN files gracefully)
echo "=== Concatenating PAF files ===" | tee filtering_stats.txt

# Concatenate Parent 1 PAF files
if [[ -s "${parent1_cen_paf}" ]]; then
    cat ${parent1_paf} ${parent1_cen_paf} | grep -v "^#" > summary_mm_sr_ARMSandCEN_against${params.parent1_name}.paf
else
    echo "Warning: CEN PAF file for ${params.parent1_name} is empty or missing, using only CHR/ARMS data" | tee -a filtering_stats.txt
    cat ${parent1_paf} | grep -v "^#" > summary_mm_sr_ARMSandCEN_against${params.parent1_name}.paf
fi

# Concatenate Parent 2 PAF files
if [[ -s "${parent2_cen_paf}" ]]; then
    cat ${parent2_paf} ${parent2_cen_paf} | grep -v "^#" > summary_mm_sr_ARMSandCEN_against${params.parent2_name}.paf
else
    echo "Warning: CEN PAF file for ${params.parent2_name} is empty or missing, using only CHR/ARMS data" | tee -a filtering_stats.txt
    cat ${parent2_paf} | grep -v "^#" > summary_mm_sr_ARMSandCEN_against${params.parent2_name}.paf
fi
    
    # Count original alignments
    parent1_original=\$(wc -l < summary_mm_sr_ARMSandCEN_against${params.parent1_name}.paf)
    parent2_original=\$(wc -l < summary_mm_sr_ARMSandCEN_against${params.parent2_name}.paf)

    echo "Original PAF alignments:" | tee -a filtering_stats.txt
    echo "  ${params.parent1_name}: \${parent1_original}" | tee -a filtering_stats.txt
    echo "  ${params.parent2_name}: \${parent2_original}" | tee -a filtering_stats.txt
    echo "" | tee -a filtering_stats.txt
    
    echo "=== Filtering parameters ===" | tee -a filtering_stats.txt
    echo "  Cluster tolerance: ${cluster_tolerance} bp" | tee -a filtering_stats.txt
    echo "  Max reads per cluster: ${max_cluster_count}" | tee -a filtering_stats.txt
    echo "" | tee -a filtering_stats.txt
    
    
   
    echo "=== Filtering ${params.parent1_name} PAF ===" | tee -a filtering_stats.txt
    awk -f ${projectDir}/scripts/13-filter_high_coverage_events.awk \
        -v tol=${cluster_tolerance} \
        -v max=${max_cluster_count} \
        -v prefix="summary_mm_sr_ARMSandCEN_against${params.parent1_name}" \
        summary_mm_sr_ARMSandCEN_against${params.parent1_name}.paf 2>> filtering_stats.txt

echo "=== Filtering ${params.parent2_name} PAF ===" | tee -a filtering_stats.txt
    awk -f ${projectDir}/scripts/13-filter_high_coverage_events.awk \
        -v tol=${cluster_tolerance} \
        -v max=${max_cluster_count} \
        -v prefix="summary_mm_sr_ARMSandCEN_against${params.parent2_name}" \
        summary_mm_sr_ARMSandCEN_against${params.parent2_name}.paf 2>> filtering_stats.txt

    echo "" | tee -a filtering_stats.txt
    echo "=== Filtering complete ===" | tee -a filtering_stats.txt
  


echo "PAF files concatenated successfully"
echo "Starting R script..."

Rscript ${projectDir}/scripts/13-plot_crossover_map.R \
    --col_bed ${parent1_bed} \
    --ler_bed ${parent2_bed} \
    --col_paf summary_mm_sr_ARMSandCEN_against${params.parent1_name}.filtered.paf \
    --ler_paf summary_mm_sr_ARMSandCEN_against${params.parent2_name}.filtered.paf \
    --sample_name ${sample_name} \
    --parent1_name ${params.parent1_name} \
    --parent2_name ${params.parent2_name}






"""
}
