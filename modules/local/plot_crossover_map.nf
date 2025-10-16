process PLOT_CROSSOVER_MAP {
    tag "crossover_plot"
    publishDir "${params.outdir}/crossover_plots", mode: 'symlink'
    
    input:
    path col_bed                    // Col-0 BED file
    path ler_bed                    // Ler-0 BED file  
    path col_arms_paf               // ARMS PAF file against Col reference
    path col_cen_paf                // CEN PAF file against Col reference
    path ler_arms_paf               // ARMS PAF file against Ler reference
    path ler_cen_paf                // CEN PAF file against Ler reference
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

    # Concatenate Col PAF files
    if [[ -s "${col_cen_paf}" ]]; then
        cat ${col_arms_paf} ${col_cen_paf} | grep -v "^#" > summary_mm_sr_ARMSandCEN_againstCol.paf
    else
        echo "Warning: CEN PAF file for Col is empty or missing, using only ARMS data" | tee -a filtering_stats.txt
        cat ${col_arms_paf} | grep -v "^#" > summary_mm_sr_ARMSandCEN_againstCol.paf
    fi

    # Concatenate Ler PAF files
    if [[ -s "${ler_cen_paf}" ]]; then
        cat ${ler_arms_paf} ${ler_cen_paf} | grep -v "^#" > summary_mm_sr_ARMSandCEN_againstLer.paf
    else
        echo "Warning: CEN PAF file for Ler is empty or missing, using only ARMS data" | tee -a filtering_stats.txt
        cat ${ler_arms_paf} | grep -v "^#" > summary_mm_sr_ARMSandCEN_againstLer.paf
    fi
    
    # Count original alignments
    col_original=\$(wc -l < summary_mm_sr_ARMSandCEN_againstCol.paf)
    ler_original=\$(wc -l < summary_mm_sr_ARMSandCEN_againstLer.paf)
    
    echo "Original PAF alignments:" | tee -a filtering_stats.txt
    echo "  Col-0: \${col_original}" | tee -a filtering_stats.txt
    echo "  Ler-0: \${ler_original}" | tee -a filtering_stats.txt
    echo "" | tee -a filtering_stats.txt
    
    echo "=== Filtering parameters ===" | tee -a filtering_stats.txt
    echo "  Cluster tolerance: ${cluster_tolerance} bp" | tee -a filtering_stats.txt
    echo "  Max reads per cluster: ${max_cluster_count}" | tee -a filtering_stats.txt
    echo "" | tee -a filtering_stats.txt
    
    
   
    echo "=== Filtering Col-0 PAF ===" | tee -a filtering_stats.txt
    awk -f ${projectDir}/scripts/13-filter_high_coverage_events.awk \
        -v tol=${cluster_tolerance} \
        -v max=${max_cluster_count} \
        -v prefix="summary_mm_sr_ARMSandCEN_againstCol" \
        summary_mm_sr_ARMSandCEN_againstCol.paf 2>> filtering_stats.txt

echo "=== Filtering Ler-0 PAF ===" | tee -a filtering_stats.txt
    awk -f ${projectDir}/scripts/13-filter_high_coverage_events.awk \
        -v tol=${cluster_tolerance} \
        -v max=${max_cluster_count} \
        -v prefix="summary_mm_sr_ARMSandCEN_againstLer" \
        summary_mm_sr_ARMSandCEN_againstLer.paf 2>> filtering_stats.txt

    echo "" | tee -a filtering_stats.txt
    echo "=== Filtering complete ===" | tee -a filtering_stats.txt
  


echo "PAF files concatenated successfully"
echo "Starting R script..."

Rscript ${projectDir}/scripts/13-plot_crossover_map.R \
    --col_bed ${col_bed} \
    --ler_bed ${ler_bed} \
    --col_paf summary_mm_sr_ARMSandCEN_againstCol.filtered.paf \
    --ler_paf summary_mm_sr_ARMSandCEN_againstLer.filtered.paf \
    --sample_name ${sample_name}






"""
}
