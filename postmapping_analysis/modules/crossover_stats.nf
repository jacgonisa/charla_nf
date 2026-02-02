process CROSSOVER_STATS {
    tag "crossover_stats_${stage}_${threshold}"
    publishDir "${params.outdir}/crossover_stats_${stage}", mode: 'symlink'

    input:
    path cowidths_files              // All files from plotting_crossover (*.cowidths and *.tsv) or empty for profiles stage
    path curated_profiles            // curated_hybrid_profiles TSV file
    val threshold                    // Threshold value
    val sample_name                  // Sample name for output files
    val stage                        // "profiles" or "mapping"
    val parent1_name                 // Name of parent 1
    val parent2_name                 // Name of parent 2

    output:
    path "*_hybrid_type_stats.tsv", emit: hybrid_type_stats
    path "*_hybrid_types_barplot.svg", emit: hybrid_types_plot
    path "*_cowidth_stats.tsv", emit: cowidth_stats, optional: true
    path "*_cowidth_histogram_*.svg", emit: cowidth_histograms, optional: true
    path "*_hybrid_type_by_parent_stats.tsv", emit: hybrid_type_by_parent_stats, optional: true
    path "*_hybrid_types_by_parent_barplot.svg", emit: hybrid_types_by_parent_plot, optional: true

script:
// Find the cowidths file from the input files
def cowidths_list = cowidths_files instanceof List ? cowidths_files : (cowidths_files ? [cowidths_files] : [])
def cowidths_file = cowidths_list.find { it.name.endsWith('.cowidths') }
def cowidths_arg = cowidths_file ? "--cowidths_file ${cowidths_file}" : "--cowidths_file NA"
"""
Rscript ${projectDir}/scripts/14-crossover_statistics.R \\
    ${cowidths_arg} \\
    --curated_profiles ${curated_profiles} \\
    --threshold ${threshold} \\
    --sample_name ${sample_name} \\
    --analysis_stage ${stage} \\
    --parent1_name ${parent1_name} \\
    --parent2_name ${parent2_name} \\
    --output_dir .
"""
}
