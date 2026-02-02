process CURATE_CONFIDENT_READS {
    tag "curate_confident_${sample_name}_thr${threshold}"
    publishDir "${params.outdir}/curated_confident_reads", mode: 'copy'

    input:
    path cowidths_file               // Cowidths file from MAPPING_ANALYSIS
    path bed_file                    // BED file from MAPPING_ANALYSIS
    path curated_profiles            // Curated hybrid profiles TSV file
    val threshold                    // Threshold value
    val sample_name                  // Sample name for output files
    val parent1_name                 // Name of parent 1
    val parent2_name                 // Name of parent 2

    output:
    path "curated_confident_reads/*_confident_*.bed", emit: confident_bed_files
    path "curated_confident_reads/*_nonconfident_*.bed", emit: nonconfident_bed_files, optional: true
    path "curated_confident_reads/*_confident_read_ids.txt", emit: confident_read_ids
    path "curated_confident_reads/*_confident_stats.tsv", emit: confident_stats
    path "curated_confident_reads/*_confident_*.svg", emit: confident_plots

    script:
    """
    Rscript ${projectDir}/scripts/15-curate_confident_reads.R \\
        --cowidths_file ${cowidths_file} \\
        --bed_file ${bed_file} \\
        --curated_profiles ${curated_profiles} \\
        --threshold ${threshold} \\
        --sample_name ${sample_name} \\
        --parent1_name ${parent1_name} \\
        --parent2_name ${parent2_name} \\
        --output_dir .
    """
}
