process PLOT_CONFIDENT_READS {
    tag "plot_confident_${sample_name}_thr${threshold}"
    publishDir "${params.outdir}/curated_confident_reads", mode: 'copy'

    input:
    path confident_bed_sco          // Confident SCO BED file
    path confident_bed_nco          // Confident NCO BED file
    path confident_bed_gcco         // Confident GC-CO BED file
    path confident_read_ids         // Confident read IDs text file
    path parent1_paf_chr            // Parent 1 CHR PAF file
    path parent2_paf_chr            // Parent 2 CHR PAF file
    path parent1_genome_fai         // Parent 1 genome FASTA index
    path parent2_genome_fai         // Parent 2 genome FASTA index
    path cowidths_file              // Cowidths file
    path curated_profiles           // Curated hybrid profiles TSV
    val threshold                   // Threshold value
    val sample_name                 // Sample name
    val parent1_name                // Name of parent 1
    val parent2_name                // Name of parent 2
    val centromere_bed              // Centromere BED file path (optional)

    output:
    path "confident_reads_plots/*_confident_crossover_plot_*.svg", emit: coverage_plots
    path "confident_reads_plots/*_confident_circos_crossovers.pdf", emit: circos_plot
    path "confident_reads_plots/*_confident_karyoplot_crossovers.pdf", emit: karyoplot

    script:
    def cen_bed_arg = centromere_bed != "NO_FILE" ? "--centromere_bed ${centromere_bed}" : "--centromere_bed NA"
    """
    Rscript ${projectDir}/scripts/18-plot_confident_reads_comprehensive.R \\
        --confident_bed_sco ${confident_bed_sco} \\
        --confident_bed_nco ${confident_bed_nco} \\
        --confident_bed_gcco ${confident_bed_gcco} \\
        --confident_read_ids ${confident_read_ids} \\
        --parent1_paf_chr ${parent1_paf_chr} \\
        --parent2_paf_chr ${parent2_paf_chr} \\
        --parent1_genome_fai ${parent1_genome_fai} \\
        --parent2_genome_fai ${parent2_genome_fai} \\
        --cowidths_file ${cowidths_file} \\
        --curated_profiles ${curated_profiles} \\
        --sample_name ${sample_name} \\
        --parent1_name ${parent1_name} \\
        --parent2_name ${parent2_name} \\
        --threshold ${threshold} \\
        ${cen_bed_arg} \\
        --output_dir .
    """
}
