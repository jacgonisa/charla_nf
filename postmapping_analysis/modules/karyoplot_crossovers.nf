process KARYOPLOT_CROSSOVERS {
    tag "karyoplot_${sample_name}_${threshold}"
    publishDir "${params.outdir}/crossover_visualizations", mode: 'symlink'

    input:
    path bed_file                    // bed_of_best.tsv file
    path parent1_genome_fai          // Parent 1 genome FASTA index
    path parent2_genome_fai          // Parent 2 genome FASTA index
    path curated_profiles            // curated_hybrid_profiles TSV file
    val threshold                    // Threshold value
    val sample_name                  // Sample name for output files
    val parent1_name                 // Name of parent 1
    val parent2_name                 // Name of parent 2

    output:
    path "*_circos_crossovers.pdf", emit: circos_plot
    path "*_karyoplot_crossovers.pdf", emit: karyoplot

    script:
    """
    Rscript ${projectDir}/scripts/16-circos_and_karyoplot_crossovers.R \\
        --bed_file ${bed_file} \\
        --parent1_genome ${parent1_genome_fai} \\
        --parent2_genome ${parent2_genome_fai} \\
        --curated_profiles ${curated_profiles} \\
        --threshold ${threshold} \\
        --sample_name ${sample_name} \\
        --parent1_name ${parent1_name} \\
        --parent2_name ${parent2_name} \\
        --output_dir .
    """
}
