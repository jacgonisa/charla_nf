nextflow.enable.dsl=2

process ANALYZE_MAFFT_ALIGNMENT {
    tag "analyze_mafft"
    label 'process_low'

    publishDir "${params.outdir}/large_indel_analysis/${params.sample_id}_large_indels/msa_results", mode: 'symlink'

    input:
    tuple val(meta), path(alignment)
    path(satellite_analysis_tsv)

    output:
    path("${params.sample_id}_satellite_analysis_updated.tsv"), emit: updated_analysis
    path("msa_visualization.png"), emit: msa_plot, optional: true
    path("versions.yml"), emit: versions, optional: true

    script:
    """
    # Run alignment analysis
    python ${projectDir}/scripts/13.4-analyze_mafft_alignment.py \\
        ${alignment} \\
        ${satellite_analysis_tsv}

    # Copy updated file with new name
    cp ${satellite_analysis_tsv} ${params.sample_id}_satellite_analysis_updated.tsv

    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
