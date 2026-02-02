nextflow.enable.dsl=2

process ANALYZE_INDEL_MONOMER_BOUNDARIES {
    tag "indel_monomer_analysis"
    label 'analysis'
    cpus params.analysis_cpus ?: 8
    memory params.analysis_memory ?: '16 GB'
    time params.analysis_time ?: '4h'

    publishDir "${params.outdir}/indel_monomer_analysis", mode: 'symlink'

    input:
    path indel_sequences        // FASTA file with indel sequences
    path indels_catalog         // TSV catalog with indel metadata
    path centromere_bed         // BED file with centromere coordinates (optional, for filtering)
    val threads

    output:
    path "indel_monomer_analysis/", emit: output_directory
    path "indel_monomer_analysis/indel_boundary_metaprofile.png", emit: metaprofile_plot, optional: true
    path "indel_monomer_analysis/inframe_analysis.png", emit: inframe_plot, optional: true
    path "indel_monomer_analysis/inframe_analysis_results.tsv", emit: inframe_results, optional: true
    path "indel_monomer_analysis/centromere_indel_coverage.png", emit: coverage_plot, optional: true
    path "indel_monomer_analysis/centromere_comparison.png", emit: comparison_plot, optional: true
    path "indel_monomer_analysis/monomer_variant_clustering.png", emit: clustering_plot, optional: true
    path "indel_monomer_analysis/monomer_variants.tsv", emit: monomer_variants, optional: true
    path "indel_monomer_analysis/indel_alignment_heatmap.png", emit: heatmap_plot, optional: true
    path "indel_monomer_analysis/*.paf", emit: alignment_files, optional: true

    script:
    bed_arg = centromere_bed.name != 'NO_FILE' ? centromere_bed : ''
    """
    # Create output directory
    mkdir -p indel_monomer_analysis

    # Check if input files exist and have content
    if [[ ! -s ${indel_sequences} ]]; then
        echo "Warning: Indel sequences file is empty or missing"
        touch indel_monomer_analysis/indel_boundary_metaprofile.png
        touch indel_monomer_analysis/inframe_analysis_results.tsv
        exit 0
    fi

    if [[ ! -s ${indels_catalog} ]]; then
        echo "Warning: Indels catalog file is empty or missing"
        touch indel_monomer_analysis/indel_boundary_metaprofile.png
        touch indel_monomer_analysis/inframe_analysis_results.tsv
        exit 0
    fi

    echo "=== Analyzing Indel Boundaries Relative to 178bp Monomer ==="

    # Run the analysis
    python ${projectDir}/scripts/15-analyze_indel_monomer_boundaries.py \\
        ${indel_sequences} \\
        ${indels_catalog} \\
        indel_monomer_analysis \\
        ${threads} \\
        ${bed_arg}

    echo "=== Indel monomer boundary analysis complete ==="

    # List generated files
    echo "Generated files:"
    ls -lh indel_monomer_analysis/
    """
}
