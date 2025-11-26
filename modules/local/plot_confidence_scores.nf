nextflow.enable.dsl=2

process PLOT_CONFIDENCE_SCORES {
    tag "confidence_plots_thr${threshold}"
    label 'plotting'
    cpus params.plot_cpus ?: 1
    memory params.plot_memory ?: '4 GB'
    time params.plot_time ?: '30m'

    publishDir "${params.outdir}/confidence_scores", mode: 'symlink'

    input:
    path confidence_scores      // TSV file with confidence scores
    val threshold              // Threshold value

    output:
    path "threshold${threshold}/*_histogram.png", emit: histogram
    path "threshold${threshold}/*_component_distributions.png", emit: component_dist
    path "threshold${threshold}/*_correlation_matrix.png", emit: correlation
    path "threshold${threshold}/*_component_scatter.png", emit: scatter
    path "threshold${threshold}/*_region_comparison.png", emit: region_comp, optional: true

    script:
    """
    threshold_dir="threshold${threshold}"
    mkdir -p \${threshold_dir}

    echo "Generating confidence score visualizations..."

    python ${projectDir}/scripts/15-plot_confidence_distribution.py \\
        --input ${confidence_scores} \\
        --output-prefix \${threshold_dir}/confidence_scores

    echo "Visualization complete!"
    ls -lh \${threshold_dir}/*.png
    """
}
