nextflow.enable.dsl=2

process VISUALIZE_SEGMENTS {
    tag "segment_visualization"
    label 'plotting'
    cpus params.plot_cpus ?: 2
    memory params.plot_memory ?: '16 GB'
    time params.plot_time ?: '4h'  // Increased for large datasets

    publishDir "${params.outdir}/segment_visualization", mode: 'symlink'

    input:
    path segments_bed                  // Segmentation BED file
    path confidence_scores, stageAs: 'confidence_scores.tsv'  // Optional confidence scores
    val max_reads

    output:
    path "*_read_structures.png", emit: structures_plot
    path "*_statistics.png", emit: statistics_plot
    path "*_summary.txt", emit: summary

    script:
    conf_arg = confidence_scores.name != 'NO_FILE' ? "--confidence-scores ${confidence_scores}" : ""
    """
    echo "=== Visualizing Read Segments ==="
    echo "NOTE: Generating visualization for all high-confidence reads (mixed types)"
    echo "For type-specific analysis, see SCO/NCO comparison plots"

    # All reads (high-confidence mixed)
    python ${projectDir}/scripts/21-visualize_read_segments.py \\
        --bed ${segments_bed} \\
        --output-prefix all_reads \\
        --max-reads ${max_reads} \\
        ${conf_arg}

    echo "=== Segment visualization complete ==="
    ls -lh *.png *.txt
    """
}
