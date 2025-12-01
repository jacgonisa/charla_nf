nextflow.enable.dsl=2

process ANALYZE_SCO_NCO {
    tag "sco_nco_analysis"
    label 'plotting'
    cpus params.plot_cpus ?: 2
    memory params.plot_memory ?: '8 GB'
    time params.plot_time ?: '1h'

    publishDir "${params.outdir}/sco_nco_analysis", mode: 'symlink'

    input:
    path confidence_scores                          // Confidence scores TSV from V3
    path raw_metrics, stageAs: 'raw_metrics.csv'    // Raw metrics (ARMS summary)
    path arms_cowidth, stageAs: 'arms.cowidths'     // ARMS cowidth file
    path cen_cowidth, stageAs: 'cen.cowidths'       // CEN cowidth file
    path segments_bed                               // Segments BED file from SEGMENT_READS

    output:
    path "sco_nco_comparison.png", emit: comparison_plot
    path "sco_nco_resolution_analysis.png", emit: resolution_plot
    path "sco_nco_summary.tsv", emit: summary_table

    script:
    """
    # Combine cowidth files
    cat arms.cowidths > combined_cowidths.txt
    tail -n +2 cen.cowidths >> combined_cowidths.txt

    # Run SCO vs NCO analysis
    python ${projectDir}/scripts/18-analyze_sco_nco_comparison.py \\
        --confidence ${confidence_scores} \\
        --raw-metrics raw_metrics.csv \\
        --cowidth combined_cowidths.txt \\
        --segments-bed ${segments_bed} \\
        --output-prefix sco_nco

    echo "=== SCO vs NCO analysis complete ==="
    ls -lh *.png *.tsv
    """
}
