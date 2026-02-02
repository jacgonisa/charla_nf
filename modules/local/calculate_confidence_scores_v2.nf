nextflow.enable.dsl=2

process CALCULATE_CONFIDENCE_SCORES_V2 {
    tag "confidence_v2_thr${threshold}"
    label 'confidence_processing'
    cpus params.confidence_cpus ?: 4
    memory params.confidence_memory ?: '8 GB'
    time params.confidence_time ?: '2h'

    publishDir "${params.outdir}/confidence_scores_v2", mode: 'symlink'

    input:
    path arms_summary, stageAs: 'summary_table_ARMS.csv'
    path cen_summary, stageAs: 'summary_table_CEN.csv'
    path arms_bed, stageAs: 'segments_ARMS.bed'
    path cen_bed, stageAs: 'segments_CEN.bed'
    path arms_mapping_results, stageAs: 'mapping_ARMS/**'  // Full ARMS mapping directory
    path cen_mapping_results, stageAs: 'mapping_CEN/**'    // Full CEN mapping directory
    path arms_cowidth, stageAs: 'cowidths_ARMS.txt'
    path cen_cowidth, stageAs: 'cowidths_CEN.txt'
    val threshold

    output:
    path "threshold${threshold}/confidence_scores.tsv", emit: confidence_scores
    path "threshold${threshold}/confidence_scores_raw_metrics.tsv", emit: raw_metrics
    path "threshold${threshold}/diagnostics/**", emit: diagnostic_plots
    path "threshold${threshold}/confidence_summary.txt", emit: summary

    script:
    """
    # Create threshold directory
    threshold_dir="threshold${threshold}"
    mkdir -p \${threshold_dir}/diagnostics

    echo "=== Improved Confidence Scoring v2 ===" > \${threshold_dir}/confidence_summary.txt
    echo "Calculating confidence scores with full diagnostic output..." >> \${threshold_dir}/confidence_summary.txt
    echo "" >> \${threshold_dir}/confidence_summary.txt

    # Run improved confidence scoring
    python ${projectDir}/scripts/14-calculate_confidence_scores_v2.py \\
        --summary-arms summary_table_ARMS.csv \\
        --summary-cen summary_table_CEN.csv \\
        --bed-arms segments_ARMS.bed \\
        --bed-cen segments_CEN.bed \\
        --mapping-base-arms mapping_ARMS \\
        --mapping-base-cen mapping_CEN \\
        --cowidth-arms cowidths_ARMS.txt \\
        --cowidth-cen cowidths_CEN.txt \\
        --output \${threshold_dir}/confidence_scores.tsv \\
        --diagnostic-dir \${threshold_dir}/diagnostics \\
        2>> \${threshold_dir}/confidence_summary.txt

    # Extract high and low confidence read lists
    if [[ -s \${threshold_dir}/confidence_scores.tsv ]]; then
        # High confidence reads (score >= 0.8)
        awk -F'\\t' 'NR>1 && \$2 >= 0.8 {print \$1}' \${threshold_dir}/confidence_scores.tsv > \${threshold_dir}/high_confidence_reads.txt

        # Low confidence reads (score < 0.5)
        awk -F'\\t' 'NR>1 && \$2 < 0.5 {print \$1}' \${threshold_dir}/confidence_scores.tsv > \${threshold_dir}/low_confidence_reads.txt

        # SCO reads
        awk -F'\\t' 'NR>1 && \$4 == "SCO" {print \$1}' \${threshold_dir}/confidence_scores.tsv > \${threshold_dir}/sco_reads.txt

        # High confidence SCO reads
        awk -F'\\t' 'NR>1 && \$4 == "SCO" && \$2 >= 0.8 {print \$1}' \${threshold_dir}/confidence_scores.tsv > \${threshold_dir}/high_confidence_sco_reads.txt

        high_count=\$(wc -l < \${threshold_dir}/high_confidence_reads.txt)
        low_count=\$(wc -l < \${threshold_dir}/low_confidence_reads.txt)
        sco_count=\$(wc -l < \${threshold_dir}/sco_reads.txt)
        high_sco_count=\$(wc -l < \${threshold_dir}/high_confidence_sco_reads.txt)

        echo "" >> \${threshold_dir}/confidence_summary.txt
        echo "Extracted read lists:" >> \${threshold_dir}/confidence_summary.txt
        echo "  High confidence reads (≥0.8): \${high_count}" >> \${threshold_dir}/confidence_summary.txt
        echo "  Low confidence reads (<0.5): \${low_count}" >> \${threshold_dir}/confidence_summary.txt
        echo "  SCO reads: \${sco_count}" >> \${threshold_dir}/confidence_summary.txt
        echo "  High confidence SCO reads: \${high_sco_count}" >> \${threshold_dir}/confidence_summary.txt
    else
        touch \${threshold_dir}/high_confidence_reads.txt
        touch \${threshold_dir}/low_confidence_reads.txt
        touch \${threshold_dir}/sco_reads.txt
        touch \${threshold_dir}/high_confidence_sco_reads.txt
        echo "Warning: No confidence scores generated" >> \${threshold_dir}/confidence_summary.txt
    fi

    echo "" >> \${threshold_dir}/confidence_summary.txt
    echo "Confidence scoring v2 completed!" >> \${threshold_dir}/confidence_summary.txt
    echo "Check diagnostic plots in: \${threshold_dir}/diagnostics/" >> \${threshold_dir}/confidence_summary.txt
    """
}
