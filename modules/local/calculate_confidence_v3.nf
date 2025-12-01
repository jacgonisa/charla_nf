nextflow.enable.dsl=2

process CALCULATE_CONFIDENCE_V3 {
    tag "confidence_v3_thr${threshold}"
    label 'confidence_processing'
    cpus params.confidence_cpus ?: 4
    memory params.confidence_memory ?: '8 GB'
    time params.confidence_time ?: '2h'

    publishDir "${params.outdir}/confidence_scores_v3", mode: 'symlink'

    input:
    path parental_proportions                      // From read structure analysis
    path curated_profiles                           // Curated hybrid profiles (SOURCE OF TRUTH for classification)
    path arms_summary, stageAs: 'arms_summary.csv'
    path cen_summary, stageAs: 'cen_summary.csv'
    path arms_bed, stageAs: 'arms.bed'
    path cen_bed, stageAs: 'cen.bed'
    val threshold

    output:
    path "threshold${threshold}/confidence_scores.tsv", emit: confidence_scores
    path "threshold${threshold}/diagnostics/**", emit: diagnostic_plots
    path "threshold${threshold}/high_confidence_reads.txt", emit: high_conf_reads
    path "threshold${threshold}/high_confidence_sco_reads.txt", emit: high_conf_sco_reads

    script:
    """
    # Create directories
    threshold_dir="threshold${threshold}"
    mkdir -p \${threshold_dir}/diagnostics

    echo "=== Confidence Scoring V3 (Redesigned Metrics) ==="

    # Run V3 scoring
    python ${projectDir}/scripts/20-calculate_confidence_v3.py \\
        --parental-proportions ${parental_proportions} \\
        --curated-profiles ${curated_profiles} \\
        --threshold ${threshold} \\
        --summary-arms arms_summary.csv \\
        --summary-cen cen_summary.csv \\
        --output \${threshold_dir}/confidence_scores.tsv \\
        --diagnostic-dir \${threshold_dir}/diagnostics

    # Extract high confidence read lists
    if [[ -s \${threshold_dir}/confidence_scores.tsv ]]; then
        # High confidence reads (score >= 0.8)
        awk -F'\\t' 'NR>1 && \$2 >= 0.8 {print \$1}' \${threshold_dir}/confidence_scores.tsv > \${threshold_dir}/high_confidence_reads.txt

        # High confidence SCO reads
        awk -F'\\t' 'NR>1 && \$8 == "SCO" && \$2 >= 0.8 {print \$1}' \${threshold_dir}/confidence_scores.tsv > \${threshold_dir}/high_confidence_sco_reads.txt

        high_count=\$(wc -l < \${threshold_dir}/high_confidence_reads.txt)
        high_sco_count=\$(wc -l < \${threshold_dir}/high_confidence_sco_reads.txt)

        echo ""
        echo "Extracted read lists:"
        echo "  High confidence reads (≥0.8): \${high_count}"
        echo "  High confidence SCO reads: \${high_sco_count}"
    else
        touch \${threshold_dir}/high_confidence_reads.txt
        touch \${threshold_dir}/high_confidence_sco_reads.txt
    fi

    echo "=== Confidence scoring V3 complete ==="
    """
}
