nextflow.enable.dsl=2

process CALCULATE_CONFIDENCE_SCORES {
    tag "confidence_scores_thr${threshold}"
    label 'confidence_processing'
    cpus params.confidence_cpus ?: 2
    memory params.confidence_memory ?: '4 GB'
    time params.confidence_time ?: '1h'

    publishDir "${params.outdir}/confidence_scores", mode: 'symlink'

    input:
    path arms_summary, stageAs: 'summary_table_ARMS.csv'      // ARMS summary_table.csv
    path cen_summary, stageAs: 'summary_table_CEN.csv'        // CEN summary_table.csv
    path arms_bed, stageAs: 'segments_ARMS.bed'               // ARMS BED file
    path cen_bed, stageAs: 'segments_CEN.bed'                 // CEN BED file
    path paf_files                                            // All PAF files from mapping
    path arms_cowidth, stageAs: 'cowidths_ARMS.txt'           // ARMS cowidths file
    path cen_cowidth, stageAs: 'cowidths_CEN.txt'             // CEN cowidths file
    val threshold                                             // Threshold value

    output:
    path "threshold${threshold}/confidence_scores.tsv", emit: confidence_scores
    path "threshold${threshold}/confidence_summary.txt", emit: summary
    path "threshold${threshold}/high_confidence_reads.txt", emit: high_confidence_reads, optional: true
    path "threshold${threshold}/low_confidence_reads.txt", emit: low_confidence_reads, optional: true

    script:
    """
    # Create threshold directory
    threshold_dir="threshold${threshold}"
    mkdir -p \${threshold_dir}

    # Create directory for PAF files
    mkdir -p \${threshold_dir}/paf_files

    # Copy PAF files to a directory
    if ls *.paf 1> /dev/null 2>&1; then
        cp *.paf \${threshold_dir}/paf_files/ || true
    fi

    echo "Calculating confidence scores for threshold ${threshold}..."

    # Run confidence scoring
    python ${projectDir}/scripts/14-calculate_confidence_scores.py \\
        --summary-arms summary_table_ARMS.csv \\
        --summary-cen summary_table_CEN.csv \\
        --bed-arms segments_ARMS.bed \\
        --bed-cen segments_CEN.bed \\
        --paf-dir \${threshold_dir}/paf_files \\
        --cowidth-arms cowidths_ARMS.txt \\
        --cowidth-cen cowidths_CEN.txt \\
        --output \${threshold_dir}/confidence_scores.tsv \\
        2> \${threshold_dir}/confidence_summary.txt

    # Extract high and low confidence read lists
    if [[ -s \${threshold_dir}/confidence_scores.tsv ]]; then
        # High confidence reads (score >= 0.8)
        awk -F'\\t' 'NR>1 && \$2 >= 0.8 {print \$1}' \${threshold_dir}/confidence_scores.tsv > \${threshold_dir}/high_confidence_reads.txt || touch \${threshold_dir}/high_confidence_reads.txt

        # Low confidence reads (score < 0.5)
        awk -F'\\t' 'NR>1 && \$2 < 0.5 {print \$1}' \${threshold_dir}/confidence_scores.tsv > \${threshold_dir}/low_confidence_reads.txt || touch \${threshold_dir}/low_confidence_reads.txt

        high_count=\$(wc -l < \${threshold_dir}/high_confidence_reads.txt)
        low_count=\$(wc -l < \${threshold_dir}/low_confidence_reads.txt)

        echo "" >> \${threshold_dir}/confidence_summary.txt
        echo "Extracted read lists:" >> \${threshold_dir}/confidence_summary.txt
        echo "  High confidence reads: \${high_count}" >> \${threshold_dir}/confidence_summary.txt
        echo "  Low confidence reads: \${low_count}" >> \${threshold_dir}/confidence_summary.txt
    else
        touch \${threshold_dir}/high_confidence_reads.txt
        touch \${threshold_dir}/low_confidence_reads.txt
        echo "Warning: No confidence scores generated" >> \${threshold_dir}/confidence_summary.txt
    fi

    echo "Confidence scoring completed!"
    """
}
