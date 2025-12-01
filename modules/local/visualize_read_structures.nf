nextflow.enable.dsl=2

process VISUALIZE_READ_STRUCTURES {
    tag "read_structures"
    label 'plotting'
    cpus params.plot_cpus ?: 4
    memory params.plot_memory ?: '16 GB'
    time params.plot_time ?: '2h'

    publishDir "${params.outdir}/read_structure_analysis", mode: 'symlink'

    input:
    path arms_mapping_dir, stageAs: 'arms_mapping/**'   // ARMS mapping results directory
    path cen_mapping_dir, stageAs: 'cen_mapping/**'     // CEN mapping results directory
    path confidence_scores                               // Optional: for filtering by type
    val max_reads

    output:
    path "output_read_structures.png", emit: structures_plot
    path "output_parental_analysis.png", emit: parental_plot
    path "output_summary.txt", emit: summary
    path "output_parental_proportions.tsv", emit: proportions_table
    path "output_SCO_read_structures.png", emit: sco_structures, optional: true
    path "output_SCO_parental_analysis.png", emit: sco_parental, optional: true
    path "output_SCO_summary.txt", emit: sco_summary, optional: true
    path "output_NCO_read_structures.png", emit: nco_structures, optional: true
    path "output_NCO_parental_analysis.png", emit: nco_parental, optional: true
    path "output_NCO_summary.txt", emit: nco_summary, optional: true

    script:
    conf_arg = confidence_scores.name != 'NO_FILE' ? "--confidence-scores ${confidence_scores}" : ""
    """
    # Create PAF directory and find all best_alignments.paf files
    mkdir -p paf_files

    echo "=== Searching for full PAF files with MAPQ data ==="

    # Find best_alignments.paf from ARMS (stageAs flattens to numbered symlinks)
    for paf in arms_mapping/*; do
        if [[ -L "\$paf" ]]; then
            target=\$(readlink "\$paf")
            if [[ "\$target" == *"best_alignments.paf" ]]; then
                # Extract mapper name from target path (e.g., summary_mm_sr)
                mapper=\$(echo "\$target" | grep -oP 'summary_[^/]+' | head -1)
                echo "Found ARMS PAF: \$mapper"
                cp "\$paf" paf_files/ARMS_\${mapper}.paf
            fi
        fi
    done

    # Find best_alignments.paf from CEN (stageAs flattens to numbered symlinks)
    for paf in cen_mapping/*; do
        if [[ -L "\$paf" ]]; then
            target=\$(readlink "\$paf")
            if [[ "\$target" == *"best_alignments.paf" ]]; then
                # Extract mapper name from target path (e.g., summary_mm_sr)
                mapper=\$(echo "\$target" | grep -oP 'summary_[^/]+' | head -1)
                echo "Found CEN PAF: \$mapper"
                cp "\$paf" paf_files/CEN_\${mapper}.paf
            fi
        fi
    done

    echo "PAF files collected:"
    ls -lh paf_files/

    echo "=== Analyzing all crossover reads ==="
    python ${projectDir}/scripts/19-visualize_read_structures.py \\
        --paf-dir paf_files \\
        --output-prefix output \\
        --max-reads ${max_reads} \\
        ${conf_arg}

    # If confidence scores provided, also create SCO and NCO specific plots
    if [ -f "${confidence_scores}" ]; then
        echo "=== Analyzing SCO reads ==="
        python ${projectDir}/scripts/19-visualize_read_structures.py \\
            --paf-dir paf_files \\
            --output-prefix output_SCO \\
            --max-reads ${max_reads} \\
            --crossover-type SCO \\
            --confidence-scores ${confidence_scores} || true

        echo "=== Analyzing NCO reads ==="
        python ${projectDir}/scripts/19-visualize_read_structures.py \\
            --paf-dir paf_files \\
            --output-prefix output_NCO \\
            --max-reads ${max_reads} \\
            --crossover-type NCO \\
            --confidence-scores ${confidence_scores} || true
    fi

    echo "=== Read structure visualization complete ==="
    ls -lh *.png *.txt *.tsv
    """
}
