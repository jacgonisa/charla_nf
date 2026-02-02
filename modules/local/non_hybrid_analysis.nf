nextflow.enable.dsl=2

process NON_HYBRID_ANALYSIS {
    tag "non_hybrid_analysis"
    label 'plotting'
    cpus params.plot_cpus ?: 2
    memory params.plot_memory ?: '8 GB'
    time params.plot_time ?: '2h'

    publishDir "${params.outdir}/non_hybrid_analysis", mode: 'symlink'

    input:
    path(mapping_directory)  // Directory from NON_HYBRID_MAPPING with BAM files

    output:
    path("${params.sample_id}_nonhybrid_analysis/"), emit: output_directory
    path("${params.sample_id}_nonhybrid_analysis/indel_analysis_results.tsv"), emit: indel_results
    path("${params.sample_id}_nonhybrid_analysis/plots/*.png"), emit: plots, optional: true

    script:
    """
    # Create output directory
    mkdir -p ${params.sample_id}_nonhybrid_analysis/plots

    # Check if BAM files exist
    if [ ! -d "${mapping_directory}/mappings" ] || [ -z "\$(find ${mapping_directory}/mappings -name '*.bam' 2>/dev/null)" ]; then
        echo "No BAM files found. Creating empty results file."
        echo -e "haplotype\\tmapping_mode\\ttotal_reads\\tmapped_reads\\tmapped_mb\\tindels\\tinsertions\\tdeletions\\tavg_indel_size\\tindel_density\\treads_with_indels\\tbackground\\tregion\\tchromosome" > ${params.sample_id}_nonhybrid_analysis/indel_analysis_results.tsv
        exit 0
    fi

    # Copy mapping directory structure to analysis directory
    # This allows the analysis script to find BAM files in expected locations
    cp -r ${mapping_directory} ${params.sample_id}_nonhybrid_temp

    # Step 3: Analyze indels from BAM files (FAST - includes plotting)
    echo "=== Analyzing indels from alignments and generating plots ==="
    python ${projectDir}/scripts/12.2-indel_analysis.py \\
        ${params.sample_id}_nonhybrid_temp \\
        ${params.sample_id}_nonhybrid_analysis/indel_analysis_results.tsv

    # Move plots to output directory if they were created
    if [ -d "${params.sample_id}_nonhybrid_temp/plots" ]; then
        mv ${params.sample_id}_nonhybrid_temp/plots/* ${params.sample_id}_nonhybrid_analysis/plots/ 2>/dev/null || true
    fi

    # Clean up temp directory
    rm -rf ${params.sample_id}_nonhybrid_temp

    echo "=== Non-hybrid reads analysis complete ==="
    """
}
