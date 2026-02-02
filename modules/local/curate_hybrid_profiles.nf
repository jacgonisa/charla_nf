nextflow.enable.dsl=2

process CURATE_HYBRID_PROFILES {
    tag "curate_hybrid_profiles"
    label 'python_processing'
    cpus params.curate_cpus ?: 8
    memory params.curate_memory ?: '8 GB'
    time params.curate_time ?: '4h'

    // PERFORMANCE NOTE: This process uses multi-threaded Rust implementation (Rayon)
    // with chunk-based processing (50k reads/chunk) to balance speed and memory usage.
    // Expected speedup: ~6-10x faster with 8 cores. Increase curate_cpus for more speed.
    // Memory usage: ~1-2GB per chunk + overhead. See docs/PERFORMANCE_OPTIMIZATION.md

    // Define output prefix based on sample
    def output_prefix = "curated_hybrid_${params.sample_id}"
    
    publishDir "${params.outdir}/hybrid_profiles", mode: 'symlink'
    
    input:
    path hybrid_profiles_file  // Output from COMBINE_HYBRID_PROFILES
    val readmer_cutoff


    output:
    path "${output_prefix}_ultracurated.tsv", emit: ultracurated_table
    path "${output_prefix}_ultracurated_cut.tsv", emit: cut_table
    path "${output_prefix}_final_with_combo.tsv", emit: final_table

   path "${output_prefix}_recombination_events.png", emit: events_plot
    path "${output_prefix}_hybrid_types_stacked.png", emit: types_plot
    path "${output_prefix}_ectopic_categories_stacked.png", emit: ectopic_plot
    path "${output_prefix}_summary_stats.txt", emit: summary_stats
    
    script:
    """
    echo "Debug: readmer_cutoff value is: '${readmer_cutoff}'"
echo "Debug: All arguments are:"
echo "Arg 1: ${hybrid_profiles_file}"
echo "Arg 2: ${output_prefix}_ultracurated.tsv"
echo "Arg 3: ${readmer_cutoff}"

    # Step 1: Create the initial curated table
#    python ${projectDir}/scripts/07-create_table_hybrid_typecrossover_ultracurated.py \\
#        ${hybrid_profiles_file} \\
#        ${output_prefix}_ultracurated.tsv ${readmer_cutoff}

##Reimplementation of this step with Rust, due to high memory Killed error

    ${projectDir}/rust/curate_profiles/target/release/curate_profiles \\
        ${hybrid_profiles_file} \\
        ${output_prefix}_ultracurated.tsv \\
        ${readmer_cutoff} \\
        ${task.cpus}


    
    # Step 2: Remove the third column
    cut -f1,2,4- ${output_prefix}_ultracurated.tsv > ${output_prefix}_ultracurated_cut.tsv
    
    # Step 3: Add summary column with event types (e.g., "CC1-LC1" or "CC1-LA2")
    python ${projectDir}/scripts/07.1-summarize_table_hybrid_addcolumncombo_largestcombo_alsonorecombinants.py \\
        ${output_prefix}_ultracurated_cut.tsv \\
        ${output_prefix}_final_with_combo.tsv

    ## Step 4: Generate summary plots from the final table
    python ${projectDir}/scripts/07.2-plot_hybrid_summary.py \\
        ${output_prefix}_final_with_combo.tsv \\
        ${output_prefix}


    """
}
