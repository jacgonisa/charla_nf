nextflow.enable.dsl=2

process COMBINE_HYBRID_PROFILES {
    tag "combine_hybrid_profiles"
    label 'hybrid_profiling'
    cpus params.hybrid_combine_cpus ?: 20
    memory params.hybrid_combine_memory ?: '100 GB'
    time params.hybrid_combine_time ?: '10h'
    
    // Define variables before they're used
    def output_prefix = "rust_hybridprofiles_${params.sample_id}" // Move this up
    def rust_binary = "${projectDir}/rust/combine_segment_fast_indixes/target/release/hybrid_profile_toolkit"
    
    input:
    val input_dir_of_txt_files // Directory containing individual .txt files from get_counts_cenhapmer_kmc
    path kmer_tsv                      // Output of process_cenhapmer_counts (the aggregated TSV)
    path qc_file                       // Output of get_counts_kmc (the readmer file)
    val min_count                      // Minimum count value

    output:
    path "${output_prefix}.hybridprofiles", emit: hybrid_profiles_file
    path 'indexes/**', emit: index_dir  // Capture all files in indexes directory
 
    publishDir "${params.outdir}/indexes", mode: 'symlink', pattern: 'indexes/**'  // Optionally publish indexes

    script:
    """
    # Debugging: Print input paths before running Rust binary
    echo "DEBUG: input_dir_of_txt_files = ${input_dir_of_txt_files}"
    echo "DEBUG: kmer_tsv = ${kmer_tsv}"
    echo "DEBUG: qc_file = ${qc_file}"
    echo "DEBUG: rust_binary = ${rust_binary}"
    
    # Ensure the Rust binary is executable
    chmod +x ${rust_binary}
    
    # Run the Rust hybrid profile toolkit
    ${rust_binary} \\
        --input-dir ${input_dir_of_txt_files} \\
        --kmer-tsv ${kmer_tsv} \\
        --qc-file ${qc_file} \\
        --output-file ${output_prefix}.hybridprofiles \\
        -v \\
        --index-dir indexes \\
        --min-count ${min_count}
    """
}
