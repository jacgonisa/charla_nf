process process_cenhapmer_counts {
    tag "process_cenhapmer_counts"
    label 'cenhapmer_postprocessing'
    cpus params.cenhapmer_postproc_cpus ?: 20
    memory params.cenhapmer_postproc_memory ?: '32 GB'
    time params.cenhapmer_postproc_time ?: '4h'
    // --- CONDA ENVIRONMENT FOR DEPENDENCIES ---
    // This will create a Conda environment with numpy and seqkit
    // Make sure Conda is enabled in your nextflow.config (which it is for 'local' profile)
    conda (params.enable_conda ? "conda-forge::numpy bioconda::seqkit" : null)
    // --- END CONDA ENVIRONMENT ---
    
    input:
    val cenhapmer_count_dir_path   // Directory path containing the count files (e.g., "results/cenhapmers")
    path reads_fasta_file          // The main FASTA file for seqkit stats
    path script_file               // The Python script file
    val output_name                // The desired output file name
    
    output:
    path "processed_cenhapmers/${output_name}", emit: kmer_matrix
    
    publishDir "${params.outdir}/processed_cenhapmers", mode: 'symlink'
    
    script:
    """
    # Extract number of reads using seqkit
    # Use the provided reads_fasta_file input
    TOTAL_READS=\$(seqkit stats -T ${reads_fasta_file} | awk 'NR==2 {print \$4}')
    echo "Total reads in input file: \$TOTAL_READS"
    
    # Create output directory
    mkdir -p processed_cenhapmers
    
    # Pass the directory path directly to the Python script
    # Use script_file directly as it will be staged into the work directory
    ${projectDir}/rust/kmer_processor/target/release/kmer_processor \\
        --base-dir ${cenhapmer_count_dir_path} \\
        --total-reads \$TOTAL_READS \\
        --output processed_cenhapmers/${output_name} \\
        --cpus ${task.cpus}
    """
}
