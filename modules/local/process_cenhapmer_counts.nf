process process_cenhapmer_counts {
    tag "process_cenhapmer_counts"
   
    publishDir "${params.outdir}/processed_cenhapmers", mode: 'symlink' // Simplified publishDir

    label 'cenhapmer_postprocessing'

    cpus params.cenhapmer_postproc_cpus ?: 20
    memory params.cenhapmer_postproc_memory ?: '8 GB'
    time params.cenhapmer_postproc_time ?: '4h'

    input:
    path cenhapmer_count_files // This will be the directory containing all collected count files
    path reads_fasta_file      // The main FASTA file for seqkit stats
    path script_file           // The Python script file (changed from 'script' to 'script_file' for clarity)
    val output_name            // The desired output file name

    output:
    path "${output_name}", emit: kmer_matrix


    script:
    """
    # Extract number of reads using seqkit
    # Use the provided reads_fasta_file input
    TOTAL_READS=\$(seqkit stats -T ${reads_fasta_file} | awk 'NR==2 {print \$4}')

    echo "Total reads in input file: \$TOTAL_READS"

    # Pass the directory containing all count files to the Python script
    # Use script_file directly as it will be staged into the work directory
    python3 ${script_file} \\
        --base_dir ${cenhapmer_count_files} \\
        --total_reads \$TOTAL_READS \\
        --output processed_cenhapmers/${output_name}
    """
}
