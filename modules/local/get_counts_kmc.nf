process get_counts_kmc {
    tag "${sample_id}"
   cpus 20
  memory '200 GB'
  time '12h'
 
    input:
    tuple path(kmc_pre), path(kmc_suf), path(fastq), val(sample_id), val(kmer_size)
    
    output:
    path "${sample_id}_k${kmer_size}_read_kmer_counts.txt"
    
    script:
    """
    get_counts_threads \\
        ${sample_id}_k${kmer_size}_kmer_counts \\
        ${fastq} \\
        20 > ${sample_id}_k${kmer_size}_read_kmer_counts.txt
    """
}
