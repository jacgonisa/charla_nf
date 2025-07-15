process generate_readmers_kmc {
    tag "${sample_id}"
    
    input:
    tuple path(fastq_file), val(sample_id), val(kmer_size)
    
    output:
    tuple path("${sample_id}_k${kmer_size}_kmer_counts.kmc_pre"), path("${sample_id}_k${kmer_size}_kmer_counts.kmc_suf"), val(sample_id), val(kmer_size)
    
    script:
    """
    mkdir -p tmp_kmc_${sample_id}
    kmc \\
        -k${kmer_size} \\
        -ci1 \\
        -cs10000 \\
        ${fastq_file} \\
        ${sample_id}_k${kmer_size}_kmer_counts \\
        tmp_kmc_${sample_id}
    """
}
