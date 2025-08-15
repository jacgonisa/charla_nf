process generate_readmers_kmc {
    tag "${sample_id}"
    publishDir "${params.outdir}/readmers", mode: 'symlink'
    
    input:
    tuple path(fasta_file), val(sample_id), val(kmer_size)
    
    output:
    tuple path("${sample_id}_k${kmer_size}_kmer_counts.kmc_pre"), path("${sample_id}_k${kmer_size}_kmer_counts.kmc_suf"), val(sample_id), val(kmer_size)
    
    script:
    """
    mkdir -p tmp_kmc_${sample_id}
    
    # Verify input is FASTA
    if ! head -n 1 ${fasta_file} | grep -q "^>"; then
        echo "ERROR: ${fasta_file} does not appear to be a valid FASTA file"
        exit 1
    fi
    
    # Count sequences for logging
    seq_count=\$(grep -c "^>" ${fasta_file})
    echo "Processing \$seq_count sequences from ${fasta_file}"
    
    # Use KMC with FASTA format flag
    kmc \\
        -k${kmer_size} \\
        -ci1 \\
        -cs10000 \\
        -fm \\
        ${fasta_file} \\
        ${sample_id}_k${kmer_size}_kmer_counts \\
        tmp_kmc_${sample_id}
    
    # Verify output files were created
    if [[ ! -f "${sample_id}_k${kmer_size}_kmer_counts.kmc_pre" || ! -f "${sample_id}_k${kmer_size}_kmer_counts.kmc_suf" ]]; then
        echo "ERROR: KMC failed to generate output files"
        exit 1
    fi
    
    echo "KMC completed successfully for ${sample_id}"
    """
}
