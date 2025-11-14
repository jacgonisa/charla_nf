process GENERATE_HISTOGRAM_KMC {
    tag "${sample_id}"

    input:
    // This input declaration MUST match the output tuple from the workflow:
    // tuple(kmc_pre, kmc_suf, sample_id, kmer_size)
    tuple path(kmc_pre_file), path(kmc_suf_file), val(sample_id), val(kmer_size)

    output:
    path "${sample_id}_k${kmer_size}_kmer_counts.histo", emit: histo

    script:
    // kmc_tools expects the base name (prefix) of the .kmc_pre/.kmc_suf files
    // The .baseName property extracts "ColLer_LEAF_k41_kmer_counts" from "ColLer_LEAF_k41_kmer_counts.kmc_pre"
    """
    kmc_tools transform \\
        ${kmc_pre_file.baseName} \\
        histogram \\
        "${sample_id}_k${kmer_size}_kmer_counts.histo"
    """

    stub:
    """
    touch "${sample_id}_k${kmer_size}_kmer_counts.histo"
    """
}
