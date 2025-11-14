process GET_COUNTS_CENHAPMER {

    tag "${acc}_${cat}_Chr${chr}"

    publishDir "${params.outdir}/cenhapmers", mode: 'copy'

// Use smaller resources now since we're running per combination
cpus params.cenhapmer_cpus ?: 1
    memory params.cenhapmer_memory ?: '5 GB'
    time '10h'


    input:
    // Correctly unpack the single incoming tuple into 5 individual variables
    tuple val(acc), val(cat), val(chr), path(fasta_file), path(cenhapmer_db_dir)

    output:
    path "${acc}_${cat}_Chr${chr}_k${params.kmer_size}.txt"

    script:
    """
    echo "DEBUG: acc=${acc}, cat=${cat}, chr=${chr}" >&2
    echo "DEBUG: fasta_file=${fasta_file}" >&2
    echo "DEBUG: cenhapmer_db_dir=${cenhapmer_db_dir}" >&2
    echo "DEBUG: params.kmer_size=${params.kmer_size}" >&2



    DB_PREFIX="${cenhapmer_db_dir}/unique_${acc}-0_${cat}_Chr${chr}_k${params.kmer_size}"

    if [[ -f "\${DB_PREFIX}.kmc_pre" && -f "\${DB_PREFIX}.kmc_suf" ]]; then
        get_counts_threads \${DB_PREFIX} ${fasta_file} 10 > ${acc}_${cat}_Chr${chr}_k${params.kmer_size}.txt
    else
        echo "Missing DB: \${DB_PREFIX}" >&2
        exit 1
    fi
    """
}

