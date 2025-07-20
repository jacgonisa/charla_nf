process get_counts_cenhapmer_kmc {

  tag "$fasta_file.simpleName"

  cpus 50
  memory '200 GB'
  time '10h'



  input:
  tuple path(fasta_file), path(cenhapmer_db_dir)

  output:
  path "cenhapmer_counts/*"

  script:
  """
  mkdir -p cenhapmer_counts

  ACCESSIONS=(Col Ler)
  CATEGORIES=(CEN noCEN)

  for ACC in \${ACCESSIONS[@]}; do
    for CAT in \${CATEGORIES[@]}; do
      for CHR in {1..5}; do
        DB_PREFIX="\${cenhapmer_db_dir}/unique_\${ACC}-0_\${CAT}_Chr\${CHR}_k${params.kmer_size}"

        if [[ -f "\${DB_PREFIX}.kmc_pre" && -f "\${DB_PREFIX}.kmc_suf" ]]; then
          OUT="cenhapmer_counts/\${ACC}_\${CAT}_Chr\${CHR}_k${params.kmer_size}.txt"
          get_counts_threads \${DB_PREFIX} ${fasta_file} 10 > \$OUT
        else
          echo "Missing DB: \${DB_PREFIX}" >&2
        fi
      done
    done
  done
  """
}

