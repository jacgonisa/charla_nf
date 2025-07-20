process FASTQ_TO_FASTA {
    tag "$meta.id"
    label 'process_low'
    
    conda "bioconda::seqtk=1.4"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.fa"), emit: fasta
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Check if input is a valid FASTQ file
    if ! head -n 4 $fastq | awk 'NR==1 && !/^@/ {exit 1} NR==3 && !/^\\+/ {exit 1}'; then
        echo "ERROR: $fastq does not appear to be a valid FASTQ file"
        exit 1
    fi
    
    # Count number of reads
    read_count=\$(wc -l < $fastq | awk '{print \$1/4}')
    echo "Processing \$read_count reads from $fastq"
    
    # Convert FASTQ to FASTA using seqtk
    seqkit fq2fa $fastq > ${prefix}.fa
    
    # Verify output
    fasta_count=\$(grep -c "^>" ${prefix}.fa)
    echo "Generated \$fasta_count sequences in ${prefix}.fa"
    
    if [ "\$read_count" != "\$fasta_count" ]; then
        echo "WARNING: Read count mismatch - FASTQ: \$read_count, FASTA: \$fasta_count"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
