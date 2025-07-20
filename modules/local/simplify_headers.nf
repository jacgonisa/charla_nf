process SIMPLIFY_HEADERS {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.simpleheader.fa"), emit: fasta
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Check if input is a valid FASTA file
    if ! head -n 1 $fasta | grep -q "^>"; then
        echo "ERROR: $fasta does not appear to be a valid FASTA file"
        exit 1
    fi
    
    # Count sequences before processing
    seq_count_before=\$(grep -c "^>" $fasta)
    echo "Processing \$seq_count_before sequences from $fasta"
    
    # Simplify headers using awk
    awk '/^>/ {sub(/ .*/, "", \$0)} {print}' $fasta > ${prefix}.simpleheader.fa
    
    # Verify output
    seq_count_after=\$(grep -c "^>" ${prefix}.simpleheader.fa)
    echo "Generated \$seq_count_after sequences in ${prefix}.simpleheader.fa"
    

# DEBUG: Attempt to delete and check result
    echo "=== DEBUG: Attempting to delete $fasta ==="
    rm -v $fasta
    exit_code=\$?
    echo "rm exit code: \$exit_code"
    
    # Check if file still exists
    if [ -f "$fasta" ]; then
        echo "WARNING: File $fasta still exists after rm command!"
        ls -la $fasta
    else
        echo "SUCCESS: File $fasta was deleted"
    fi
    
    # List all files in current directory
    echo "=== DEBUG: All files in work directory ==="
    ls -la




    if [ "\$seq_count_before" != "\$seq_count_after" ]; then
        echo "ERROR: Sequence count mismatch - Input: \$seq_count_before, Output: \$seq_count_after"
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n1 | sed 's/^.*awk //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.simpleheader.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n1 | sed 's/^.*awk //; s/ .*\$//')
    END_VERSIONS
    """
}
