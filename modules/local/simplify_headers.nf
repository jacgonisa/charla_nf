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
    awk '/^>/ {sub(/ .*/, "", \$0)} {print}' $fasta > ${prefix}.simpleheader.fa

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
