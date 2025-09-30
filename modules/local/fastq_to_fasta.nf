process FASTQ_TO_FASTA {
    tag "$meta.id"
    label 'process_low'

    // Software environment (if needed)
    // conda "bioconda::seqkit=1.4 bioconda::biopython=1.81"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.fa"), emit: fasta
    tuple val(meta), path("*_quality_masked.fa"), optional: true, emit: masked_fasta
    path "versions.yml", emit: versions

    script:
    // Put Groovy vars here so Nextflow parses them correctly
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def min_quality = params.min_quality ?: 14
    def do_masking  = params.do_quality_masking ?: false

    """
    # Convert original FASTQ → non-masked FASTA
    seqkit fq2fa $fastq > ${prefix}.fa

    # Optionally produce a masked FASTA
    ${do_masking ? """
        echo "=== Creating quality-masked FASTA (Q < ${min_quality} → N) ==="
        python3 ${projectDir}/scripts/quality_mask_fastq.py \
            $fastq \
            ${prefix}_quality_masked.fq \
            ${min_quality}
        seqkit fq2fa ${prefix}_quality_masked.fq > ${prefix}_quality_masked.fa
        rm ${prefix}_quality_masked.fq
    """ : ""}

    # Versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(echo \$(seqkit 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}

