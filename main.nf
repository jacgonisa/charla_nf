#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHARLA Mapping Pipeline
    Simplified mapping-only pipeline based on CHARLA_NF (without recombination analysis)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/jacgonisa/charla_nf
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Print help message
def helpMessage() {
    log.info"""
    ========================================
     CHARLA Mapping Pipeline
    ========================================

    Usage:
      nextflow run main.nf --input <input.fq.gz> --sample_id <sample_name> [options]

    Required arguments:
      --input               Path to input FASTQ or FASTA file (long-read sequencing data)
      --sample_id           Sample identifier (e.g., 'sample1', 'mutant_rep1')

    Optional arguments:
      --outdir              Output directory (default: 'results')
      --reference           Reference genome FASTA (default: './Col-0_chr1-5_renamed.fa')
      --reference_dir       Reference genome directory for mapping (default: './index')
      --kmer_db             K-mer database directory (default: './merylDB')
      --kmer_size           K-mer size for analysis (default: 31)
      --input_format        Input format: 'fastq' or 'fasta' (default: 'fastq')

    Resource options:
      --max_cpus            Maximum CPUs (default: 20)
      --max_memory          Maximum memory (default: '100.GB')
      --max_time            Maximum time (default: '48.h')

    Examples:
      # Run with FASTQ input
      nextflow run main.nf --input data/sample1_long_reads.fq.gz --sample_id sample1

      # Run with custom reference
      nextflow run main.nf --input data/mutant_reads.fq.gz \\
                           --sample_id mutant1 \\
                           --reference /path/to/custom_ref.fa
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Validate required parameters
if (!params.input) {
    log.error "ERROR: --input parameter is required"
    helpMessage()
    exit 1
}

if (!params.sample_id) {
    log.error "ERROR: --sample_id parameter is required"
    helpMessage()
    exit 1
}

log.info """
========================================
 CHARLA Mapping Pipeline
========================================
Input file       : ${params.input}
Sample ID        : ${params.sample_id}
Input format     : ${params.input_format}
Reference genome : ${params.reference}
Reference dir    : ${params.reference_dir}
K-mer database   : ${params.kmer_db}
K-mer size       : ${params.kmer_size}
Output directory : ${params.outdir}
========================================
"""

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CHARLA_MAPPING } from './workflows/charla_mapping'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    // Create input channel
    ch_input = Channel.fromPath(params.input, checkIfExists: true)
        .map { file ->
            def meta = [id: params.sample_id]
            tuple(meta, file)
        }

    // Create reference channels
    ch_reference = Channel.fromPath(params.reference, checkIfExists: true)
    ch_reference_dir = Channel.fromPath(params.reference_dir, checkIfExists: true, type: 'dir')
    ch_kmer_db = Channel.fromPath(params.kmer_db, checkIfExists: true, type: 'dir')

    // Run mapping workflow
    CHARLA_MAPPING(
        ch_input,
        ch_reference,
        ch_reference_dir,
        ch_kmer_db
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
