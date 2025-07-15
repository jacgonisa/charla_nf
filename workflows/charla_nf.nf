nextflow.enable.dsl=2
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Keep commented out as before
// include { paramsSummaryMap } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
// Keep commented out as before
// include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_charla_nf_pipeline'

// New modules for preprocessing
include { FASTQ_TO_FASTA } from '../modules/local/fastq_to_fasta.nf'
include { SIMPLIFY_HEADERS } from '../modules/local/simplify_headers.nf'

// Existing modules
include { generate_readmers_kmc } from '../modules/local/generate_readmers_kmc.nf'
include { generate_histogram_kmc } from '../modules/local/generate_histogram_kmc.nf'
include { get_counts_kmc } from '../modules/local/get_counts_kmc.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow CHARLA_NF {
    take:
    input_tuple // input channel
    
    main:
    // Transform input tuple to match FASTQ_TO_FASTA expected format
    // Input: [fastq_file, sample_id, kmer_size]
    // Output: [meta, fastq_file]
    fastq_input = input_tuple.map { fastq_file, sample_id, kmer_size ->
        tuple([id: sample_id], fastq_file)
    }
    
    // Step 1: Convert FASTQ to FASTA
    fastq_input | FASTQ_TO_FASTA
    
    // Step 2: Simplify headers
    FASTQ_TO_FASTA.out.fasta | SIMPLIFY_HEADERS
    
    // Update the input for downstream processes to use the processed FASTA
    // We need to recreate the original tuple format: [fasta_file, sample_id, kmer_size]
    processed_input = SIMPLIFY_HEADERS.out.fasta
        .cross(input_tuple)
        .filter { processed, original -> 
            processed[0].id == original[1] }  // match sample_id
        .map { processed, original -> 
            tuple(processed[1], original[1], original[2]) }  // [fasta_file, sample_id, kmer_size]
    
    // Generate k-mer database files (now using FASTA instead of FASTQ)
    processed_input | generate_readmers_kmc
    
    // Generate histogram from k-mer database
    generate_readmers_kmc.out | generate_histogram_kmc
    
    // Prepare input for get_counts_kmc by adding the processed fasta file
    // Combine the kmc output with the processed input to get all required data
    counts_input = generate_readmers_kmc.out
        .cross(processed_input)
        .filter { kmc_output, processed_input -> 
            kmc_output[2] == processed_input[1] && kmc_output[3] == processed_input[2] }
        .map { kmc_output, processed_input -> 
            tuple(kmc_output[0], kmc_output[1], processed_input[0], kmc_output[2], kmc_output[3]) }
    
    // Debug: Check what goes into get_counts_kmc
    counts_input.view { "Counts input: $it" }
 
    // Get k-mer counts for reads (now using processed FASTA)
    counts_input | get_counts_kmc
    
    emit:
    fasta = SIMPLIFY_HEADERS.out.fasta
    kmc_pre = generate_readmers_kmc.out.map { kmc_pre, kmc_suf, sample_id, kmer_size -> kmc_pre }
    kmc_suf = generate_readmers_kmc.out.map { kmc_pre, kmc_suf, sample_id, kmer_size -> kmc_suf }
    histo = generate_histogram_kmc.out.histo
    counts = get_counts_kmc.out
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
*/
