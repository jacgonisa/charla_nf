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



include { get_counts_cenhapmer_kmc } from '../modules/local/get_counts_cenhapmer_kmc.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow CHARLA_NF {
    take:
    input_tuple // [fastq_file, sample_id, kmer_size]
    cenhapmer_db_dir     // directory where unique kmc DBs are
    
    main:
    // Debug: Check initial input
    input_tuple.view { "Initial input: $it" }
    
    // Transform input tuple to match FASTQ_TO_FASTA expected format
    // Keep original metadata alongside
    fastq_input = input_tuple.map { fastq_file, sample_id, kmer_size ->
        tuple([id: sample_id], fastq_file, sample_id, kmer_size)
    }
    
    // Debug: Check transformed input
    fastq_input.view { "FASTQ input: $it" }
    
    // Step 1: Convert FASTQ to FASTA
    fastq_for_conversion = fastq_input.map { meta, fastq, sample_id, kmer_size -> 
        tuple(meta, fastq) 
    }
    
    fastq_for_conversion | FASTQ_TO_FASTA
    
    // Step 2: Simplify headers
    FASTQ_TO_FASTA.out.fasta | SIMPLIFY_HEADERS
    
    // Debug: Check SIMPLIFY_HEADERS output
    SIMPLIFY_HEADERS.out.fasta.view { "SIMPLIFY_HEADERS output: $it" }
    
    // Reconstruct the input for KMC processes using join
    // This is much cleaner than cross/filter
    metadata_for_join = fastq_input.map { meta, fastq, sample_id, kmer_size -> 
        tuple(meta, sample_id, kmer_size) 
    }
    
    processed_input = SIMPLIFY_HEADERS.out.fasta
        .join(metadata_for_join)
        .map { meta, fasta_file, sample_id, kmer_size -> 
            tuple(fasta_file, sample_id, kmer_size) 
        }
    // Extract FASTA from SIMPLIFY_HEADERS and create tuple with cenhapmer DB dir
   cenhapmer_input = SIMPLIFY_HEADERS.out.fasta.map { fasta_file ->
           tuple(fasta_file, params.cenhapmer_db_dir)
     }

    // Debug: Check processed_input
    processed_input.view { "Processed input for KMC: $it" }



   // ====================================================================
    // Branch A: Standard KMC Profiling
    // ====================================================================


    
    // Generate k-mer database files (now using FASTA)
    processed_input | generate_readmers_kmc
    
    // Debug: Check KMC output
    generate_readmers_kmc.out.view { "KMC output: $it" }
    
    // Generate histogram from k-mer database
    generate_readmers_kmc.out | generate_histogram_kmc
    
    // Debug: Check histogram output
    generate_histogram_kmc.out.view { "Histogram output: $it" }
    
    // Prepare input for get_counts_kmc
    // We need to combine KMC output with the original FASTA file
    fasta_for_counts = processed_input.map { fasta, sample_id, kmer_size -> 
        tuple(sample_id, fasta) 
    }
    
    counts_input = generate_readmers_kmc.out
        .map { kmc_pre, kmc_suf, sample_id, kmer_size -> 
            tuple(sample_id, kmc_pre, kmc_suf, kmer_size) 
        }
        .join(fasta_for_counts)
        .map { sample_id, kmc_pre, kmc_suf, kmer_size, fasta -> 
            tuple(kmc_pre, kmc_suf, fasta, sample_id, kmer_size) 
        }
    
    // Debug: Check counts input
    counts_input.view { "Counts input: $it" }
 
    // Get k-mer counts for reads (now using processed FASTA)
    counts_input | get_counts_kmc
    
    // Debug: Check final output
    get_counts_kmc.out.view { "Final counts output: $it" }
    

 /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Branch B: CENHAPMER PROFILING SECTION (Parallel Branch)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // 1. Validate cenhapmer_db_dir parameter
    if (!cenhapmer_db_dir) {
        error "ERROR: The --cenhapmer_db_dir parameter must be provided for cenhapmer analysis."
    }


    cenhapmer_input | get_counts_cenhapmer_kmc

    get_counts_cenhapmer_kmc.out.view { "Cenhapmer output: $it" }
    


    emit:
    fasta = SIMPLIFY_HEADERS.out.fasta
    kmc_pre = generate_readmers_kmc.out.map { kmc_pre, kmc_suf, sample_id, kmer_size -> kmc_pre }
    kmc_suf = generate_readmers_kmc.out.map { kmc_pre, kmc_suf, sample_id, kmer_size -> kmc_suf }
    histo = generate_histogram_kmc.out.histo
    counts = get_counts_kmc.out
    cenhapmer_counts = get_counts_cenhapmer_kmc.out

}




/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
*/
