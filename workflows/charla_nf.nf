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
include { process_cenhapmer_counts } from '../modules/local/process_cenhapmer_counts.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CHARLA_NF {
    take:
    input_tuple         // [fastq_file, sample_id, kmer_size]
    cenhapmer_db_dir    // directory where unique kmc DBs are
    ch_kmc_combinations // [acc, cat, chr] - Add this to the take block

    main:
    // Debug: Check initial input
    input_tuple.view { "Initial input: $it" }
    
    // Determine if conversion is needed based on params.input_format
    input_data = input_tuple.map { file_path, sample_id, kmer_size ->
        def meta = [id: sample_id]
        tuple(meta, file_path, sample_id, kmer_size)
    }
    
    // Branch depending on FASTQ or FASTA
    processed_data = Channel.empty()
    processed_fasta_files = Channel.empty()
    
    if (params.input_format == 'fastq') {
        log.info "[CHARLA_NF] Detected input format: FASTQ → Will convert to FASTA"
        input_data
            .map { meta, fastq, sample_id, kmer_size -> tuple(meta, fastq) }
            | FASTQ_TO_FASTA
        
        FASTQ_TO_FASTA.out.fasta | SIMPLIFY_HEADERS
        
        // Store processed FASTA files for emit
        processed_fasta_files = SIMPLIFY_HEADERS.out.fasta
        
        // Reconstruct full tuple after processing
        processed_data = SIMPLIFY_HEADERS.out.fasta
            .combine(input_data.map { meta, fastq, sample_id, kmer_size -> tuple(meta.id, sample_id, kmer_size) })
            .filter { fasta_meta, fasta_file, sample_id, orig_sample_id, kmer_size -> 
                fasta_meta.id == orig_sample_id 
            }
            .map { fasta_meta, fasta_file, sample_id, orig_sample_id, kmer_size -> 
                tuple(fasta_file, sample_id, kmer_size) 
            }
            
    } else if (params.input_format == 'fasta') {
        log.info "[CHARLA_NF] Detected input format: FASTA → Will skip conversion"
        processed_data = input_data.map { meta, fasta, sample_id, kmer_size -> 
            tuple(fasta, sample_id, kmer_size) 
        }
        // Store FASTA files for emit (no processing needed)
        processed_fasta_files = input_data.map { meta, fasta, sample_id, kmer_size -> 
            tuple(meta, fasta) 
        }
    } else {
        error "Invalid value for --input_format: '${params.input_format}'. Must be 'fastq' or 'fasta'"
    }

    // Debug: Check processed input
    processed_data.view { "Processed input for KMC: $it" }

    // ====================================================================
    // Branch A: Standard KMC Profiling
    // ====================================================================
    
    // Generate k-mer database files (now using FASTA)
    processed_data | generate_readmers_kmc
    
    // Debug: Check KMC output
    generate_readmers_kmc.out.view { "KMC output: $it" }
    
    // Generate histogram from k-mer database
    generate_readmers_kmc.out | generate_histogram_kmc
    
    // Debug: Check histogram output
    generate_histogram_kmc.out.view { "Histogram output: $it" }
    
    // Prepare input for get_counts_kmc
    // We need to combine KMC output with the original FASTA file
    fasta_for_counts = processed_data.map { fasta, sample_id, kmer_size -> 
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

    // ~~~ Parallelized CENHAPMER section ~~~
    // Combine each (acc, cat, chr) from ch_kmc_combinations
    // with the single item from processed_data (which contains fasta_file)
    // Create the channel in a simpler way
    // ~~~ THE CORRECT, NON-BLOCKING WAY ~~~
    cenhapmer_parallel_input = ch_kmc_combinations
        .combine(processed_data)
        .map { acc, cat, chr, fasta_file, sample_id, kmer_size ->
            // The process needs the path to the directory, not the directory channel itself
            def db_dir = file(params.cenhapmer_db_dir)
            return tuple(acc, cat, chr, fasta_file, db_dir)
        }

    // Add debug to see what's being passed to cenhapmer
    cenhapmer_parallel_input.view { 
        "Task: ${it[0]}_${it[1]}_${it[2]} - Fasta: ${it[3].name} - DB: ${it[4]}" 
    }

    // Run the process in parallel and assign its output to the variable
    cenhapmer_parallel_input | get_counts_cenhapmer_kmc

    // View output
    get_counts_cenhapmer_kmc.out.view { "Cenhapmer output: $it" }

    // --- NEW MODULE CALL: process_cenhapmer_counts ---
    // FIXED: Pass the absolute directory path instead of collecting files
    ch_cenhapmer_dir_path = Channel.value("${workflow.launchDir}/${params.outdir}/cenhapmers")
    
    // Create a channel for the script path, assuming it's in the 'scripts' directory relative to the pipeline
    // Nextflow will stage this script into the work directory for the process
    ch_script_file = Channel.fromPath("${baseDir}/scripts/process_cenhapmer_script.py")
    
    // Get the single fasta file for seqkit stats
    ch_fasta_for_postproc = processed_fasta_files.map { meta, fasta -> fasta }.first()
    
    // Now, call process_cenhapmer_counts with the directory path
    process_cenhapmer_counts_output = process_cenhapmer_counts(
        ch_cenhapmer_dir_path,   // Directory path (e.g., "results/cenhapmers")
        ch_fasta_for_postproc,   // The main FASTA file for seqkit stats
        ch_script_file,          // The Python script file
        params.sample_id         // Use sample_id as the output_name for the final matrix
    )
    // --- END NEW MODULE CALL ---

    emit:
    fasta = processed_fasta_files // Use consistent channel that works for both FASTQ and FASTA
    kmc_pre = generate_readmers_kmc.out.map { kmc_pre, kmc_suf, sample_id, kmer_size -> kmc_pre }
    kmc_suf = generate_readmers_kmc.out.map { kmc_pre, kmc_suf, sample_id, kmer_size -> kmc_suf }
    histo = generate_histogram_kmc.out.histo
    counts = get_counts_kmc.out
    cenhapmer_counts = get_counts_cenhapmer_kmc.out
    processed_cenhapmer_matrix = process_cenhapmer_counts_output.kmer_matrix // Emit the final processed matrix
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
*/
