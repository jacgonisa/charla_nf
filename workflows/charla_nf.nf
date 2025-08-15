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
include { COMBINE_HYBRID_PROFILES } from '../modules/local/combine_hybrid_profiles.nf'
include { CURATE_HYBRID_PROFILES } from '../modules/local/curate_hybrid_profiles.nf'
include { SEGMENT_READS } from '../modules/local/segment_reads.nf'
include { MAPPING_ANALYSIS } from '../modules/local/mapping_analysis.nf'

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
 //   cenhapmer_parallel_input | get_counts_cenhapmer_kmc

   get_counts_cenhapmer_kmc_output = get_counts_cenhapmer_kmc(cenhapmer_parallel_input)

    // View output
   // get_counts_cenhapmer_kmc.out.view { "Cenhapmer output: $it" }

   get_counts_cenhapmer_kmc_output.view { "Cenhapmer output: $it" }

  

   // WAIT for all cenhapmer processes to complete before proceeding
    all_cenhapmer_done = get_counts_cenhapmer_kmc_output.collect()


   // Debug what files are collected
all_cenhapmer_done.view { files ->
    println "=== ALL CENHAPMER FILES COLLECTED ==="
    files.each { file ->
        println "File: ${file.name}"
    }
    println "Total files: ${files.size()}"
    return files
}

    // --- NEW MODULE CALL: process_cenhapmer_counts ---
    // FIXED: Pass the absolute directory path instead of collecting files
   // ch_cenhapmer_dir_path = Channel.value("${workflow.launchDir}/${params.outdir}/cenhapmers")
 

// Create dependency barrier - this ensures process_cenhapmer_counts waits
ch_cenhapmer_dir_path = all_cenhapmer_done
    .map { files -> "${workflow.launchDir}/${params.outdir}/cenhapmers" }

   
    // Create a channel for the script path, assuming it's in the 'scripts' directory relative to the pipeline
    // Nextflow will stage this script into the work directory for the process
    ch_script_file = Channel.fromPath("${baseDir}/scripts/process_cenhapmer_script.py")
    
    // Get the single fasta file for seqkit stats
    ch_fasta_for_postproc = processed_fasta_files.map { meta, fasta -> fasta }.first()
     

    // --- END NEW MODULE CALL ---


    // Now, call process_cenhapmer_counts with the directory path
    process_cenhapmer_counts_output = process_cenhapmer_counts(
        ch_cenhapmer_dir_path,   // Directory path (e.g., "results/cenhapmers")
        ch_fasta_for_postproc,   // The main FASTA file for seqkit stats
        ch_script_file,          // The Python script file
        params.sample_id         // Use sample_id as the output_name for the final matrix
    ) 
  

  // --- END NEW MODULE CALL ---

   /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Branch C: COMBINE HYBRID PROFILES
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    // 1. Validate index_dir parameter
    if (!params.index_dir) {
        error "ERROR: The --index_dir parameter must be provided for hybrid profile combination."
    }

    // Prepare inputs for COMBINE_HYBRID_PROFILES
    // input_dir_of_txt_files: This is the directory where get_counts_cenhapmer_kmc publishes its individual .txt files.
    // The 'type: 'dir'' in the module definition will ensure this channel is staged as a directory.
    ch_input_dir_for_hybrid = Channel.value("${workflow.launchDir}/${params.outdir}/cenhapmers")


    // kmer_tsv: This is the output of process_cenhapmer_counts (the aggregated TSV)
    ch_kmer_tsv = process_cenhapmer_counts_output.kmer_matrix

    // qc_file: This is the output of get_counts_kmc (the readmer file)
    // The get_counts_kmc.out channel emits a single file.
    ch_qc_file = get_counts_kmc.out

    // Ensure the Rust binary path is treated as a file for staging
    ch_rust_binary = Channel.fromPath("${baseDir}/rust/combine_segment_fast_indixes/target/release/hybrid_profile_toolkit")
    
    // Call the new module
    combine_hybrid_profiles_output = COMBINE_HYBRID_PROFILES(
        ch_input_dir_for_hybrid,
        ch_kmer_tsv,
        ch_qc_file,
        params.min_count                
    )

   /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Branch D: CURATE HYBRID PROFILES
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
     def readmer_cutoff = params.readmer_cutoff ? params.readmer_cutoff.toInteger() : 10

    // Process the hybrid profiles through curation steps
   curate_hybrid_profiles_output = CURATE_HYBRID_PROFILES(
    combine_hybrid_profiles_output.hybrid_profiles_file.collect(),
    readmer_cutoff  // Add readmer_cutoff parameter with default value
)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Branch E: SEGMENT READS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

 // Debug the channels
 // Debug the channels
    curate_hybrid_profiles_output.final_table.view { "final_table: $it" }
    combine_hybrid_profiles_output.hybrid_profiles_file.view { "hybrid_profiles_file: $it" }
    processed_fasta_files.view { "processed_fasta_files: $it" }
    
    // Segment reads based on threshold and generate BED files
    segment_reads_output = SEGMENT_READS(
        curate_hybrid_profiles_output.final_table.collect(),           
        combine_hybrid_profiles_output.hybrid_profiles_file.collect(), 
        processed_fasta_files.map { meta, file -> file },  // Extract file from [meta, file]
        params.threshold ?: 50,
	params.kmer_size                               
    )


/*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Branch F: MAPPING ANALYSIS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
// Validate reference genomes directory
    if (!params.reference_genomes_dir) {
        error "ERROR: The --reference_genomes_dir parameter must be provided for mapping analysis."
    } 
    // Create channel for reference genomes directory
 //   ch_reference_genomes = Channel.fromPath(params.reference_genomes_dir, type: 'dir')
     ch_reference_genomes = Channel.fromPath(params.reference_genomes_dir, type: 'dir')
    .ifEmpty { error "❌ No reference genome directories found in: ${params.reference_genomes_dir}" }


     // Debug: Check what SEGMENT_READS produced
    segment_reads_output.segmented_fasta.view { "Available segmented_fasta: $it" }
    segment_reads_output.bed_files.view { "Available bed_files: $it" }
   

//    segment_reads_output.segmented_fasta.view { "Segmented FASTA entry: ${it} (class: ${it.getClass()})" }
    
// Flatten the segmented_fasta channel to emit individual files
// Then, filter for the specific files we need based on the actual file names


arms_fasta_files = segment_reads_output.segmented_fasta
    .flatten() // Crucial: Each file from the list becomes a separate item
    .filter {
        // Now 'it' will be a Path object for each individual file
        it.name.contains('ARMS.segments.fa') || it.name.contains('ARMS')
    }
    .ifEmpty {
        // Emit an error if no ARMS FASTA files are found after filtering
        error "❌ No ARMS FASTA segment files (e.g., 'ARMS.segments.fa') found in the segmented FASTA output."
    }

// Apply the same flatten and filter logic for CEN FASTA files
cen_fasta_files = segment_reads_output.segmented_fasta
    .flatten() // Crucial: Each file from the list becomes a separate item
    .filter {
        it.name.contains('CEN.segments.fa') || it.name.contains('CEN')
    }
    .ifEmpty {
        error "❌ No CEN FASTA segment files (e.g., 'CEN.segments.fa') found in the segmented FASTA output."
    }

// Apply the same flatten and filter logic for ARMS BED files
arms_bed_files = segment_reads_output.bed_files
    .flatten() // Assuming bed_files also emits a list
    .filter {
        it.name.contains('ARMS.bed') || it.name.contains('ARMS')
    }
    .ifEmpty {
        error "❌ No ARMS BED files (e.g., 'ARMS.bed') found in the BED files output."
    }

// Apply the same flatten and filter logic for CEN BED files
cen_bed_files = segment_reads_output.bed_files
    .flatten() // Assuming bed_files also emits a list
    .filter {
        it.name.contains('CEN.bed') || it.name.contains('CEN')
    }
    .ifEmpty {
        error "❌ No CEN BED files (e.g., 'CEN.bed') found in the BED files output."
    }


 
    // Debug what we found
    arms_fasta_files.view { "Found ARMS fasta: $it" }
    cen_fasta_files.view { "Found CEN fasta: $it" }
    arms_bed_files.view { "Found ARMS bed: $it" }
    cen_bed_files.view { "Found CEN bed: $it" }
    
    // Run MAPPING_ANALYSIS with the first available file of each type
    mapping_analysis_output = MAPPING_ANALYSIS(
        arms_fasta_files.first(),   // ARMS segments FASTA
        cen_fasta_files.first(),    // CEN segments FASTA  
        arms_bed_files.first(),     // ARMS BED
        cen_bed_files.first(),      // CEN BED
        ch_reference_genomes,       // Reference genomes
        params.threshold ?: 50      // Threshold
    )





    emit:
    fasta = processed_fasta_files // Use consistent channel that works for both FASTQ and FASTA
    kmc_pre = generate_readmers_kmc.out.map { kmc_pre, kmc_suf, sample_id, kmer_size -> kmc_pre }
    kmc_suf = generate_readmers_kmc.out.map { kmc_pre, kmc_suf, sample_id, kmer_size -> kmc_suf }
    histo = generate_histogram_kmc.out.histo
    counts = get_counts_kmc.out
    cenhapmer_counts = get_counts_cenhapmer_kmc.out
    processed_cenhapmer_matrix = process_cenhapmer_counts_output.kmer_matrix // Emit the final processed matrix
    hybrid_profiles = combine_hybrid_profiles_output.hybrid_profiles_file // Emit the final hybrid profiles
    curated_profiles = curate_hybrid_profiles_output.final_table // Emit the final curated table
    curated_intermediate = curate_hybrid_profiles_output.ultracurated_table // Emit intermediate for debugging
    curated_cut = curate_hybrid_profiles_output.cut_table // Emit cut version for debuggingig

        // Segmentation outputs
    filtered_reads = segment_reads_output.filtered_table // Threshold-filtered reads
    nonectopic_candidates = segment_reads_output.nonectopic_candidates // Non-ectopic candidates
    arms_candidates = segment_reads_output.arms_candidates // ARMS candidates
    cen_candidates = segment_reads_output.cen_candidates // CEN candidates
    arms_fasta = segment_reads_output.arms_fasta // ARMS FASTA sequences
    cen_fasta = segment_reads_output.cen_fasta // CEN FASTA sequences
    bed_files = segment_reads_output.bed_files // All BED files
    segmented_fasta = segment_reads_output.segmented_fasta // Final segmented FASTA files



	// Mapping outputs
    
    arms_mapping = mapping_analysis_output.arms_mapping_results
    cen_mapping = mapping_analysis_output.cen_mapping_results
    plotting_results = mapping_analysis_output.plotting_results
    paf_files = mapping_analysis_output.paf_files
    cowidth_files = mapping_analysis_output.cowidth_files


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
*/
