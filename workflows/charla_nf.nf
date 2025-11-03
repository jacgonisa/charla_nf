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
include { SIMPLIFY_HEADERS as SIMPLIFY_HEADERS_READMER } from '../modules/local/simplify_headers.nf'
include { SIMPLIFY_HEADERS as SIMPLIFY_HEADERS_CENHAPMER } from '../modules/local/simplify_headers.nf'


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
include { PLOT_CROSSOVER_MAP } from '../modules/local/plot_crossover_map.nf'
include { NON_HYBRID_READS_ANALYSIS } from '../modules/local/non_hybrid_reads_analysis.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CHARLA_NF {
    take:
    input_tuple           // [fastq_file, sample_id, kmer_size]
    cenhapmer_db_dir      // directory where unique kmc DBs are
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
    
   if (params.input_format == 'fastq') {
    log.info "[CHARLA_NF] Detected input format: FASTQ → Will convert to FASTA"

 // Map input tuple for FASTQ_TO_FASTA
        fastq_ch = input_data.map { meta, fastq, sample_id, kmer_size -> tuple(meta, fastq) }

        // Run FASTQ_TO_FASTA
        fastq_to_fasta_ch = fastq_ch | FASTQ_TO_FASTA

        // Separate the multi-output channel into individual channels
        // FASTQ_TO_FASTA emits a single non-masked FASTA file
        fasta_ch = fastq_to_fasta_ch.fasta


    masked_fastq_ch = fastq_to_fasta_ch.masked_fasta   // quality-masked FASTQ

 // Debug outputs
    fasta_ch.view { meta, fasta_file -> "[DEBUG] FASTQ_TO_FASTA → ${meta.id}: $fasta_file" }
    masked_fastq_ch.view { meta, masked_file -> "[DEBUG] FASTQ_TO_FASTA MASKED → ${meta.id}: $masked_file" }

    log.info "[CHARLA_NF] Deleteme"
readmer_input_ch = fasta_ch.map { meta, fasta_file -> tuple(meta, fasta_file) }
   

readmer_simplified_ch = SIMPLIFY_HEADERS_READMER(readmer_input_ch)



// Simplify headers for Cenhapmer (masked FASTA if available)

//cenhapmer_input_ch = fasta_ch.map { meta, fasta_file -> tuple(meta, fasta_file) }

cenhapmer_input_ch = masked_fastq_ch.map { meta, masked_file -> tuple(meta, masked_file) }

cenhapmer_simplified_ch = SIMPLIFY_HEADERS_CENHAPMER(cenhapmer_input_ch)


  // Define the processed data channels
    processed_data_readmer_ch = readmer_simplified_ch.fasta
    processed_data_cenhapmer_ch = cenhapmer_simplified_ch.fasta

 // Reconstruct tuple for processed_data (if needed for backward compatibility)
    processed_data = processed_data_readmer_ch.combine(
        input_data.map { meta, fastq, sample_id, kmer_size -> 
            tuple(meta.id, sample_id, kmer_size) 
        }
    ).filter { fasta_meta, fasta_file, sample_id, orig_sample_id, kmer_size -> 
        fasta_meta.id == orig_sample_id 
    }.map { fasta_meta, fasta_file, sample_id, orig_sample_id, kmer_size -> 
        tuple(fasta_file, sample_id, kmer_size) 
    }



} else if (params.input_format == 'fasta') {
    log.info "[CHARLA_NF] Detected input format: FASTA → Will simplify headers"

 // Map input for SIMPLIFY_HEADERS
    fasta_input_ch = input_data.map { meta, fasta, sample_id, kmer_size -> tuple(meta, fasta) }
    
    // Use same process for both readmer and cenhapmer
    simplified_ch = SIMPLIFY_HEADERS(fasta_input_ch)

    // Both branches use the same simplified FASTA
    processed_data_readmer_ch = simplified_ch.fasta
    processed_data_cenhapmer_ch = simplified_ch.fasta

    // Debug
    processed_data_readmer_ch.view { "Processed input for readmer profiling: $it" }
    processed_data_cenhapmer_ch.view { "Processed input for cenhapmer profiling: $it" }
 
} else {
    error "Invalid value for --input_format: '${params.input_format}'. Must be 'fastq' or 'fasta'"
}
 
    // Create a channel for the final simplified FASTA file
    
processed_fasta_files = processed_data_cenhapmer_ch // Or processed_data_cenhapmer_ch
    
// Debug
processed_fasta_files.view { "Processed input for final processing: $it" }


// ====================================================================
// READMER BRANCH (always non-masked fasta)
// ====================================================================


// Prepare readmer input channel with correct tuple structure
readmer_input_channel = processed_data_readmer_ch
    .map { meta, fasta_file ->
        tuple(fasta_file, meta.id, params.kmer_size)
    }

readmer_input_channel.view { "Input to readmer branch processes: $it" }

readmer_input_channel | generate_readmers_kmc

// Normalize generate_readmers_kmc output
kmc_out = generate_readmers_kmc.out.map { kmc_pre, kmc_suf, sample_id, kmer_size ->
    tuple(sample_id, kmer_size, kmc_pre, kmc_suf)
}

kmc_out.view { "kmc_out structure: $it" }


generate_readmers_kmc.out.view { "KMC output (readmer): $it" }

// Histogram
generate_readmers_kmc.out | generate_histogram_kmc
generate_histogram_kmc.out.view { "Histogram output (readmer): $it" }



// Ensure readmer_input_channel also starts with (sample_id, kmer_size)
readmer_for_join = readmer_input_channel.map { fasta_file, sample_id, kmer_size ->
    tuple(sample_id, kmer_size, fasta_file)
}
readmer_for_join.view { "readmer_for_join structure: $it" }



// Now join cleanly
counts_input = kmc_out.join(readmer_for_join, by: [0, 1])  // Join by sample_id and kmer_size
    .map { sample_id, kmer_size, kmc_pre, kmc_suf, fasta_file ->
        tuple(kmc_pre, kmc_suf, fasta_file, sample_id, kmer_size)
    }


// Prepare counts input
//counts_input = generate_readmers_kmc.out
//    .join(readmer_input_channel, by: [1, 2]) // Join on sample_id and kmer_size
//    .map { sample_id, kmer_size, kmc_pre, kmc_suf, fasta_file ->
//        tuple(kmc_pre, kmc_suf, fasta_file, sample_id, kmer_size)
//    }



counts_input.view { "Counts input (readmer): $it" }
counts_input | get_counts_kmc
get_counts_kmc.out.view { "Final counts output (readmer): $it" }


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Branch B: CENHAPMER PROFILING SECTION (Parallel Branch)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // 1. Validate cenhapmer_db_dir parameter
    if (!cenhapmer_db_dir) {
        error "ERROR: The --cenhapmer_db_dir parameter must be provided for cenhapmer analysis."
    }

// Prepare cenhapmer input channel
cenhapmer_input_channel = processed_data_cenhapmer_ch.map { meta, fasta_file ->
    tuple(fasta_file, meta.id, params.kmer_size)
}

cenhapmer_input_channel.view { "Cenhapmer input channel: $it" }

// Create parallel cenhapmer input
cenhapmer_parallel_input = ch_kmc_combinations
    .combine(cenhapmer_input_channel)
    .map { acc, cat, chr, fasta_file, sample_id, kmer_size ->
        def db_dir = file(params.cenhapmer_db_dir)
        return tuple(acc, cat, chr, fasta_file, db_dir)
    }


cenhapmer_parallel_input.view { "Cenhapmer parallel input: $it" }

get_counts_cenhapmer_kmc_output = get_counts_cenhapmer_kmc(cenhapmer_parallel_input)


//    get_counts_cenhapmer_kmc_output = get_counts_cenhapmer_kmc(cenhapmer_parallel_input)
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

    // Create dependency barrier
    ch_cenhapmer_dir_path = all_cenhapmer_done
        .map { files ->
            def outdir = params.outdir.replaceAll('/$', '')  // Remove trailing slash
            "${workflow.launchDir}/${outdir}/cenhapmers"
        }

    ch_script_file = Channel.fromPath("${baseDir}/rust/kmer_processor/target/release/kmer_processor")
    ch_fasta_for_postproc = processed_fasta_files.map { meta, fasta -> fasta }.first()
    
    // Now, call process_cenhapmer_counts with the directory path
    process_cenhapmer_counts_output = process_cenhapmer_counts(
        ch_cenhapmer_dir_path,    // Directory path (e.g., "results/cenhapmers")
        ch_fasta_for_postproc,    // The main FASTA file for seqkit stats
        ch_script_file,           // The Python script file
        params.sample_id          // Use sample_id as the output_name for the final matrix
    )  

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Branch C: COMBINE HYBRID PROFILES
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    // 1. Validate index_dir parameter
    if (!params.index_dir) {
        error "ERROR: The --index_dir parameter must be provided for hybrid profile combination."
    }

    def outdir_normalized = params.outdir.replaceAll('/$', '')  // Remove trailing slash
    ch_input_dir_for_hybrid = Channel.value("${workflow.launchDir}/${outdir_normalized}/cenhapmers")
    ch_kmer_tsv = process_cenhapmer_counts_output.kmer_matrix
    ch_qc_file = get_counts_kmc.out
    ch_rust_binary = Channel.fromPath("${baseDir}/rust/combine_segment_fast_indixes/target/release/hybrid_profile_toolkit")
    
    combine_hybrid_profiles_output = COMBINE_HYBRID_PROFILES(
        ch_input_dir_for_hybrid,
        ch_kmer_tsv,
        ch_qc_file,
        params.min_count              
    )


    // ADD THE CLEANUP RIGHT HERE - after COMBINE_HYBRID_PROFILES completes
    combine_hybrid_profiles_output.hybrid_profiles_file.subscribe { hybrid_file ->
        def outdir_clean = params.outdir.replaceAll('/$', '')  // Remove trailing slash
        def cenhapmer_dir = file("${workflow.launchDir}/${outdir_clean}/cenhapmers")
        if (cenhapmer_dir.exists()) {
            log.info "🧹 Cleaning up entire cenhapmers directory to save disk space..."
            if (cenhapmer_dir.deleteDir()) {
                log.info "✅ Deleted cenhapmers directory"
            }
        }
    }


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Branch D: CURATE HYBRID PROFILES
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    def readmer_cutoff = params.readmer_cutoff ? params.readmer_cutoff.toInteger() : 10
    curate_hybrid_profiles_output = CURATE_HYBRID_PROFILES(
        combine_hybrid_profiles_output.hybrid_profiles_file.collect(),
        readmer_cutoff
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Branch E: SEGMENT READS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Debug the channels
    curate_hybrid_profiles_output.final_table.view { "final_table: $it" }
    combine_hybrid_profiles_output.hybrid_profiles_file.view { "hybrid_profiles_file: $it" }
    processed_fasta_files.view { "processed_fasta_files: $it" }
    
    // Segment reads based on threshold and generate BED files
    segment_reads_output = SEGMENT_READS(
        curate_hybrid_profiles_output.final_table.collect(),            
        combine_hybrid_profiles_output.hybrid_profiles_file.collect(),
        processed_data_readmer_ch.map { meta, file -> file }, // Extract file from [meta, file] and Use readmer (non-masked) data
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
    
    ch_reference_genomes = Channel.fromPath(params.reference_genomes_dir, type: 'dir')
        .ifEmpty { error "❌ No reference genome directories found in: ${params.reference_genomes_dir}" }

    // Debug: Check what SEGMENT_READS produced
    segment_reads_output.segmented_fasta.view { "Available segmented_fasta: $it" }
    segment_reads_output.bed_files.view { "Available bed_files: $it" }
    
    arms_fasta_files = segment_reads_output.segmented_fasta
        .flatten()
        .filter { it.name.contains('ARMS.segments.fa') || it.name.contains('ARMS') }
    
    cen_fasta_files = segment_reads_output.segmented_fasta
        .flatten()
        .filter { it.name.contains('CEN.segments.fa') || it.name.contains('CEN') }
    
    arms_bed_files = segment_reads_output.bed_files
        .flatten()
        .filter { it.name.contains('ARMS.bed') || it.name.contains('ARMS') }
    
    cen_bed_files = segment_reads_output.bed_files
        .flatten()
        .filter { it.name.contains('CEN.bed') || it.name.contains('CEN') }

    // Run MAPPING_ANALYSIS with the first available file of each type
    mapping_analysis_output = MAPPING_ANALYSIS(
        arms_fasta_files.first(),   // ARMS segments FASTA
        cen_fasta_files.first(),    // CEN segments FASTA  
        arms_bed_files.first(),     // ARMS BED
        cen_bed_files.first(),      // CEN BED
        ch_reference_genomes,       // Reference genomes
        params.threshold ?: 50      // Threshold
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        CROSSOVER PLOTTING ANALYSIS (OPTIONAL - CENTROMERE-AWARE MODE)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Check if running in centromere-aware mode (BED files provided)
    def centromere_aware_mode = params.parent1_bed && params.parent2_bed

    if (centromere_aware_mode) {
        log.info "[CHARLA_NF] Running in CENTROMERE-AWARE mode - crossover plotting enabled"

        // Create channels for BED files (now using generic parent names)
        ch_parent1_bed = Channel.fromPath(params.parent1_bed)
            .ifEmpty { error "❌ ${params.parent1_name} BED file not found: ${params.parent1_bed}" }

        ch_parent2_bed = Channel.fromPath(params.parent2_bed)
            .ifEmpty { error "❌ ${params.parent2_name} BED file not found: ${params.parent2_bed}" }

        // Filter PAF files for different types and references (now using dynamic parent names)
        parent1_arms_paf_files = mapping_analysis_output.paf_files
            .flatten()
            .filter { it.name.contains("ARMS_against${params.parent1_name}.paf") }

        parent1_cen_paf_files = mapping_analysis_output.paf_files
            .flatten()
            .filter { it.name.contains("CEN_against${params.parent1_name}.paf") }

        parent2_arms_paf_files = mapping_analysis_output.paf_files
            .flatten()
            .filter { it.name.contains("ARMS_against${params.parent2_name}.paf") }

        parent2_cen_paf_files = mapping_analysis_output.paf_files
            .flatten()
            .filter { it.name.contains("CEN_against${params.parent2_name}.paf") }

        // Debug what we found
        parent1_arms_paf_files.view { "Found ${params.parent1_name} ARMS PAF: $it" }
        parent1_cen_paf_files.view { "Found ${params.parent1_name} CEN PAF: $it" }
        parent2_arms_paf_files.view { "Found ${params.parent2_name} ARMS PAF: $it" }
        parent2_cen_paf_files.view { "Found ${params.parent2_name} CEN PAF: $it" }

        // Run crossover plotting with the filtered channels (using generic parent references)
        crossover_plot_output = PLOT_CROSSOVER_MAP(
            ch_parent1_bed.first(),
            ch_parent2_bed.first(),
            parent1_arms_paf_files.first(),
            parent1_cen_paf_files.first(),
            parent2_arms_paf_files.first(),
            parent2_cen_paf_files.first(),
            params.threshold ?: 50,
            params.sample_id
        )
    } else {
        log.info "[CHARLA_NF] Running in CENTROMERE-UNAWARE mode - crossover plotting SKIPPED"
        log.info "[CHARLA_NF] To enable crossover plotting, provide --parent1_bed and --parent2_bed parameters"
        crossover_plot_output = null
    }


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NON-HYBRID READS ANALYSIS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Run non-hybrid reads analysis using the kmer matrix from process_cenhapmer_counts
non_hybrid_analysis_output = NON_HYBRID_READS_ANALYSIS(
    process_cenhapmer_counts_output.kmer_matrix,              // TSV file with kmer counts per read (from process_cenhapmer_counts)
    processed_data_readmer_ch.map { meta, file -> file }.first(),  // Main FASTA file (non-masked)
    ch_reference_genomes.first()                               // Reference genomes directory (reuse from mapping analysis)
)



    emit:
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        EMIT BLOCK
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    // Pre-processing outputs
    fasta = processed_fasta_files 
    
    // K-mer profiling outputs
    kmc_pre = generate_readmers_kmc.out.map { kmc_pre, kmc_suf, sample_id, kmer_size -> kmc_pre }
    kmc_suf = generate_readmers_kmc.out.map { kmc_pre, kmc_suf, sample_id, kmer_size -> kmc_suf }
    histo = generate_histogram_kmc.out.histo
    counts = get_counts_kmc.out
    cenhapmer_counts = get_counts_cenhapmer_kmc_output
    processed_cenhapmer_matrix = process_cenhapmer_counts_output.kmer_matrix
    
    // Hybrid profiling outputs
    hybrid_profiles = combine_hybrid_profiles_output.hybrid_profiles_file
    curated_profiles = curate_hybrid_profiles_output.final_table
    curated_intermediate = curate_hybrid_profiles_output.ultracurated_table
    curated_cut = curate_hybrid_profiles_output.cut_table

    // Segmentation outputs
    filtered_reads = segment_reads_output.filtered_table
    nonectopic_candidates = segment_reads_output.nonectopic_candidates
    arms_candidates = segment_reads_output.arms_candidates
    cen_candidates = segment_reads_output.cen_candidates
    arms_fasta = segment_reads_output.arms_fasta
    cen_fasta = segment_reads_output.cen_fasta
    bed_files = segment_reads_output.bed_files
    segmented_fasta = segment_reads_output.segmented_fasta

    // Mapping outputs
    arms_mapping = mapping_analysis_output.arms_mapping_results
    cen_mapping = mapping_analysis_output.cen_mapping_results
    plotting_results = mapping_analysis_output.plotting_results
    paf_files = mapping_analysis_output.paf_files
    cowidth_files = mapping_analysis_output.cowidth_files

    // Plotting output (only if centromere-aware mode)
    crossover_plots = centromere_aware_mode ? crossover_plot_output.crossover_plots : Channel.empty()
    coverage_data = centromere_aware_mode ? crossover_plot_output.coverage_data : Channel.empty()
 //   plot_log_files = crossover_plot_output.log_files

// Non-hybrid analysis outputs
    nonhybrid_directory = non_hybrid_analysis_output.output_directory
    nonhybrid_read_lists = non_hybrid_analysis_output.read_lists
    nonhybrid_sequences = non_hybrid_analysis_output.extracted_sequences
    nonhybrid_bam_files = non_hybrid_analysis_output.bam_files
    nonhybrid_bai_files = non_hybrid_analysis_output.bai_files
    indel_analysis = non_hybrid_analysis_output.indel_results

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
*/
