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
include { PLOT_CROSSOVER_MAP } from '../modules/local/plot_crossover_map.nf'

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
        input_data
            .map { meta, fastq, sample_id, kmer_size -> tuple(meta, fastq) }
            | FASTQ_TO_FASTA
        
        // Pass the output of FASTQ_TO_FASTA to SIMPLIFY_HEADERS
        FASTQ_TO_FASTA.out.fasta | SIMPLIFY_HEADERS
        
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
        log.info "[CHARLA_NF] Detected input format: FASTA → Will simplify headers"
        
        // Pass the raw fasta input to SIMPLIFY_HEADERS
        input_data
            .map { meta, fasta, sample_id, kmer_size -> tuple(meta, fasta) }
            | SIMPLIFY_HEADERS

        // Reconstruct full tuple after processing
        processed_data = SIMPLIFY_HEADERS.out.fasta
            .combine(input_data.map { meta, fasta, sample_id, kmer_size -> tuple(meta.id, sample_id, kmer_size) })
            .filter { fasta_meta, fasta_file, sample_id, orig_sample_id, kmer_size ->
                fasta_meta.id == orig_sample_id
            }
            .map { fasta_meta, fasta_file, sample_id, orig_sample_id, kmer_size ->
                tuple(fasta_file, sample_id, kmer_size)
            }

    } else {
        error "Invalid value for --input_format: '${params.input_format}'. Must be 'fastq' or 'fasta'"
    }
    
    // Create a channel for the final simplified FASTA file
    processed_fasta_files = SIMPLIFY_HEADERS.out.fasta
    
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
    cenhapmer_parallel_input = ch_kmc_combinations
        .combine(processed_data)
        .map { acc, cat, chr, fasta_file, sample_id, kmer_size ->
            def db_dir = file(params.cenhapmer_db_dir)
            return tuple(acc, cat, chr, fasta_file, db_dir)
        }

    // Add debug to see what's being passed to cenhapmer
    cenhapmer_parallel_input.view { 
        "Task: ${it[0]}_${it[1]}_${it[2]} - Fasta: ${it[3].name} - DB: ${it[4]}" 
    }

    get_counts_cenhapmer_kmc_output = get_counts_cenhapmer_kmc(cenhapmer_parallel_input)
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
        .map { files -> "${workflow.launchDir}/${params.outdir}/cenhapmers" }

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

    ch_input_dir_for_hybrid = Channel.value("${workflow.launchDir}/${params.outdir}/cenhapmers")
    ch_kmer_tsv = process_cenhapmer_counts_output.kmer_matrix
    ch_qc_file = get_counts_kmc.out
    ch_rust_binary = Channel.fromPath("${baseDir}/rust/combine_segment_fast_indixes/target/release/hybrid_profile_toolkit")
    
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
        processed_fasta_files.map { meta, file -> file }, // Extract file from [meta, file]
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
        CROSSOVER PLOTTING ANALYSIS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Validate input BED files directory
    if (!params.col_bed_file) {
        error "ERROR: The --col_bed_file parameter must be provided for crossover plotting."
    }

    if (!params.ler_bed_file) {
        error "ERROR: The --ler_bed_file parameter must be provided for crossover plotting."
    }

    // Create channels for BED files
    ch_col_bed = Channel.fromPath(params.col_bed_file)
        .ifEmpty { error "❌ Col BED file not found: ${params.col_bed_file}" }

    ch_ler_bed = Channel.fromPath(params.ler_bed_file)
        .ifEmpty { error "❌ Ler BED file not found: ${params.ler_bed_file}" }

    // Filter PAF files for different types and references
    col_arms_paf_files = mapping_analysis_output.paf_files
        .flatten()
        .filter { it.name.contains('ARMS_againstCol.paf') }

    col_cen_paf_files = mapping_analysis_output.paf_files
        .flatten()
        .filter { it.name.contains('CEN_againstCol.paf') }

    ler_arms_paf_files = mapping_analysis_output.paf_files
        .flatten()
        .filter { it.name.contains('ARMS_againstLer.paf') }

    ler_cen_paf_files = mapping_analysis_output.paf_files
        .flatten()
        .filter { it.name.contains('CEN_againstLer.paf') }

    // Debug what we found
    col_arms_paf_files.view { "Found Col ARMS PAF: $it" }
    col_cen_paf_files.view { "Found Col CEN PAF: $it" }
    ler_arms_paf_files.view { "Found Ler ARMS PAF: $it" }
    ler_cen_paf_files.view { "Found Ler CEN PAF: $it" }
    
    // Now, run crossover plotting with the filtered channels.
   // crossover_plot_output = PLOT_CROSSOVER_MAP(
   //     col_arms_paf_files.first(),
   //     col_cen_paf_files.first(),
   //     ler_arms_paf_files.first(),
   //     ler_cen_paf_files.first(),
   //     ch_col_bed.first(),
   //     ch_ler_bed.first(),
   //     params.threshold ?: 50,
   //     params.sample_id
   // )


    // Now, run crossover plotting with the filtered channels.
    crossover_plot_output = PLOT_CROSSOVER_MAP(
    ch_col_bed.first(),
    ch_ler_bed.first(),
    col_arms_paf_files.first(),
    col_cen_paf_files.first(),
    ler_arms_paf_files.first(),
    ler_cen_paf_files.first(),
    params.threshold ?: 50,
    params.sample_id
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

    // Plotting output
    crossover_plots = crossover_plot_output.crossover_plots
    coverage_data = crossover_plot_output.coverage_data
 //   plot_log_files = crossover_plot_output.log_files
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
*/
