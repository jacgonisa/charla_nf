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

// Post-mapping analysis modules (INTEGRATED into pipeline)
// These are in postmapping_analysis/ folder but integrated here
include { CROSSOVER_STATS as CROSSOVER_STATS_PROFILES } from '../postmapping_analysis/modules/crossover_stats.nf'
include { CROSSOVER_STATS as CROSSOVER_STATS_MAPPING } from '../postmapping_analysis/modules/crossover_stats.nf'
include { CURATE_CONFIDENT_READS } from '../postmapping_analysis/modules/curate_confident_reads.nf'
include { CREATE_COMPREHENSIVE_ALIGNMENT_OUTPUT } from '../postmapping_analysis/modules/create_comprehensive_alignment_output.nf'

// Post-mapping analysis modules (NOT YET INTEGRATED - available for manual use)
// Uncomment these when ready to integrate visualization into pipeline:
// include { PLOT_CONFIDENT_READS } from '../postmapping_analysis/modules/plot_confident_reads.nf'
// include { KARYOPLOT_CROSSOVERS } from '../postmapping_analysis/modules/karyoplot_crossovers.nf'

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

// Detect centromere-aware vs centromere-unaware mode
def centromere_aware_mode = params.parent1_bed && params.parent2_bed
log.info "[CHARLA_NF] Mode: ${centromere_aware_mode ? 'CENTROMERE-AWARE (CEN/ARMS)' : 'CENTROMERE-UNAWARE (CHR only)'}"

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
            // Check if outdir is already an absolute path
            def cenhapmer_dir = outdir.startsWith('/') ?
                "${outdir}/cenhapmers" :
                "${workflow.launchDir}/${outdir}/cenhapmers"
            cenhapmer_dir
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
    // Check if outdir is already an absolute path
    def input_dir = outdir_normalized.startsWith('/') ?
        "${outdir_normalized}/cenhapmers" :
        "${workflow.launchDir}/${outdir_normalized}/cenhapmers"
    ch_input_dir_for_hybrid = Channel.value(input_dir)
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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CROSSOVER STATISTICS - STAGE 1: K-MER PROFILE ANALYSIS (Before Mapping)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    log.info "[CHARLA_NF] Generating crossover statistics (k-mer profile stage)"

    // Run crossover statistics on k-mer profiles (without cowidths)
    crossover_stats_profiles_output = CROSSOVER_STATS_PROFILES(
        [],  // No cowidths/bed files at this stage
        curate_hybrid_profiles_output.final_table,
        params.threshold ?: 50,
        params.sample_id,
        "profiles",
        params.parent1_name,
        params.parent2_name
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
        params.kmer_size,
        centromere_aware_mode  // Pass mode to handle CHR vs CEN/ARMS splitting
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
    
    // Verify the reference genomes directory exists
    def ref_dir = file(params.reference_genomes_dir)
    if (!ref_dir.exists() || !ref_dir.isDirectory()) {
        error "❌ Reference genomes directory does not exist or is not a directory: ${params.reference_genomes_dir}"
    }

    // Create a channel with the directory path
    ch_reference_genomes = Channel.value(ref_dir)

    // Debug: Check what SEGMENT_READS produced
    segment_reads_output.segmented_fasta.view { "Available segmented_fasta: $it" }
    segment_reads_output.bed_files.view { "Available bed_files: $it" }

    // Handle file filtering based on mode
    if (centromere_aware_mode) {
        log.info "[CHARLA_NF] Filtering for ARMS/CEN files (centromere-aware mode)"
        arms_fasta_files = segment_reads_output.segmented_fasta
            .flatten()
            .filter { it.name.contains('ARMS.segments.fa') }

        cen_fasta_files = segment_reads_output.segmented_fasta
            .flatten()
            .filter { it.name.contains('CEN.segments.fa') }

        arms_bed_files = segment_reads_output.bed_files
            .flatten()
            .filter { it.name.contains('ARMS.bed') }

        cen_bed_files = segment_reads_output.bed_files
            .flatten()
            .filter { it.name.contains('CEN.bed') }
    } else {
        log.info "[CHARLA_NF] Filtering for CHR files (centromere-unaware mode)"
        // In CHR mode, use CHR files for "arms" input and create empty channel for "cen"
        arms_fasta_files = segment_reads_output.segmented_fasta
            .flatten()
            .filter { it.name.contains('CHR.segments.fa') }

        cen_fasta_files = Channel.empty()

        arms_bed_files = segment_reads_output.bed_files
            .flatten()
            .filter { it.name.contains('CHR.bed') }

        cen_bed_files = Channel.empty()
    }

    // Run MAPPING_ANALYSIS with the first available file of each type
    // Create unique placeholders for empty CEN files in centromere-unaware mode
    def empty_cen_fasta = file("${workDir}/NO_CEN_FASTA")
    def empty_cen_bed = file("${workDir}/NO_CEN_BED")

    mapping_analysis_output = MAPPING_ANALYSIS(
        arms_fasta_files.first(),   // ARMS/CHR segments FASTA
        cen_fasta_files.ifEmpty { empty_cen_fasta }.first(),    // CEN segments FASTA (or placeholder)
        arms_bed_files.first(),     // ARMS/CHR BED
        cen_bed_files.ifEmpty { empty_cen_bed }.first(),      // CEN BED (or placeholder)
        ch_reference_genomes,       // Reference genomes
        params.threshold ?: 50      // Threshold
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        CROSSOVER PLOTTING ANALYSIS (OPTIONAL - CENTROMERE-AWARE MODE)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Check if running in centromere-aware mode (already defined at workflow start)
    if (centromere_aware_mode) {
        log.info "[CHARLA_NF] Running in CENTROMERE-AWARE mode - crossover plotting with centromere annotations"

        // Filter PAF files for different types and references
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

        // Run crossover plotting with centromere annotations
        crossover_plot_output = PLOT_CROSSOVER_MAP(
            params.parent1_bed,
            params.parent2_bed,
            parent1_arms_paf_files.first(),
            parent1_cen_paf_files.first(),
            parent2_arms_paf_files.first(),
            parent2_cen_paf_files.first(),
            params.threshold ?: 50,
            params.sample_id
        )
    } else {
        log.info "[CHARLA_NF] Running in CENTROMERE-UNAWARE mode - chromosome-level crossover plotting"

        // Filter PAF files for CHR (chromosome-level)
        parent1_chr_paf_files = mapping_analysis_output.paf_files
            .flatten()
            .filter { it.name.contains("CHR_against${params.parent1_name}.paf") }

        parent2_chr_paf_files = mapping_analysis_output.paf_files
            .flatten()
            .filter { it.name.contains("CHR_against${params.parent2_name}.paf") }

        // Get empty CEN PAF files created by MAPPING_ANALYSIS as placeholders
        parent1_cen_empty = mapping_analysis_output.paf_files
            .flatten()
            .filter { it.name.contains("CEN_against${params.parent1_name}.paf") }

        parent2_cen_empty = mapping_analysis_output.paf_files
            .flatten()
            .filter { it.name.contains("CEN_against${params.parent2_name}.paf") }

        // Debug what we found
        parent1_chr_paf_files.view { "Found ${params.parent1_name} CHR PAF: $it" }
        parent2_chr_paf_files.view { "Found ${params.parent2_name} CHR PAF: $it" }
        parent1_cen_empty.view { "Found ${params.parent1_name} empty CEN PAF: $it" }
        parent2_cen_empty.view { "Found ${params.parent2_name} empty CEN PAF: $it" }

        // Run crossover plotting without centromere annotations
        crossover_plot_output = PLOT_CROSSOVER_MAP(
            "NA",  // No BED file for parent 1
            "NA",  // No BED file for parent 2
            parent1_chr_paf_files.first(),
            parent1_cen_empty.first(),  // Empty placeholder for parent1 CEN
            parent2_chr_paf_files.first(),
            parent2_cen_empty.first(),  // Empty placeholder for parent2 CEN
            params.threshold ?: 50,
            params.sample_id
        )
    }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CROSSOVER STATISTICS - STAGE 2: MAPPING ANALYSIS (After Mapping)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    log.info "[CHARLA_NF] Generating crossover statistics with width analysis (mapping stage)"

    // Get all plotting crossover files (cowidths and bed_of_best.tsv)
    plotting_files = mapping_analysis_output.plotting_results
        .flatten()
        .filter { it.name.endsWith('.cowidths') || it.name.endsWith('.bed_of_best.tsv') }
        .collect()

    // Get the curated profiles file from curate_hybrid
    curated_profiles_file = curate_hybrid_profiles_output.final_table

    // Run crossover statistics with cowidths and bed files
    crossover_stats_mapping_output = CROSSOVER_STATS_MAPPING(
        plotting_files,
        curated_profiles_file,
        params.threshold ?: 50,
        params.sample_id,
        "mapping",
        params.parent1_name,
        params.parent2_name
    )

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CURATE CONFIDENT READS (After Crossover Statistics)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    log.info "[CHARLA_NF] Curating confident recombinant reads with segment count filtering"

    // Get specific files from mapping analysis
    cowidths = mapping_analysis_output.cowidth_files
    bed_of_best = plotting_files.flatten().filter { it.name.contains("bed_of_best.tsv") }

    // Run curation module (PAF filtering disabled to avoid file staging issues)
    curated_confident_reads_output = CURATE_CONFIDENT_READS(
        cowidths,
        bed_of_best,
        curated_profiles_file,
        params.threshold ?: 50,
        params.sample_id,
        params.parent1_name,
        params.parent2_name
    )

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CREATE COMPREHENSIVE ALIGNMENT OUTPUT (Addressing Collaborator Feedback)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    log.info "[CHARLA_NF] Creating comprehensive alignment output with PAF headers and coordinate data"

    // Get segment PAF files from mapping analysis
    segment_paf_b73 = mapping_analysis_output.segment_paf_parent1.first()
    segment_paf_mo17 = mapping_analysis_output.segment_paf_parent2.first()

    // Get curated confident reads directory
    curated_dir = curated_confident_reads_output.confident_bed_files.first().parent

    // Create comprehensive alignment output
    comprehensive_alignment_output = CREATE_COMPREHENSIVE_ALIGNMENT_OUTPUT(
        curated_dir,
        segment_paf_b73,
        segment_paf_mo17,
        params.sample_id
    )

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PLOT CONFIDENT READS (Comprehensive Visualization)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    log.info "[CHARLA_NF] Generating comprehensive plots for confident reads"

    // Get confident BED files from curation output
    confident_bed_sco = curated_confident_reads_output.confident_bed_files
        .flatten()
        .filter { it.name.contains("_confident_SCO.bed") }
        .first()

    confident_bed_nco = curated_confident_reads_output.confident_bed_files
        .flatten()
        .filter { it.name.contains("_confident_NCO.bed") }
        .first()

    confident_bed_gcco = curated_confident_reads_output.confident_bed_files
        .flatten()
        .filter { it.name.contains("_confident_GC-CO.bed") }
        .first()

    // Get PAF files from mapping analysis
    parent1_paf_chr = plotting_files
        .flatten()
        .filter { it.name.contains("summary_mm_sr_CHR_against${params.parent1_name}.paf") }
        .first()

    parent2_paf_chr = plotting_files
        .flatten()
        .filter { it.name.contains("summary_mm_sr_CHR_against${params.parent2_name}.paf") }
        .first()

    // Create genome index channels (use different variable names to avoid collision)
    def parent1_fai_plot_path = file("${params.reference_genomes_dir}/${params.parent1_name}.fa.fai")
    def parent2_fai_plot_path = file("${params.reference_genomes_dir}/${params.parent2_name}.fa.fai")

    // Check if files exist and log warnings, but run module anyway
    if (!parent1_fai_plot_path.exists()) {
        log.warn "⚠️  Parent1 genome index not found: ${params.reference_genomes_dir}/${params.parent1_name}.fa.fai"
    }
    if (!parent2_fai_plot_path.exists()) {
        log.warn "⚠️  Parent2 genome index not found: ${params.reference_genomes_dir}/${params.parent2_name}.fa.fai"
    }

    // Create channels
    parent1_fai_plot = Channel.fromPath(parent1_fai_plot_path)
    parent2_fai_plot = Channel.fromPath(parent2_fai_plot_path)

    // Centromere BED file (optional)
    def cen_bed_val = (params.parent1_bed && file(params.parent1_bed).exists()) ? file(params.parent1_bed) : "NO_FILE"

    // TODO: PLOT_CONFIDENT_READS module not yet integrated
    // Uncomment when ready to integrate visualization
    /*
    confident_plots_output = PLOT_CONFIDENT_READS(
        confident_bed_sco,
        confident_bed_nco,
        confident_bed_gcco,
        curated_confident_reads_output.confident_read_ids,
        parent1_paf_chr,
        parent2_paf_chr,
        parent1_fai_plot,
        parent2_fai_plot,
        cowidths,
        curated_profiles_file,
        params.threshold ?: 50,
        params.sample_id,
        params.parent1_name,
        params.parent2_name,
        cen_bed_val
    )
    */
    log.info "⏩ Skipping PLOT_CONFIDENT_READS (not yet integrated - available in postmapping_analysis/)"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    KARYOPLOT VISUALIZATION (After Mapping)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    log.info "[CHARLA_NF] Generating karyotype plot of crossovers"

    // NOTE: This requires genome FASTA index (.fai) files for both parents
    // These should be located in the reference_genomes_dir with the naming convention:
    //   ${parent1_name}.fa.fai and ${parent2_name}.fa.fai
    // Generate .fai files using: samtools faidx <genome.fa>

    // Get the bed_of_best.tsv file from plotting_files
    bed_file = plotting_files
        .flatten()
        .filter { it.name.endsWith('.bed_of_best.tsv') }
        .first()

    // Create genome index file channels
    // Check if .fai files exist before creating channels
    def parent1_fai_karyo_path = file("${params.reference_genomes_dir}/${params.parent1_name}.fa.fai")
    def parent2_fai_karyo_path = file("${params.reference_genomes_dir}/${params.parent2_name}.fa.fai")

    def both_fai_exist_karyo = parent1_fai_karyo_path.exists() && parent2_fai_karyo_path.exists()

    if (!parent1_fai_karyo_path.exists()) {
        log.warn "⚠️  Parent1 genome index not found: ${params.reference_genomes_dir}/${params.parent1_name}.fa.fai"
        log.warn "⚠️  Skipping karyoplot generation. Generate with: samtools faidx ${params.parent1_name}.fa"
    }

    if (!parent2_fai_karyo_path.exists()) {
        log.warn "⚠️  Parent2 genome index not found: ${params.reference_genomes_dir}/${params.parent2_name}.fa.fai"
        log.warn "⚠️  Skipping karyoplot generation. Generate with: samtools faidx ${params.parent2_name}.fa"
    }

    // TODO: KARYOPLOT_CROSSOVERS module not yet integrated
    // Uncomment when ready to integrate karyoplot visualization
    /*
    if (both_fai_exist_karyo) {
        parent1_fai = Channel.fromPath(parent1_fai_karyo_path)
        parent2_fai = Channel.fromPath(parent2_fai_karyo_path)

        karyoplot_output = KARYOPLOT_CROSSOVERS(
            bed_file,
            parent1_fai,
            parent2_fai,
            curated_profiles_file,
            params.threshold ?: 50,
            params.sample_id,
            params.parent1_name,
            params.parent2_name
        )
    } else {
        log.info "⏩ Skipping karyoplot generation (genome index files not found)"
    }
    */
    log.info "⏩ Skipping KARYOPLOT_CROSSOVERS (not yet integrated - available in postmapping_analysis/)"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NON-HYBRID READS ANALYSIS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Run non-hybrid reads analysis only in centromere-aware mode (requires CEN/ARMS regions)
// In centromere-unaware mode (CHR only), skip this analysis as it's hardcoded for CEN/ARMS categories
if (centromere_aware_mode) {
    non_hybrid_analysis_output = NON_HYBRID_READS_ANALYSIS(
        process_cenhapmer_counts_output.kmer_matrix,              // TSV file with kmer counts per read (from process_cenhapmer_counts)
        processed_data_readmer_ch.map { meta, file -> file }.first(),  // Main FASTA file (non-masked)
        ch_reference_genomes.first()                               // Reference genomes directory (reuse from mapping analysis)
    )
} else {
    // Create empty channels for centromere-unaware mode
    log.info "⏩ Skipping non-hybrid reads analysis (only available in centromere-aware mode)"
    non_hybrid_analysis_output = [
        output_directory: Channel.empty(),
        read_lists: Channel.empty(),
        extracted_sequences: Channel.empty(),
        bam_files: Channel.empty(),
        bai_files: Channel.empty(),
        indel_results: Channel.empty()
    ]
}



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
    chr_candidates = segment_reads_output.chr_candidates
    arms_fasta = segment_reads_output.arms_fasta
    cen_fasta = segment_reads_output.cen_fasta
    chr_fasta = segment_reads_output.chr_fasta
    bed_files = segment_reads_output.bed_files
    segmented_fasta = segment_reads_output.segmented_fasta

    // Mapping outputs
    arms_mapping = mapping_analysis_output.arms_mapping_results
    cen_mapping = mapping_analysis_output.cen_mapping_results
    chr_mapping = mapping_analysis_output.chr_mapping_results
    plotting_results = mapping_analysis_output.plotting_results
    paf_files = mapping_analysis_output.paf_files
    cowidth_files = mapping_analysis_output.cowidth_files

    // Plotting output (available in both centromere-aware and centromere-unaware modes)
    crossover_plots = crossover_plot_output.crossover_plots
    coverage_data = crossover_plot_output.coverage_data
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
