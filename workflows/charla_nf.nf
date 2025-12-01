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
// REMOVED: V2 confidence scoring - replaced by V3
// include { CALCULATE_CONFIDENCE_SCORES_V2 } from '../modules/local/calculate_confidence_scores_v2.nf'
include { CALCULATE_CONFIDENCE_V3 } from '../modules/local/calculate_confidence_v3.nf'
include { ANALYZE_SCO_NCO } from '../modules/local/analyze_sco_nco.nf'
include { VISUALIZE_READ_STRUCTURES } from '../modules/local/visualize_read_structures.nf'
include { VISUALIZE_SEGMENTS } from '../modules/local/visualize_segments.nf'
// REMOVED: Old confidence plotting module (not needed - v2 generates its own diagnostic plots)
// include { PLOT_CONFIDENCE_SCORES } from '../modules/local/plot_confidence_scores.nf'
include { PLOT_CROSSOVER_MAP } from '../modules/local/plot_crossover_map.nf'
include { NON_HYBRID_MAPPING } from '../modules/local/non_hybrid_mapping.nf'
include { NON_HYBRID_ANALYSIS } from '../modules/local/non_hybrid_analysis.nf'
include { LARGE_INDEL_ANALYSIS } from '../modules/local/large_indel_analysis.nf'
include { ANALYZE_INDEL_MONOMER_BOUNDARIES } from '../modules/local/analyze_indel_monomer_boundaries.nf'
// REMOVED: MAFFT modules (replaced with indel monomer analysis)
// include { MAFFT_ALIGN } from '../modules/local/mafft/align.nf'
// include { ANALYZE_MAFFT_ALIGNMENT } from '../modules/local/analyze_mafft_alignment.nf'

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
  fasta_ch = fastq_to_fasta_ch.fasta       // main FASTA output

    masked_fastq_ch = fastq_to_fasta_ch.masked_fasta   // quality-masked FASTQ


//fasta_ch = fastq_to_fasta_ch.fasta.map { meta, files ->
    // if multiple .fa files are emitted, keep only the non-masked one
  //  def non_masked = files.find { !it.name.endsWith("_quality_masked.fa") }
  //     println "[DEBUG] non_masked = $non_masked" 
  // tuple(meta, non_masked)
//}


//    masked_fastq_ch = fastq_to_fasta_ch.masked_fasta   // quality-masked FASTQ

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

  //  ch_script_file = Channel.fromPath("${baseDir}/rust/kmer_processor/target/release/kmer_processor")
    ch_fasta_for_postproc = processed_fasta_files.map { meta, fasta -> fasta }.first()
    
    // Now, call process_cenhapmer_counts with the directory path
    process_cenhapmer_counts_output = process_cenhapmer_counts(
        all_cenhapmer_done,    // Directory path (e.g., "results/cenhapmers")
        ch_fasta_for_postproc,    // The main FASTA file for seqkit stats
     //   ch_script_file,          // The Python script file
        params.sample_id          // Use sample_id as the output_name for the final matrix
    )  

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Branch C: COMBINE HYBRID PROFILES
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Branch C: COMBINE HYBRID PROFILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// 1. Validate index_dir parameter
if (!params.index_dir) {
    error "ERROR: The --index_dir parameter must be provided for hybrid profile combination."
}

// Remove trailing slash from outdir
def outdir_clean = params.outdir.replaceAll('/$', '')

// Build correct path to cenhapmers directory
// If outdir is absolute, use it directly.
// If relative, prepend workflow.launchDir.
def cenhapmer_dir_path = params.outdir.startsWith('/') ?
    "${outdir_clean}/cenhapmers" :
    "${workflow.launchDir}/${outdir_clean}/cenhapmers"

log.info "📁 Cenhapmer input directory resolved to: ${cenhapmer_dir_path}"

// Create channels
ch_input_dir_for_hybrid = Channel.value(cenhapmer_dir_path)
ch_kmer_tsv             = process_cenhapmer_counts_output.kmer_matrix
ch_qc_file              = get_counts_kmc.out
ch_rust_binary          = Channel.fromPath("${baseDir}/rust/combine_segment_fast_indixes/target/release/hybrid_profile_toolkit")

// Run process
combine_hybrid_profiles_output = COMBINE_HYBRID_PROFILES(
    ch_input_dir_for_hybrid,
    ch_kmer_tsv,
    ch_qc_file,
    params.min_count
)

// Cleanup after completion
combine_hybrid_profiles_output.hybrid_profiles_file.subscribe { hybrid_file ->
    def cleanup_dir = file(cenhapmer_dir_path)

    if (cleanup_dir.exists()) {
        log.info "🧹 Cleaning cenhapmers directory: ${cleanup_dir}"
        if (cleanup_dir.deleteDir()) {
            log.info "✅ Deleted cenhapmers directory"
        } else {
            log.warn "⚠ Failed to delete cenhapmers directory"
        }
    } else {
        log.warn "⚠ Cleanup skipped: cenhapmers directory does not exist at ${cleanup_dir}"
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
        params.kmer_size,
        params.min_segment_length ?: 50,  // Minimum segment length in bp
        params.min_token_count ?: 50      // Minimum k-mer token count
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
        CONFIDENCE SCORING
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Step 1: Visualize read structures with parental haplotypes (doesn't need confidence scores)
    read_structure_output = VISUALIZE_READ_STRUCTURES(
        mapping_analysis_output.arms_mapping_results.first(),  // ARMS mapping directory (contains best_alignments.paf)
        mapping_analysis_output.cen_mapping_results.first(),   // CEN mapping directory (contains best_alignments.paf)
        file('NO_FILE'),                                       // No confidence scores needed for initial run
        params.max_structure_reads ?: 200
    )

    // Step 2: Calculate confidence scores V3 (uses parental proportions + curated profiles for TRUE classification)
    confidence_v3_output = CALCULATE_CONFIDENCE_V3(
        read_structure_output.proportions_table.first(),
        curate_hybrid_profiles_output.cut_table.first(),  // SOURCE OF TRUTH for classification (ultracurated_cut)
        mapping_analysis_output.arms_summary.first(),
        mapping_analysis_output.cen_summary.first(),
        arms_bed_files.first(),
        cen_bed_files.first(),
        params.threshold ?: 50
    )

    // Step 3: Filter segments BED file for SCO/NCO analysis
    segments_bed_for_analysis = segment_reads_output.bed_files
        .flatten()
        .filter { it.name.contains('segments_') && it.name.contains('_filtered.bed') }
        .first()

    // Step 3: Analyze SCO vs NCO comparison (uses V3 confidence scores)
    sco_nco_analysis_output = ANALYZE_SCO_NCO(
        confidence_v3_output.confidence_scores.first(),
        mapping_analysis_output.arms_summary.first(),  // Use summary as raw metrics
        mapping_analysis_output.cowidth_files.flatten()
            .filter { it.name.contains('ARMS') }
            .first(),
        mapping_analysis_output.cowidth_files.flatten()
            .filter { it.name.contains('CEN') }
            .first(),
        segments_bed_for_analysis  // Add BED file input
    )

    // Step 4: Visualize segment structures from BED files (shows unassigned regions)
    // Filter for the main segments BED file (not ARMS/CEN specific ones)
    segments_bed = segment_reads_output.bed_files
        .flatten()
        .filter { it.name.contains('segments_') && it.name.contains('_filtered.bed') }
        .first()

    segment_visualization_output = VISUALIZE_SEGMENTS(
        segments_bed,
        confidence_v3_output.confidence_scores.first(),
        params.max_structure_reads ?: 200
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


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NON-HYBRID READS MAPPING (Step 1: Time-consuming alignment)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Run non-hybrid reads mapping - classify, extract, and map sequences
// This step is time-consuming and will be cached by Nextflow
non_hybrid_mapping_output = NON_HYBRID_MAPPING(
    process_cenhapmer_counts_output.kmer_matrix,              // TSV file with kmer counts per read (from process_cenhapmer_counts)
    processed_data_readmer_ch.map { meta, file -> file }.first(),  // Main FASTA file (non-masked)
    ch_reference_genomes.first()                               // Reference genomes directory (reuse from mapping analysis)
)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NON-HYBRID READS ANALYSIS (Step 2: Fast indel analysis and plotting)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Run non-hybrid reads analysis - analyze indels and create plots
// This step is fast and can be re-run independently when plots are modified
non_hybrid_analysis_output = NON_HYBRID_ANALYSIS(
    non_hybrid_mapping_output.output_directory                 // Directory with BAM files from mapping step
)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LARGE INDEL ANALYSIS (>50bp) - Satellite Structure Investigation (OPTIONAL)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.run_large_indel_analysis) {
    // Run large indel analysis to investigate 178bp satellite patterns
    large_indel_analysis_output = LARGE_INDEL_ANALYSIS(
        non_hybrid_mapping_output.output_directory,                // Directory from NON_HYBRID_MAPPING (contains BAM files)
        processed_data_readmer_ch.map { meta, file -> file }.first(),  // Main FASTA file (non-masked)
        ch_reference_genomes.first(),                              // Reference genomes directory
        ch_col_bed.first(),                                    // Parent 1 BED file for chromosome reconstruction
        ch_ler_bed.first()                                     // Parent 2 BED file for chromosome reconstruction
    )

    // Analyze large indel boundaries relative to the 178bp CEN178 satellite monomer
    // This replaces the MAFFT alignment which was found to be less useful
    indel_monomer_analysis_output = ANALYZE_INDEL_MONOMER_BOUNDARIES(
        large_indel_analysis_output.indel_sequences,       // FASTA of indel sequences
        large_indel_analysis_output.large_indels_catalog,  // TSV catalog
        ch_col_bed.first(),                                 // BED with centromere coords (using Col as reference)
        params.threads ?: 8                                 // Number of threads
    )
} else {
    // Create empty channels when indel analysis is disabled
    large_indel_analysis_output = [
        output_directory: Channel.empty(),
        large_indels_catalog: Channel.empty(),
        indel_sequences: Channel.empty(),
        mapping_files: Channel.empty(),
        genomic_plots: Channel.empty(),
        satellite_blocks: Channel.empty(),
        satellite_analysis: Channel.empty(),
        structure_plots: Channel.empty()
    ]
    indel_monomer_analysis_output = [
        output_directory: Channel.empty(),
        metaprofile_plot: Channel.empty(),
        inframe_plot: Channel.empty(),
        inframe_results: Channel.empty(),
        coverage_plot: Channel.empty(),
        comparison_plot: Channel.empty(),
        clustering_plot: Channel.empty(),
        monomer_variants: Channel.empty(),
        heatmap_plot: Channel.empty(),
        alignment_files: Channel.empty()
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

    // Confidence scoring outputs (V3 - Redesigned with proper metrics)
    confidence_scores = confidence_v3_output.confidence_scores
    confidence_diagnostic_plots = confidence_v3_output.diagnostic_plots
    confidence_high_conf_reads = confidence_v3_output.high_conf_reads
    confidence_high_sco_reads = confidence_v3_output.high_conf_sco_reads

    // SCO vs NCO analysis outputs
    sco_nco_comparison_plot = sco_nco_analysis_output.comparison_plot
    sco_nco_resolution_plot = sco_nco_analysis_output.resolution_plot
    sco_nco_summary = sco_nco_analysis_output.summary_table

    // Read structure visualization outputs
    read_structures_plot = read_structure_output.structures_plot
    read_parental_plot = read_structure_output.parental_plot
    read_structure_summary = read_structure_output.summary
    parental_proportions = read_structure_output.proportions_table
    read_structures_sco = read_structure_output.sco_structures
    read_structures_nco = read_structure_output.nco_structures


// Non-hybrid mapping outputs (from mapping step)
    nonhybrid_mapping_directory = non_hybrid_mapping_output.output_directory
    nonhybrid_read_lists = non_hybrid_mapping_output.read_lists
    nonhybrid_sequences = non_hybrid_mapping_output.extracted_sequences

// Non-hybrid analysis outputs (from analysis step)
    nonhybrid_analysis_directory = non_hybrid_analysis_output.output_directory
    indel_analysis = non_hybrid_analysis_output.indel_results
    indel_plots = non_hybrid_analysis_output.plots

    // Large indel analysis outputs
    large_indel_directory = large_indel_analysis_output.output_directory
    large_indels_catalog = large_indel_analysis_output.large_indels_catalog
    large_indel_sequences = large_indel_analysis_output.indel_sequences
    large_indel_mappings = large_indel_analysis_output.mapping_files
    large_indel_genomic_plots = large_indel_analysis_output.genomic_plots
    large_indel_satellite_blocks = large_indel_analysis_output.satellite_blocks
    large_indel_satellite_analysis = large_indel_analysis_output.satellite_analysis
    large_indel_structure_plots = large_indel_analysis_output.structure_plots

    // Indel monomer boundaries analysis outputs (replaces MAFFT)
    indel_monomer_directory = indel_monomer_analysis_output.output_directory
    indel_monomer_metaprofile = indel_monomer_analysis_output.metaprofile_plot
    indel_monomer_inframe_plot = indel_monomer_analysis_output.inframe_plot
    indel_monomer_inframe_results = indel_monomer_analysis_output.inframe_results
    indel_monomer_coverage_plot = indel_monomer_analysis_output.coverage_plot
    indel_monomer_comparison_plot = indel_monomer_analysis_output.comparison_plot
    indel_monomer_clustering_plot = indel_monomer_analysis_output.clustering_plot
    indel_monomer_variants = indel_monomer_analysis_output.monomer_variants
    indel_monomer_heatmap = indel_monomer_analysis_output.heatmap_plot
    indel_monomer_alignments = indel_monomer_analysis_output.alignment_files

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
*/
