/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHARLA Mapping Workflow
    Based on CHARLA_NF modules for mapping and indel analysis
    (without recombination analysis)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQ_TO_FASTA         } from '../modules/local/fastq_to_fasta'
include { SIMPLIFY_HEADERS       } from '../modules/local/simplify_headers'
include { GENERATE_READMERS_KMC  } from '../modules/local/generate_readmers_kmc'
include { GENERATE_HISTOGRAM_KMC } from '../modules/local/generate_histogram_kmc'
include { GET_COUNTS_KMC         } from '../modules/local/get_counts_kmc'
include { GET_COUNTS_CENHAPMER   } from '../modules/local/get_counts_cenhapmer_kmc'
include { PROCESS_CENHAPMER_COUNTS } from '../modules/local/process_cenhapmer_counts'
include { NON_HYBRID_READS_ANALYSIS } from '../modules/local/non_hybrid_reads_analysis'
include { LARGE_INDEL_ANALYSIS   } from '../modules/local/large_indel_analysis'
include { MAFFT_ALIGN            } from '../modules/local/mafft/align'
include { ANALYZE_MAFFT_ALIGNMENT } from '../modules/local/analyze_mafft_alignment'

workflow CHARLA_MAPPING {

    take:
    ch_input          // tuple val(meta), path(input_file)
    ch_reference      // path(reference_genome)
    ch_reference_dir  // path(reference_directory)
    ch_kmer_db        // path(kmer_database)

    main:

    //
    // Step 1: Convert FASTQ to FASTA if needed
    //
    if (params.input_format == 'fastq') {
        FASTQ_TO_FASTA(ch_input)
        ch_fasta = FASTQ_TO_FASTA.out.fasta
    } else {
        ch_fasta = ch_input
    }

    //
    // Step 2: Simplify FASTA headers for downstream processing
    //
    SIMPLIFY_HEADERS(ch_fasta)
    ch_simplified = SIMPLIFY_HEADERS.out.fasta

    //
    // Step 3: Generate readmers (k-mers from reads) using KMC
    //
    // Prepare input channel: tuple(fasta_file, sample_id, kmer_size)
    ch_readmer_input = ch_simplified.map { meta, fasta ->
        tuple(fasta, params.sample_id, params.kmer_size)
    }

    GENERATE_READMERS_KMC(ch_readmer_input)

    //
    // Step 4: Generate histogram from readmers
    //
    GENERATE_HISTOGRAM_KMC(GENERATE_READMERS_KMC.out)

    //
    // Step 5: Get k-mer counts per read (readmer profiling)
    //
    // Prepare input: tuple(kmc_pre, kmc_suf, fastq, sample_id, kmer_size)
    ch_get_counts_input = GENERATE_READMERS_KMC.out
        .combine(ch_simplified.map { meta, fasta -> fasta })

    GET_COUNTS_KMC(ch_get_counts_input)

    //
    // Step 6: Get cenhapmer k-mer counts per read
    // Create channel for each accession/category/chromosome combination
    //
    // Define accessions and categories based on parent names
    def accessions = [params.parent1_name]
    if (params.parent2_name) {
        accessions << params.parent2_name
    }
    def categories = ['ARMS', 'CEN']
    def chromosomes = 1..params.num_chromosomes

    // Create combinations channel
    ch_cenhapmer_combinations = Channel.of(accessions)
        .flatMap { acc -> categories.collect { cat -> [acc, cat] } }
        .flatMap { acc, cat -> chromosomes.collect { chr -> [acc, cat, chr] } }
        .combine(ch_simplified.map { meta, fasta -> fasta })
        .combine(ch_kmer_db)
        .map { acc, cat, chr, fasta, db -> tuple(acc, cat, chr, fasta, db) }

    GET_COUNTS_CENHAPMER(ch_cenhapmer_combinations)

    //
    // Step 7: Process cenhapmer k-mer counts into matrix format
    //
    // Create channels for PROCESS_CENHAPMER_COUNTS
    ch_cenhapmer_dir_path = Channel.value("${params.outdir}/cenhapmers")
    ch_fasta_for_postproc = ch_simplified.map { meta, fasta -> fasta }.first()
    ch_script_file = Channel.fromPath("${projectDir}/rust/kmer_processor/target/release/kmer_processor")

    // Wait for all cenhapmer counts to complete before processing
    GET_COUNTS_CENHAPMER.out.collect().map { it -> "done" }

    PROCESS_CENHAPMER_COUNTS(
        ch_cenhapmer_dir_path,
        ch_fasta_for_postproc,
        ch_script_file,
        params.sample_id
    )

    //
    // Step 8: Analyze non-hybrid reads and perform mapping
    //
    NON_HYBRID_READS_ANALYSIS(
        PROCESS_CENHAPMER_COUNTS.out.kmer_matrix,
        ch_simplified.map { meta, file -> file }.first(),
        ch_reference_dir.first()
    )

    //
    // Step 9: Analyze large indels (>50bp) and extract satellite structures
    //
    // Note: BED files for chromosomal features are optional.
    // If not provided, analysis will proceed without genomic context visualization
    //
    ch_parent1_bed = params.parent1_bed ?
        Channel.fromPath(params.parent1_bed, checkIfExists: true) :
        Channel.empty()
    ch_parent2_bed = params.parent2_bed ?
        Channel.fromPath(params.parent2_bed, checkIfExists: true) :
        Channel.empty()

    // If BED files are not provided, create empty channels
    if (!params.parent1_bed) {
        ch_parent1_bed = Channel.value(file('NO_FILE'))
    }
    if (!params.parent2_bed) {
        ch_parent2_bed = Channel.value(file('NO_FILE'))
    }

    LARGE_INDEL_ANALYSIS(
        NON_HYBRID_READS_ANALYSIS.out.output_directory,
        ch_simplified.map { meta, file -> file }.first(),
        ch_reference_dir.first(),
        ch_parent1_bed.first(),
        ch_parent2_bed.first()
    )

    //
    // Step 10: Multiple sequence alignment of satellite blocks with MAFFT
    //
    ch_satellite_blocks = LARGE_INDEL_ANALYSIS.out.satellite_blocks
        .map { blocks_file ->
            def meta = [id: "${params.sample_id}_satellite_blocks"]
            tuple(meta, blocks_file)
        }

    MAFFT_ALIGN(ch_satellite_blocks)

    //
    // Step 11: Analyze MAFFT alignment and update satellite analysis
    //
    ANALYZE_MAFFT_ALIGNMENT(
        MAFFT_ALIGN.out.fas,
        LARGE_INDEL_ANALYSIS.out.satellite_analysis
    )

    emit:
    // FASTA outputs
    simplified_fasta = SIMPLIFY_HEADERS.out.fasta

    // Readmer analysis outputs
    readmer_kmc_output = GENERATE_READMERS_KMC.out
    readmer_histogram = GENERATE_HISTOGRAM_KMC.out.histo
    readmer_counts = GET_COUNTS_KMC.out

    // Cenhapmer analysis outputs
    cenhapmer_counts = GET_COUNTS_CENHAPMER.out
    kmer_matrix = PROCESS_CENHAPMER_COUNTS.out.kmer_matrix

    // Non-hybrid reads analysis
    nonhybrid_directory = NON_HYBRID_READS_ANALYSIS.out.output_directory
    nonhybrid_read_lists = NON_HYBRID_READS_ANALYSIS.out.read_lists
    nonhybrid_sequences = NON_HYBRID_READS_ANALYSIS.out.extracted_sequences
    nonhybrid_bam_files = NON_HYBRID_READS_ANALYSIS.out.bam_files
    nonhybrid_bai_files = NON_HYBRID_READS_ANALYSIS.out.bai_files
    nonhybrid_indel_results = NON_HYBRID_READS_ANALYSIS.out.indel_results

    // Large indel analysis
    large_indel_directory = LARGE_INDEL_ANALYSIS.out.output_directory
    large_indels_catalog = LARGE_INDEL_ANALYSIS.out.large_indels_catalog
    large_indel_sequences = LARGE_INDEL_ANALYSIS.out.indel_sequences
    large_indel_mapping = LARGE_INDEL_ANALYSIS.out.mapping_files
    large_indel_satellite_blocks = LARGE_INDEL_ANALYSIS.out.satellite_blocks
    large_indel_satellite_analysis = LARGE_INDEL_ANALYSIS.out.satellite_analysis

    // MAFFT alignment
    mafft_alignment = MAFFT_ALIGN.out.fas
    mafft_analysis_updated = ANALYZE_MAFFT_ALIGNMENT.out.updated_analysis
    mafft_visualization = ANALYZE_MAFFT_ALIGNMENT.out.msa_plot
}
