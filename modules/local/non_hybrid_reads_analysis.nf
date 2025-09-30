nextflow.enable.dsl=2

process NON_HYBRID_READS_ANALYSIS {
    tag "non_hybrid_analysis"
    label 'mapping_analysis'
    cpus params.nonhybrid_cpus ?: 8
    memory params.nonhybrid_memory ?: '16 GB'
    time params.nonhybrid_time ?: '4h'
    
    publishDir "${params.outdir}/non_hybrid_analysis", mode: 'symlink'
    
    input:
    path(kmer_matrix_file)        // TSV file with kmer counts per read (from process_cenhapmer_counts)
    path(main_fasta_file)         // Main FASTA file with all reads
    path(reference_genomes_dir)   // Directory with reference genomes
 
    output:
    path("${params.sample_id}_nonhybrid/"), emit: output_directory
    path("${params.sample_id}_nonhybrid/*_reads.txt"), emit: read_lists
    path("${params.sample_id}_nonhybrid/*_sequences.fa"), emit: extracted_sequences, optional: true
    path("${params.sample_id}_nonhybrid/mappings/*.bam"), emit: bam_files, optional: true
    path("${params.sample_id}_nonhybrid/mappings/*.bam.bai"), emit: bai_files, optional: true
    path("${params.sample_id}_nonhybrid/indel_analysis_results.tsv"), emit: indel_results
    
    script:
    """
    # Create necessary directories using sample ID


    
    # Step 1: Classify and extract reads
    echo "=== Step 1: Classifying and extracting non-hybrid reads ==="
    python ${projectDir}/scripts/12-classify_and_extract_reads.py \\
        ${kmer_matrix_file} \\
        ${main_fasta_file} \\
        ${params.sample_id}_nonhybrid \\
        ${params.kmer_size}
    
    # Check if any sequences were extracted
    if [ ! -d "${params.sample_id}_nonhybrid" ] || [ -z "\$(find ${params.sample_id}_nonhybrid -name '*_sequences.fa' -size +0c 2>/dev/null)" ]; then
        echo "No sequences extracted. Creating empty results file."
        echo -e "haplotype\\tmapping_mode\\ttotal_reads\\tmapped_reads\\tindels\\tavg_indel_size\\tindel_density" > ${params.sample_id}/indel_analysis_results.tsv
        exit 0
    fi
    
    # Step 2: Map sequences to reference genomes
    echo "=== Step 2: Mapping sequences to reference genomes ==="
    python ${projectDir}/scripts/12.1-mapping_analysis.py \\
        ${params.sample_id}_nonhybrid \\
        ${task.cpus}
    
    # Step 3: Analyze indels from BAM files
    echo "=== Step 3: Analyzing indels from alignments ==="
    python ${projectDir}/scripts/12.2-indel_analysis.py \\
        ${params.sample_id}_nonhybrid \\
        ${params.sample_id}_nonhybrid/indel_analysis_results.tsv
    
    echo "=== Non-hybrid reads analysis complete ==="
    """
}
