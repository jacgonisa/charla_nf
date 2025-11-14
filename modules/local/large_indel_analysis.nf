nextflow.enable.dsl=2

process LARGE_INDEL_ANALYSIS {
    tag "large_indel_analysis"
    label 'indel_analysis'
    cpus params.large_indel_cpus ?: 8
    memory params.large_indel_memory ?: '16 GB'
    time params.large_indel_time ?: '4h'

    publishDir "${params.outdir}/large_indel_analysis", mode: 'symlink'

    input:
    path(nonhybrid_dir)              // Directory from NON_HYBRID_READS_ANALYSIS
    path(main_fasta_file)            // Main FASTA file with all reads
    path(reference_genomes_dir)      // Directory with reference genomes
    path(parent1_bed)                // Parent 1 BED file for chromosome structure
    path(parent2_bed)                // Parent 2 BED file for chromosome structure

    output:
    path("${params.sample_id}_large_indels/"), emit: output_directory
    path("${params.sample_id}_large_indels/large_indels_catalog.tsv"), emit: large_indels_catalog
    path("${params.sample_id}_large_indels/indel_sequences.fa"), emit: indel_sequences, optional: true
    path("${params.sample_id}_large_indels/indel_mappings/*.paf"), emit: mapping_files, optional: true
    path("${params.sample_id}_large_indels/indel_mappings/*.png"), emit: genomic_plots, optional: true
    path("${params.sample_id}_large_indels/msa_results/satellite_blocks.fa"), emit: satellite_blocks, optional: true
    path("${params.sample_id}_large_indels/satellite_analysis.tsv"), emit: satellite_analysis
    path("${params.sample_id}_large_indels/satellite_structure_analysis.png"), emit: structure_plots, optional: true

    script:
    """
    # Create necessary directories
    mkdir -p ${params.sample_id}_large_indels/indel_mappings
    mkdir -p ${params.sample_id}_large_indels/msa_results

    # Step 1: Extract large indels (>50bp) and their sequences
    echo "=== Step 1: Extracting large indels (>50bp) and their sequences ==="
    python ${projectDir}/scripts/13.1-extract_large_indels.py \\
        ${nonhybrid_dir} \\
        ${main_fasta_file} \\
        ${params.sample_id}_large_indels \\
        50

    # Check if any large indels were found
    if [ ! -f "${params.sample_id}_large_indels/indel_sequences.fa" ] || [ ! -s "${params.sample_id}_large_indels/indel_sequences.fa" ]; then
        echo "No large indels found. Creating empty results files."
        touch ${params.sample_id}_large_indels/large_indels_catalog.tsv
        echo -e "indel_id\\tsize\\tmultiple_of_178\\tsatellite_blocks\\talignment_quality" > ${params.sample_id}_large_indels/satellite_analysis.tsv
        exit 0
    fi

    # Step 2: Map indel sequences back to reference genomes
    echo "=== Step 2: Mapping indel sequences to reference genomes ==="
    python ${projectDir}/scripts/13.2-map_indels_to_genome.py \\
        ${params.sample_id}_large_indels/indel_sequences.fa \\
        ${params.sample_id}_large_indels/large_indels_catalog.tsv \\
        ${nonhybrid_dir} \\
        ${reference_genomes_dir} \\
        ${params.sample_id}_large_indels/indel_mappings \\
        ${task.cpus} \\
        ${parent1_bed} \\
        ${parent2_bed}

    # Step 3: Analyze indel sequences for 178bp satellite blocks and prepare blocks for MSA
    echo "=== Step 3: Analyzing indel sequences for 178bp satellite blocks ==="
    python ${projectDir}/scripts/13.3-analyze_satellite_structure.py \\
        ${params.sample_id}_large_indels/indel_sequences.fa \\
        ${params.sample_id}_large_indels/large_indels_catalog.tsv \\
        ${params.sample_id}_large_indels/msa_results \\
        ${params.sample_id}_large_indels/satellite_analysis.tsv \\
        ${task.cpus}

    echo "=== Large indel analysis complete ==="
    echo "Note: MAFFT alignment will be performed as a separate module in the workflow"
    """
}
