nextflow.enable.dsl=2

process SEGMENT_READS {
    tag "segment_reads"
    label 'segment_processing'
    cpus params.segment_cpus ?: 4
    memory params.segment_memory ?: '16 GB'
    time params.segment_time ?: '4h'
    
    publishDir "${params.outdir}/read_segments", mode: 'copy'
    
    input:
    path curated_table           // Final curated table from CURATE_HYBRID_PROFILES
    path hybrid_profiles_file    // Original hybrid profiles file
    path input_fasta            // Original FASTA file
    val threshold               // Threshold value (e.g., 50)
    val kmer_size    

    output:
    path "threshold*/segments_*_filtered.tsv", emit: filtered_table
    path "threshold*/CANDIDATE_reads_nonectopic", emit: nonectopic_candidates
    path "threshold*/CANDIDATE_reads_nonectopic_ARMS", emit: arms_candidates
    path "threshold*/CANDIDATE_reads_nonectopic_CEN", emit: cen_candidates
    path "threshold*/CANDIDATE_reads_nonectopic_ARMS.txt", emit: arms_readids
    path "threshold*/CANDIDATE_reads_nonectopic_CEN.txt", emit: cen_readids
    path "threshold*/CANDIDATE_reads_nonectopic_ARMS.fa", emit: arms_fasta
    path "threshold*/CANDIDATE_reads_nonectopic_CEN.fa", emit: cen_fasta
    path "threshold*/*.bed", emit: bed_files
    path "threshold*/*.segments.fa", emit: segmented_fasta
    
    script:
    """
    # Define variables inside script block
    output_prefix="segments_${params.sample_id}_thr${threshold}"
    threshold_dir="threshold${threshold}"
    
    # Create threshold directory
    mkdir -p \${threshold_dir}
    cd \${threshold_dir}
    
    # Step 8.1: Filter reads by threshold and criteria
    echo "Step 8.1: Filtering reads with threshold ${threshold}..."
    awk -F'\\t' '\$2 == ${threshold} && \$5 == "yes|2"' ../${curated_table} > \${output_prefix}_filtered.tsv
    
    # Step 8.2: Get non-ectopic candidates
    echo "Step 8.2: Getting non-ectopic candidates..."
    awk -F'\\t' '\$7 == "no"' \${output_prefix}_filtered.tsv > CANDIDATE_reads_nonectopic
    
    # Separate ARMS and CEN candidates
    grep CA CANDIDATE_reads_nonectopic > CANDIDATE_reads_nonectopic_ARMS || touch CANDIDATE_reads_nonectopic_ARMS
    grep CC CANDIDATE_reads_nonectopic > CANDIDATE_reads_nonectopic_CEN || touch CANDIDATE_reads_nonectopic_CEN
    
    # Extract read IDs
    cut -f1 CANDIDATE_reads_nonectopic_CEN > CANDIDATE_reads_nonectopic_CEN.txt
    cut -f1 CANDIDATE_reads_nonectopic_ARMS > CANDIDATE_reads_nonectopic_ARMS.txt
    
    # Extract FASTA sequences using seqkit
    if [[ -s CANDIDATE_reads_nonectopic_CEN.txt ]]; then
    seqkit grep -f CANDIDATE_reads_nonectopic_CEN.txt ../${input_fasta} | seqkit seq -w 0 > CANDIDATE_reads_nonectopic_CEN.fa
else
    touch CANDIDATE_reads_nonectopic_CEN.fa
fi

if [[ -s CANDIDATE_reads_nonectopic_ARMS.txt ]]; then
    seqkit grep -f CANDIDATE_reads_nonectopic_ARMS.txt ../${input_fasta} | seqkit seq -w 0 > CANDIDATE_reads_nonectopic_ARMS.fa
else
    touch CANDIDATE_reads_nonectopic_ARMS.fa
fi
     


 
    # Step 8.3: Generate BED files
    echo "Step 8.3: Generating BED files..."
    
    # BED file for all filtered reads
 #   python ${projectDir}/scripts/08-getting_bed_read_merged_nosmallsegments_allmerged_thr${threshold}.py \\
 #       ../${hybrid_profiles_file} \\
 #       \${output_prefix}_filtered.tsv \\
 #       \${output_prefix}_filtered.bed ${kmer_size}
   
   ${projectDir}/rust/bed_processor/target/release/bed_processor \\
      ../${hybrid_profiles_file} \\
       \${output_prefix}_filtered.tsv \\
       \${output_prefix}_filtered.bed ${kmer_size}

    
  
    # BED file for CEN candidates (if file exists and not empty)
    if [[ -s CANDIDATE_reads_nonectopic_CEN ]]; then
        #python ${projectDir}/scripts/08-getting_bed_read_merged_nosmallsegments_allmerged_thr${threshold}.py \\
        #    ../${hybrid_profiles_file} \\
        #    CANDIDATE_reads_nonectopic_CEN \\
        #    CANDIDATE_reads_nonectopic_CEN.bed ${kmer_size}

         ${projectDir}/rust/bed_processor/target/release/bed_processor \\
      	 	../${hybrid_profiles_file} \\
		CANDIDATE_reads_nonectopic_CEN \\
                CANDIDATE_reads_nonectopic_CEN.bed ${kmer_size}
    

    else
        touch CANDIDATE_reads_nonectopic_CEN.bed
    fi
    
    # BED file for ARMS candidates (if file exists and not empty)
    if [[ -s CANDIDATE_reads_nonectopic_ARMS ]]; then
        #python ${projectDir}/scripts/08-getting_bed_read_merged_nosmallsegments_allmerged_thr${threshold}.py \\
        #    ../${hybrid_profiles_file} \\
        #    CANDIDATE_reads_nonectopic_ARMS \\
        #    CANDIDATE_reads_nonectopic_ARMS.bed ${kmer_size}

     ${projectDir}/rust/bed_processor/target/release/bed_processor \\
                ../${hybrid_profiles_file} \\
                CANDIDATE_reads_nonectopic_ARMS \\
                CANDIDATE_reads_nonectopic_ARMS.bed ${kmer_size}

    else
        touch CANDIDATE_reads_nonectopic_ARMS.bed
    fi
    
    # Step 8.4: Segment reads using BED files
    echo "Step 8.4: Segmenting reads..."
    
    # Segment CEN reads (if FASTA and BED exist and are not empty)
    if [[ -s CANDIDATE_reads_nonectopic_CEN.fa && -s CANDIDATE_reads_nonectopic_CEN.bed ]]; then
        python ${projectDir}/scripts/08-segment_bed_into_fasta.py \\
            CANDIDATE_reads_nonectopic_CEN.fa \\
            CANDIDATE_reads_nonectopic_CEN.bed \\
            CANDIDATE_reads_nonectopic_CEN.segments.fa
    else
        touch CANDIDATE_reads_nonectopic_CEN.segments.fa
    fi
    
    # Segment ARMS reads (if FASTA and BED exist and are not empty)
    if [[ -s CANDIDATE_reads_nonectopic_ARMS.fa && -s CANDIDATE_reads_nonectopic_ARMS.bed ]]; then
        python ${projectDir}/scripts/08-segment_bed_into_fasta.py \\
            CANDIDATE_reads_nonectopic_ARMS.fa \\
            CANDIDATE_reads_nonectopic_ARMS.bed \\
            CANDIDATE_reads_nonectopic_ARMS.segments.fa
    else
        touch CANDIDATE_reads_nonectopic_ARMS.segments.fa
    fi
    
    echo "Read segmentation completed!"
    echo "Summary:"
    echo "- Filtered reads: \$(wc -l < \${output_prefix}_filtered.tsv)"
    echo "- Non-ectopic candidates: \$(wc -l < CANDIDATE_reads_nonectopic)"
    echo "- ARMS candidates: \$(wc -l < CANDIDATE_reads_nonectopic_ARMS)"
    echo "- CEN candidates: \$(wc -l < CANDIDATE_reads_nonectopic_CEN)"
    """
}
