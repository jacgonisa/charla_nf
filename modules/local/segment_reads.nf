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
    val centromere_aware        // Boolean: true for CEN/ARMS mode, false for CHR-only mode    

    output:
    path "threshold*/segments_*_filtered.tsv", emit: filtered_table
    path "threshold*/CANDIDATE_reads_nonectopic", emit: nonectopic_candidates
    path "threshold*/CANDIDATE_reads_nonectopic_ARMS", emit: arms_candidates, optional: true
    path "threshold*/CANDIDATE_reads_nonectopic_CEN", emit: cen_candidates, optional: true
    path "threshold*/CANDIDATE_reads_nonectopic_CHR", emit: chr_candidates, optional: true
    path "threshold*/CANDIDATE_reads_nonectopic_ARMS.txt", emit: arms_readids, optional: true
    path "threshold*/CANDIDATE_reads_nonectopic_CEN.txt", emit: cen_readids, optional: true
    path "threshold*/CANDIDATE_reads_nonectopic_CHR.txt", emit: chr_readids, optional: true
    path "threshold*/CANDIDATE_reads_nonectopic_ARMS.fa", emit: arms_fasta, optional: true
    path "threshold*/CANDIDATE_reads_nonectopic_CEN.fa", emit: cen_fasta, optional: true
    path "threshold*/CANDIDATE_reads_nonectopic_CHR.fa", emit: chr_fasta, optional: true
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

    # Mode detection: centromere-aware (CEN/ARMS) vs centromere-unaware (CHR)
    if [[ "${centromere_aware}" == "true" ]]; then
        echo "Mode: Centromere-aware (CEN/ARMS splitting)"
        # Separate ARMS and CEN candidates (pattern: ?A?=ARMS, ?C?=CEN)
        grep -E '[A-Z]A[0-9]' CANDIDATE_reads_nonectopic > CANDIDATE_reads_nonectopic_ARMS || touch CANDIDATE_reads_nonectopic_ARMS
        grep -E '[A-Z]C[0-9]' CANDIDATE_reads_nonectopic > CANDIDATE_reads_nonectopic_CEN || touch CANDIDATE_reads_nonectopic_CEN

        # Extract read IDs
        cut -f1 CANDIDATE_reads_nonectopic_CEN > CANDIDATE_reads_nonectopic_CEN.txt || touch CANDIDATE_reads_nonectopic_CEN.txt
        cut -f1 CANDIDATE_reads_nonectopic_ARMS > CANDIDATE_reads_nonectopic_ARMS.txt || touch CANDIDATE_reads_nonectopic_ARMS.txt

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
    else
        echo "Mode: Centromere-unaware (CHR-only, no splitting)"
        # In CHR mode, all candidates are CHR regions
        cp CANDIDATE_reads_nonectopic CANDIDATE_reads_nonectopic_CHR

        # Extract read IDs
        cut -f1 CANDIDATE_reads_nonectopic_CHR > CANDIDATE_reads_nonectopic_CHR.txt || touch CANDIDATE_reads_nonectopic_CHR.txt

        # Extract FASTA sequences using seqkit
        if [[ -s CANDIDATE_reads_nonectopic_CHR.txt ]]; then
            seqkit grep -f CANDIDATE_reads_nonectopic_CHR.txt ../${input_fasta} | seqkit seq -w 0 > CANDIDATE_reads_nonectopic_CHR.fa
        else
            touch CANDIDATE_reads_nonectopic_CHR.fa
        fi
    fi
     


 
    # Step 8.3: Generate BED files
    echo "Step 8.3: Generating BED files..."

    # BED file for all filtered reads
    ${projectDir}/rust/bed_processor/target/release/bed_processor \\
        ../${hybrid_profiles_file} \\
        \${output_prefix}_filtered.tsv \\
        \${output_prefix}_filtered.bed ${kmer_size}

    # Generate BED files based on mode
    if [[ "${centromere_aware}" == "true" ]]; then
        # BED file for CEN candidates (if file exists and not empty)
        if [[ -s CANDIDATE_reads_nonectopic_CEN ]]; then
            ${projectDir}/rust/bed_processor/target/release/bed_processor \\
                ../${hybrid_profiles_file} \\
                CANDIDATE_reads_nonectopic_CEN \\
                CANDIDATE_reads_nonectopic_CEN.bed ${kmer_size}
        else
            touch CANDIDATE_reads_nonectopic_CEN.bed
        fi

        # BED file for ARMS candidates (if file exists and not empty)
        if [[ -s CANDIDATE_reads_nonectopic_ARMS ]]; then
            ${projectDir}/rust/bed_processor/target/release/bed_processor \\
                ../${hybrid_profiles_file} \\
                CANDIDATE_reads_nonectopic_ARMS \\
                CANDIDATE_reads_nonectopic_ARMS.bed ${kmer_size}
        else
            touch CANDIDATE_reads_nonectopic_ARMS.bed
        fi
    else
        # BED file for CHR candidates
        if [[ -s CANDIDATE_reads_nonectopic_CHR ]]; then
            ${projectDir}/rust/bed_processor/target/release/bed_processor \\
                ../${hybrid_profiles_file} \\
                CANDIDATE_reads_nonectopic_CHR \\
                CANDIDATE_reads_nonectopic_CHR.bed ${kmer_size}
        else
            touch CANDIDATE_reads_nonectopic_CHR.bed
        fi
    fi
    
    # Step 8.4: Segment reads using BED files
    echo "Step 8.4: Segmenting reads..."

    if [[ "${centromere_aware}" == "true" ]]; then
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
    else
        # Segment CHR reads (if FASTA and BED exist and are not empty)
        if [[ -s CANDIDATE_reads_nonectopic_CHR.fa && -s CANDIDATE_reads_nonectopic_CHR.bed ]]; then
            python ${projectDir}/scripts/08-segment_bed_into_fasta.py \\
                CANDIDATE_reads_nonectopic_CHR.fa \\
                CANDIDATE_reads_nonectopic_CHR.bed \\
                CANDIDATE_reads_nonectopic_CHR.segments.fa
        else
            touch CANDIDATE_reads_nonectopic_CHR.segments.fa
        fi

        echo "Read segmentation completed!"
        echo "Summary:"
        echo "- Filtered reads: \$(wc -l < \${output_prefix}_filtered.tsv)"
        echo "- Non-ectopic candidates: \$(wc -l < CANDIDATE_reads_nonectopic)"
        echo "- CHR candidates: \$(wc -l < CANDIDATE_reads_nonectopic_CHR)"
    fi
    """
}
