nextflow.enable.dsl=2

process MAPPING_ANALYSIS {
    tag "mapping_analysis_thr${threshold}"
    label 'mapping_processing'
    cpus params.mapping_cpus ?: 16
    memory params.mapping_memory ?: '32 GB'
    time params.mapping_time ?: '8h'
    
    publishDir "${params.outdir}/mapping_analysis", mode: 'copy'
    
    input:
    path arms_segments_fasta     // ARMS segments FASTA from SEGMENT_READS
    path cen_segments_fasta      // CEN segments FASTA from SEGMENT_READS
    path arms_bed_file          // ARMS BED file from SEGMENT_READS
    path cen_bed_file           // CEN BED file from SEGMENT_READS
    path reference_genomes_dir   // Directory containing reference genomes
    val threshold               // Threshold value
    
    output:
    path "threshold${threshold}/mapping_full_read/${params.sample_id}/nonectopic/ARMS_segments_alnto${params.parent1_name}_alnto${params.parent2_name}/**", emit: arms_mapping_results, optional: true
    path "threshold${threshold}/mapping_full_read/${params.sample_id}/nonectopic/CEN_segments_alnto${params.parent1_name}_alnto${params.parent2_name}/**", emit: cen_mapping_results, optional: true
    path "threshold${threshold}/mapping_full_read/${params.sample_id}/nonectopic/CHR_segments_alnto${params.parent1_name}_alnto${params.parent2_name}/**", emit: chr_mapping_results, optional: true
    path "threshold${threshold}/plotting_crossover/${params.sample_id}/**", emit: plotting_results
    path "threshold${threshold}/mapping_full_read/${params.sample_id}/nonectopic/*_segments_alnto${params.parent1_name}_alnto${params.parent2_name}/summary_table.csv", emit: summary_tables, optional: true
    path "threshold${threshold}/plotting_crossover/${params.sample_id}/*.paf", emit: paf_files
    path "threshold${threshold}/plotting_crossover/${params.sample_id}/*.cowidths", emit: cowidth_files
    
    script:
    """
    # Create directory structure
    threshold_dir="threshold${threshold}"
    mapping_base="\${threshold_dir}/mapping_full_read/${params.sample_id}/nonectopic"
    plotting_base="\${threshold_dir}/plotting_crossover/${params.sample_id}"

    # Detect region types from input files
    if [[ "${arms_segments_fasta}" == *"CHR"* ]]; then
        REGION_TYPE="CHR"
        echo "Detected centromere-unaware mode (CHR regions)"
    else
        REGION_TYPE="ARMS"
        echo "Detected centromere-aware mode (ARMS/CEN regions)"
    fi

    # Create dynamic directory names
    ARMS_DIR="\${mapping_base}/\${REGION_TYPE}_segments_alnto${params.parent1_name}_alnto${params.parent2_name}"
    CEN_DIR="\${mapping_base}/CEN_segments_alnto${params.parent1_name}_alnto${params.parent2_name}"

    mkdir -p "\${ARMS_DIR}"
    if [[ "\${REGION_TYPE}" == "ARMS" ]]; then
        mkdir -p "\${CEN_DIR}"
    fi
    mkdir -p "\${plotting_base}"

    echo "Step 9.1: Performing mapping analysis..."

    # Step 9.1: Mapping using the bash script
    # Check if segments files exist and are not empty before mapping
    if [[ -s ${arms_segments_fasta} ]]; then
        echo "Mapping \${REGION_TYPE} segments..."
        bash ${projectDir}/scripts/09-mapping_240225.sh \\
            ${arms_segments_fasta} \\
            1 \\
            1 \\
            "\${ARMS_DIR}/" \\
            "${params.parent1_name}" \\
            "${params.parent2_name}" \\
            ${reference_genomes_dir}
    else
        echo "\${REGION_TYPE} segments file is empty, creating placeholder directories..."
        mkdir -p "\${ARMS_DIR}/summary_mm_sr"
        touch "\${ARMS_DIR}/summary_mm_sr/best_alignments.paf"
    fi

    if [[ "\${REGION_TYPE}" == "ARMS" ]] && [[ -s ${cen_segments_fasta} ]]; then
        echo "Mapping CEN segments..."
        bash ${projectDir}/scripts/09-mapping_240225.sh \\
            ${cen_segments_fasta} \\
            1 \\
            1 \\
            "\${CEN_DIR}/" \\
            "${params.parent1_name}" \\
            "${params.parent2_name}" \\
            ${reference_genomes_dir}
    elif [[ "\${REGION_TYPE}" == "ARMS" ]]; then
        echo "CEN segments file is empty, creating placeholder directories..."
        mkdir -p "\${CEN_DIR}/summary_mm_sr"
        touch "\${CEN_DIR}/summary_mm_sr/best_alignments.paf"
    fi
   


 
    echo "Step 9.2: Detecting inversions, chimeras and crossovers..."
    
    # Step 9.2: Detect inversions, chimeras for ARMS/CHR
    if [[ -d "\${ARMS_DIR}" && -s ${arms_segments_fasta} ]]; then
        cd "\${ARMS_DIR}"
        bash ${projectDir}/scripts/09.2-has_inversion_chimera_oneNM_tie_differentchromosomes_wrongecotype.sh "${projectDir}/scripts/09.2-has_inversion_chimera_oneNM_tie_differentchromosomes_wrongecotype.py"  ./
        cd - > /dev/null
    fi

    # Step 9.2: Detect inversions, chimeras for CEN (only in centromere-aware mode)
    if [[ "\${REGION_TYPE}" == "ARMS" ]] && [[ -d "\${CEN_DIR}" && -s ${cen_segments_fasta} ]]; then
        cd "\${CEN_DIR}"
        bash ${projectDir}/scripts/09.2-has_inversion_chimera_oneNM_tie_differentchromosomes_wrongecotype.sh "${projectDir}/scripts/09.2-has_inversion_chimera_oneNM_tie_differentchromosomes_wrongecotype.py" ./
        cd - > /dev/null
    fi
    
    echo "Step 9.3: Making summary tables and intersection plots..."
    
    # Step 9.3: Make summary tables for \${REGION_TYPE}
    if [[ -d "\${ARMS_DIR}" ]]; then
        cd "\${ARMS_DIR}"
        # Check if summary directories exist (not files)
        if ls -d summary_*/ 1> /dev/null 2>&1; then
            python ${projectDir}/scripts/09.3_makesummary.py summary_*/
            if [[ -f summary_table.csv ]]; then
                python ${projectDir}/scripts/09.3_plotintersection.py summary_table.csv
            fi
        else
            # Create empty summary table if no summary folders exist
            echo "read_id,status" > summary_table.csv
        fi
        cd - > /dev/null
    fi

    # Step 9.3: Make summary tables for CEN (only in centromere-aware mode)
    if [[ "\${REGION_TYPE}" == "ARMS" ]] && [[ -d "\${CEN_DIR}" ]]; then
        cd "\${CEN_DIR}"
        # Check if summary directories exist (not files)
        if ls -d summary_*/ 1> /dev/null 2>&1; then
            python ${projectDir}/scripts/09.3_makesummary.py summary_*/
            if [[ -f summary_table.csv ]]; then
                python ${projectDir}/scripts/09.3_plotintersection.py summary_table.csv
            fi
        else
            # Create empty summary table if no summary folders exist
            echo "read_id,status" > summary_table.csv
        fi
        cd - > /dev/null
    fi
    
    echo "Step 10: Getting best read pairs and plotting..."
    
    # Step 10: Get best read pairs for \${REGION_TYPE}
    if [[ -f "\${ARMS_DIR}/summary_table.csv" ]]; then
        python ${projectDir}/scripts/10-get_best_read_pairs.py \\
            "\${ARMS_DIR}/summary_table.csv" \\
            "\${ARMS_DIR}/best_best_other_reads.txt"
    else
        touch "\${ARMS_DIR}/best_best_other_reads.txt"
    fi

    # Step 10: Get best read pairs for CEN (only in centromere-aware mode)
    if [[ "\${REGION_TYPE}" == "ARMS" ]]; then
        CSV_FILE="\${CEN_DIR}/summary_table.csv"
        if [[ -f "\$CSV_FILE" ]] && \
           awk -F, 'NR==1 {exit !(/summary_lrhq/ && /summary_lrhqae/ && /summary_mm_ont/ && /summary_mm_sr/ && /summary_wm/)}' "\$CSV_FILE"; then

            echo "Valid summary table with all required columns found for CEN. Getting best read pairs..."
            python ${projectDir}/scripts/10-get_best_read_pairs.py \\
                "\$CSV_FILE" \\
                "\${CEN_DIR}/best_best_other_reads.txt"
        else
            echo "CEN summary table is empty or missing required columns. Creating placeholder and continuing..."
            touch "\${CEN_DIR}/best_best_other_reads.txt"
        fi
    fi 
    echo "Extracting PAF files for plotting..."

    # Extract PAF files for \${REGION_TYPE} - Against Parent1
    if [[ -s "\${ARMS_DIR}/best_best_other_reads.txt" && \\
          -f "\${ARMS_DIR}/summary_mm_sr/best_alignments.paf" ]]; then
        grep -F -f "\${ARMS_DIR}/best_best_other_reads.txt" \\
             "\${ARMS_DIR}/summary_mm_sr/best_alignments.paf" | \\
             grep -E "_${params.parent1_name[0]}[HCA][0-9]+" | cut -f1,6-9  > "\${plotting_base}/summary_mm_sr_\${REGION_TYPE}_against${params.parent1_name}.paf" || \\
             touch "\${plotting_base}/summary_mm_sr_\${REGION_TYPE}_against${params.parent1_name}.paf"
    else
        touch "\${plotting_base}/summary_mm_sr_\${REGION_TYPE}_against${params.parent1_name}.paf"
    fi

    # Extract PAF files for CEN - Against Parent1 (only in centromere-aware mode)
    if [[ "\${REGION_TYPE}" == "ARMS" ]] && [[ -s "\${CEN_DIR}/best_best_other_reads.txt" && \\
          -f "\${CEN_DIR}/summary_mm_sr/best_alignments.paf" ]]; then
        grep -F -f "\${CEN_DIR}/best_best_other_reads.txt" \\
             "\${CEN_DIR}/summary_mm_sr/best_alignments.paf" | \\
             grep -E "_${params.parent1_name[0]}[HCA][0-9]+" | cut -f1,6-9 > "\${plotting_base}/summary_mm_sr_CEN_against${params.parent1_name}.paf" || \\
             touch "\${plotting_base}/summary_mm_sr_CEN_against${params.parent1_name}.paf"
    elif [[ "\${REGION_TYPE}" == "ARMS" ]]; then
        touch "\${plotting_base}/summary_mm_sr_CEN_against${params.parent1_name}.paf"
    fi

    # Extract PAF files for \${REGION_TYPE} - Against Parent2
    if [[ -s "\${ARMS_DIR}/best_best_other_reads.txt" && \\
          -f "\${ARMS_DIR}/summary_mm_sr/best_alignments.paf" ]]; then
        grep -F -f "\${ARMS_DIR}/best_best_other_reads.txt" \\
             "\${ARMS_DIR}/summary_mm_sr/best_alignments.paf" | \\
             grep -E "_${params.parent2_name[0]}[HCA][0-9]+" | cut -f1,6-9 > "\${plotting_base}/summary_mm_sr_\${REGION_TYPE}_against${params.parent2_name}.paf" || \\
             touch "\${plotting_base}/summary_mm_sr_\${REGION_TYPE}_against${params.parent2_name}.paf"
    else
        touch "\${plotting_base}/summary_mm_sr_\${REGION_TYPE}_against${params.parent2_name}.paf"
    fi

    # Extract PAF files for CEN - Against Parent2 (only in centromere-aware mode)
    if [[ "\${REGION_TYPE}" == "ARMS" ]] && [[ -s "\${CEN_DIR}/best_best_other_reads.txt" && \\
          -f "\${CEN_DIR}/summary_mm_sr/best_alignments.paf" ]]; then
        grep -F -f "\${CEN_DIR}/best_best_other_reads.txt" \\
             "\${CEN_DIR}/summary_mm_sr/best_alignments.paf" | \\
             grep -E "_${params.parent2_name[0]}[HCA][0-9]+" | cut -f1,6-9 > "\${plotting_base}/summary_mm_sr_CEN_against${params.parent2_name}.paf" || \\
             touch "\${plotting_base}/summary_mm_sr_CEN_against${params.parent2_name}.paf"
    elif [[ "\${REGION_TYPE}" == "ARMS" ]]; then
        touch "\${plotting_base}/summary_mm_sr_CEN_against${params.parent2_name}.paf"
    fi
    
    echo "Calculating crossover widths..."
    
    # Calculate CO widths for \${REGION_TYPE}
    if [[ -s "\${ARMS_DIR}/best_best_other_reads.txt" && -s ${arms_bed_file} ]]; then
        grep -F -f "\${ARMS_DIR}/best_best_other_reads.txt" ${arms_bed_file} > \\
            "\${plotting_base}/CANDIDATE_reads_nonectopic_\${REGION_TYPE}.bed_of_best.tsv" || \\
            touch "\${plotting_base}/CANDIDATE_reads_nonectopic_\${REGION_TYPE}.bed_of_best.tsv"

        if [[ -s "\${plotting_base}/CANDIDATE_reads_nonectopic_\${REGION_TYPE}.bed_of_best.tsv" ]]; then
            python ${projectDir}/scripts/11-get_widths_co.py \\
                "\${plotting_base}/CANDIDATE_reads_nonectopic_\${REGION_TYPE}.bed_of_best.tsv" \\
                "\${plotting_base}/CANDIDATE_reads_nonectopic_\${REGION_TYPE}.bed_of_best.cowidths"
        else
            touch "\${plotting_base}/CANDIDATE_reads_nonectopic_\${REGION_TYPE}.bed_of_best.cowidths"
        fi
    else
        touch "\${plotting_base}/CANDIDATE_reads_nonectopic_\${REGION_TYPE}.bed_of_best.tsv"
        touch "\${plotting_base}/CANDIDATE_reads_nonectopic_\${REGION_TYPE}.bed_of_best.cowidths"
    fi

    # Calculate CO widths for CEN (only in centromere-aware mode)
    if [[ "\${REGION_TYPE}" == "ARMS" ]] && [[ -s "\${CEN_DIR}/best_best_other_reads.txt" && -s ${cen_bed_file} ]]; then
        grep -F -f "\${CEN_DIR}/best_best_other_reads.txt" ${cen_bed_file} > \\
            "\${plotting_base}/CANDIDATE_reads_nonectopic_CEN.bed_of_best.tsv" || \\
            touch "\${plotting_base}/CANDIDATE_reads_nonectopic_CEN.bed_of_best.tsv"

        if [[ -s "\${plotting_base}/CANDIDATE_reads_nonectopic_CEN.bed_of_best.tsv" ]]; then
            python ${projectDir}/scripts/11-get_widths_co.py \\
                "\${plotting_base}/CANDIDATE_reads_nonectopic_CEN.bed_of_best.tsv" \\
                "\${plotting_base}/CANDIDATE_reads_nonectopic_CEN.bed_of_best.cowidths"
        else
            touch "\${plotting_base}/CANDIDATE_reads_nonectopic_CEN.bed_of_best.cowidths"
        fi
    elif [[ "\${REGION_TYPE}" == "ARMS" ]]; then
        touch "\${plotting_base}/CANDIDATE_reads_nonectopic_CEN.bed_of_best.tsv"
        touch "\${plotting_base}/CANDIDATE_reads_nonectopic_CEN.bed_of_best.cowidths"
    fi

    # In CHR mode, create empty CEN PAF files as placeholders for downstream processes
    if [[ "\${REGION_TYPE}" == "CHR" ]]; then
        touch "\${plotting_base}/summary_mm_sr_CEN_against${params.parent1_name}.paf"
        touch "\${plotting_base}/summary_mm_sr_CEN_against${params.parent2_name}.paf"
    fi
    
    echo "Mapping analysis completed!"
    echo "Results saved in:"
    echo "- Mapping results: \${mapping_base}/"
    echo "- Plotting files: \${plotting_base}/"
    """
}
