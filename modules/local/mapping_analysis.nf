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
    path "threshold${threshold}/mapping_full_read/${params.sample_id}/nonectopic/ARMS_segments_alntoCol_alntoLer/**", emit: arms_mapping_results
    path "threshold${threshold}/mapping_full_read/${params.sample_id}/nonectopic/CEN_segments_alntoCol_alntoLer/**", emit: cen_mapping_results
    path "threshold${threshold}/plotting_crossover/${params.sample_id}/**", emit: plotting_results
    path "threshold${threshold}/mapping_full_read/${params.sample_id}/nonectopic/ARMS_segments_alntoCol_alntoLer/summary_table.csv", emit: arms_summary
    path "threshold${threshold}/mapping_full_read/${params.sample_id}/nonectopic/CEN_segments_alntoCol_alntoLer/summary_table.csv", emit: cen_summary
    path "threshold${threshold}/plotting_crossover/${params.sample_id}/*.paf", emit: paf_files
    path "threshold${threshold}/plotting_crossover/${params.sample_id}/*.cowidths", emit: cowidth_files
    
    script:
    """
    # Create directory structure
    threshold_dir="threshold${threshold}"
    mapping_base="\${threshold_dir}/mapping_full_read/${params.sample_id}/nonectopic"
    plotting_base="\${threshold_dir}/plotting_crossover/${params.sample_id}"
    
    mkdir -p "\${mapping_base}/CEN_segments_alntoCol_alntoLer"
    mkdir -p "\${mapping_base}/ARMS_segments_alntoCol_alntoLer"
    mkdir -p "\${plotting_base}"
    
    echo "Step 9.1: Performing mapping analysis..."
    
    # Step 9.1: Mapping using the bash script
    # Check if segments files exist and are not empty before mapping
    if [[ -s ${arms_segments_fasta} ]]; then
        echo "Mapping ARMS segments..."
        bash ${projectDir}/scripts/09-mapping_240225.sh \\
            ${arms_segments_fasta} \\
            1 \\
            1 \\
            \${mapping_base}/ARMS_segments_alntoCol_alntoLer/
    else
        echo "ARMS segments file is empty, creating placeholder directories..."
        mkdir -p "\${mapping_base}/ARMS_segments_alntoCol_alntoLer/summary_mm_sr"
        touch "\${mapping_base}/ARMS_segments_alntoCol_alntoLer/summary_mm_sr/best_alignments.paf"
    fi
    
    if [[ -s ${cen_segments_fasta} ]]; then
        echo "Mapping CEN segments..."
        bash ${projectDir}/scripts/09-mapping_240225.sh \\
            ${cen_segments_fasta} \\
            1 \\
            1 \\
            \${mapping_base}/CEN_segments_alntoCol_alntoLer/
    else
        echo "CEN segments file is empty, creating placeholder directories..."
        mkdir -p "\${mapping_base}/CEN_segments_alntoCol_alntoLer/summary_mm_sr"
        touch "\${mapping_base}/CEN_segments_alntoCol_alntoLer/summary_mm_sr/best_alignments.paf"
    fi
   


 
    echo "Step 9.2: Detecting inversions, chimeras and crossovers..."
    
    # Step 9.2: Detect inversions, chimeras for ARMS
    if [[ -d "\${mapping_base}/ARMS_segments_alntoCol_alntoLer" && -s ${arms_segments_fasta} ]]; then
        cd "\${mapping_base}/ARMS_segments_alntoCol_alntoLer"
        bash ${projectDir}/scripts/09.2-has_inversion_chimera_oneNM_tie_differentchromosomes_wrongecotype.sh "${projectDir}/scripts/09.2-has_inversion_chimera_oneNM_tie_differentchromosomes_wrongecotype.py"  ./
        cd - > /dev/null
    fi
    
    # Step 9.2: Detect inversions, chimeras for CEN
    if [[ -d "\${mapping_base}/CEN_segments_alntoCol_alntoLer" && -s ${cen_segments_fasta} ]]; then
        cd "\${mapping_base}/CEN_segments_alntoCol_alntoLer"
        bash ${projectDir}/scripts/09.2-has_inversion_chimera_oneNM_tie_differentchromosomes_wrongecotype.sh "${projectDir}/scripts/09.2-has_inversion_chimera_oneNM_tie_differentchromosomes_wrongecotype.py" ./
        cd - > /dev/null
    fi
    
    echo "Step 9.3: Making summary tables and intersection plots..."
    
    # Step 9.3: Make summary tables for ARMS
    if [[ -d "\${mapping_base}/ARMS_segments_alntoCol_alntoLer" ]]; then
        cd "\${mapping_base}/ARMS_segments_alntoCol_alntoLer"
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

    # Step 9.3: Make summary tables for CEN
    if [[ -d "\${mapping_base}/CEN_segments_alntoCol_alntoLer" ]]; then
        cd "\${mapping_base}/CEN_segments_alntoCol_alntoLer"
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
    
    # Step 10: Get best read pairs for ARMS
    if [[ -f "\${mapping_base}/ARMS_segments_alntoCol_alntoLer/summary_table.csv" ]]; then
        python ${projectDir}/scripts/10-get_best_read_pairs.py \\
            "\${mapping_base}/ARMS_segments_alntoCol_alntoLer/summary_table.csv" \\
            "\${mapping_base}/ARMS_segments_alntoCol_alntoLer/best_best_other_reads.txt"
    else
        touch "\${mapping_base}/ARMS_segments_alntoCol_alntoLer/best_best_other_reads.txt"
    fi
    
    # Step 10: Get best read pairs for CEN
  #  if [[ -f "\${mapping_base}/CEN_segments_alntoCol_alntoLer/summary_table.csv" ]]; then
  #      python ${projectDir}/scripts/10-get_best_read_pairs.py \\
  #          "\${mapping_base}/CEN_segments_alntoCol_alntoLer/summary_table.csv" \\
  #          "\${mapping_base}/CEN_segments_alntoCol_alntoLer/best_best_other_reads.txt"
  #  else
  #      touch "\${mapping_base}/CEN_segments_alntoCol_alntoLer/best_best_other_reads.txt"
  #  fi
   
## Step 10: Get best read pairs for CEN
CSV_FILE="\${mapping_base}/CEN_segments_alntoCol_alntoLer/summary_table.csv"
if [[ -f "\$CSV_FILE" ]] && \
   awk -F, 'NR==1 {exit !(/summary_lrhq/ && /summary_lrhqae/ && /summary_mm_ont/ && /summary_mm_sr/ && /summary_wm/)}' "\$CSV_FILE"; then
    
    echo "Valid summary table with all required columns found for CEN. Getting best read pairs..."
    python ${projectDir}/scripts/10-get_best_read_pairs.py \\
        "\$CSV_FILE" \\
        "\${mapping_base}/CEN_segments_alntoCol_alntoLer/best_best_other_reads.txt"
else
    echo "CEN summary table is empty or missing required columns. Creating placeholder and continuing..."
    touch "\${mapping_base}/CEN_segments_alntoCol_alntoLer/best_best_other_reads.txt"
fi 
    echo "Extracting PAF files for plotting..."
    
    # Extract PAF files for ARMS - Against Col
    if [[ -s "\${mapping_base}/ARMS_segments_alntoCol_alntoLer/best_best_other_reads.txt" && \\
          -f "\${mapping_base}/ARMS_segments_alntoCol_alntoLer/summary_mm_sr/best_alignments.paf" ]]; then
        grep -F -f "\${mapping_base}/ARMS_segments_alntoCol_alntoLer/best_best_other_reads.txt" \\
             "\${mapping_base}/ARMS_segments_alntoCol_alntoLer/summary_mm_sr/best_alignments.paf" | \\
             grep Col | cut -f1,6-9  > "\${plotting_base}/summary_mm_sr_ARMS_againstCol.paf" || \\
             touch "\${plotting_base}/summary_mm_sr_ARMS_againstCol.paf"
    else
        touch "\${plotting_base}/summary_mm_sr_ARMS_againstCol.paf"
    fi
    
    # Extract PAF files for CEN - Against Col
    if [[ -s "\${mapping_base}/CEN_segments_alntoCol_alntoLer/best_best_other_reads.txt" && \\
          -f "\${mapping_base}/CEN_segments_alntoCol_alntoLer/summary_mm_sr/best_alignments.paf" ]]; then
        grep -F -f "\${mapping_base}/CEN_segments_alntoCol_alntoLer/best_best_other_reads.txt" \\
             "\${mapping_base}/CEN_segments_alntoCol_alntoLer/summary_mm_sr/best_alignments.paf" | \\
             grep Col | cut -f1,6-9 > "\${plotting_base}/summary_mm_sr_CEN_againstCol.paf" || \\
             touch "\${plotting_base}/summary_mm_sr_CEN_againstCol.paf"
    else
        touch "\${plotting_base}/summary_mm_sr_CEN_againstCol.paf"
    fi
    
    # Extract PAF files for ARMS - Against Ler
    if [[ -s "\${mapping_base}/ARMS_segments_alntoCol_alntoLer/best_best_other_reads.txt" && \\
          -f "\${mapping_base}/ARMS_segments_alntoCol_alntoLer/summary_mm_sr/best_alignments.paf" ]]; then
        grep -F -f "\${mapping_base}/ARMS_segments_alntoCol_alntoLer/best_best_other_reads.txt" \\
             "\${mapping_base}/ARMS_segments_alntoCol_alntoLer/summary_mm_sr/best_alignments.paf" | \\
             grep Ler | cut -f1,6-9 > "\${plotting_base}/summary_mm_sr_ARMS_againstLer.paf" || \\
             touch "\${plotting_base}/summary_mm_sr_ARMS_againstLer.paf"
    else
        touch "\${plotting_base}/summary_mm_sr_ARMS_againstLer.paf"
    fi
    
    # Extract PAF files for CEN - Against Ler
    if [[ -s "\${mapping_base}/CEN_segments_alntoCol_alntoLer/best_best_other_reads.txt" && \\
          -f "\${mapping_base}/CEN_segments_alntoCol_alntoLer/summary_mm_sr/best_alignments.paf" ]]; then
        grep -F -f "\${mapping_base}/CEN_segments_alntoCol_alntoLer/best_best_other_reads.txt" \\
             "\${mapping_base}/CEN_segments_alntoCol_alntoLer/summary_mm_sr/best_alignments.paf" | \\
             grep Ler | cut -f1,6-9 > "\${plotting_base}/summary_mm_sr_CEN_againstLer.paf" || \\
             touch "\${plotting_base}/summary_mm_sr_CEN_againstLer.paf"
    else
        touch "\${plotting_base}/summary_mm_sr_CEN_againstLer.paf"
    fi
    
    echo "Calculating crossover widths..."
    
    # Calculate CO widths for ARMS
    if [[ -s "\${mapping_base}/ARMS_segments_alntoCol_alntoLer/best_best_other_reads.txt" && -s ${arms_bed_file} ]]; then
        grep -F -f "\${mapping_base}/ARMS_segments_alntoCol_alntoLer/best_best_other_reads.txt" ${arms_bed_file} > \\
            "\${plotting_base}/CANDIDATE_reads_nonectopic_ARMS.bed_of_best.tsv" || \\
            touch "\${plotting_base}/CANDIDATE_reads_nonectopic_ARMS.bed_of_best.tsv"
        
        if [[ -s "\${plotting_base}/CANDIDATE_reads_nonectopic_ARMS.bed_of_best.tsv" ]]; then
            python ${projectDir}/scripts/11-get_widths_co.py \\
                "\${plotting_base}/CANDIDATE_reads_nonectopic_ARMS.bed_of_best.tsv" \\
                "\${plotting_base}/CANDIDATE_reads_nonectopic_ARMS.bed_of_best.cowidths"
        else
            touch "\${plotting_base}/CANDIDATE_reads_nonectopic_ARMS.bed_of_best.cowidths"
        fi
    else
        touch "\${plotting_base}/CANDIDATE_reads_nonectopic_ARMS.bed_of_best.tsv"
        touch "\${plotting_base}/CANDIDATE_reads_nonectopic_ARMS.bed_of_best.cowidths"
    fi
    
    # Calculate CO widths for CEN
    if [[ -s "\${mapping_base}/CEN_segments_alntoCol_alntoLer/best_best_other_reads.txt" && -s ${cen_bed_file} ]]; then
        grep -F -f "\${mapping_base}/CEN_segments_alntoCol_alntoLer/best_best_other_reads.txt" ${cen_bed_file} > \\
            "\${plotting_base}/CANDIDATE_reads_nonectopic_CEN.bed_of_best.tsv" || \\
            touch "\${plotting_base}/CANDIDATE_reads_nonectopic_CEN.bed_of_best.tsv"
        
        if [[ -s "\${plotting_base}/CANDIDATE_reads_nonectopic_CEN.bed_of_best.tsv" ]]; then
            python ${projectDir}/scripts/11-get_widths_co.py \\
                "\${plotting_base}/CANDIDATE_reads_nonectopic_CEN.bed_of_best.tsv" \\
                "\${plotting_base}/CANDIDATE_reads_nonectopic_CEN.bed_of_best.cowidths"
        else
            touch "\${plotting_base}/CANDIDATE_reads_nonectopic_CEN.bed_of_best.cowidths"
        fi
    else
        touch "\${plotting_base}/CANDIDATE_reads_nonectopic_CEN.bed_of_best.tsv"
        touch "\${plotting_base}/CANDIDATE_reads_nonectopic_CEN.bed_of_best.cowidths"
    fi
    
    echo "Mapping analysis completed!"
    echo "Results saved in:"
    echo "- Mapping results: \${mapping_base}/"
    echo "- Plotting files: \${plotting_base}/"
    """
}
