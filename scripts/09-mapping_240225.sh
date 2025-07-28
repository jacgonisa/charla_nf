#!/bin/bash
#
# Usage:
#   ./align_to_genomes.sh reads.fa MIN_MAPQ MAX_ALNS
#
# Aligns reads against Col-0 and Ler-0 genomes using Winnowmap and Minimap2.
#
# Outputs:
#   - PAF files with CIGAR (-c option) in "alignment_results/"
#   - Filters out alignments with MAPQ < MIN_MAPQ
#   - Skips secondary alignments (--secondary=no)
#

# Check for required arguments
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 reads.fa MIN_MAPQ MAX_ALNS OUTDIR"
    exit 1
fi

# Input parameters
READS="$1"         # e.g. reads.fa
MIN_MAPQ="$2"      # e.g. 20
MAX_ALNS="$3"      # e.g. 5
OUTDIR="$4"

# Genome names (fixed)
#GENOMES=("ColLer_renamed")
GENOMES=("Col-0_renamed" "Ler-0_renamed")

# Define alignment commands with -c (PAF with CIGAR) and --secondary=no
declare -A MODE_CMD
MODE_CMD["wm_mapont"]="winnowmap -c --secondary=no -W index/\${genome}_repetitive_k15.txt -x map-ont -t 10 -p 1.0 -N ${MAX_ALNS} index/\${genome}.fa \${READS}"
MODE_CMD["mm_lr_hq"]="minimap2 -c --secondary=no -x lr:hq -t 10 -p 1.0 -N ${MAX_ALNS} index/\${genome}.fa \${READS}"
MODE_CMD["mm_ont"]="minimap2 -c --secondary=no -x map-ont -t 10 -p 1.0 -N ${MAX_ALNS} index/\${genome}.fa \${READS}"
MODE_CMD["mm_lr_hqae"]="minimap2 -c --secondary=no -x lr:hqae -t 10 -p 1.0 -N ${MAX_ALNS} index/\${genome}.fa \${READS}"
MODE_CMD["mm_sr"]="minimap2 -c --secondary=no -x sr -t 10 -p 1.0 -N ${MAX_ALNS} index/\${genome}.fa \${READS}"

# Array of modes to loop over
modes=("wm_mapont" "mm_ont" "mm_lr_hq" "mm_lr_hqae" "mm_sr")

# Create a directory for output
#OUTDIR="alignment_results"
mkdir -p "${OUTDIR}"

# Function to process a PAF file and compute simple metrics: number of alignments and average MAPQ.
process_paf() {
    local paf_file="$1"
    local num_alns=$(wc -l < "${paf_file}")
    local avg_mapq=$(awk -v minq="${MIN_MAPQ}" '{
         if ($12 >= minq) { sum+=$12; count++ }
     } END { if(count>0) printf "%.2f", sum/count; else print "0" }' "${paf_file}")
    echo "${num_alns},${avg_mapq}"
}


# Loop over each alignment mode and genome
for mode in "${modes[@]}"; do
    echo "----------------------------------------------"
    echo "Running mode: ${mode}"
    for genome in "${GENOMES[@]}"; do
        echo "Aligning to genome: ${genome}"
        base_name=$(basename "${READS}" .fa)
        paf_out="${OUTDIR}/${base_name}_${mode}_alnTo_${genome}.paf"
        
        # Build and execute the alignment command
        align_cmd=${MODE_CMD[${mode}]}
        echo "Running: ${align_cmd}"
        eval "${align_cmd}" | awk -v minq="${MIN_MAPQ}" '($12 >= minq)' > "${paf_out}"
        
        # Remove empty output files
        if [ ! -s "${paf_out}" ]; then
            echo "No alignments found for ${mode} on ${genome}; removing empty file."
            rm -f "${paf_out}"
        else
            # Compute alignment stats
            metrics=$(process_paf "${paf_out}")
            num_alns=$(echo "${metrics}" | cut -d, -f1)
            avg_mapq=$(echo "${metrics}" | cut -d, -f2)
            echo "Mode ${mode}, Genome ${genome}: Alignments = ${num_alns}, Average MAPQ = ${avg_mapq}"
        fi
        
        echo "Finished processing ${mode} on genome ${genome}."
        echo
    done
done

echo "All alignments complete. Check '${OUTDIR}' for PAF files."

