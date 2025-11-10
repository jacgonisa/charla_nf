#!/bin/bash
#
# Usage:
#   ./align_to_genomes.sh reads.fa MIN_MAPQ MAX_ALNS OUTDIR PARENT1 PARENT2 GENOME_DIR
#
# Aligns reads against parent genomes using Winnowmap and Minimap2.
#
# Outputs:
#   - PAF files with CIGAR (-c option) in OUTDIR
#   - Filters out alignments with MAPQ < MIN_MAPQ
#   - Skips secondary alignments (--secondary=no)
#

# Check for required arguments
if [ "$#" -ne 7 ]; then
    echo "Usage: $0 reads.fa MIN_MAPQ MAX_ALNS OUTDIR PARENT1 PARENT2 GENOME_DIR"
    exit 1
fi

# Input parameters
READS="$1"         # e.g. reads.fa
MIN_MAPQ="$2"      # e.g. 20
MAX_ALNS="$3"      # e.g. 5
OUTDIR="$4"
PARENT1="$5"       # e.g. C57BL6J or Col-0
PARENT2="$6"       # e.g. DBA2J or Ler-0
GENOME_DIR="$7"    # e.g. /path/to/genomes

# Genome names (dynamic - use provided parent names, no suffix)
GENOMES=("${PARENT1}" "${PARENT2}")

# Function to find genome file with any extension
find_genome_file() {
    local genome_name="$1"
    local genome_dir="$2"

    # Try common extensions
    for ext in .fa .fasta .fna; do
        if [[ -f "${genome_dir}/${genome_name}${ext}" ]]; then
            echo "${genome_dir}/${genome_name}${ext}"
            return 0
        fi
    done

    # If not found, return error
    echo "ERROR: Could not find genome file for ${genome_name} in ${genome_dir}" >&2
    return 1
}

# Function to generate repetitive k-mer file using meryl
generate_repetitive_kmers() {
    local genome_file="$1"
    local genome_name="$2"
    local genome_dir="$3"

    local meryl_db="${genome_dir}/${genome_name}_merylDB"
    local repetitive_txt="${genome_dir}/${genome_name}_repetitive_k15.txt"

    # Check if repetitive k-mer file already exists
    if [[ -f "${repetitive_txt}" ]]; then
        echo "Repetitive k-mer file already exists: ${repetitive_txt}"
        echo "${repetitive_txt}"
        return 0
    fi

    echo "Generating repetitive k-mers for ${genome_name}..."

    # Generate meryl database
    echo "Running: meryl count k=15 output ${meryl_db} ${genome_file}"
    meryl count k=15 output "${meryl_db}" "${genome_file}"

    if [[ $? -ne 0 ]]; then
        echo "ERROR: meryl count failed for ${genome_name}" >&2
        return 1
    fi

    # Extract repetitive k-mers (those appearing in >99.98% of positions)
    echo "Running: meryl print greater-than distinct=0.9998 ${meryl_db} > ${repetitive_txt}"
    meryl print greater-than distinct=0.9998 "${meryl_db}" > "${repetitive_txt}"

    if [[ $? -ne 0 ]]; then
        echo "ERROR: meryl print failed for ${genome_name}" >&2
        return 1
    fi

    echo "Repetitive k-mer file generated: ${repetitive_txt}"
    echo "${repetitive_txt}"
    return 0
}

# Array of modes to loop over (including winnowmap)
modes=("wm_mapont" "mm_ont" "mm_lr_hq" "mm_lr_hqae" "mm_sr")

# Create a directory for output
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

        # Find the genome file
        genome_file=$(find_genome_file "${genome}" "${GENOME_DIR}")
        if [[ $? -ne 0 ]]; then
            echo "ERROR: Genome file not found for ${genome}, skipping..."
            continue
        fi

        echo "Using genome file: ${genome_file}"

        base_name=$(basename "${READS}" .fa)
        paf_out="${OUTDIR}/${base_name}_${mode}_alnTo_${genome}.paf"

        # Build alignment command dynamically based on mode
        case "${mode}" in
            wm_mapont)
                # Generate repetitive k-mer file for winnowmap
                repetitive_kmers=$(generate_repetitive_kmers "${genome_file}" "${genome}" "${GENOME_DIR}")
                if [[ $? -ne 0 ]]; then
                    echo "ERROR: Failed to generate repetitive k-mers for ${genome}, skipping winnowmap..."
                    continue
                fi
                align_cmd="winnowmap -c --secondary=no -W \"${repetitive_kmers}\" -x map-ont -t 10 -p 1.0 -N ${MAX_ALNS} \"${genome_file}\" \"${READS}\""
                ;;
            mm_ont)
                align_cmd="minimap2 -c --secondary=no -x map-ont -t 10 -p 1.0 -N ${MAX_ALNS} \"${genome_file}\" \"${READS}\""
                ;;
            mm_lr_hq)
                align_cmd="minimap2 -c --secondary=no -x lr:hq -t 10 -p 1.0 -N ${MAX_ALNS} \"${genome_file}\" \"${READS}\""
                ;;
            mm_lr_hqae)
                align_cmd="minimap2 -c --secondary=no -x lr:hqae -t 10 -p 1.0 -N ${MAX_ALNS} \"${genome_file}\" \"${READS}\""
                ;;
            mm_sr)
                align_cmd="minimap2 -c --secondary=no -x sr -t 10 -p 1.0 -N ${MAX_ALNS} \"${genome_file}\" \"${READS}\""
                ;;
            *)
                echo "ERROR: Unknown mode ${mode}"
                continue
                ;;
        esac

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

