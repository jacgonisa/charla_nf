#!/bin/bash

###############################################################################
# Script: 02-generate_kmers.sh
# Description: Generate k-mer databases for each accession, chromosome, and
#              region (CEN/ARMS) using KMC
#
# Usage: ./02-generate_kmers.sh [options]
#
# Options:
#   --input_dir DIR       Directory with split genomes (default: genomes_split)
#   --output_dir DIR      Output directory for k-mer databases (default: 01-all_kmers_per_Chr_CEN_and_accesion)
#   --parent1 NAME        Name of parent 1 accession (default: Col-0)
#   --parent2 NAME        Name of parent 2 accession (default: Ler-0)
#   --num_chr N           Number of chromosomes (default: 5)
#   --chr_prefix PREFIX   Chromosome prefix (default: Chr)
#   --kmer_sizes SIZES    Comma-separated k-mer sizes (default: 21,31,41)
#   --min_count N         Minimum k-mer count (default: 1)
#   --max_count N         Maximum k-mer count (default: 100000000)
#   --threads N           Number of threads for KMC (default: 4)
#   --help                Show this help message
#
# Requirements:
#   - KMC (k-mer Counter) tool
#   - Split genome FASTA files from step 01
#
# Output:
#   - ${output_dir}/${ACCESSION}/${CHR}/{CEN,ARMS}/${ACCESSION}_{CEN|ARMS}_${CHR}_k${K}.kmc_{pre,suf}
###############################################################################

# Default parameters
INPUT_DIR="genomes_split"
OUTPUT_DIR="01-all_kmers_per_Chr_CEN_and_accesion"
PARENT1="Col-0"
PARENT2="Ler-0"
NUM_CHR=5
CHR_PREFIX="Chr"
KMER_SIZES="21,31,41"
MIN_COUNT=1
MAX_COUNT=100000000
THREADS=4
SKIP_CEN=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --input_dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        --output_dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --parent1)
            PARENT1="$2"
            shift 2
            ;;
        --parent2)
            PARENT2="$2"
            shift 2
            ;;
        --num_chr)
            NUM_CHR="$2"
            shift 2
            ;;
        --chr_prefix)
            CHR_PREFIX="$2"
            shift 2
            ;;
        --kmer_sizes)
            KMER_SIZES="$2"
            shift 2
            ;;
        --min_count)
            MIN_COUNT="$2"
            shift 2
            ;;
        --max_count)
            MAX_COUNT="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --skip-cen)
            SKIP_CEN=true
            shift
            ;;
        --help)
            head -n 35 "$0" | grep "^#" | sed 's/^# //'
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Parse k-mer sizes into array
IFS=',' read -ra K_SIZES_ARRAY <<< "$KMER_SIZES"

# Generate chromosome array
CHROMOSOMES=()
for i in $(seq 1 $NUM_CHR); do
    CHROMOSOMES+=("${CHR_PREFIX}${i}")
done

echo "=== K-mer Generation Configuration ==="
echo "Input directory:    $INPUT_DIR"
echo "Output directory:   $OUTPUT_DIR"
echo "Parent 1:           $PARENT1"
echo "Parent 2:           $PARENT2"
echo "Chromosomes:        ${CHROMOSOMES[@]}"
echo "K-mer sizes:        ${K_SIZES_ARRAY[@]}"
echo "Min count:          $MIN_COUNT"
echo "Max count:          $MAX_COUNT"
echo "Threads:            $THREADS"
echo "Skip centromeres:   $SKIP_CEN"
echo "======================================="
echo ""

# Create output and temp directories
mkdir -p "${OUTPUT_DIR}/tmp"

# Loop over accessions
for ACCESSION in "$PARENT1" "$PARENT2"; do
    echo "🧬 Processing accession: $ACCESSION"

    ACCESSION_DIR="${OUTPUT_DIR}/${ACCESSION}"
    mkdir -p "$ACCESSION_DIR"

    for CHR in "${CHROMOSOMES[@]}"; do
        echo "  🧬 Processing chromosome: $CHR"

        CHR_DIR="${ACCESSION_DIR}/${CHR}"

        if [[ "$SKIP_CEN" == "true" ]]; then
            # Skip-CEN mode: Process whole chromosomes
            mkdir -p "${CHR_DIR}/CHR"

            CHR_FILE="${INPUT_DIR}/${ACCESSION}/CHR/${ACCESSION}_CHR_${CHR}.fa"
            echo "    📂 CHR file:  $CHR_FILE"

            # Validate input file
            if [[ ! -f "$CHR_FILE" ]]; then
                echo "    ⚠️  Warning: CHR file not found, skipping"
                continue
            fi

            # Check if file is empty
            if [[ ! -s "$CHR_FILE" ]]; then
                echo "    ⚠️  Warning: CHR file is empty, skipping"
                continue
            fi

            # Process each k-mer size
            for K in "${K_SIZES_ARRAY[@]}"; do
                echo "      🔢 Processing k-mer size: $K"

                CHR_OUTPUT="${CHR_DIR}/CHR/${ACCESSION}_CHR_${CHR}_k${K}"

                # Run KMC on whole chromosome
                echo "      ⚙️  Running KMC on chromosome..."
                kmc -fa -k${K} -ci${MIN_COUNT} -cs${MAX_COUNT} -t${THREADS} \
                    "$CHR_FILE" "$CHR_OUTPUT" "${OUTPUT_DIR}/tmp"

                if [[ $? -eq 0 ]]; then
                    echo "      ✅ CHR k-mer generation successful"
                else
                    echo "      ❌ CHR k-mer generation failed"
                fi
            done
        else
            # Original mode: Process CEN and ARMS separately
            mkdir -p "${CHR_DIR}/CEN"
            mkdir -p "${CHR_DIR}/ARMS"

            # Input FASTA paths
            CEN_FILE="${INPUT_DIR}/${ACCESSION}/CEN/${ACCESSION}_CEN_${CHR}.fa"
            ARMS_FILE="${INPUT_DIR}/${ACCESSION}/ARMS/${ACCESSION}_ARMS_${CHR}.fa"

            echo "    📂 CEN file:  $CEN_FILE"
            echo "    📂 ARMS file: $ARMS_FILE"

            # Validate input files
            if [[ ! -f "$CEN_FILE" ]]; then
                echo "    ⚠️  Warning: CEN file not found, skipping"
                continue
            fi
            if [[ ! -f "$ARMS_FILE" ]]; then
                echo "    ⚠️  Warning: ARMS file not found, skipping"
                continue
            fi

            # Check if files are empty
            if [[ ! -s "$CEN_FILE" ]]; then
                echo "    ⚠️  Warning: CEN file is empty, skipping"
            fi
            if [[ ! -s "$ARMS_FILE" ]]; then
                echo "    ⚠️  Warning: ARMS file is empty, skipping"
            fi

            # Process each k-mer size
            for K in "${K_SIZES_ARRAY[@]}"; do
                echo "      🔢 Processing k-mer size: $K"

                CEN_OUTPUT="${CHR_DIR}/CEN/${ACCESSION}_CEN_${CHR}_k${K}"
                ARMS_OUTPUT="${CHR_DIR}/ARMS/${ACCESSION}_ARMS_${CHR}_k${K}"

                # Run KMC on centromere (only if file has content)
                if [[ -s "$CEN_FILE" ]]; then
                    echo "      ⚙️  Running KMC on centromere..."
                    kmc -fa -k${K} -ci${MIN_COUNT} -cs${MAX_COUNT} -t${THREADS} \
                        "$CEN_FILE" "$CEN_OUTPUT" "${OUTPUT_DIR}/tmp"

                    if [[ $? -eq 0 ]]; then
                        echo "      ✅ CEN k-mer generation successful"
                    else
                        echo "      ❌ CEN k-mer generation failed"
                    fi
                else
                    echo "      ⏭️  Skipping empty CEN file"
                fi

                # Run KMC on chromosome arms
                if [[ -s "$ARMS_FILE" ]]; then
                    echo "      ⚙️  Running KMC on chromosome arms..."
                    kmc -fa -k${K} -ci${MIN_COUNT} -cs${MAX_COUNT} -t${THREADS} \
                        "$ARMS_FILE" "$ARMS_OUTPUT" "${OUTPUT_DIR}/tmp"

                    if [[ $? -eq 0 ]]; then
                        echo "      ✅ ARMS k-mer generation successful"
                    else
                        echo "      ❌ ARMS k-mer generation failed"
                    fi
                else
                    echo "      ⏭️  Skipping empty ARMS file"
                fi
            done
        fi
    done

    echo ""
done

# Clean up temporary directory
echo "🧹 Cleaning up temporary files..."
rm -rf "${OUTPUT_DIR}/tmp"

echo "✅ All k-mer processing completed!"
echo "Output directory: $OUTPUT_DIR"
