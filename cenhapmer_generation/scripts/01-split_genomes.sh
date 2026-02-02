#!/bin/bash

###############################################################################
# Script: 01-split_genomes.sh
# Description: Split genome FASTAs into centromere (CEN) and chromosome arms
#              (ARMS) regions using BED files
#
# Usage: ./01-split_genomes.sh [options]
#
# Options:
#   --fasta_dir DIR       Directory containing genome FASTA files (default: genomes)
#   --fasta_suffix SUFFIX FASTA file suffix (default: .ragtag_scaffolds.fa)
#   --bed_dir DIR         Directory containing BED annotation files (default: annotations)
#   --out_dir DIR         Output directory for split genomes (default: genomes_split)
#   --parent1 NAME        Name of parent 1 accession (default: Col-0)
#   --parent2 NAME        Name of parent 2 accession (default: Ler-0)
#   --num_chr N           Number of chromosomes (default: 5)
#   --chr_prefix PREFIX   Chromosome prefix in BED files (default: Chr)
#   --help                Show this help message
#
# Requirements:
#   - bedtools (for getfasta command)
#   - Genome FASTA files named: ${ACCESSION}.ragtag_scaffolds.fa
#   - BED files for CEN regions: ${ACCESSION}.ragtag_scaffolds_CEN.bed
#   - BED files for ARMS regions: ${ACCESSION}.ragtag_scaffolds_ARMS.bed
#
# Output:
#   - genomes_split/${ACCESSION}/CEN/${ACCESSION}_CEN_${CHR}.fa
#   - genomes_split/${ACCESSION}/ARMS/${ACCESSION}_ARMS_${CHR}.fa
###############################################################################

# Default parameters
FASTA_DIR="genomes"
BED_DIR="annotations"
OUT_DIR="genomes_split"
PARENT1="Col-0"
PARENT2="Ler-0"
NUM_CHR=5
CHR_PREFIX="Chr"
SKIP_CEN=false
FASTA_SUFFIX=".ragtag_scaffolds.fa"
BED_CEN_SUFFIX=".ragtag_scaffolds_CEN.bed"
BED_ARMS_SUFFIX=".ragtag_scaffolds_ARMS.bed"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --fasta_dir)
            FASTA_DIR="$2"
            shift 2
            ;;
        --fasta_suffix)
            FASTA_SUFFIX="$2"
            shift 2
            ;;
        --bed_dir)
            BED_DIR="$2"
            shift 2
            ;;
        --out_dir)
            OUT_DIR="$2"
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
        --skip-cen)
            SKIP_CEN=true
            shift
            ;;
        --help)
            head -n 30 "$0" | grep "^#" | sed 's/^# //'
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Generate chromosome array
CHROMOSOMES=()
for i in $(seq 1 $NUM_CHR); do
    CHROMOSOMES+=("${CHR_PREFIX}${i}")
done

echo "=== Genome Splitting Configuration ==="
echo "FASTA directory:    $FASTA_DIR"
echo "BED directory:      $BED_DIR"
echo "Output directory:   $OUT_DIR"
echo "Parent 1:           $PARENT1"
echo "Parent 2:           $PARENT2"
echo "Chromosomes:        ${CHROMOSOMES[@]}"
echo "Skip centromeres:   $SKIP_CEN"
echo "======================================"
echo ""

# Create output directory
mkdir -p "$OUT_DIR"

# Process each accession
for ACCESSION in "$PARENT1" "$PARENT2"; do
    echo "[INFO] Processing $ACCESSION..."

    # Input files
    FASTA="${FASTA_DIR}/${ACCESSION}${FASTA_SUFFIX}"

    # Validate input files
    if [[ ! -f "$FASTA" ]]; then
        echo "[ERROR] Missing FASTA: $FASTA"
        echo "[ERROR] Skipping $ACCESSION"
        continue
    fi

    if [[ "$SKIP_CEN" == "true" ]]; then
        # Mode: Chromosome-level hapmers (no CEN/ARMS split)
        echo "[INFO] Skip-CEN mode: extracting whole chromosomes"

        # Output folder for whole chromosomes
        CHR_OUT="${OUT_DIR}/${ACCESSION}/CHR"
        mkdir -p "$CHR_OUT"

        # Extract each chromosome as a whole
        for CHR in "${CHROMOSOMES[@]}"; do
            echo "  [CHR] Extracting ${CHR}..."

            # Extract chromosome and convert to single-line FASTA format (required for KMC)
            if command -v samtools &> /dev/null; then
                # Use samtools to extract, then convert to single-line format
                samtools faidx "$FASTA" "$CHR" 2>/dev/null | \
                    awk '/^>/ {if (NR > 1) printf "\n"; print; next} {printf "%s", $0} END {printf "\n"}' \
                    > "${CHR_OUT}/${ACCESSION}_CHR_${CHR}.fa" || \
                    echo "  [WARNING] Could not extract ${CHR} from $ACCESSION"
            else
                # Fallback: use awk to extract from FASTA (works with single-line FASTA)
                # This extracts the chromosome and ensures single-line format
                awk -v chr="$CHR" '
                BEGIN { p=0; seq="" }
                /^>/ {
                    # If we were printing a previous chromosome, output it
                    if (p && seq != "") {
                        printf "%s\n", seq
                        seq = ""
                    }
                    p = 0
                    # Check if this is our target chromosome
                    # Match exact: >chr10 or >chr10<space> but not >chr1
                    if ($0 == ">"chr || $0 ~ "^>"chr"[ \t]" || $0 ~ "^>"chr"$") {
                        print $0
                        p = 1
                    }
                    next
                }
                p { seq = seq $0 }
                END {
                    if (p && seq != "") printf "%s\n", seq
                }
                ' "$FASTA" > "${CHR_OUT}/${ACCESSION}_CHR_${CHR}.fa"

                if [[ ! -s "${CHR_OUT}/${ACCESSION}_CHR_${CHR}.fa" ]]; then
                    echo "  [WARNING] No sequence found for ${CHR} in $ACCESSION"
                fi
            fi
        done
    else
        # Original mode: CEN/ARMS split using BED files
        CEN_BED="${BED_DIR}/${ACCESSION}${BED_CEN_SUFFIX}"
        ARMS_BED="${BED_DIR}/${ACCESSION}${BED_ARMS_SUFFIX}"

        if [[ ! -f "$CEN_BED" ]]; then
            echo "[ERROR] Missing CEN BED file: $CEN_BED"
            echo "[ERROR] Skipping $ACCESSION"
            continue
        fi
        if [[ ! -f "$ARMS_BED" ]]; then
            echo "[ERROR] Missing ARMS BED file: $ARMS_BED"
            echo "[ERROR] Skipping $ACCESSION"
            continue
        fi

        # Output folders
        CEN_OUT="${OUT_DIR}/${ACCESSION}/CEN"
        ARMS_OUT="${OUT_DIR}/${ACCESSION}/ARMS"
        mkdir -p "$CEN_OUT" "$ARMS_OUT"

        # Process each chromosome
        for CHR in "${CHROMOSOMES[@]}"; do
            echo "  [CEN] Extracting ${CHR}..."
            awk -v chr="$CHR" '$1 == chr' "$CEN_BED" > tmp_cen_${ACCESSION}_${CHR}.bed

            if [[ -s tmp_cen_${ACCESSION}_${CHR}.bed ]]; then
                bedtools getfasta -fi "$FASTA" -bed tmp_cen_${ACCESSION}_${CHR}.bed \
                    -fo "${CEN_OUT}/${ACCESSION}_CEN_${CHR}.fa" -name
            else
                echo "  [WARNING] No CEN regions found for ${CHR} in $ACCESSION"
                touch "${CEN_OUT}/${ACCESSION}_CEN_${CHR}.fa"
            fi

            echo "  [ARMS] Extracting ${CHR}..."
            awk -v chr="$CHR" '$1 == chr' "$ARMS_BED" > tmp_arms_${ACCESSION}_${CHR}.bed

            if [[ -s tmp_arms_${ACCESSION}_${CHR}.bed ]]; then
                bedtools getfasta -fi "$FASTA" -bed tmp_arms_${ACCESSION}_${CHR}.bed \
                    -fo "${ARMS_OUT}/${ACCESSION}_ARMS_${CHR}.fa" -name
            else
                echo "  [WARNING] No ARMS regions found for ${CHR} in $ACCESSION"
                touch "${ARMS_OUT}/${ACCESSION}_ARMS_${CHR}.fa"
            fi

            rm -f tmp_cen_${ACCESSION}_${CHR}.bed tmp_arms_${ACCESSION}_${CHR}.bed
        done
    fi

    echo "[DONE] Finished $ACCESSION"
    echo ""
done

echo "✅ Genome splitting complete!"
echo "Output directory: $OUT_DIR"
