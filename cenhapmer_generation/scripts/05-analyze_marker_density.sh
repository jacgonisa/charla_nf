#!/bin/bash

###############################################################################
# Script: 05-analyze_marker_density.sh
# Description: Analyze marker (cenhapmer) density across chromosomes
#              Compares RAW k-mers (step 01) vs UNIQUE hapmers (step 03)
#              Provides both absolute counts and normalized to genome size
#
# Usage: ./05-analyze_marker_density.sh [options]
#
# Options:
#   --raw_kmer_dir DIR     Directory with raw k-mers from step 01 (default: 01-all_kmers_per_Chr_CEN_and_accesion)
#   --cenhapmer_dir DIR    Directory with unique hapmers from step 03 (default: 03-cenhapmers)
#   --genome_dir DIR       Directory with genome FASTAs (default: genomes)
#   --output_dir DIR       Output directory for analysis (default: 04-marker_density)
#   --parent1 NAME         Name of parent 1 accession (default: Col-0)
#   --parent2 NAME         Name of parent 2 accession (default: Ler-0)
#   --kmer_sizes SIZES     Comma-separated k-mer sizes (default: 21,31,41)
#   --threads N            Number of threads for KMC (default: 20)
#   --fasta_suffix SUFFIX  FASTA file suffix (default: .ragtag_scaffolds.fa)
#   --skip-cen             Skip centromere mode (chromosome-level hapmers)
#   --help                 Show this help message
#
# Requirements:
#   - KMC tools (kmc_dump, kmc_tools)
#   - Raw k-mer databases from step 01
#   - Unique hapmer databases from step 03
#   - Parent genome FASTA files
#
# Output:
#   - ${output_dir}/k${K}/raw_kmers_summary.csv
#   - ${output_dir}/k${K}/unique_hapmers_summary.csv
#   - ${output_dir}/k${K}/comparison_summary.csv
###############################################################################

set -euo pipefail

# Default parameters
RAW_KMER_DIR="01-all_kmers_per_Chr_CEN_and_accesion"
CENHAPMER_DIR="03-cenhapmers"
GENOME_DIR="genomes"
OUTPUT_DIR="04-marker_density"
PARENT1="Col-0"
PARENT2="Ler-0"
KMER_SIZES="21,31,41"
THREADS=20
FASTA_SUFFIX=".ragtag_scaffolds.fa"
SKIP_CEN=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --raw_kmer_dir)
            RAW_KMER_DIR="$2"
            shift 2
            ;;
        --cenhapmer_dir)
            CENHAPMER_DIR="$2"
            shift 2
            ;;
        --genome_dir)
            GENOME_DIR="$2"
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
        --kmer_sizes)
            KMER_SIZES="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --fasta_suffix)
            FASTA_SUFFIX="$2"
            shift 2
            ;;
        --skip-cen)
            SKIP_CEN=true
            shift
            ;;
        --help)
            head -n 34 "$0" | grep "^#" | sed 's/^# //'
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

echo "╔════════════════════════════════════════════════════════════╗"
echo "║        Marker Density Analysis (RAW vs UNIQUE)            ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
echo "Configuration:"
echo "  Raw k-mer directory:   $RAW_KMER_DIR"
echo "  Unique hapmer directory: $CENHAPMER_DIR"
echo "  Genome directory:      $GENOME_DIR"
echo "  Output directory:      $OUTPUT_DIR"
echo "  Parent 1:              $PARENT1"
echo "  Parent 2:              $PARENT2"
echo "  K-mer sizes:           ${K_SIZES_ARRAY[@]}"
echo "  Threads:               $THREADS"
echo "  Skip centromeres:      $SKIP_CEN"
echo "════════════════════════════════════════════════════════════"
echo ""

# Check dependencies
echo "Checking dependencies..."
if ! command -v kmc_tools &> /dev/null; then
    echo "❌ Error: kmc_tools is not installed or not in PATH"
    exit 1
fi
if ! command -v kmc_dump &> /dev/null; then
    echo "❌ Error: kmc_dump is not installed or not in PATH"
    exit 1
fi
echo "  ✅ kmc_tools found"
echo "  ✅ kmc_dump found"
echo ""

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Function to get genome size
get_genome_size() {
    local fasta_file=$1
    # Count all non-header characters
    grep -v "^>" "$fasta_file" | tr -d '\n' | wc -c
}

# Function to analyze k-mer database
analyze_kmer_db() {
    local db_path=$1
    local db_name=$2
    local output_dir=$3
    local genome_size=$4

    # Use kmc_dump to count unique k-mers
    local unique_kmers=$(kmc_dump "$db_path" /dev/stdout 2>/dev/null | wc -l)

    # Use kmc_tools transform to get histogram
    kmc_tools transform "$db_path" histogram "${output_dir}/${db_name}.hist" 2>/dev/null || {
        echo "0,0,0,0.0"
        return
    }

    # Calculate total k-mers (sum of all occurrences)
    local total_kmers=$(awk '{sum+=$1*$2} END {print sum}' "${output_dir}/${db_name}.hist" 2>/dev/null || echo "0")

    # Calculate normalized values (per Mbp)
    local genome_mbp=$(awk "BEGIN {printf \"%.2f\", $genome_size/1000000}")
    local unique_per_mbp=$(awk "BEGIN {printf \"%.2f\", $unique_kmers/$genome_mbp}")
    local total_per_mbp=$(awk "BEGIN {printf \"%.2f\", $total_kmers/$genome_mbp}")

    # Return: unique_kmers,total_kmers,unique_per_mbp,total_per_mbp
    echo "$unique_kmers,$total_kmers,$unique_per_mbp,$total_per_mbp"
}

# Process each k-mer size
for K in "${K_SIZES_ARRAY[@]}"; do
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║ Processing k-mer size: $K"
    echo "╚════════════════════════════════════════════════════════════╝"
    echo ""

    RAW_K_DIR="${RAW_KMER_DIR}/k${K}"
    UNIQUE_K_DIR="${CENHAPMER_DIR}/k${K}"
    OUTPUT_K_DIR="${OUTPUT_DIR}/k${K}"
    mkdir -p "$OUTPUT_K_DIR"

    # Create CSV headers
    echo "database,parent,region,chromosome,unique_kmers,total_kmers,unique_per_mbp,total_per_mbp,genome_size_bp" > "${OUTPUT_K_DIR}/raw_kmers_summary.csv"
    echo "database,parent,region,chromosome,unique_kmers,total_kmers,unique_per_mbp,total_per_mbp,genome_size_bp" > "${OUTPUT_K_DIR}/unique_hapmers_summary.csv"

    # Process each parent
    for PARENT in "$PARENT1" "$PARENT2"; do
        echo "  📊 Processing $PARENT..."

        # Genome FASTA file
        GENOME_FASTA="${GENOME_DIR}/${PARENT}${FASTA_SUFFIX}"

        if [[ ! -f "$GENOME_FASTA" ]]; then
            echo "    ⚠️  Warning: Genome FASTA not found: $GENOME_FASTA"
            echo "    Skipping $PARENT for k=$K"
            continue
        fi

        # Calculate genome size once per parent
        echo "    📏 Calculating genome size..."
        GENOME_SIZE=$(get_genome_size "$GENOME_FASTA")
        GENOME_MBP=$(awk "BEGIN {printf \"%.2f\", $GENOME_SIZE/1000000}")
        echo "       Genome size: ${GENOME_SIZE} bp (${GENOME_MBP} Mbp)"

        # ═══════════════════════════════════════════════════════════
        # PART 1: Analyze RAW k-mers from step 01
        # ═══════════════════════════════════════════════════════════
        echo ""
        echo "    🧬 Analyzing RAW k-mers (step 01)..."

        # Raw k-mers are organized as: {RAW_KMER_DIR}/{PARENT}/{CHR}/{REGION}/{PARENT}_{REGION}_{CHR}_k{K}.kmc_*
        PARENT_RAW_DIR="${RAW_KMER_DIR}/${PARENT}"

        if [[ -d "$PARENT_RAW_DIR" ]]; then
            shopt -s nullglob
            # Find all .kmc_pre files for this parent and k-mer size
            RAW_FILES=("${PARENT_RAW_DIR}"/*/*/*_k${K}.kmc_pre)
            shopt -u nullglob

            if [[ ${#RAW_FILES[@]} -gt 0 ]]; then
                echo "       Found ${#RAW_FILES[@]} raw k-mer databases"

                for PRE_FILE in "${RAW_FILES[@]}"; do
                    BASE=$(basename "$PRE_FILE" .kmc_pre)
                    DIR=$(dirname "$PRE_FILE")
                    FULLPATH="${DIR}/${BASE}"

                    # Extract region and chromosome from filename
                    # Format: {PARENT}_{REGION}_{CHR}_k{K}
                    REGION=$(echo "$BASE" | cut -d'_' -f2)
                    CHR=$(echo "$BASE" | cut -d'_' -f3)

                    echo "         Analyzing ${BASE}..."

                    # Get statistics
                    STATS=$(analyze_kmer_db "$FULLPATH" "${BASE}_raw" "$OUTPUT_K_DIR" "$GENOME_SIZE")
                    IFS=',' read -r UNIQUE TOTAL UNIQUE_NORM TOTAL_NORM <<< "$STATS"

                    echo "           Unique: ${UNIQUE} (${UNIQUE_NORM}/Mbp) | Total: ${TOTAL} (${TOTAL_NORM}/Mbp)"

                    # Save to CSV
                    echo "${BASE},${PARENT},${REGION},${CHR},${UNIQUE},${TOTAL},${UNIQUE_NORM},${TOTAL_NORM},${GENOME_SIZE}" >> "${OUTPUT_K_DIR}/raw_kmers_summary.csv"
                done
            else
                echo "       ⚠️  No raw k-mer files found for k=${K}"
            fi
        else
            echo "       ⚠️  Raw k-mer directory not found: $PARENT_RAW_DIR"
        fi

        # ═══════════════════════════════════════════════════════════
        # PART 2: Analyze UNIQUE hapmers from step 03
        # ═══════════════════════════════════════════════════════════
        echo ""
        echo "    ✨ Analyzing UNIQUE hapmers (step 03)..."

        if [[ -d "$UNIQUE_K_DIR" ]]; then
            shopt -s nullglob
            UNIQUE_FILES=("${UNIQUE_K_DIR}"/unique_${PARENT}_*_k${K}.kmc_pre)
            shopt -u nullglob

            if [[ ${#UNIQUE_FILES[@]} -gt 0 ]]; then
                echo "       Found ${#UNIQUE_FILES[@]} unique hapmer databases"

                for PRE_FILE in "${UNIQUE_FILES[@]}"; do
                    BASE=$(basename "$PRE_FILE" .kmc_pre)
                    FULLPATH="${UNIQUE_K_DIR}/${BASE}"

                    # Extract region and chromosome from filename
                    # Format: unique_{PARENT}_{REGION}_{CHR}_k{K}
                    REGION=$(echo "$BASE" | cut -d'_' -f3)
                    CHR=$(echo "$BASE" | cut -d'_' -f4)

                    echo "         Analyzing ${BASE}..."

                    # Get statistics
                    STATS=$(analyze_kmer_db "$FULLPATH" "${BASE}_unique" "$OUTPUT_K_DIR" "$GENOME_SIZE")
                    IFS=',' read -r UNIQUE TOTAL UNIQUE_NORM TOTAL_NORM <<< "$STATS"

                    echo "           Unique: ${UNIQUE} (${UNIQUE_NORM}/Mbp) | Total: ${TOTAL} (${TOTAL_NORM}/Mbp)"

                    # Save to CSV
                    echo "${BASE},${PARENT},${REGION},${CHR},${UNIQUE},${TOTAL},${UNIQUE_NORM},${TOTAL_NORM},${GENOME_SIZE}" >> "${OUTPUT_K_DIR}/unique_hapmers_summary.csv"
                done
            else
                echo "       ⚠️  No unique hapmer files found"
            fi
        else
            echo "       ⚠️  Unique hapmer directory not found: $UNIQUE_K_DIR"
        fi

        echo ""
    done

    echo "  ✅ k=${K} analysis complete!"
    echo ""
done

echo "╔════════════════════════════════════════════════════════════╗"
echo "║             Analysis Complete!                             ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
echo "Results saved to: $OUTPUT_DIR"
echo ""
echo "Generated files:"
for K in "${K_SIZES_ARRAY[@]}"; do
    echo "  k=${K}:"
    if [[ -f "${OUTPUT_DIR}/k${K}/raw_kmers_summary.csv" ]]; then
        echo "    ✅ raw_kmers_summary.csv"
    fi
    if [[ -f "${OUTPUT_DIR}/k${K}/unique_hapmers_summary.csv" ]]; then
        echo "    ✅ unique_hapmers_summary.csv"
    fi
done
echo ""
echo "Next step: Generate comparison plots"
echo "  python3 scripts/07-plot_marker_comparison.py \\"
echo "    --marker_dir $OUTPUT_DIR \\"
echo "    --output_dir 05-marker_plots \\"
echo "    --parent1 $PARENT1 \\"
echo "    --parent2 $PARENT2 \\"
echo "    --kmer_sizes $KMER_SIZES"
echo ""
