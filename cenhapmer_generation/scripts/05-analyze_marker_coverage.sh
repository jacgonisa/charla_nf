#!/bin/bash

###############################################################################
# Script: 05-analyze_marker_coverage.sh
# Description: Analyze per-base marker (cenhapmer) coverage across chromosomes
#              Uses get_counts_threads to map k-mers back to genome positions
#
# Usage: ./05-analyze_marker_coverage.sh [options]
#
# Options:
#   --cenhapmer_dir DIR    Directory with cenhapmers (default: 03-cenhapmers)
#   --genome_dir DIR       Directory with genome FASTAs (default: genomes)
#   --output_dir DIR       Output directory for analysis (default: 04-marker_coverage)
#   --parent1 NAME         Name of parent 1 accession (default: Col-0)
#   --parent2 NAME         Name of parent 2 accession (default: Ler-0)
#   --kmer_sizes SIZES     Comma-separated k-mer sizes (default: 21,31,41)
#   --threads N            Number of threads for get_counts_threads (default: 20)
#   --fasta_suffix SUFFIX  FASTA file suffix (default: .ragtag_scaffolds.fa)
#   --skip-cen             Skip centromere mode (chromosome-level hapmers)
#   --help                 Show this help message
#
# Requirements:
#   - get_counts_threads (from KMC tools)
#   - Cenhapmer databases from step 04
#   - Parent genome FASTA files
#
# Output:
#   - ${output_dir}/k${K}/${DATABASE_NAME}.counts
#   - Per-base marker coverage files (comma-separated values)
###############################################################################

set -euo pipefail

# Default parameters
CENHAPMER_DIR="03-cenhapmers"
GENOME_DIR="genomes"
OUTPUT_DIR="04-marker_coverage"
PARENT1="Col-0"
PARENT2="Ler-0"
KMER_SIZES="21,31,41"
THREADS=20
FASTA_SUFFIX=".ragtag_scaffolds.fa"
SKIP_CEN=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
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
            head -n 32 "$0" | grep "^#" | sed 's/^# //'
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
echo "║        Per-Base Marker Coverage Analysis                  ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
echo "Configuration:"
echo "  Cenhapmer directory: $CENHAPMER_DIR"
echo "  Genome directory:    $GENOME_DIR"
echo "  Output directory:    $OUTPUT_DIR"
echo "  Parent 1:            $PARENT1"
echo "  Parent 2:            $PARENT2"
echo "  K-mer sizes:         ${K_SIZES_ARRAY[@]}"
echo "  Threads:             $THREADS"
echo "  Skip centromeres:    $SKIP_CEN"
echo "════════════════════════════════════════════════════════════"
echo ""

# Check dependencies
echo "Checking dependencies..."
if ! command -v get_counts_threads &> /dev/null; then
    echo "❌ Error: get_counts_threads is not installed or not in PATH"
    echo "   This tool is part of KMC (k-mer counter) suite"
    echo "   Make sure KMC is properly installed and in your PATH"
    exit 1
fi
echo "  ✅ get_counts_threads found"
echo ""

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Process each k-mer size
for K in "${K_SIZES_ARRAY[@]}"; do
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║ Processing k-mer size: $K"
    echo "╚════════════════════════════════════════════════════════════╝"
    echo ""

    K_DIR="${CENHAPMER_DIR}/k${K}"
    OUTPUT_K_DIR="${OUTPUT_DIR}/k${K}"
    mkdir -p "$OUTPUT_K_DIR"

    if [[ ! -d "$K_DIR" ]]; then
        echo "⚠️  Warning: Directory not found: $K_DIR"
        echo "   Skipping k=$K"
        continue
    fi

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

        # Find all cenhapmer databases for this parent
        shopt -s nullglob
        CENHAPMER_FILES=("${K_DIR}"/unique_${PARENT}_*_k${K}.kmc_pre)
        shopt -u nullglob

        if [[ ${#CENHAPMER_FILES[@]} -eq 0 ]]; then
            echo "    ⚠️  Warning: No cenhapmer files found for $PARENT"
            continue
        fi

        echo "    Found ${#CENHAPMER_FILES[@]} cenhapmer databases"

        # Process each cenhapmer database
        for PRE_FILE in "${CENHAPMER_FILES[@]}"; do
            # Remove .kmc_pre suffix to get base name
            BASE=$(basename "$PRE_FILE" .kmc_pre)
            FULLPATH="${K_DIR}/${BASE}"

            echo "      🔍 Analyzing $BASE..."

            OUTPUT_FILE="${OUTPUT_K_DIR}/${BASE}.counts"

            # Run get_counts_threads to map k-mers to genome positions
            # This creates a file with comma-separated counts for each base position
            if get_counts_threads "$FULLPATH" "$GENOME_FASTA" "$THREADS" > "$OUTPUT_FILE" 2>/dev/null; then
                # Get file size for verification
                FILE_SIZE=$(stat -f%z "$OUTPUT_FILE" 2>/dev/null || stat -c%s "$OUTPUT_FILE" 2>/dev/null || echo "unknown")
                echo "      ✅ Coverage file created: ${BASE}.counts (${FILE_SIZE} bytes)"
            else
                echo "      ⚠️  Warning: get_counts_threads failed for $BASE"
                rm -f "$OUTPUT_FILE"
                continue
            fi
        done

        echo ""
    done

    echo ""
done

echo "╔════════════════════════════════════════════════════════════╗"
echo "║             Per-Base Coverage Analysis Complete!          ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
echo "Results saved to: $OUTPUT_DIR"
echo ""
echo "Next step: Run binning script to generate 1kb window statistics"
echo "  python3 scripts/06-bin_marker_coverage.py --coverage_dir $OUTPUT_DIR"
echo ""
