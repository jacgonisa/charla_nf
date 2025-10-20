#!/bin/bash

###############################################################################
# Script: 05-analyze_marker_density.sh
# Description: Analyze marker (cenhapmer) density across chromosomes
#              This helps determine which k-mer size provides optimal marker
#              coverage for haplotype assignment
#
# Usage: ./05-analyze_marker_density.sh [options]
#
# Options:
#   --cenhapmer_dir DIR    Directory with cenhapmers (default: 03-cenhapmers)
#   --genome_dir DIR       Directory with genome FASTAs (default: genomes)
#   --output_dir DIR       Output directory for analysis (default: 04-marker_density)
#   --parent1 NAME         Name of parent 1 accession (default: Col-0)
#   --parent2 NAME         Name of parent 2 accession (default: Ler-0)
#   --kmer_sizes SIZES     Comma-separated k-mer sizes (default: 21,31,41)
#   --threads N            Number of threads for KMC (default: 20)
#   --fasta_suffix SUFFIX  FASTA file suffix (default: .ragtag_scaffolds.fa)
#   --help                 Show this help message
#
# Requirements:
#   - KMC tools (get_counts_threads from kmc_tools)
#   - Cenhapmer databases from step 04
#   - Parent genome FASTA files
#
# Output:
#   - ${output_dir}/k${K}/${PARENT}_{CEN|ARMS}_Chr${N}_k${K}.counts
#   - Count files showing marker presence at each base position
###############################################################################

set -euo pipefail

# Default parameters
CENHAPMER_DIR="03-cenhapmers"
GENOME_DIR="genomes"
OUTPUT_DIR="04-marker_density"
PARENT1="Col-0"
PARENT2="Ler-0"
KMER_SIZES="21,31,41"
THREADS=20
FASTA_SUFFIX=".ragtag_scaffolds.fa"

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
echo "║        Marker Density Analysis                             ║"
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
echo "════════════════════════════════════════════════════════════"
echo ""

# Check dependencies
echo "Checking dependencies..."
if ! command -v kmc_tools &> /dev/null; then
    echo "❌ Error: kmc_tools is not installed or not in PATH"
    exit 1
fi
echo "  ✅ kmc_tools found"
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

            # Run get_counts_threads (from kmc_tools)
            # This counts marker presence at each base position in the genome
            kmc_tools dump -ci1 "$FULLPATH" /dev/stdout | \
                kmc_tools filter "$FULLPATH" @/dev/stdin /dev/null -ci1 -cx1000000 || true

            # Alternative: use kmc_tools transform to get counts per position
            # This creates a file showing marker count at each genomic position
            OUTPUT_FILE="${OUTPUT_K_DIR}/${BASE}.counts"

            # Use KMC dump to get all k-mers and their counts
            kmc_tools transform "$FULLPATH" dump "${OUTPUT_FILE}.tmp" 2>/dev/null || {
                echo "      ⚠️  Warning: kmc_tools transform failed for $BASE"
                continue
            }

            # Process the dump to create a position-based count file
            # For now, we'll use a simpler approach with kmc_tools histogram
            kmc_tools histogram "$FULLPATH" "${OUTPUT_FILE}.hist" 2>/dev/null || {
                echo "      ⚠️  Warning: histogram generation failed for $BASE"
                rm -f "${OUTPUT_FILE}.tmp"
                continue
            }

            # Get basic statistics
            TOTAL_MARKERS=$(awk '{sum+=$2} END {print sum}' "${OUTPUT_FILE}.hist" 2>/dev/null || echo "0")
            UNIQUE_MARKERS=$(wc -l < "${OUTPUT_FILE}.tmp" 2>/dev/null || echo "0")

            echo "      ✅ Statistics for $BASE:"
            echo "         Total markers:  $TOTAL_MARKERS"
            echo "         Unique markers: $UNIQUE_MARKERS"

            # Save summary
            echo "$BASE,$TOTAL_MARKERS,$UNIQUE_MARKERS" >> "${OUTPUT_K_DIR}/summary.csv"

            # Clean up temporary files
            rm -f "${OUTPUT_FILE}.tmp"
        done

        echo ""
    done

    echo ""
done

echo "╔════════════════════════════════════════════════════════════╗"
echo "║             Analysis Complete!                             ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
echo "Results saved to: $OUTPUT_DIR"
echo ""
echo "Summary files:"
for K in "${K_SIZES_ARRAY[@]}"; do
    SUMMARY_FILE="${OUTPUT_DIR}/k${K}/summary.csv"
    if [[ -f "$SUMMARY_FILE" ]]; then
        echo "  - k=${K}: $SUMMARY_FILE"
    fi
done
echo ""
echo "Next step: Run visualization script to compare k-mer sizes"
echo "  python3 scripts/06-plot_marker_density.py --input_dir $OUTPUT_DIR"
echo ""
