#!/bin/bash

###############################################################################
# Master Script: generate_cenhapmers.sh
# Description: Complete pipeline for generating cenhapmer databases from
#              parent genome assemblies
#
# This script coordinates all four steps of cenhapmer generation:
#   1. Split genomes into CEN and ARMS regions
#   2. Generate k-mer databases for each region
#   3. Create KMC operation files
#   4. Execute operations to find unique k-mers (cenhapmers)
#
# Usage: ./generate_cenhapmers.sh [options]
#
# Options:
#   --parent1 NAME         Name of parent 1 accession (default: Col-0)
#   --parent2 NAME         Name of parent 2 accession (default: Ler-0)
#   --num_chr N            Number of chromosomes (default: 5)
#   --chr_prefix PREFIX    Chromosome prefix in BED files (default: Chr)
#   --kmer_sizes SIZES     Comma-separated k-mer sizes (default: 21,31,41)
#   --fasta_dir DIR        Directory with genome FASTAs (default: genomes)
#   --bed_dir DIR          Directory with BED files (default: annotations)
#   --threads N            Number of threads for KMC (default: 4)
#   --parallel N           Parallel kmc_tools processes (default: 1)
#   --clean                Remove intermediate files after completion
#   --skip-step N          Skip step N (1-4)
#   --help                 Show this help message
#
# Required Input Files:
#   genomes/${PARENT}.ragtag_scaffolds.fa
#   annotations/${PARENT}.ragtag_scaffolds_CEN.bed
#   annotations/${PARENT}.ragtag_scaffolds_ARMS.bed
#
# Output:
#   03-cenhapmers/k{21,31,41}/unique_${PARENT}_{CEN|ARMS}_Chr${N}_k${K}.kmc_{pre,suf}
###############################################################################

set -euo pipefail

# Default parameters
PARENT1="Col-0"
PARENT2="Ler-0"
NUM_CHR=5
CHR_PREFIX="Chr"
KMER_SIZES="21,31,41"
FASTA_DIR="genomes"
BED_DIR="annotations"
THREADS=4
PARALLEL=1
CLEAN=false
SKIP_STEPS=()

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/scripts"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
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
        --fasta_dir)
            FASTA_DIR="$2"
            shift 2
            ;;
        --bed_dir)
            BED_DIR="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --parallel)
            PARALLEL="$2"
            shift 2
            ;;
        --clean)
            CLEAN=true
            shift
            ;;
        --skip-step)
            SKIP_STEPS+=("$2")
            shift 2
            ;;
        --help)
            head -n 40 "$0" | grep "^#" | sed 's/^# //'
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Function to check if step should be skipped
should_skip() {
    local step=$1
    for skip in "${SKIP_STEPS[@]}"; do
        if [[ "$skip" == "$step" ]]; then
            return 0
        fi
    done
    return 1
}

# Print configuration
echo ""
echo "╔════════════════════════════════════════════════════════════╗"
echo "║        CHARLA Cenhapmer Generation Pipeline               ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
echo "Configuration:"
echo "  Parent 1:           $PARENT1"
echo "  Parent 2:           $PARENT2"
echo "  Chromosomes:        $NUM_CHR (prefix: $CHR_PREFIX)"
echo "  K-mer sizes:        $KMER_SIZES"
echo "  FASTA directory:    $FASTA_DIR"
echo "  BED directory:      $BED_DIR"
echo "  KMC threads:        $THREADS"
echo "  Parallel jobs:      $PARALLEL"
echo "  Clean intermediate: $CLEAN"
echo ""
echo "────────────────────────────────────────────────────────────"
echo ""

# Check required dependencies
echo "Checking dependencies..."
for cmd in kmc kmc_tools bedtools python3; do
    if ! command -v $cmd &> /dev/null; then
        echo "❌ Error: $cmd is not installed or not in PATH"
        exit 1
    fi
    echo "  ✅ $cmd found"
done
echo ""

START_TIME=$(date +%s)

# ============================================================================
# Step 1: Split genomes into CEN and ARMS regions
# ============================================================================
if should_skip 1; then
    echo "⏭️  Skipping Step 1: Split genomes"
else
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║ Step 1/4: Split genomes into CEN and ARMS regions         ║"
    echo "╚════════════════════════════════════════════════════════════╝"
    echo ""

    bash "${SCRIPT_DIR}/01-split_genomes.sh" \
        --fasta_dir "$FASTA_DIR" \
        --bed_dir "$BED_DIR" \
        --out_dir "genomes_split" \
        --parent1 "$PARENT1" \
        --parent2 "$PARENT2" \
        --num_chr "$NUM_CHR" \
        --chr_prefix "$CHR_PREFIX"

    echo ""
    echo "✅ Step 1 complete!"
    echo ""
fi

# ============================================================================
# Step 2: Generate k-mer databases
# ============================================================================
if should_skip 2; then
    echo "⏭️  Skipping Step 2: Generate k-mer databases"
else
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║ Step 2/4: Generate k-mer databases with KMC               ║"
    echo "╚════════════════════════════════════════════════════════════╝"
    echo ""

    bash "${SCRIPT_DIR}/02-generate_kmers.sh" \
        --input_dir "genomes_split" \
        --output_dir "01-all_kmers_per_Chr_CEN_and_accesion" \
        --parent1 "$PARENT1" \
        --parent2 "$PARENT2" \
        --num_chr "$NUM_CHR" \
        --chr_prefix "$CHR_PREFIX" \
        --kmer_sizes "$KMER_SIZES" \
        --threads "$THREADS"

    echo ""
    echo "✅ Step 2 complete!"
    echo ""
fi

# ============================================================================
# Step 3: Generate KMC operation files
# ============================================================================
if should_skip 3; then
    echo "⏭️  Skipping Step 3: Generate KMC operation files"
else
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║ Step 3/4: Generate KMC operation files                    ║"
    echo "╚════════════════════════════════════════════════════════════╝"
    echo ""

    python3 "${SCRIPT_DIR}/03-generate_kmc_operations.py" \
        --input_dir "01-all_kmers_per_Chr_CEN_and_accesion" \
        --output_dir "02-kmc_complex_op_files" \
        --parent1 "$PARENT1" \
        --parent2 "$PARENT2" \
        --num_chr "$NUM_CHR" \
        --chr_prefix "$CHR_PREFIX" \
        --kmer_sizes "$KMER_SIZES"

    echo ""
    echo "✅ Step 3 complete!"
    echo ""
fi

# ============================================================================
# Step 4: Execute KMC operations to generate cenhapmers
# ============================================================================
if should_skip 4; then
    echo "⏭️  Skipping Step 4: Generate cenhapmers"
else
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║ Step 4/4: Execute KMC operations (generate cenhapmers)    ║"
    echo "╚════════════════════════════════════════════════════════════╝"
    echo ""

    python3 "${SCRIPT_DIR}/04-run_kmc_operations.py" \
        --input_dir "02-kmc_complex_op_files" \
        --output_dir "03-cenhapmers" \
        --kmer_sizes "$KMER_SIZES" \
        --parallel "$PARALLEL"

    echo ""
    echo "✅ Step 4 complete!"
    echo ""
fi

# ============================================================================
# Optional: Clean up intermediate files
# ============================================================================
if [[ "$CLEAN" == "true" ]]; then
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║ Cleaning up intermediate files                            ║"
    echo "╚════════════════════════════════════════════════════════════╝"
    echo ""

    echo "Removing intermediate directories..."
    rm -rf genomes_split
    rm -rf 01-all_kmers_per_Chr_CEN_and_accesion
    rm -rf 02-kmc_complex_op_files

    echo "✅ Cleanup complete!"
    echo ""
fi

# ============================================================================
# Final summary
# ============================================================================
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
HOURS=$((DURATION / 3600))
MINUTES=$(((DURATION % 3600) / 60))
SECONDS=$((DURATION % 60))

echo "╔════════════════════════════════════════════════════════════╗"
echo "║                 PIPELINE COMPLETE! 🎉                      ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
echo "Total runtime: ${HOURS}h ${MINUTES}m ${SECONDS}s"
echo ""
echo "Cenhapmer databases generated:"
echo "  Directory: 03-cenhapmers/"
echo "  K-mer sizes: $KMER_SIZES"
echo "  Accessions: $PARENT1, $PARENT2"
echo "  Chromosomes: ${NUM_CHR} (${CHR_PREFIX}1-${CHR_PREFIX}${NUM_CHR})"
echo ""
echo "Next steps:"
echo "  1. Verify cenhapmer databases:"
echo "     ls -lh 03-cenhapmers/k*/*.kmc_*"
echo ""
echo "  2. Use these cenhapmers with CHARLA pipeline:"
echo "     nextflow run jacgonisa/charla_nf \\"
echo "       --reads your_reads.fq \\"
echo "       --cenhapmer_db_dir 03-cenhapmers/k21 \\"
echo "       --parent1_name $PARENT1 \\"
echo "       --parent2_name $PARENT2 \\"
echo "       ..."
echo ""
echo "For more information, see: cenhapmer_generation/docs/README.md"
echo ""
