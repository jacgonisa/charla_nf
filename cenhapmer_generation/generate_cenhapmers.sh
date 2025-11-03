#!/bin/bash

###############################################################################
# Master Script: generate_cenhapmers.sh
# Description: Complete pipeline for generating cenhapmer databases from
#              parent genome assemblies
#
# This script coordinates cenhapmer generation with optional analysis:
#   1. Split genomes into CEN and ARMS regions
#   2. Generate k-mer databases for each region
#   3. Create KMC operation files
#   4. Execute operations to find unique k-mers (cenhapmers)
#   5. (Optional) Analyze marker density statistics
#   6. (Optional) Generate genome-wide marker coverage visualizations
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
#   --fasta_suffix SUFFIX  FASTA file suffix (default: .ragtag_scaffolds.fa)
#   --bed_dir DIR          Directory with BED files (default: annotations)
#   --threads N            Number of threads for KMC (default: 4)
#   --parallel N           Parallel kmc_tools processes (default: 1)
#   --analyze-markers      Run marker density analysis (optional, recommended)
#   --visualize-coverage   Run marker coverage visualization (optional, creates plots)
#   --skip-cen             Skip centromere annotation (generate chromosome-level hapmers only)
#   --clean                Remove intermediate files after completion
#   --skip-step N          Skip step N (1-6)
#   --help                 Show this help message
#
# Required Input Files:
#   genomes/${PARENT}.ragtag_scaffolds.fa
#   annotations/${PARENT}.ragtag_scaffolds_CEN.bed  (if --skip-cen not used)
#   annotations/${PARENT}.ragtag_scaffolds_ARMS.bed (if --skip-cen not used)
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
FASTA_SUFFIX=".ragtag_scaffolds.fa"
BED_DIR="annotations"
THREADS=4
PARALLEL=1
ANALYZE_MARKERS=false
VISUALIZE_COVERAGE=false
CLEAN=false
SKIP_CEN=false
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
        --fasta_suffix)
            FASTA_SUFFIX="$2"
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
        --analyze-markers)
            ANALYZE_MARKERS=true
            shift
            ;;
        --visualize-coverage)
            VISUALIZE_COVERAGE=true
            shift
            ;;
        --skip-cen)
            SKIP_CEN=true
            shift
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
echo "  Skip centromeres:   $SKIP_CEN"
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

    SPLIT_CMD="${SCRIPT_DIR}/01-split_genomes.sh \
        --fasta_dir $FASTA_DIR \
        --fasta_suffix $FASTA_SUFFIX \
        --bed_dir $BED_DIR \
        --out_dir genomes_split \
        --parent1 $PARENT1 \
        --parent2 $PARENT2 \
        --num_chr $NUM_CHR \
        --chr_prefix $CHR_PREFIX"

    if [[ "$SKIP_CEN" == "true" ]]; then
        SPLIT_CMD="$SPLIT_CMD --skip-cen"
    fi

    bash $SPLIT_CMD

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

    KMER_CMD="${SCRIPT_DIR}/02-generate_kmers.sh \
        --input_dir genomes_split \
        --output_dir 01-all_kmers_per_Chr_CEN_and_accesion \
        --parent1 $PARENT1 \
        --parent2 $PARENT2 \
        --num_chr $NUM_CHR \
        --chr_prefix $CHR_PREFIX \
        --kmer_sizes $KMER_SIZES \
        --threads $THREADS"

    if [[ "$SKIP_CEN" == "true" ]]; then
        KMER_CMD="$KMER_CMD --skip-cen"
    fi

    bash $KMER_CMD

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

    OP_CMD="python3 ${SCRIPT_DIR}/03-generate_kmc_operations.py \
        --input_dir 01-all_kmers_per_Chr_CEN_and_accesion \
        --output_dir 02-kmc_complex_op_files \
        --parent1 $PARENT1 \
        --parent2 $PARENT2 \
        --num_chr $NUM_CHR \
        --chr_prefix $CHR_PREFIX \
        --kmer_sizes $KMER_SIZES"

    if [[ "$SKIP_CEN" == "true" ]]; then
        OP_CMD="$OP_CMD --skip-cen"
    fi

    $OP_CMD

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
# Step 5 (Optional): Analyze marker density
# ============================================================================
if [[ "$ANALYZE_MARKERS" == "true" ]]; then
    if should_skip 5; then
        echo "⏭️  Skipping Step 5: Marker density analysis"
    else
        echo "╔════════════════════════════════════════════════════════════╗"
        echo "║ Step 5/5: Analyze marker density (optional)               ║"
        echo "╚════════════════════════════════════════════════════════════╝"
        echo ""

        ANALYZE_CMD="bash ${SCRIPT_DIR}/05-analyze_marker_density.sh \
            --cenhapmer_dir 03-cenhapmers \
            --genome_dir $FASTA_DIR \
            --output_dir 04-marker_density \
            --parent1 $PARENT1 \
            --parent2 $PARENT2 \
            --kmer_sizes $KMER_SIZES \
            --threads $THREADS \
            --fasta_suffix $FASTA_SUFFIX"

        PLOT_CMD="python3 ${SCRIPT_DIR}/06-plot_marker_density.py \
            --marker_dir 04-marker_density \
            --output_dir 05-marker_plots \
            --parent1 $PARENT1 \
            --parent2 $PARENT2 \
            --kmer_sizes $KMER_SIZES \
            --num_chr $NUM_CHR \
            --chr_prefix $CHR_PREFIX"

        if [[ "$SKIP_CEN" == "true" ]]; then
            ANALYZE_CMD="$ANALYZE_CMD --skip-cen"
            PLOT_CMD="$PLOT_CMD --skip-cen"
        fi

        $ANALYZE_CMD

        echo ""
        echo "Generating visualizations..."
        $PLOT_CMD

        echo ""
        echo "✅ Step 5 complete!"
        echo ""
        echo "📊 Marker density analysis results:"
        echo "   Analysis data: 04-marker_density/"
        echo "   Plots:         05-marker_plots/"
        echo ""
    fi
fi

# ============================================================================
# Step 6 (Optional): Marker Coverage Visualization
# ============================================================================
if [[ "$VISUALIZE_COVERAGE" == "true" ]]; then
    if should_skip 6; then
        echo "⏭️  Skipping Step 6: Marker coverage visualization"
    else
        echo "╔════════════════════════════════════════════════════════════╗"
        echo "║ Step 6/6: Marker coverage visualization (optional)        ║"
        echo "╚════════════════════════════════════════════════════════════╝"
        echo ""

        # Step 6a: Analyze marker coverage (map hapmers to genome positions)
        COVERAGE_CMD="bash ${SCRIPT_DIR}/05-analyze_marker_coverage.sh \
            --cenhapmer_dir 03-cenhapmers \
            --genome_dir $FASTA_DIR \
            --output_dir 04-marker_coverage \
            --parent1 $PARENT1 \
            --parent2 $PARENT2 \
            --kmer_sizes $KMER_SIZES \
            --threads $THREADS \
            --fasta_suffix $FASTA_SUFFIX"

        # Step 6b: Bin coverage into windows
        BIN_CMD="python3 ${SCRIPT_DIR}/06-bin_marker_coverage.py \
            --coverage_dir 04-marker_coverage \
            --output_dir 04-marker_coverage_plots \
            --window_size 1000 \
            --kmer_sizes $KMER_SIZES"

        # Step 6c: Generate plots with improved multi-row layout
        PLOT_CMD="python3 ${SCRIPT_DIR}/07-plot_marker_coverage.py"

        if [[ "$SKIP_CEN" == "true" ]]; then
            COVERAGE_CMD="$COVERAGE_CMD --skip-cen"
            BIN_CMD="$BIN_CMD --skip-cen"
        fi

        echo "Step 6a: Analyzing marker coverage..."
        $COVERAGE_CMD

        echo ""
        echo "Step 6b: Binning coverage into windows..."
        $BIN_CMD

        echo ""
        echo "Step 6c: Generating visualizations..."

        # Parse k-mer sizes into array
        IFS=',' read -ra K_SIZES_ARRAY <<< "$KMER_SIZES"

        # Generate plots for each parent and k-mer size
        for K in "${K_SIZES_ARRAY[@]}"; do
            for PARENT in "$PARENT1" "$PARENT2"; do
                echo "  📊 Plotting $PARENT k=$K..."

                CSV_FILE="04-marker_coverage_plots/k${K}_marker_counts_1000bp.csv"

                if [[ ! -f "$CSV_FILE" ]]; then
                    echo "    ⚠️  Warning: CSV file not found: $CSV_FILE"
                    continue
                fi

                # Determine genome file
                GENOME_FILE="${FASTA_DIR}/${PARENT}${FASTA_SUFFIX}"

                if [[ ! -f "$GENOME_FILE" ]]; then
                    echo "    ⚠️  Warning: Genome file not found: $GENOME_FILE"
                    continue
                fi

                # Run plotting with SVG/PDF output and multi-row layout
                python3 "${SCRIPT_DIR}/07-plot_marker_coverage.py" \
                    -i "$CSV_FILE" \
                    -s "$PARENT" \
                    --kmer "$K" \
                    -g "$GENOME_FILE" \
                    -w 1000 \
                    -o 04-marker_coverage_plots \
                    --format png svg pdf \
                    --ncols 5 \
                    --log-scale \
                    --max-points 1000 || echo "    ⚠️  Warning: Plotting failed for $PARENT k=$K"
            done
        done

        echo ""
        echo "✅ Step 6 complete!"
        echo ""
        echo "📊 Marker coverage visualization results:"
        echo "   Coverage data:  04-marker_coverage/"
        echo "   Binned data:    04-marker_coverage_plots/*.csv"
        echo "   Plots (PNG/SVG/PDF): 04-marker_coverage_plots/"
        echo ""
    fi
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
