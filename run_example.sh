#!/bin/bash
# Example script for running CHARLA mapping pipeline
# Modify the parameters below for your specific dataset

set -e

# Configuration - MODIFY THESE
INPUT_FILE="../data/your_sample.fq.gz"
SAMPLE_ID="your_sample_name"
OUTPUT_DIR="results_${SAMPLE_ID}"

# Reference files (modify if using different references)
REFERENCE="./Col-0_chr1-5_renamed.fa"
REFERENCE_DIR="./index"
KMER_DB="./merylDB"

# Optional BED files for chromosomal context (set to empty string if not needed)
PARENT1_BED=""
PARENT2_BED=""

echo "========================================="
echo "Running CHARLA mapping pipeline"
echo "Sample: ${SAMPLE_ID}"
echo "========================================="

# Build the command
CMD="nextflow run main.nf \
    --input ${INPUT_FILE} \
    --sample_id ${SAMPLE_ID} \
    --outdir ${OUTPUT_DIR} \
    --reference ${REFERENCE} \
    --reference_dir ${REFERENCE_DIR} \
    --kmer_db ${KMER_DB} \
    --max_cpus 20 \
    --max_memory '100.GB' \
    -profile singularity \
    -resume"

# Add BED files if provided
if [ -n "${PARENT1_BED}" ]; then
    CMD="${CMD} --parent1_bed ${PARENT1_BED}"
fi

if [ -n "${PARENT2_BED}" ]; then
    CMD="${CMD} --parent2_bed ${PARENT2_BED}"
fi

# Run the pipeline
eval ${CMD}

echo ""
echo "========================================="
echo "Pipeline completed!"
echo "Results in: ${OUTPUT_DIR}/"
echo "========================================="
