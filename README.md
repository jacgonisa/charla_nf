# CHARLA Mapping Pipeline

A simplified mapping-only pipeline based on [CHARLA_NF](https://github.com/jacgonisa/charla_nf) for analyzing long-read sequencing data with **single-parent k-mer classification** but without recombination analysis.

## Overview

This pipeline runs four key analysis modules from CHARLA_NF with **single-parent mode**:

1. **K-mer Profiling & Classification** - Classify reads based on k-mer profiles against a single reference
2. **NON_HYBRID_READS_ANALYSIS** - Classify and map reads based on k-mer dominance
3. **LARGE_INDEL_ANALYSIS** - Extract and analyze large indels (>50bp), identify satellite structures
4. **MAFFT_ALIGN** - Multiple sequence alignment of satellite blocks
5. **ANALYZE_MAFFT_ALIGNMENT** - Analyze alignment results and update satellite analysis

Perfect for analyzing structural variations, indels, and satellite DNA in long-read sequencing data using **k-mer-based read classification** without the overhead of crossover/recombination detection between two parents.

## Pipeline Flow

```
Long-read FASTQ/FASTA Input
    ↓
FASTQ → FASTA conversion (if needed)
    ↓
Header simplification
    ↓
K-mer generation (KMC)
    └─ Generate k-mers from each read
    ↓
Cenhapmer k-mer counting (against single-parent merylDB)
    └─ Count k-mers matching reference ARMS/CEN regions
    ↓
K-mer matrix processing
    └─ Create TSV matrix: read_id × k-mer_counts
    ↓
NON_HYBRID_READS_ANALYSIS
    ├─ Read classification (based on k-mer dominance)
    │   └─ Assign reads: "Col_ARMS_Chr1", "Col_CEN_Chr2", etc.
    ├─ Sequence extraction (seqkit)
    ├─ Mapping to reference (minimap2)
    └─ Indel analysis (from BAM files)
    ↓
LARGE_INDEL_ANALYSIS
    ├─ Extract large indels (>50bp)
    ├─ Map indels to genome (PAF format)
    ├─ Analyze satellite structures (178bp repeats)
    └─ Extract satellite blocks (FASTA)
    ↓
MAFFT_ALIGN
    └─ Multiple sequence alignment of satellite blocks
    ↓
ANALYZE_MAFFT_ALIGNMENT
    └─ Update satellite analysis with alignment metrics
```

## Single-Parent K-mer Classification

This pipeline performs **k-mer-based read classification** against a **single reference genome** (e.g., Col-0):

1. **K-mer counting**: For each read, count k-mers matching reference ARMS/CEN regions
2. **Classification**: Assign reads to chromosome regions based on k-mer dominance
3. **Mapping**: Map classified reads to reference genome with minimap2
4. **Structural analysis**: Identify indels and satellite structures

Unlike full CHARLA_NF (which detects crossovers between two parents), this pipeline:
- ✅ Keeps k-mer-based read classification
- ✅ Analyzes structural variants and indels
- ✅ Studies satellite DNA organization
- ❌ Does NOT detect recombination/crossovers between parents
- ❌ Does NOT require two-parent k-mer databases

## Requirements

- Nextflow >= 23.04.0
- Singularity or Docker
- Reference genome (FASTA)
- K-mer database (merylDB format)
- Reference index directory

## Quick Start

### 1. Setup

Make sure you have your reference genome, index, and k-mer database ready:

```bash
cd /path/to/charla_mapping

# Link or copy your reference files
ln -s /path/to/reference.fa ./
ln -s /path/to/index ./
ln -s /path/to/merylDB ./
```

### 2. Run the pipeline

**Basic usage:**
```bash
nextflow run main.nf \
    --input data/long_reads.fq.gz \
    --sample_id sample1 \
    --outdir results \
    -profile singularity
```

**With custom reference:**
```bash
nextflow run main.nf \
    --input data/mutant_reads.fq.gz \
    --sample_id mutant1 \
    --reference /path/to/my_reference.fa \
    --reference_dir /path/to/my_index \
    --kmer_db /path/to/my_merylDB \
    --outdir results_mutant1 \
    -profile singularity
```

### 3. Resume a failed run

```bash
nextflow run main.nf \
    --input data/long_reads.fq.gz \
    --sample_id sample1 \
    -profile singularity \
    -resume
```

## Parameters

### Required
- `--input` - Path to input FASTQ or FASTA file (long-read sequencing data)
- `--sample_id` - Sample identifier (e.g., 'sample1', 'mutant_rep1')

### Optional
- `--outdir` - Output directory (default: 'results')
- `--reference` - Reference genome FASTA (default: './Col-0_chr1-5_renamed.fa')
- `--reference_dir` - Reference directory for mapping (default: './index')
- `--kmer_db` - K-mer database directory (default: './merylDB')
- `--kmer_size` - K-mer size (default: 31)
- `--input_format` - 'fastq' or 'fasta' (default: 'fastq')
- `--parent1_bed` - BED file for chromosomal features (optional)
- `--parent2_bed` - BED file for second parent/reference (optional)

### Resource Options
- `--max_cpus` - Maximum CPUs (default: 20)
- `--max_memory` - Maximum memory (default: '100.GB')
- `--max_time` - Maximum time (default: '48.h')

Module-specific resources can be adjusted in `nextflow.config`:
- `--readmer_cpus`, `--readmer_memory`, `--readmer_time`
- `--cenhapmer_cpus`, `--cenhapmer_memory`, `--cenhapmer_time`
- `--nonhybrid_cpus`, `--nonhybrid_memory`, `--nonhybrid_time`
- `--large_indel_cpus`, `--large_indel_memory`, `--large_indel_time`

## Output Structure

```
results/
├── fastq_to_fasta/           # FASTA conversion (if input was FASTQ)
├── simplify_headers/         # Simplified FASTA
├── generate_readmers_kmc/    # K-mer generation
├── get_counts_cenhapmer/     # K-mer counts
├── process_cenhapmer_counts/ # K-mer matrix
├── non_hybrid_analysis/      # Non-hybrid reads analysis
│   └── <sample_id>_nonhybrid/
│       ├── read_lists/       # Read classification lists
│       ├── sequences/        # Extracted FASTA sequences
│       ├── bam/             # Mapped BAM files
│       └── indel_analysis.tsv
├── large_indel_analysis/     # Large indel analysis
│   └── <sample_id>_large_indels/
│       ├── large_indels_catalog.tsv
│       ├── indel_sequences.fa
│       ├── indel_mappings.paf
│       ├── genomic_plots/
│       ├── satellite_blocks.fa
│       ├── satellite_analysis.tsv
│       └── msa_results/
│           ├── <sample_id>_satellite_blocks.aln
│           ├── <sample_id>_satellite_analysis_updated.tsv
│           └── msa_visualization.png
└── pipeline_info/            # Execution reports
    ├── execution_timeline.html
    ├── execution_report.html
    ├── execution_trace.txt
    └── pipeline_dag.svg
```

## Key Output Files

### Non-hybrid Reads Analysis
- `read_lists/*.txt` - Lists of reads for each haplotype/category
- `sequences/*.fa` - Extracted sequences per category
- `bam/*.bam` - Mapped reads (BAM format)
- `indel_analysis.tsv` - Summary statistics of indels

### Large Indel Analysis
- `large_indels_catalog.tsv` - Catalog of all indels >50bp
  - Columns: indel_id, type, size, chromosome, position, sequence, etc.
- `indel_sequences.fa` - FASTA file with indel sequences
- `indel_mappings.paf` - Genomic positions of indels (PAF format)
- `satellite_blocks.fa` - Extracted satellite block sequences
- `satellite_analysis.tsv` - Satellite structure statistics

### MAFFT Alignment
- `msa_results/*.aln` - Multiple sequence alignment of satellite blocks
- `*_satellite_analysis_updated.tsv` - Satellite analysis enriched with alignment metrics
- `msa_visualization.png` - Visual representation of alignment (optional)

## Use Cases

This pipeline is suitable for:

- Analyzing structural variations in mutant vs wild-type comparisons
- Detecting and characterizing large indels (>50bp) in long-read data
- Studying satellite DNA organization and variation
- Mapping long reads without crossover detection overhead
- General long-read sequencing analysis focusing on structural variants

## Notes

### BED Files (Optional)

BED files for chromosomal features are optional. If provided, they enable:
- Better genomic context visualization
- Enhanced plots showing indel positions on chromosomes

If not provided, analysis proceeds without genomic context plots.

### K-mer Database (cenhapmer_db_dir)

**IMPORTANT:** This pipeline uses **single-parent k-mer classification**. The k-mer database directory must contain **KMC format** k-mer databases for:

**Required KMC database files (for k=31, accession=Col-0):**
```
unique_Col-0_ARMS_Chr1_k31.kmc_pre
unique_Col-0_ARMS_Chr1_k31.kmc_suf
unique_Col-0_ARMS_Chr2_k31.kmc_pre
unique_Col-0_ARMS_Chr2_k31.kmc_suf
... (through Chr5)
unique_Col-0_CEN_Chr1_k31.kmc_pre
unique_Col-0_CEN_Chr1_k31.kmc_suf
unique_Col-0_CEN_Chr2_k31.kmc_pre
unique_Col-0_CEN_Chr2_k31.kmc_suf
... (through Chr5)
```

Replace "Col-0" with your accession name and adjust k-mer size if using different parameters.

**K-mer Classification Logic:**
1. Reads are classified based on k-mer counts matching your reference
2. A read is assigned to a haplotype (e.g., Col_ARMS_Chr1) if:
   - It has >50 k-mers matching one haplotype
   - AND other haplotype k-mer counts sum to ≤10
3. Reads not matching these criteria are excluded (no clear classification)

**For single-parent analysis (e.g., Col-0 only):**
- Generate KMC k-mer databases for ARMS and CEN regions of your reference genome
- Reads will be classified as "matching Col ARMS/CEN regions" or "not matching"
- No second parent k-mers needed (unlike full CHARLA_NF two-parent pipeline)

## Troubleshooting

### Container Issues
If using Singularity, ensure it's available:
```bash
singularity --version
```

For Docker:
```bash
nextflow run main.nf --input data.fq.gz --sample_id test -profile docker
```

### Memory Issues
For large datasets, increase memory allocation:
```bash
nextflow run main.nf \
    --input data.fq.gz \
    --sample_id test \
    --readmer_memory '150 GB' \
    --cenhapmer_memory '80 GB'
```

### Check Pipeline DAG
Generate the workflow diagram:
```bash
nextflow run main.nf --help
# DAG will be in results/pipeline_info/pipeline_dag.svg
```

## Example: Arabidopsis Data

For Arabidopsis samples with Col-0 reference:

```bash
# Run analysis
nextflow run main.nf \
    --input data/sample_long_reads.fq.gz \
    --sample_id sample1 \
    --reference Col-0_chr1-5_renamed.fa \
    --reference_dir index/ \
    --kmer_db merylDB/ \
    -profile singularity
```

## Generating K-mer Database for Single-Parent Analysis

To generate a cenhapmer k-mer database for single-parent analysis, you need **KMC format** k-mer databases for:
1. **ARMS regions** (chromosome arms) for each chromosome
2. **CEN regions** (centromeres) for each chromosome

### Example: Creating Col-0 k-mer database with KMC

```bash
# Prerequisites: KMC installed, BED files defining ARMS and CEN regions
REFERENCE="Col-0_chr1-5.fa"
ACCESSION="Col-0"
KMER_SIZE=31
OUTPUT_DIR="cenhapmer_db"

mkdir -p ${OUTPUT_DIR}

# For each chromosome (1-5)
for CHR in {1..5}; do
    # Extract ARMS region sequence
    bedtools getfasta -fi ${REFERENCE} -bed ${ACCESSION}_ARMS_Chr${CHR}.bed \
        -fo ${ACCESSION}_ARMS_Chr${CHR}.fa

    # Extract CEN region sequence
    bedtools getfasta -fi ${REFERENCE} -bed ${ACCESSION}_CEN_Chr${CHR}.bed \
        -fo ${ACCESSION}_CEN_Chr${CHR}.fa

    # Generate k-mer databases with KMC
    # ARMS regions
    kmc -k${KMER_SIZE} -m16 -ci1 -cs1000000 -t4 \
        ${ACCESSION}_ARMS_Chr${CHR}.fa \
        ${OUTPUT_DIR}/unique_${ACCESSION}_ARMS_Chr${CHR}_k${KMER_SIZE} \
        tmp_kmc/

    # CEN regions
    kmc -k${KMER_SIZE} -m16 -ci1 -cs1000000 -t4 \
        ${ACCESSION}_CEN_Chr${CHR}.fa \
        ${OUTPUT_DIR}/unique_${ACCESSION}_CEN_Chr${CHR}_k${KMER_SIZE} \
        tmp_kmc/
done

# Clean up temporary files
rm -rf tmp_kmc/
rm ${ACCESSION}_*.fa

echo "K-mer database created in ${OUTPUT_DIR}/"
```

**Expected output files:**
```
cenhapmer_db/
├── unique_Col-0_ARMS_Chr1_k31.kmc_pre
├── unique_Col-0_ARMS_Chr1_k31.kmc_suf
├── unique_Col-0_ARMS_Chr2_k31.kmc_pre
├── unique_Col-0_ARMS_Chr2_k31.kmc_suf
... (through Chr5)
├── unique_Col-0_CEN_Chr1_k31.kmc_pre
├── unique_Col-0_CEN_Chr1_k31.kmc_suf
... (through Chr5)
```

**Note:**
- You may already have these k-mer databases if you've run CHARLA_NF before
- The pipeline expects files named exactly as: `unique_{accession}_{region}_Chr{N}_k{size}.kmc_*`
- For single-parent mode, you only need one accession (e.g., Col-0), not two

## Credits

This pipeline is based on [CHARLA_NF](https://github.com/jacgonisa/charla_nf) by jacgonisa.

Original publication: [Add citation when available]

## License

Same as CHARLA_NF
