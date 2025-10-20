# Cenhapmer Generation Pipeline

This directory contains a complete pipeline for generating **cenhapmers** (centromeric haplotype k-mers) from parent genome assemblies. Cenhapmers are k-mers unique to specific chromosome regions (centromeres or arms) of specific parental accessions, and are essential for the CHARLA pipeline to identify parental haplotypes in F1 hybrid reads.

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Quick Start](#quick-start)
- [Input Files](#input-files)
- [Pipeline Steps](#pipeline-steps)
- [Usage Examples](#usage-examples)
- [Output](#output)
- [Troubleshooting](#troubleshooting)

## Overview

The cenhapmer generation pipeline identifies k-mers that are unique to specific chromosome regions of each parent genome. This is accomplished through a four-step process, with an optional fifth step for quality assessment:

1. **Split genomes**: Separate each chromosome into centromere (CEN) and chromosome arm (ARMS) regions
2. **Generate k-mers**: Use KMC to count all k-mers in each region
3. **Create operations**: Generate KMC operation files to find unique k-mers
4. **Execute operations**: Run KMC complex operations to generate final cenhapmer databases
5. **Analyze markers** (optional): Assess marker density to determine optimal k-mer size

### What are cenhapmers?

Cenhapmers are k-mers that appear in one specific chromosome region of one parent but not in any other chromosome region of either parent. For example, a `Col-0_CEN_Chr1` cenhapmer appears in the centromere of Chr1 in Col-0, but not in:
- Chr1 arms of Col-0
- Any region of Chr2-5 in Col-0
- Any region of any chromosome in Ler-0 (or other parent)

This uniqueness property allows CHARLA to unambiguously assign F1 hybrid reads to their parental origin.

## Requirements

### Software Dependencies

- **KMC** (v3.2.1 or later): K-mer Counter for efficient k-mer counting
  - `kmc`: Main k-mer counting tool
  - `kmc_tools`: Tool for set operations on k-mer databases
  - Installation: https://github.com/refresh-bio/KMC

- **bedtools** (v2.30.0 or later): For extracting FASTA sequences from BED regions
  - Installation: https://bedtools.readthedocs.io/

- **Python 3** (v3.6 or later): For generating operation files

- **Bash** (v4.0 or later): For running shell scripts

### Hardware Requirements

- **Disk space**: ~50-100 GB per parent genome (varies with k-mer size and genome complexity)
- **Memory**: At least 16 GB RAM (32+ GB recommended for large genomes)
- **CPU**: Multi-core processor recommended (pipeline can use 4-20 threads)

### Input Data Requirements

For each parent genome, you need:

1. **Genome assembly** (FASTA format)
2. **Centromere annotations** (BED format)
3. **Chromosome arm annotations** (BED format)

## Quick Start

```bash
# Navigate to your working directory
cd /path/to/your/project

# Copy the cenhapmer generation pipeline
git clone https://github.com/jacgonisa/charla_nf.git
cd charla_nf/cenhapmer_generation

# Make scripts executable
chmod +x generate_cenhapmers.sh
chmod +x scripts/*.sh
chmod +x scripts/*.py

# Prepare your input files (see Input Files section)
# Then run the complete pipeline
./generate_cenhapmers.sh \\
    --parent1 Col-0 \\
    --parent2 Ler-0 \\
    --num_chr 5 \\
    --kmer_sizes 21,31,41 \\
    --threads 8 \\
    --parallel 4
```

## Input Files

### Directory Structure

Create the following directory structure:

```
working_directory/
├── genomes/
│   ├── Col-0.ragtag_scaffolds.fa
│   └── Ler-0.ragtag_scaffolds.fa
├── annotations/
│   ├── Col-0.ragtag_scaffolds_CEN.bed
│   ├── Col-0.ragtag_scaffolds_ARMS.bed
│   ├── Ler-0.ragtag_scaffolds_CEN.bed
│   └── Ler-0.ragtag_scaffolds_ARMS.bed
└── cenhapmer_generation/  # This pipeline
```

### Genome FASTA Files

Genome assemblies should be in FASTA format. Chromosome names in the FASTA should match those in the BED files.

**Example:**
```
>Chr1
ATCGATCGATCG...
>Chr2
GCTAGCTAGCTA...
...
```

### BED Annotation Files

Two BED files are required per parent:

1. **CEN.bed**: Centromere regions
2. **ARMS.bed**: Chromosome arm regions (non-centromeric)

**BED Format** (tab-separated):
```
Chr1    14840000    15540000    centromere_1
Chr2    3670000     4430000     centromere_2
...
```

**Important**:
- Chromosome names (column 1) must match FASTA headers
- CEN and ARMS regions should be non-overlapping and cover the entire chromosome
- Coordinates are 0-based (BED standard)

#### How to create BED files

If you have centromere coordinates, you can generate the BED files:

**Example for Arabidopsis Col-0:**
```bash
# Create CEN.bed
cat > Col-0.ragtag_scaffolds_CEN.bed << EOF
Chr1    14840000    15540000    CEN1
Chr2    3670000     4430000     CEN2
Chr3    13590000    14290000    CEN3
Chr4    4450000     5150000     CEN4
Chr5    11790000    12490000    CEN5
EOF

# Create ARMS.bed (complement of CEN regions)
# This requires knowing chromosome lengths
bedtools complement -i Col-0.ragtag_scaffolds_CEN.bed \\
    -g <(samtools faidx Col-0.ragtag_scaffolds.fa -i chromsizes) \\
    > Col-0.ragtag_scaffolds_ARMS.bed
```

## Pipeline Steps

### Step 1: Split Genomes

**Script**: `scripts/01-split_genomes.sh`

Extracts CEN and ARMS regions for each chromosome using bedtools.

**Output**:
```
genomes_split/
├── Col-0/
│   ├── CEN/
│   │   ├── Col-0_CEN_Chr1.fa
│   │   ├── Col-0_CEN_Chr2.fa
│   │   └── ...
│   └── ARMS/
│       ├── Col-0_ARMS_Chr1.fa
│       └── ...
└── Ler-0/
    └── ...
```

### Step 2: Generate K-mers

**Script**: `scripts/02-generate_kmers.sh`

Runs KMC on each FASTA file to count k-mers.

**Output**:
```
01-all_kmers_per_Chr_CEN_and_accesion/
├── Col-0/
│   ├── Chr1/
│   │   ├── CEN/
│   │   │   ├── Col-0_CEN_Chr1_k21.kmc_pre
│   │   │   ├── Col-0_CEN_Chr1_k21.kmc_suf
│   │   │   └── ...
│   │   └── ARMS/
│   │       └── ...
│   └── ...
└── Ler-0/
    └── ...
```

### Step 3: Generate KMC Operations

**Script**: `scripts/03-generate_kmc_operations.py`

Creates operation files that define how to find unique k-mers.

**Example operation file** (`unique_Col-0_CEN_Chr1_k21.op`):
```
INPUT:
set1 = 01-all_kmers_per_Chr_CEN_and_accesion/Col-0/Chr1/CEN/Col-0_CEN_Chr1_k21
set2 = 01-all_kmers_per_Chr_CEN_and_accesion/Col-0/Chr1/ARMS/Col-0_ARMS_Chr1_k21
set3 = 01-all_kmers_per_Chr_CEN_and_accesion/Col-0/Chr2/CEN/Col-0_CEN_Chr2_k21
...

OUTPUT:
unique_Col-0_CEN_Chr1_k21 = set1 - (set2 + set3 + ... + set20)
```

**Output**:
```
02-kmc_complex_op_files/
├── k21/
│   ├── unique_Col-0_CEN_Chr1_k21.op
│   ├── unique_Col-0_ARMS_Chr1_k21.op
│   └── ...
├── k31/
│   └── ...
└── k41/
    └── ...
```

### Step 4: Execute KMC Operations

**Script**: `scripts/04-run_kmc_operations.py`

Runs `kmc_tools complex` on each operation file to generate unique k-mer databases.

**Output** (final cenhapmers):
```
03-cenhapmers/
├── k21/
│   ├── unique_Col-0_CEN_Chr1_k21.kmc_pre
│   ├── unique_Col-0_CEN_Chr1_k21.kmc_suf
│   ├── unique_Col-0_ARMS_Chr1_k21.kmc_pre
│   ├── unique_Col-0_ARMS_Chr1_k21.kmc_suf
│   ├── unique_Ler-0_CEN_Chr1_k21.kmc_pre
│   ├── unique_Ler-0_CEN_Chr1_k21.kmc_suf
│   └── ...
├── k31/
│   └── ...
└── k41/
    └── ...
```

## Usage Examples

### Example 1: Arabidopsis Col-0 × Ler-0 (default)

```bash
./generate_cenhapmers.sh \\
    --parent1 Col-0 \\
    --parent2 Ler-0 \\
    --num_chr 5 \\
    --kmer_sizes 21,31,41 \\
    --threads 8
```

### Example 2: Different Species with 7 Chromosomes

```bash
./generate_cenhapmers.sh \\
    --parent1 StrainA \\
    --parent2 StrainB \\
    --num_chr 7 \\
    --chr_prefix Chr \\
    --kmer_sizes 21,31 \\
    --fasta_dir my_genomes \\
    --bed_dir my_annotations \\
    --threads 16 \\
    --parallel 4
```

### Example 3: Using Custom Directory Names

```bash
./generate_cenhapmers.sh \\
    --parent1 Cvi \\
    --parent2 Est \\
    --num_chr 5 \\
    --kmer_sizes 21 \\
    --fasta_dir reference_genomes \\
    --bed_dir genome_annotations
```

### Example 4: With Marker Density Analysis (Recommended)

```bash
./generate_cenhapmers.sh \\
    --parent1 Col-0 \\
    --parent2 Ler-0 \\
    --num_chr 5 \\
    --kmer_sizes 21,31,41 \\
    --threads 8 \\
    --parallel 4 \\
    --analyze-markers
```

This will generate additional plots and statistics to help you choose the best k-mer size.

### Example 5: Parallel Processing with Cleanup

```bash
./generate_cenhapmers.sh \\
    --parent1 Col-0 \\
    --parent2 Ler-0 \\
    --num_chr 5 \\
    --kmer_sizes 21,31,41 \\
    --threads 20 \\
    --parallel 8 \\
    --clean
```

### Example 6: Resume from Specific Step

If a step failed, you can skip completed steps:

```bash
# Skip steps 1 and 2, resume from step 3
./generate_cenhapmers.sh \\
    --parent1 Col-0 \\
    --parent2 Ler-0 \\
    --skip-step 1 \\
    --skip-step 2
```

### Example 7: Running Individual Steps

You can also run steps individually:

```bash
# Step 1 only
bash scripts/01-split_genomes.sh \\
    --parent1 Col-0 \\
    --parent2 Ler-0 \\
    --num_chr 5

# Step 2 only
bash scripts/02-generate_kmers.sh \\
    --parent1 Col-0 \\
    --parent2 Ler-0 \\
    --num_chr 5 \\
    --kmer_sizes 21 \\
    --threads 8

# Step 3 only
python3 scripts/03-generate_kmc_operations.py \\
    --parent1 Col-0 \\
    --parent2 Ler-0 \\
    --kmer_sizes 21

# Step 4 only
python3 scripts/04-run_kmc_operations.py \\
    --kmer_sizes 21 \\
    --parallel 4
```

## Output

### Cenhapmer Databases

The final output is in `03-cenhapmers/`:

```
03-cenhapmers/
├── k21/
│   ├── unique_Col-0_CEN_Chr1_k21.kmc_pre
│   ├── unique_Col-0_CEN_Chr1_k21.kmc_suf
│   ├── unique_Col-0_ARMS_Chr1_k21.kmc_pre
│   ├── unique_Col-0_ARMS_Chr1_k21.kmc_suf
│   └── ... (40 files total: 2 parents × 5 chr × 2 regions × 2 file types)
├── k31/
│   └── ... (same structure)
└── k41/
    └── ... (same structure)
```

Each cenhapmer database consists of two files:
- `.kmc_pre`: Index file
- `.kmc_suf`: Data file

### Using Cenhapmers with CHARLA

Once generated, use these cenhapmers with the CHARLA pipeline:

```bash
nextflow run jacgonisa/charla_nf \\
    --reads your_F1_pollen.fq \\
    --sample_id my_sample \\
    --cenhapmer_db_dir /path/to/03-cenhapmers/k21 \\
    --parent1_name Col-0 \\
    --parent2_name Ler-0 \\
    --parent1_bed annotations/Col-0_renamed.bed \\
    --parent2_bed annotations/Ler-0_renamed.bed \\
    --num_chromosomes 5 \\
    --kmer_size 21 \\
    --outdir results \\
    -profile singularity
```

### Verifying Output

Check that all expected files were created:

```bash
# Count cenhapmer files for k=21
ls 03-cenhapmers/k21/*.kmc_pre | wc -l
# Expected: 20 files (2 parents × 5 chr × 2 regions)

# Check file sizes (should not be empty)
ls -lh 03-cenhapmers/k21/

# Inspect k-mer counts using kmc_tools
kmc_tools info 03-cenhapmers/k21/unique_Col-0_CEN_Chr1_k21
```

## Troubleshooting

### Common Issues

#### 1. "kmc: command not found"

**Solution**: Install KMC and add it to your PATH:
```bash
# Download KMC
wget https://github.com/refresh-bio/KMC/releases/download/v3.2.1/KMC3.2.1.linux.tar.gz
tar -xzf KMC3.2.1.linux.tar.gz
sudo mv kmc* /usr/local/bin/
```

#### 2. "bedtools: command not found"

**Solution**: Install bedtools:
```bash
# Ubuntu/Debian
sudo apt-get install bedtools

# macOS
brew install bedtools

# Conda
conda install -c bioconda bedtools
```

#### 3. Empty cenhapmer files

**Cause**: No unique k-mers found for that region.

**Solution**: This can happen with:
- Very similar parent genomes
- Small chromosome regions
- High k-mer size relative to region size

Try using a smaller k-mer size or check that your BED files are correct.

#### 4. "Out of memory" errors during KMC

**Solution**: Reduce the number of threads or use disk-based mode:
```bash
# Reduce threads
./generate_cenhapmers.sh --threads 2

# Or modify Step 2 script to add -m parameter for max memory
kmc -fa -k21 -ci1 -cs100000000 -m8 -t4 ...
```

#### 5. Pipeline runs slowly

**Optimization tips**:
- Increase `--threads` for KMC (Step 2)
- Increase `--parallel` for kmc_tools (Step 4)
- Use SSD storage for intermediate files
- Process different k-mer sizes in separate runs

### Getting Help

If you encounter issues:

1. Check that all dependencies are installed correctly
2. Verify input file formats and paths
3. Review error messages carefully
4. Check disk space and memory availability
5. Try running individual steps to isolate the problem

For further assistance:
- GitHub Issues: https://github.com/jacgonisa/charla_nf/issues
- CHARLA Documentation: https://github.com/jacgonisa/charla_nf/docs

## Marker Density Analysis (Step 5)

The optional marker density analysis step helps you choose the optimal k-mer size for your CHARLA analysis by visualizing how many unique markers (cenhapmers) are distributed across chromosomes.

### Why Analyze Marker Density?

Different k-mer sizes involve important tradeoffs:
- **Smaller k (e.g., k=21)**: Fewer unique markers, less specific, but **more resilient to sequencing errors**
- **Larger k (e.g., k=41)**: More unique markers, higher specificity, but **less tolerant of read errors**

**Key insight**: A single sequencing error in a read will break a k=41 marker but may still preserve several k=21 markers from the same region. For noisy long reads, smaller k-mer sizes are often more robust.

The marker density analysis provides:
1. **Visual comparison** of marker distribution across chromosomes for different k-mer sizes
2. **Quantitative statistics** on total markers per region (CEN vs ARMS)
3. **Recommendations** for optimal k-mer size based on your data

### Running Marker Density Analysis

Add the `--analyze-markers` flag to your pipeline run:

```bash
./generate_cenhapmers.sh \\
    --parent1 Col-0 \\
    --parent2 Ler-0 \\
    --kmer_sizes 21,31,41 \\
    --analyze-markers
```

### Output Files

The analysis generates:

```
05-marker_plots/
├── marker_counts_per_chromosome.png      # Bar charts comparing k-mer sizes
├── marker_counts_cen_vs_arms.png         # CEN vs ARMS distribution
├── kmer_size_comparison.png              # Overall comparison
├── marker_density_summary.txt            # Text summary with recommendations
└── marker_density_summary.csv            # Detailed statistics (CSV format)
```

### Interpreting Results

1. **Total Markers**: Higher is generally better (more coverage)
2. **CEN vs ARMS Balance**: Both regions should have adequate markers
3. **Chromosome Distribution**: Markers should be distributed across all chromosomes

**Example Interpretation**:
```
K-mer size: 21
  Total markers:   9,123,456
  Mean per region: 456,173
  Target density:  42 markers/kb (25th percentile)

K-mer size: 31
  Total markers:   12,456,789
  Mean per region: 622,839
  Target density:  120 markers/kb (25th percentile)

K-mer size: 41
  Total markers:   15,234,567
  Mean per region: 763,728
  Target density:  180 markers/kb (25th percentile)

Note: Higher k-mer sizes generate MORE unique markers but are less resilient to errors
Recommendation: k=21 or k=31 for Arabidopsis
  - k=21: Best for noisy long reads (ONT, PacBio CLR)
  - k=31: Good for high-quality reads (PacBio HiFi, Illumina)
```

### When to Use Which K-mer Size

- **k=21**:
  - ✅ Best for noisy long reads (ONT, older PacBio)
  - ✅ More resilient to sequencing errors
  - ✅ Faster processing
  - ⚠️ Fewer unique markers, less specific

- **k=31**:
  - ✅ Good balance of error tolerance and specificity
  - ✅ Recommended for high-quality long reads (PacBio HiFi)
  - ✅ More markers than k=21

- **k=41**:
  - ✅ Maximum unique markers and specificity
  - ✅ Best discrimination for complex genomes
  - ⚠️ Requires high-quality reads (any error breaks the k-mer)
  - ⚠️ Not recommended for noisy data

## Advanced Topics

### Custom K-mer Sizes

The default k-mer sizes (21, 31, 41) work well for Arabidopsis, but you may need different sizes for:
- **Larger genomes**: Use larger k (31-51) for more specificity
- **Smaller genomes**: Smaller k (17-25) may provide better coverage
- **High heterozygosity**: Larger k may help distinguish haplotypes

### Filtering Parameters

You can modify the KMC filtering parameters in `scripts/02-generate_kmers.sh`:
- `-ci`: Minimum k-mer count (default: 1)
- `-cs`: Maximum k-mer count (default: 100000000)
- `-m`: Maximum RAM usage in GB

### Alternative Naming Schemes

If your files use different naming conventions, modify the suffix variables in scripts:
```bash
# In 01-split_genomes.sh
FASTA_SUFFIX=".fa"  # Instead of .ragtag_scaffolds.fa
BED_CEN_SUFFIX="_CEN.bed"  # Instead of .ragtag_scaffolds_CEN.bed
```

## Citation

If you use this pipeline, please cite:

```
González-Isa, J. et al. (2025). CHARLA: Centromeric Haplotype Analysis
of Recombination using Long-reads of Arabidopsis. [In preparation]
```

## License

This pipeline is part of the CHARLA project and is distributed under the MIT License.

## Contact

For questions or suggestions:
- Jacob González Isa <jg2070@cam.ac.uk>
- GitHub: https://github.com/jacgonisa/charla_nf
