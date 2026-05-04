# CHARLA: Centromeric Haplotype Analysis of Recombination using Long-reads of Arabidopsis

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A524.10.5-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.3.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.3.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/jacgonisa/charla_nf)

![][logo]

[logo]: https://github.com/jacgonisa/charla_nf/blob/master/images/CHARLA_image.png

## Introduction

**CHARLA** (**C**entromeric **H**aplotype **A**nalysis of **R**ecombination using **L**ong-reads of **A**rabidopsis) is a bioinformatics pipeline designed to detect and analyze crossover events in *Arabidopsis thaliana* hybrid pollen using long-read sequencing data (Oxford Nanopore or PacBio). The pipeline performs k-mer-based haplotype profiling to identify recombination events, including crossovers in centromeric regions, and generates comprehensive visualizations and analysis reports.

### Pipeline Overview

CHARLA uses a k-mer profiling approach to distinguish parental haplotypes in F1 hybrid reads. The pipeline:

1. **Processes sequencing data** - Converts FASTQ to FASTA, simplifies headers, and optionally performs quality-based masking
2. **Profiles read k-mers (readmers)** - Generates k-mer databases from input reads using KMC3
3. **Profiles haplotype-specific k-mers (cenhapmers)** - Counts parental-specific k-mers in each read using pre-computed databases
4. **Identifies hybrid reads** - Combines profiles to detect reads with both parental haplotypes (potential crossover reads)
5. **Segments reads** - Identifies crossover breakpoints and segments reads into arm and centromeric regions
6. **Maps and validates** - Aligns segmented reads to reference genomes using minimap2
7. **Visualizes crossovers** - Generates genome-wide crossover maps and coverage plots
8. **Analyzes non-hybrid reads** - Performs indel analysis on parental-only reads

### Workflow Diagram

```
FASTQ/FASTA Input
      │
      ├─→ [Optional] FASTQ_TO_FASTA (with quality masking)
      │
      ├─→ SIMPLIFY_HEADERS
      │
      ├────────────────────────────────┬──────────────────────────────────┐
      │                                │                                  │
      ▼                                ▼                                  ▼
Branch A: READMER                Branch B: CENHAPMER            Branch C: NON-HYBRID
Profiling                        Profiling                       Analysis
      │                                │                                  │
      ├─→ generate_readmers_kmc        ├─→ get_counts_cenhapmer_kmc     │
      ├─→ generate_histogram_kmc       ├─→ process_cenhapmer_counts     │
      └─→ get_counts_kmc               │                                 │
                                       │                                 │
                    ┌──────────────────┴──────────────────┐              │
                    ▼                                      ▼              │
            COMBINE_HYBRID_PROFILES              NON_HYBRID_READS_ANALYSIS
                    │                                                     │
                    ▼                                                     ▼
            CURATE_HYBRID_PROFILES                        Indel Analysis & Mapping
                    │
                    ▼
              SEGMENT_READS
                    │
                    ├─→ MAPPING_ANALYSIS (minimap2)
                    │
                    └─→ PLOT_CROSSOVER_MAP
                              │
                              ▼
                    Crossover Maps & Coverage Plots
```

## Quick Start

### Prerequisites

- **Nextflow** >= 24.10.5 ([Installation guide](https://www.nextflow.io/docs/latest/getstarted.html))
- **Container system**: Singularity (default), Docker, Conda, or Apptainer
- **Pre-computed cenhapmer databases**: K-mer databases for parental haplotypes (Col-0 and Ler-0)
- **Reference genomes**: Indexed reference genomes for both parental accessions
- **BED files**: Genomic feature annotations for Col-0 and Ler-0

### Test the Pipeline

```bash
nextflow run jacgonisa/charla_nf -profile test,singularity
```

## Usage

### Basic Command

```bash
nextflow run jacgonisa/charla_nf \
    --reads <input.fastq|input.fasta> \
    --sample_id <sample_name> \
    --input_format <fastq|fasta> \
    --cenhapmer_db_dir <path/to/kmer_databases> \
    --reference_genomes_dir <path/to/references> \
    --col_bed_file <path/to/Col-0.bed> \
    --ler_bed_file <path/to/Ler-0.bed> \
    --outdir <output_directory> \
    --kmer_size <21|31> \
    -profile singularity
```

### Required Parameters

| Parameter | Description |
|-----------|-------------|
| `--reads` | Path to input FASTQ or FASTA file containing sequencing reads |
| `--sample_id` | Unique identifier for the sample |
| `--input_format` | Input file format: `fastq` or `fasta` |
| `--cenhapmer_db_dir` | Directory containing pre-computed KMC databases for parental haplotypes |
| `--reference_genomes_dir` | Directory containing reference genome FASTA files for both parents |
| `--col_bed_file` | BED file with genomic annotations for Col-0 reference |
| `--ler_bed_file` | BED file with genomic annotations for Ler-0 reference |
| `--outdir` | Directory for output files |
| `--kmer_size` | K-mer size for analysis (typically 21 or 31) |

### Optional Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--do_quality_masking` | `false` | Enable quality-based masking of low-quality bases in FASTQ files |
| `--readmer_cutoff` | `10` | Minimum k-mer count threshold for readmer profiling |
| `--min_count` | `25` | Minimum k-mer count for hybrid profile combination |
| `--threshold` | `50` | Threshold for read segmentation and filtering |
| `--readmer_cpus` | `20` | Number of CPUs for readmer profiling |
| `--readmer_memory` | `100 GB` | Memory allocation for readmer profiling |
| `--cenhapmer_cpus` | `1` | Number of CPUs per cenhapmer task (parallel tasks) |
| `--cenhapmer_memory` | `50 GB` | Memory allocation for cenhapmer profiling |

### Example Commands

#### Example 1: FASTQ Input with Quality Masking

```bash
nextflow run jacgonisa/charla_nf \
    --reads data/ColLer_F1_pollen_minq20.fq \
    --sample_id ColLer_F1_pollen_minq20 \
    --input_format fastq \
    --cenhapmer_db_dir /path/to/cenhapmers/k21 \
    --reference_genomes_dir /path/to/references \
    --col_bed_file annotations/Col-0_renamed.bed \
    --ler_bed_file annotations/Ler-0_renamed.bed \
    --outdir results/pollen_analysis \
    --kmer_size 21 \
    --do_quality_masking true \
    --readmer_cutoff 10 \
    --min_count 25 \
    -profile singularity \
    -with-trace trace.txt \
    -with-report report.html \
    -with-timeline timeline.html \
    -resume
```

#### Example 2: FASTA Input (Pre-processed Reads)

```bash
nextflow run jacgonisa/charla_nf \
    --reads data/simulated_reads.fasta \
    --sample_id simulation_k31 \
    --input_format fasta \
    --cenhapmer_db_dir /path/to/cenhapmers/k31 \
    --reference_genomes_dir /path/to/references \
    --col_bed_file annotations/Col-0_renamed.bed \
    --ler_bed_file annotations/Ler-0_renamed.bed \
    --outdir results/simulation \
    --kmer_size 31 \
    -profile singularity \
    -resume
```

### Running with Different Profiles

CHARLA supports multiple execution profiles:

```bash
# Singularity (default, recommended for HPC)
-profile singularity

# Docker (for local machines)
-profile docker

# Conda (no containers)
-profile conda

# Local execution (for testing)
-profile local

# HPC cluster with SLURM
-profile cluster
```

## Output Structure

The pipeline generates organized output in the specified `--outdir`:

```
results/
├── data/                          # Pre-processed FASTA files
├── readmers/                      # Read k-mer profiling results
│   ├── *.kmc_pre                  # KMC database files
│   ├── *.kmc_suf
│   ├── *.histo                    # K-mer frequency histograms
│   └── *.counts.txt               # Per-read k-mer counts
├── processed_cenhapmers/          # Processed cenhapmer matrices
│   └── *.cenhapmer_matrix.tsv     # K-mer count matrix per read
├── hybrid_profiles/               # Hybrid read identification
│   ├── *.hybrid_profiles.txt      # Combined profiles
│   ├── *.final_table.txt          # Curated hybrid reads
│   └── *.ultracurated_table.txt   # High-confidence hybrids
├── read_segments/                 # Segmented reads and BED files
│   ├── *.ARMS.segments.fa         # Arm region segments
│   ├── *.CEN.segments.fa          # Centromeric segments
│   ├── *.ARMS.bed                 # Arm coordinates
│   └── *.CEN.bed                  # Centromere coordinates
├── mapping_analysis/              # Minimap2 alignment results
│   ├── *.paf                      # PAF alignment files
│   ├── *.cowidth.txt              # Crossover width analysis
│   └── plots/                     # Mapping visualizations
├── crossover_plots/               # Final crossover maps
│   ├── *.crossover_map.pdf        # Genome-wide crossover visualization
│   └── coverage_data/             # Coverage statistics
├── nonhybrid_analysis/            # Non-hybrid read analysis
│   ├── read_lists/                # Parental read identifiers
│   ├── sequences/                 # Extracted FASTA sequences
│   ├── bam_files/                 # Aligned BAM files
│   └── indel_analysis/            # Indel detection results
└── pipeline_info/                 # Execution reports
    ├── execution_report.html
    ├── execution_timeline.html
    └── execution_trace.txt
```

## Key Features

### 1. Quality-Based Masking
When `--do_quality_masking true` is specified, low-quality bases in FASTQ files are masked (converted to 'N') before k-mer profiling, improving the accuracy of haplotype assignment.

### 2. Parallel Cenhapmer Processing
The pipeline efficiently processes multiple cenhapmer databases in parallel, significantly reducing runtime for large-scale analyses.

### 3. Hybrid Read Detection
Uses a sophisticated k-mer profiling approach to identify reads containing both parental haplotypes, indicating potential crossover events.

### 4. Crossover Localization
Segments hybrid reads to precisely identify crossover breakpoints at k-mer resolution, distinguishing arm and centromeric crossovers.

### 5. Comprehensive Validation
Validates putative crossovers through minimap2 alignment to parental reference genomes, filtering out false positives.

### 6. Non-Hybrid Analysis
Performs separate analysis of parental-only reads, including indel detection and structural variant analysis.

## Performance Optimization

### Resource Configuration

The pipeline automatically configures resources based on the specified parameters. For optimal performance:

**For FASTQ input with masking:**
- `--readmer_cpus 20` (for k-mer database generation)
- `--readmer_memory 100GB`
- `--cenhapmer_cpus 1` with high parallelization

**For large datasets (>10 million reads):**
- Increase `--cenhapmer_memory` to 100GB
- Use `--readmer_cutoff` to filter low-abundance k-mers
- Consider running on HPC with `-profile cluster`

### Resuming Failed Runs

Nextflow's `-resume` flag allows you to continue from the last successful step:

```bash
nextflow run jacgonisa/charla_nf ... -resume
```

## Troubleshooting

### Common Issues

**Issue:** Out of memory errors during readmer profiling
- **Solution:** Increase `--readmer_memory` or use `--readmer_cutoff` to reduce k-mer space

**Issue:** No hybrid reads detected
- **Solution:** Check k-mer database quality, verify input is from F1 hybrid, adjust `--min_count` threshold

**Issue:** Container/Singularity errors
- **Solution:** Ensure container cache directory has sufficient space: `export NXF_SINGULARITY_CACHEDIR=/path/to/cache`

### Debug Mode

Enable detailed logging:

```bash
nextflow run jacgonisa/charla_nf ... -profile debug
```

## Credits

CHARLA was originally developed by **Jacob González Isa** at the University of Cambridge.

### Contributors

We thank the following people for their extensive assistance in the development of this pipeline:

- **Namil Son** - Algorithm development and testing
- **Matthew Naish** - Biological validation and experimental design

## Citations

If you use CHARLA for your research, please cite:

> González Isa, J., Son, N., Naish, M., et al. (2025). CHARLA: A Nextflow pipeline for centromeric haplotype analysis of recombination using long-read sequencing. *In preparation*.

### Tools Used

This pipeline uses the following tools:

- **Nextflow** - Workflow management (Di Tommaso et al., 2017)
- **KMC3** - K-mer counting (Kokot et al., 2017)
- **minimap2** - Read alignment (Li, 2018)
- **seqkit** - Sequence manipulation (Shen et al., 2016)
- Custom Rust implementations for performance-critical steps

### nf-core

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For support and questions:
- **Issues:** [GitHub Issues](https://github.com/jacgonisa/charla_nf/issues)
- **Email:** jg2070@cam.ac.uk

## License

This project is licensed under the MIT License - see the LICENSE file for details.
a