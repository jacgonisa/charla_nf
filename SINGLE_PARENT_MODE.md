# Single-Parent K-mer Classification Mode

## Overview

This pipeline implements **single-parent k-mer classification** - it classifies long reads based on their k-mer profiles against a **single reference genome** (e.g., Col-0), without detecting recombination between two parents.

## What This Means

### ✅ What the Pipeline DOES:
1. **K-mer-based read classification** - Classifies reads based on k-mer counts matching ARMS/CEN regions
2. **Structural variant analysis** - Maps reads and identifies indels
3. **Large indel detection** - Extracts and analyzes indels >50bp
4. **Satellite DNA analysis** - Studies satellite structure organization
5. **Multiple sequence alignment** - Aligns satellite blocks with MAFFT

### ❌ What the Pipeline DOES NOT DO:
1. **Crossover detection** - No recombination analysis between two parents
2. **Two-parent classification** - Only classifies against ONE reference
3. **Hybrid read profiling** - No detection of reads with mixed parental k-mers

## Key Difference from Full CHARLA_NF

| Feature | Full CHARLA_NF | This Pipeline (Single-Parent) |
|---------|----------------|-------------------------------|
| Parents | 2 (e.g., Col + Ler) | 1 (e.g., Col only) |
| K-mer DB | Both parents' k-mers | One parent's k-mers only |
| Read Classification | "Col", "Ler", "Hybrid" | "Col_ARMS_Chr1", "Col_CEN_Chr2", etc. |
| Crossover Detection | YES | NO |
| Recombination Analysis | YES | NO |
| Structural Variants | YES | YES |
| Indel Analysis | YES | YES |
| Satellite Analysis | YES | YES |

## How Read Classification Works

### Single-Parent K-mer Classification Logic

For each read:
1. Count k-mers matching each chromosome region (ARMS/CEN) of the reference
2. Assign read to region if:
   - It has **>50 k-mers** matching that region
   - AND other regions have **≤10 k-mers total**
3. Reads not meeting criteria are excluded (no clear classification)

### Example Classification

Given a Col-0 reference with ARMS and CEN regions:

```
Read_001:
  - Col_ARMS_Chr1: 85 k-mers → ✅ ASSIGNED
  - Col_CEN_Chr1: 3 k-mers
  - All others: 2 k-mers total

Read_002:
  - Col_ARMS_Chr1: 45 k-mers → ❌ NOT ASSIGNED (< 50)
  - Col_CEN_Chr2: 30 k-mers

Read_003:
  - Col_ARMS_Chr1: 60 k-mers → ❌ NOT ASSIGNED
  - Col_ARMS_Chr2: 55 k-mers → (total others > 10)
```

## Required K-mer Database

### Directory Structure

Your `cenhapmer_db` directory must contain **KMC format** k-mer databases:

```
cenhapmer_db/
├── unique_Col-0_ARMS_Chr1_k31.kmc_pre
├── unique_Col-0_ARMS_Chr1_k31.kmc_suf
├── unique_Col-0_ARMS_Chr2_k31.kmc_pre
├── unique_Col-0_ARMS_Chr2_k31.kmc_suf
├── unique_Col-0_ARMS_Chr3_k31.kmc_pre
├── unique_Col-0_ARMS_Chr3_k31.kmc_suf
├── unique_Col-0_ARMS_Chr4_k31.kmc_pre
├── unique_Col-0_ARMS_Chr4_k31.kmc_suf
├── unique_Col-0_ARMS_Chr5_k31.kmc_pre
├── unique_Col-0_ARMS_Chr5_k31.kmc_suf
├── unique_Col-0_CEN_Chr1_k31.kmc_pre
├── unique_Col-0_CEN_Chr1_k31.kmc_suf
├── unique_Col-0_CEN_Chr2_k31.kmc_pre
├── unique_Col-0_CEN_Chr2_k31.kmc_suf
├── unique_Col-0_CEN_Chr3_k31.kmc_pre
├── unique_Col-0_CEN_Chr3_k31.kmc_suf
├── unique_Col-0_CEN_Chr4_k31.kmc_pre
├── unique_Col-0_CEN_Chr4_k31.kmc_suf
├── unique_Col-0_CEN_Chr5_k31.kmc_pre
└── unique_Col-0_CEN_Chr5_k31.kmc_suf
```

**Naming Convention:**
```
unique_{ACCESSION}_{REGION}_Chr{N}_k{SIZE}.kmc_pre
unique_{ACCESSION}_{REGION}_Chr{N}_k{SIZE}.kmc_suf
```

Where:
- `ACCESSION` = Your reference name (e.g., "Col-0")
- `REGION` = "ARMS" or "CEN"
- `N` = Chromosome number (1-5 for Arabidopsis)
- `SIZE` = K-mer size (default: 31)

### Important Note

**You do NOT need k-mers from a second parent!**

In full CHARLA_NF, you would need both Col AND Ler k-mers. For single-parent mode, you only need Col (or your reference of choice).

## Use Cases

This single-parent mode is perfect for:

1. **Mutant vs Wild-Type Analysis**
   - Compare ATXR5/6 mutant reads to Col-0 reference
   - Identify structural changes without hybrid/crossover context

2. **Structural Variant Detection**
   - Focus on indels and structural changes
   - No need for two-parent k-mer databases

3. **Satellite DNA Studies**
   - Analyze satellite organization in single genotype
   - Study 178bp repeat patterns

4. **Long-Read Classification**
   - Classify reads by chromosome region using k-mers
   - Map and analyze without recombination overhead

## Configuration Parameters

### Key Parameters for Single-Parent Mode

```bash
nextflow run main.nf \
    --input sample.fq.gz \
    --sample_id sample1 \
    --reference Col-0.fa \
    --reference_dir ./index \
    --kmer_db ./cenhapmer_db \    # Only Col-0 k-mers needed!
    --kmer_size 31 \
    --parent1_name 'Col' \         # Single parent
    --parent2_name null \          # No second parent!
    -profile singularity
```

### Pipeline Config (`nextflow.config`)

```groovy
params {
    // Single parent mode
    parent1_name = 'Col'     // Your reference
    parent2_name = null      // Not used in single-parent mode
    parent1_bed  = null      // Optional
    parent2_bed  = null      // Not needed

    kmer_db = './cenhapmer_db'  // Only ONE parent's k-mers
}
```

## Output Interpretation

### Read Lists

Classified reads will be saved as:
- `CA1_reads.txt` - Reads matching Col ARMS Chr1
- `CA2_reads.txt` - Reads matching Col ARMS Chr2
- `CA3_reads.txt` - Reads matching Col ARMS Chr3
- `CA4_reads.txt` - Reads matching Col ARMS Chr4
- `CA5_reads.txt` - Reads matching Col ARMS Chr5
- `CC1_reads.txt` - Reads matching Col CEN Chr1
- `CC2_reads.txt` - Reads matching Col CEN Chr2
- `CC3_reads.txt` - Reads matching Col CEN Chr3
- `CC4_reads.txt` - Reads matching Col CEN Chr4
- `CC5_reads.txt` - Reads matching Col CEN Chr5

**Abbreviation Key:**
- `CA` = Col ARMS
- `CC` = Col CEN
- Number = Chromosome

### Mapped Files

For each classified read group:
- `{group}_sequences.fa` - Extracted sequences
- `bam/{group}.bam` - Mapped reads
- Indel analysis results

## Questions?

**Q: Why do I need k-mer classification if I'm only using one reference?**

A: K-mer classification allows you to:
1. Identify reads that clearly match your reference (high-quality matches)
2. Separate reads by chromosome/region (ARMS vs CEN)
3. Filter out ambiguous or low-quality reads
4. Analyze structural variants in specific genomic contexts

**Q: Can I skip the k-mer classification?**

A: The k-mer classification is integral to this pipeline's workflow. If you want simpler mapping without k-mer classification, consider using minimap2 directly.

**Q: What if I want to add a second parent later?**

A: Generate k-mers for the second parent using the same naming convention (e.g., `unique_Ler-0_ARMS_Chr1_k31.*`) and place them in the same directory. The pipeline will automatically detect both parents' k-mers.

## Summary

This is **CHARLA with single-parent k-mer classification** - you get all the powerful k-mer-based read classification and structural variant analysis, without the complexity of two-parent recombination detection. Perfect for analyzing mutants, variants, or any long-read data against a single reference genome!
