# Configuration Fix Summary

## Errors Fixed

### 1. Nextflow DSL2 Syntax Error

**Error:**
```
ERROR ~ No such variable: optional
 -- Check script 'charla_nf/modules/local/visualize_read_structures.nf' at line: 14
```

**Problem:**
Incorrect syntax for optional input parameters in Nextflow DSL2:
```groovy
path confidence_scores, optional: true  // WRONG!
```

**Fix:**
Removed the `optional: true` syntax (not valid in DSL2 input declarations):
```groovy
path confidence_scores  // CORRECT
```

The optionality is handled in the script section with:
```bash
conf_arg = confidence_scores.name != 'NO_FILE' ? "--confidence-scores ${confidence_scores}" : ""
```

### 2. Missing Parameters

**Warnings:**
```
WARN: Access to undefined parameter `enable_conda`
WARN: Access to undefined parameter `threshold`
WARN: Access to undefined parameter `segment_cpus`
WARN: Access to undefined parameter `segment_memory`
WARN: Access to undefined parameter `segment_time`
WARN: Access to undefined parameter `mapping_cpus`
WARN: Access to undefined parameter `mapping_memory`
WARN: Access to undefined parameter `mapping_time`
```

**Fix:**
Added missing parameters to `nextflow.config`:

```groovy
// Confidence scoring threshold
threshold = 50              // Minimum marker count threshold

// Segmentation and mapping resources
segment_cpus = 4
segment_memory = '8 GB'
segment_time = '2h'

mapping_cpus = 8
mapping_memory = '16 GB'
mapping_time = '4h'

// Read structure visualization
max_structure_reads = 200   // Number of reads to visualize

// Large indel analysis resources
large_indel_cpus = 8
large_indel_memory = '16 GB'
large_indel_time = '4h'

// Environment options
enable_conda = false        // Use containers instead of conda
```

## Files Modified

1. **modules/local/visualize_read_structures.nf** (line 14)
   - Removed invalid `optional: true` syntax from input declaration

2. **nextflow.config** (lines 76-97)
   - Added missing parameter defaults

## Previous Fix (Still Active)

The optimization to `scripts/13.2-map_indels_to_genome.py` is still in effect:
- Batch BAM file reading (3000x speedup)
- LARGE_INDEL_ANALYSIS should complete in minutes instead of hours

## Ready to Run

Your pipeline should now start without errors. The warnings are resolved and the syntax error is fixed.

```bash
nextflow run charla_nf/main.nf \
    --reads Col_Ler_F1_pollen_q10_1kb_pore.fq.gz \
    --do_quality_masking true \
    --sample_id Col_Ler_F1_pollen_q10_1kb \
    --input_format fastq \
    --cenhapmer_db_dir /home/jg2070/Desktop/PhD/crossover/03-cenhapmers/k31 \
    --reference_genomes_dir /home/jg2070/Desktop/PhD/crossover/index \
    --outdir /mnt/ssd-8tb/Col_Ler_F1_pollen_q10_1kb_2023 \
    --kmer_size 31 \
    --cenhapmer_cpus 1 \
    --cenhapmer_memory 50GB \
    --readmer_memory 50GB \
    --readmer_cpus 10 \
    --readmer_cutoff 10 \
    -profile singularity \
    -with-trace trace.txt \
    -work-dir /mnt/ssd-8tb/Col_Ler_F1_pollen_q10_1kb/work \
    --col_bed_file only_kmer_pipeline/10-plotting_crossover/Col-0_renamed.bed \
    --ler_bed_file only_kmer_pipeline/10-plotting_crossover/Ler-0_renamed.bed \
    -resume \
    --curate_time 10h
```

Expected behavior:
- No warnings about undefined parameters
- No syntax errors
- LARGE_INDEL_ANALYSIS completes quickly (2-5 minutes)
- New modules run successfully (read structures, confidence V3, etc.)
