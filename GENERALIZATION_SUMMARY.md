# CHARLA Pipeline Generalization - Summary

## What Was Done

Successfully generalized the CHARLA pipeline to support **any F1 hybrid cross**, not just Col×Ler! The pipeline can now process pollen and sperm sequencing data from different Arabidopsis accessions or even potentially other species with different chromosome numbers.

## Key Changes

### 1. **New Generic Parameters** (nextflow.config:45-54)

Replace accession-specific parameters with flexible, generic ones:

```groovy
// NEW GENERIC PARAMETERS
parent1_name = 'Col'            // Name of parent 1 accession (e.g., 'Col', 'Cvi', 'Ws')
parent2_name = 'Ler'            // Name of parent 2 accession (e.g., 'Ler', 'Est', 'Sha')
parent1_bed = null              // Path to parent 1 BED file
parent2_bed = null              // Path to parent 2 BED file
num_chromosomes = 5             // Number of chromosomes (default: 5 for Arabidopsis)
```

### 2. **Dynamic K-mer Combination Generation** (main.nf:60-73)

Automatically creates k-mer combinations based on your specified parents:

```groovy
ch_kmc_combinations = Channel
    .from([params.parent1_name, params.parent2_name])
    .flatMap { acc ->
        ['CEN', 'ARMS'].collectMany { cat ->
            (1..params.num_chromosomes).collect { chr ->
                [acc, cat, chr]
            }
        }
    }
```

This replaces the hardcoded `['Col', 'Ler']` and chromosome range `(1..5)`.

### 3. **Backward Compatibility Layer** (main.nf:32-42)

Ensures old commands still work with a deprecation warning:

```groovy
// Backward compatibility: convert old col/ler params to new parent params
if (params.col_bed_file && !params.parent1_bed) {
    params.parent1_bed = params.col_bed_file
    params.parent1_name = 'Col'
    log.warn "⚠️  DEPRECATED: --col_bed_file is deprecated. Use --parent1_bed and --parent1_name instead."
}
```

### 4. **Dynamic PAF File Filtering** (workflows/charla_nf.nf:407-428)

PAF files are now filtered dynamically based on parent names:

```groovy
// OLD (hardcoded):
col_arms_paf_files = ... .filter { it.name.contains('ARMS_againstCol.paf') }

// NEW (dynamic):
parent1_arms_paf_files = ... .filter { it.name.contains("ARMS_against${params.parent1_name}.paf") }
```

## How to Use the Generalized Pipeline

### Example 1: Col×Ler Cross (Default - Backward Compatible)

```bash
nextflow run jacgonisa/charla_nf \
    --reads data/ColLer_F1_pollen.fq.gz \
    --sample_id ColLer_F1_pollen \
    --input_format fastq \
    --cenhapmer_db_dir /path/to/cenhapmers_ColLer \
    --reference_genomes_dir /path/to/references \
    --parent1_bed annotations/Col-0.bed \
    --parent2_bed annotations/Ler-0.bed \
    --outdir results/ColLer \
    -profile singularity
```

**OR** using old parameter names (with deprecation warning):
```bash
nextflow run jacgonisa/charla_nf \
    --reads data/ColLer_F1_pollen.fq.gz \
    --col_bed_file annotations/Col-0.bed \
    --ler_bed_file annotations/Ler-0.bed \
    --outdir results/ColLer \
    -profile singularity
```

### Example 2: Cvi×Est Cross (New Hybrid)

```bash
nextflow run jacgonisa/charla_nf \
    --reads data/CviEst_F1_pollen.fq.gz \
    --sample_id CviEst_F1_pollen \
    --input_format fastq \
    --parent1_name Cvi \
    --parent2_name Est \
    --cenhapmer_db_dir /path/to/cenhapmers_CviEst \
    --reference_genomes_dir /path/to/references \
    --parent1_bed annotations/Cvi.bed \
    --parent2_bed annotations/Est.bed \
    --outdir results/CviEst \
    -profile singularity
```

### Example 3: Different Species (e.g., with 7 chromosomes)

```bash
nextflow run jacgonisa/charla_nf \
    --reads data/species_F1.fq.gz \
    --sample_id species_F1 \
    --parent1_name StrainA \
    --parent2_name StrainB \
    --num_chromosomes 7 \
    --cenhapmer_db_dir /path/to/cenhapmers \
    --parent1_bed annotations/StrainA.bed \
    --parent2_bed annotations/StrainB.bed \
    --outdir results/species \
    -profile singularity
```

## What Changed Under the Hood

### File Structure
```
charla_nf/
├── GENERALIZATION_PLAN.md         # NEW: Detailed implementation plan
├── GENERALIZATION_SUMMARY.md      # NEW: This file - user guide
├── main.nf                         # UPDATED: Dynamic parent names + backward compat
├── nextflow.config                 # UPDATED: Generic parameters
└── workflows/charla_nf.nf         # UPDATED: Dynamic PAF filtering
```

### Key Code Sections Updated

| File | Line Range | What Changed |
|------|------------|--------------|
| `nextflow.config` | 14-15 | Made `reads` and `sample_id` parameters generic (removed Col/Ler defaults) |
| `nextflow.config` | 45-54 | Added new parent1/parent2 parameters and `num_chromosomes` |
| `nextflow.config` | 393 | Updated help command to show new parameters |
| `main.nf` | 32-42 | Added backward compatibility layer |
| `main.nf` | 60-73 | Changed from hardcoded `['Col', 'Ler']` to `[params.parent1_name, params.parent2_name]` |
| `workflows/charla_nf.nf` | 392-405 | Changed BED file validation to use generic parent names |
| `workflows/charla_nf.nf` | 407-428 | Updated PAF file filtering with dynamic parent names |
| `workflows/charla_nf.nf` | 444-453 | Updated PLOT_CROSSOVER_MAP call with generic parent references |

## Benefits of This Generalization

1. **Flexibility**: Works with ANY F1 hybrid cross, not just Col×Ler
2. **Extensibility**: Easy to add new accessions - just change parameter values
3. **Species-agnostic**: Can work with organisms that have different chromosome numbers
4. **Backward compatible**: Old commands still work (with warnings)
5. **Future-proof**: Easy to extend for batch processing (sample sheets)
6. **Clear naming**: `parent1`/`parent2` is more intuitive than `col`/`ler`

## What Still Needs to Be Done (Optional)

### 1. Sample Sheet for Batch Processing
Create an input CSV schema to process multiple samples at once:

```csv
sample_id,reads,input_format,parent1,parent2,parent1_bed,parent2_bed
ColLer_pollen,data/pollen1.fq.gz,fastq,Col,Ler,annotations/Col.bed,annotations/Ler.bed
CviEst_sperm,data/sperm1.fq.gz,fastq,Cvi,Est,annotations/Cvi.bed,annotations/Est.bed
```

### 2. Update README
Add examples showing:
- How to use new parameters
- Examples with different crosses
- Migration guide from old to new parameter names

### 3. Test with Real Data
Test the generalized pipeline with a non-Col×Ler cross to verify:
- K-mer database compatibility
- PAF file naming conventions
- Plotting scripts handle dynamic names

## Testing the Changes

The changes are backward-compatible, so your existing Col×Ler data should work identically:

```bash
# Test with your existing data (should work exactly as before)
nextflow run . \
    --reads /path/to/your/ColLer_data.fq.gz \
    --sample_id ColLer_test \
    --parent1_bed /path/to/Col-0.bed \
    --parent2_bed /path/to/Ler-0.bed \
    --cenhapmer_db_dir /path/to/your/cenhapmers \
    --reference_genomes_dir /path/to/references \
    --outdir test_generalization \
    -profile local
```

## Branches

- `feature/comprehensive-test-suite` - Comprehensive nf-test suite (17/17 tests passing ✅)
- `feature/generalize-accessions` - **Current branch** - Generalization work
- `master` - Original code

## Next Steps

1. **Test the generalization** with your existing Col×Ler data
2. **Try a different cross** (if you have data for Cvi, Est, Sha, Ws, etc.)
3. **Merge to master** when confident everything works
4. **Optional**: Implement sample sheet for batch processing

## Questions to Consider

1. Do your cenhapmer databases follow a naming convention that includes parent names (e.g., `Col_CEN_chr1.kmc_pre`)?
2. Do your reference genomes follow a similar naming convention?
3. Do you want to process multiple samples in a single run (batch mode)?
4. Do you have pollen/sperm data from different crosses you want to test?

Let me know if you'd like to proceed with sample sheet implementation or README updates!
