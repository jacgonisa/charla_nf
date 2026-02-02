# Integration Status - Post-Mapping Analysis

This document clearly shows which parts of the post-mapping analysis module are **INTEGRATED** into the main CHARLA Nextflow pipeline and which are **STANDALONE** (manual use only).

Last updated: 2025-02-02

---

## ✅ INTEGRATED Modules

These modules are **already working** in the main CHARLA pipeline (`workflows/charla_nf.nf`):

### 1. CROSSOVER_STATS
**Module**: `postmapping_analysis/modules/crossover_stats.nf`
**Script**: `postmapping_analysis/scripts/14-crossover_statistics.R`
**Status**: ✅ **INTEGRATED**
**When it runs**: Automatically after mapping analysis
**Output**: Crossover statistics (event types, segment sizes, etc.)

### 2. CURATE_CONFIDENT_READS
**Module**: `postmapping_analysis/modules/curate_confident_reads.nf`
**Script**: `postmapping_analysis/scripts/15-curate_confident_reads.R`
**Status**: ✅ **INTEGRATED**
**When it runs**: Automatically after mapping analysis
**Output**:
- Curated confident read BED files (SCO, NCO, GC-CO)
- Pre-mapping coordinates (read positions)
- Post-mapping coordinates (genome positions)
- Read ID lists

### 3. CREATE_COMPREHENSIVE_ALIGNMENT_OUTPUT
**Module**: `postmapping_analysis/modules/create_comprehensive_alignment_output.nf`
**Script**: `postmapping_analysis/scripts/19-create_comprehensive_alignment_output.R`
**Status**: ✅ **INTEGRATED**
**When it runs**: Automatically after confident reads curation
**Output**: Single-row comprehensive TSV with:
- Read positions from k-mer analysis
- Reference genome positions from PAF alignments
- Segments ordered by reference position
- Complete strand information

**This is the CRITICAL module** that provides the collaborator-compatible format!

---

## ⚠️ NOT YET INTEGRATED (Available for Manual Use)

These modules are **available** in `postmapping_analysis/` but **NOT automatically run** by the pipeline. You can use them manually:

### 4. PLOT_CONFIDENT_READS
**Module**: `postmapping_analysis/modules/plot_confident_reads.nf`
**Script**: `postmapping_analysis/scripts/18-plot_confident_reads_comprehensive.R`
**Status**: ⚠️ **NOT INTEGRATED**
**Why**: Optional visualization, requires additional dependencies
**How to use manually**:
```bash
Rscript postmapping_analysis/scripts/18-plot_confident_reads_comprehensive.R \
  --confident_bed_sco path/to/SCO.bed \
  --confident_bed_nco path/to/NCO.bed \
  --sample_id sample_name \
  --output_dir output/
```

**To integrate**: Uncomment lines 666-684 in `workflows/charla_nf.nf`

---

### 5. KARYOPLOT_CROSSOVERS
**Module**: `postmapping_analysis/modules/karyoplot_crossovers.nf`
**Script**: `postmapping_analysis/scripts/15-karyoplot_crossovers.R`
**Status**: ⚠️ **NOT INTEGRATED**
**Why**: Requires karyoploteR R package and BSgenome objects
**How to use manually**:
```bash
Rscript postmapping_analysis/scripts/15-karyoplot_crossovers.R \
  --bed_file path/to/crossovers.bed \
  --parent1_fai parent1.fa.fai \
  --parent2_fai parent2.fa.fai \
  --sample_id sample_name \
  --output_dir output/
```

**To integrate**: Uncomment lines 726-742 in `workflows/charla_nf.nf`

---

### 6. Single Chromosome Plot (Script 17)
**Script**: `postmapping_analysis/scripts/17-single_chr_crossover_plot.R`
**Module**: ❌ **NO MODULE YET**
**Status**: ⚠️ **STANDALONE ONLY**
**How to use manually**:
```bash
Rscript postmapping_analysis/scripts/17-single_chr_crossover_plot.R \
  --input comprehensive_alignments.tsv \
  --chr chr5 \
  --output output/chr5_plot.pdf
```

**To integrate**: Create Nextflow module first

---

### 7. Circos Plot (Script 16)
**Script**: `postmapping_analysis/scripts/16-circos_and_karyoplot_crossovers.R`
**Module**: ❌ **NO MODULE YET**
**Status**: ⚠️ **STANDALONE ONLY**
**How to use manually**:
```bash
Rscript postmapping_analysis/scripts/16-circos_and_karyoplot_crossovers.R \
  --input comprehensive_alignments.tsv \
  --output output/circos_plot.pdf
```

**To integrate**: Create Nextflow module first

---

## Summary Table

| Module/Script | Integration Status | Auto-Run | Manual Use | Notes |
|--------------|-------------------|----------|------------|-------|
| **CROSSOVER_STATS** (14) | ✅ Integrated | Yes | Yes | Works automatically |
| **CURATE_CONFIDENT_READS** (15) | ✅ Integrated | Yes | Yes | Works automatically |
| **CREATE_COMPREHENSIVE_OUTPUT** (19) | ✅ Integrated | Yes | Yes | **Critical - provides main output** |
| **PLOT_CONFIDENT_READS** (18) | ⚠️ Commented out | No | Yes | Easy to integrate (uncomment) |
| **KARYOPLOT_CROSSOVERS** (15) | ⚠️ Commented out | No | Yes | Easy to integrate (uncomment) |
| **Single Chr Plot** (17) | ❌ Not integrated | No | Yes | Needs module creation |
| **Circos Plot** (16) | ❌ Not integrated | No | Yes | Needs module creation |

---

## What Gets Generated Automatically?

When you run the CHARLA pipeline now, you **automatically get**:

### From Integrated Modules

1. **Comprehensive Alignment TSV** (script 19)
   ```
   sample_id_comprehensive_alignments.tsv
   ```
   - Single-row per read format
   - Both read and genome coordinates
   - Ordered by reference position
   - **Ready for analysis!**

2. **Curated Confident Reads** (script 15)
   ```
   sample_id_confident_SCO_premapping_readcoords.bed
   sample_id_confident_SCO_postmapping_refcoords.bed
   sample_id_confident_NCO_premapping_readcoords.bed
   sample_id_confident_NCO_postmapping_refcoords.bed
   sample_id_confident_GC-CO_premapping_readcoords.bed
   sample_id_confident_GC-CO_postmapping_refcoords.bed
   sample_id_confident_read_ids.txt
   ```

3. **Crossover Statistics** (script 14)
   ```
   sample_id_crossover_statistics.tsv
   sample_id_segment_size_distribution.tsv
   ```

### NOT Generated Automatically (manual only)

4. **Visualization plots** (scripts 15, 16, 17, 18)
   - Karyoplots
   - Circos plots
   - Single chromosome plots
   - Comprehensive read plots

---

## How to Enable Visualization in Pipeline

### Option 1: Quick Enable (Uncommenting)

Edit `workflows/charla_nf.nf`:

1. **Uncomment include statements** (lines 44-45):
   ```groovy
   include { PLOT_CONFIDENT_READS } from '../postmapping_analysis/modules/plot_confident_reads.nf'
   include { KARYOPLOT_CROSSOVERS } from '../postmapping_analysis/modules/karyoplot_crossovers.nf'
   ```

2. **Uncomment PLOT_CONFIDENT_READS call** (lines 666-684):
   ```groovy
   confident_plots_output = PLOT_CONFIDENT_READS(...)
   ```

3. **Uncomment KARYOPLOT_CROSSOVERS call** (lines 726-742):
   ```groovy
   karyoplot_output = KARYOPLOT_CROSSOVERS(...)
   ```

### Option 2: Parameter-Based Enable (Better)

Add to `nextflow.config`:
```groovy
params {
    // ... existing parameters ...

    // Post-mapping visualization (optional)
    run_visualization = false  // Set to true to enable
}
```

Then in workflow, wrap calls in:
```groovy
if (params.run_visualization) {
    confident_plots_output = PLOT_CONFIDENT_READS(...)
    karyoplot_output = KARYOPLOT_CROSSOVERS(...)
}
```

---

## Dependencies

### Always Required (Integrated Modules)
```r
dplyr
tidyr
optparse
```

### Only for Visualization (Manual Use)
```r
ggplot2
karyoploteR    # For karyoplots
BSgenome       # For karyoplots
regioneR       # For karyoplots
RColorBrewer   # For plotting
viridis        # For plotting
```

---

## Testing Integration

To verify the integrated modules work:

```bash
nextflow run main.nf \
  --input samples.csv \
  --reference_genomes_dir genomes/ \
  --parent1_name B73 \
  --parent2_name Mo17 \
  --outdir results/
```

**Expected outputs** in `results/curated_confident_reads/`:
- ✅ `*_comprehensive_alignments.tsv` (from script 19)
- ✅ `*_confident_*_premapping_readcoords.bed` (from script 15)
- ✅ `*_confident_*_postmapping_refcoords.bed` (from script 15)
- ✅ `*_crossover_statistics.tsv` (from script 14)
- ❌ Plots (not integrated yet)

---

## Next Steps for Full Integration

1. ✅ **Core modules integrated** (scripts 14, 15, 19)
2. ⚠️ **Add parameter for visualization** (optional)
3. ⚠️ **Test with real data**
4. ⚠️ **Create modules for scripts 16, 17**
5. ⚠️ **Add containerization** (Docker/Singularity)

---

## Questions?

- **"Why isn't visualization integrated?"** → Optional, requires extra dependencies
- **"Can I use scripts manually?"** → Yes! All scripts work standalone
- **"Which output should I use for analysis?"** → `*_comprehensive_alignments.tsv` from script 19
- **"How do I add karyoplots?"** → Uncomment lines in workflow or run script 15 manually

---

*Last updated: 2025-02-02 - Integration of core modules (14, 15, 19) complete*
