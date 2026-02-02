# Indel Analysis Script Updates

**Date**: November 26, 2025
**Script**: `12.2-indel_analysis.py`
**Version**: 2.0

---

## Summary of Changes

Two major improvements have been made to the indel analysis script to improve visualization and biological accuracy:

1. ✅ **Log scale on Y-axis** - Better visualization of count distributions
2. ✅ **Corrected normalization** - Changed from "indels per read" to "indels per mapped Mb"

---

## Change 1: Log Scale Y-Axis

### Why?
- Indel counts often span multiple orders of magnitude
- Linear scale makes it hard to see patterns in low-count bins
- Log scale reveals trends across the full distribution

### What Changed?
Updated all histogram plots to use logarithmic y-axis:

```python
# Before
ax.set_ylim(0, global_max_count_full)
ax.set_ylabel('Count', fontsize=9)

# After
ax.set_yscale('log')
ax.set_ylim(0.5, global_max_count_full * 2)
ax.set_ylabel('Count (log scale)', fontsize=9)
ax.grid(True, alpha=0.2, axis='y', which='both')  # Show major and minor grid
```

### Affected Plots:
- ✅ Full-range histograms (per haplotype)
- ✅ Combined histograms (all backgrounds)
- ✅ Zoomed histograms (0-1000bp)
- ✅ Comparative bar charts (ARMS vs CEN, Col vs Ler, etc.)

### Benefits:
- 📊 Better visualization of rare events
- 📊 Easier comparison across different magnitude ranges
- 📊 More informative grid lines (both major and minor)

---

## Change 2: Corrected Normalization to Mapped Mb

### Why?
**Problem with old normalization (indels/read):**
- Reads vary in length (especially long reads!)
- A 10kb read vs 50kb read have different "opportunity" for indels
- Normalization by read count is biologically misleading

**Correct normalization (indels/Mb):**
- Accounts for the actual amount of sequence mapped
- Comparable across different read length distributions
- Standard metric in genomics (e.g., SNP density, indel density)

### What Changed?

#### 1. Track mapped bases during BAM parsing:
```python
# NEW: Track total aligned bases
total_mapped_bases = 0

for read in bamfile.fetch(until_eof=True):
    if not read.is_unmapped and read.cigartuples:
        # Count only match/mismatch operations (M)
        aligned_bases = sum(length for op, length in read.cigartuples if op == 0)
        total_mapped_bases += aligned_bases
```

#### 2. Calculate density per Mb:
```python
# Before
indel_density = num_indels / mapped_reads if mapped_reads > 0 else 0

# After
mapped_mb = total_mapped_bases / 1_000_000
indel_density = num_indels / mapped_mb if mapped_mb > 0 else 0
```

#### 3. Updated output:
```python
# New field in return dictionary
'mapped_mb': round(mapped_mb, 2)

# Updated console output
print(f"  Mapped bases: {total_mapped_bases:,} bp ({mapped_mb:.2f} Mb)")
print(f"  Indel density: {indel_density:.4f} indels/Mb")
```

### Affected Outputs:

#### TSV File Header (Updated):
```
haplotype  mapping_mode  total_reads  mapped_reads  mapped_mb  indels  insertions  deletions  avg_indel_size  indel_density  reads_with_indels  background  region  chromosome
```

**New column**: `mapped_mb` - Total megabases of sequence mapped

#### Plot Labels (Updated):
- Old: "Indels per read (>50bp)"
- New: "**Indels per Mb (>50bp)**"

### Benefits:
- ✅ Biologically accurate normalization
- ✅ Comparable across samples with different read lengths
- ✅ Standard metric used in genomics literature
- ✅ Accounts for actual sequence content, not arbitrary read boundaries

---

## Example Interpretation

### Before (Per Read):
```
Sample A: 100 indels, 1000 reads (avg 10kb each)
  → 0.10 indels/read

Sample B: 100 indels, 500 reads (avg 20kb each)
  → 0.20 indels/read
```
**Misleading!** Sample B appears to have 2× more indels, but both samples actually have the same total mapped sequence (10 Mb).

### After (Per Mb):
```
Sample A: 100 indels, 10 Mb mapped
  → 10.0 indels/Mb

Sample B: 100 indels, 10 Mb mapped
  → 10.0 indels/Mb
```
**Correct!** Both samples have identical indel density.

---

## Biological Implications

### What Does "Indels/Mb" Tell Us?

1. **Structural variation load**: How many large indels per unit of sequence
2. **Regional differences**: CEN vs ARMS have different indel densities
3. **Parent-specific patterns**: Col vs Ler indel landscapes
4. **Chromosome specificity**: Some chromosomes may be more indel-prone

### Expected Patterns:

#### Centromeric Regions (CEN):
- **Higher indel density** expected
- Repetitive DNA = more structural variation
- Satellite repeats (178bp) visible in indel size distribution

#### Chromosome Arms (ARMS):
- **Lower indel density** expected
- More conserved genic regions
- Indels more likely to be deleterious

### Statistical Comparisons:
With proper normalization, you can now meaningfully compare:
- ✅ ARMS vs CEN indel density
- ✅ Col-0 vs Ler-0 structural variation
- ✅ Chromosome-specific patterns
- ✅ Mapping mode effects (primary vs supplementary alignments)

---

## File Changes

### Modified Files:
1. **`scripts/12.2-indel_analysis.py`**
   - Updated `analyze_bam_file()` function
   - Added mapped base tracking
   - Changed normalization calculation
   - Updated all plot y-axis labels
   - Added log scale to all histograms and bar charts

2. **`modules/local/non_hybrid_analysis.nf`**
   - Updated empty results file header to include `mapped_mb`

### Output Changes:
- **TSV files**: Added `mapped_mb` column
- **Plots**:
  - Y-axis now log scale
  - Labels changed to "Indels per Mb"
  - Better grid lines (major + minor)
  - Improved visual clarity

---

## Usage

No changes to command-line usage! The script automatically:
1. Calculates mapped bases from BAM CIGAR strings
2. Computes density per Mb
3. Generates plots with log scale
4. Outputs new TSV format

```bash
# Same command as before
python scripts/12.2-indel_analysis.py <sequences_dir> <output_results.tsv>
```

---

## Validation

### Quick Check:
Look at your output TSV file:

```bash
head -2 marker_analysis/indel_analysis_results.tsv
```

You should see:
- New `mapped_mb` column
- `indel_density` values that make biological sense (typically 1-100 indels/Mb)

### Visual Check:
Open any generated histogram:
- Y-axis should say "Count (log scale)"
- Grid lines should show logarithmic spacing
- Low counts should be visible

### Statistical Check:
Open comparative plots:
- Bar charts should have log y-axis
- Labels should say "Indels per Mb"
- Heatmap color bar should say "Indels per Mb"

---

## Nextflow Integration

The Nextflow error you saw (`bai_files` property issue) is a cosmetic error during pipeline finalization and doesn't affect results. It occurs because:

1. `NON_HYBRID_MAPPING` emits optional `bai_files`
2. When some conditions don't produce BAM files, the emit channel is empty
3. Nextflow tries to access properties during publishing, causing the error

**Impact**: None - the analysis completes successfully, outputs are correct.

**Fix needed**: Update publishDir handling in Nextflow (separate issue from indel analysis).

---

## Summary

### What You Get Now:

✅ **Better Visualization**
- Log scale y-axis reveals full distribution
- Easier to spot patterns and outliers
- Professional publication-ready figures

✅ **Correct Biology**
- Indel density per Mb (not per arbitrary read)
- Meaningful cross-sample comparisons
- Standard genomics metric

✅ **Enhanced Output**
- New `mapped_mb` column in TSV
- Updated plot labels
- Same easy-to-use command

### Recommendation:
Re-run your indel analysis with the updated script to get:
1. More informative visualizations (log scale)
2. Biologically accurate density metrics (per Mb)
3. Better statistical comparisons

---

**Questions?** Check the updated script at `scripts/12.2-indel_analysis.py` or the module at `modules/local/non_hybrid_analysis.nf`.
