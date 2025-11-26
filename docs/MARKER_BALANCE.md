# Marker Balance Analysis and Normalization

## Table of Contents
1. [Overview](#overview)
2. [Why Marker Balance Matters](#why-marker-balance-matters)
3. [Analyzing Marker Balance](#analyzing-marker-balance)
4. [Normalization Strategies](#normalization-strategies)
5. [Integration with Pipeline](#integration-with-pipeline)
6. [Troubleshooting](#troubleshooting)

---

## Overview

Cenhapmers (centromeric haplotype-specific k-mers) are unique markers that distinguish between parental haplotypes. However, the marker database may have **imbalances** that can bias crossover detection:

- **Parent imbalance**: Col-0 may have more unique markers than Ler-0 (or vice versa)
- **Region imbalance**: ARMS regions typically have more markers than CEN (centromeric) regions
- **Chromosome imbalance**: Some chromosomes may have more divergence and thus more markers

These imbalances can lead to:
- ✗ False crossover calls (preferentially switching to marker-rich parent)
- ✗ Missed crossovers in marker-poor regions
- ✗ Biased haplotype assignment

---

## Why Marker Balance Matters

### Problem 1: Unequal Parent Markers

**Scenario**: Col-0 has 1.5× more markers than Ler-0

**Impact**:
```
Read with true genotype: Col-Ler-Col (two crossovers)

Without normalization:
- Col segments: High k-mer counts (many markers available)
- Ler segments: Low k-mer counts (fewer markers available)
- Risk: Ler segment may not pass threshold → missed crossover!

With normalization:
- Counts converted to probabilities (count / total_markers)
- Both parents now comparable
- Crossover correctly detected
```

### Problem 2: Region-Specific Density

**Scenario**: ARMS have 3× more markers than CEN

**Impact**:
- Crossovers in centromeric regions may be less confident
- Segmentation may fail in low-marker CEN regions
- Crossover resolution is poorer in CEN

**Mitigation**: Region-aware confidence scoring (already in pipeline!)

### Problem 3: Chromosome Variability

**Scenario**: Chr1 has 2× more markers than Chr5

**Impact**:
- Inconsistent detection sensitivity across chromosomes
- May miss crossovers on marker-poor chromosomes

**Mitigation**: Chromosome-specific thresholding or global normalization

---

## Analyzing Marker Balance

### Step 1: Run Marker Balance Analysis

```bash
python scripts/16-analyze_marker_balance.py \
    --cenhapmer-dir /path/to/cenhapmers_k31 \
    --output-prefix results/marker_analysis/k31_balance
```

**Required Input**:
- Cenhapmer KMC databases (`.kmc_pre` / `.kmc_suf` files)
- Expected structure:
  ```
  cenhapmers_k31/
  ├── Col-0_CA1.kmc_pre/suf  # Col chromosome 1 arms
  ├── Col-0_CC1.kmc_pre/suf  # Col chromosome 1 centromere
  ├── Ler-0_LA1.kmc_pre/suf  # Ler chromosome 1 arms
  ├── Ler-0_LC1.kmc_pre/suf  # Ler chromosome 1 centromere
  └── ...
  ```

**Outputs Generated**:

1. **`k31_balance_marker_counts.tsv`** - Raw marker counts per database
   ```
   database    parent  region_type  chromosome  kmer_count
   Col-0_CA1   Col-0   ARMS         1           125432
   Col-0_CC1   Col-0   CEN          1           45678
   Ler-0_LA1   Ler-0   ARMS         1           118234
   ...
   ```

2. **`k31_balance_balance_summary.txt`** - Statistical summary
   ```
   MARKER BALANCE SUMMARY
   ======================

   Parent Totals:
     Col-0: 1,234,567
     Ler-0: 1,098,234

   Imbalance Ratio: 1.12:1
   Imbalance Percentage: 12.4%

   ⚡ CAUTION: >10% imbalance detected
   ```

3. **`k31_balance_marker_distribution.png`** - Visualizations
   - Total markers per parent (bar chart)
   - Markers by region type (CEN vs ARMS)
   - Markers per chromosome
   - Heatmap: parent × region × chromosome

4. **`k31_balance_per_database.png`** - Detailed per-database bar chart

5. **`k31_balance_normalization_recommendations.txt`** - Actionable guidance
   ```
   STATUS: ⚠️  HIGH IMBALANCE (>20%)

   RECOMMENDED STRATEGIES:
   1. REGENERATE CENHAPMER DATABASE (BEST):
      - Use stricter filtering for over-represented parent
      - Aim for <10% imbalance

   2. NORMALIZE k-mer counts (current approach):
      - Use probability scores: count / total_markers
      - Apply weighted normalization
   ...
   ```

### Step 2: Interpret Results

**Imbalance Levels**:
- **<10%**: ✓ Well-balanced, no action needed
- **10-20%**: ⚡ Moderate, consider normalization for precision work
- **>20%**: ⚠️ High imbalance, normalization strongly recommended

**Region Imbalance (ARMS:CEN ratio)**:
- **Expected**: ARMS have 2-4× more markers than CEN (centromeres are repetitive)
- **Impact**: Use region-aware confidence scoring (already integrated!)

**Chromosome Imbalance (CV >30%)**:
- **Impact**: Detection sensitivity varies across chromosomes
- **Mitigation**: Global normalization or chromosome-specific thresholds

---

## Normalization Strategies

### Strategy 1: Probability Normalization (BEST for high imbalance)

**Method**: Convert raw counts to probabilities
```
Normalized_score = (raw_count / total_markers_in_region) × scale_factor
```

**When to use**: Imbalance >20%

**Command**:
```bash
python scripts/17-normalize_kmer_counts.py \
    --input results/processed_cenhapmers/kmer_matrix.tsv \
    --output results/processed_cenhapmers/kmer_matrix_normalized.tsv \
    --marker-stats results/marker_analysis/k31_balance_marker_counts.tsv \
    --method probability \
    --scale-factor 10000
```

**Pros**:
- Fully corrects for marker availability
- Directly comparable scores across parents
- Works even with extreme imbalances

**Cons**:
- Scores are scaled probabilities, not raw counts
- May need to adjust downstream thresholds

### Strategy 2: Weighted Normalization (GOOD for moderate imbalance)

**Method**: Weight by inverse frequency relative to median
```
Weight = median_markers / markers_in_region
Normalized_count = raw_count × weight
```

**When to use**: Imbalance 10-30%

**Command**:
```bash
python scripts/17-normalize_kmer_counts.py \
    --input results/processed_cenhapmers/kmer_matrix.tsv \
    --output results/processed_cenhapmers/kmer_matrix_normalized.tsv \
    --marker-stats results/marker_analysis/k31_balance_marker_counts.tsv \
    --method weighted
```

**Pros**:
- Balances marker availability
- Preserves relative count magnitudes
- Works with existing thresholds

**Cons**:
- Less effective for extreme imbalances

### Strategy 3: Global Mean Normalization (SIMPLE)

**Method**: Normalize by global mean across all regions
```
Weight = global_mean / markers_in_region
Normalized_count = raw_count × weight
```

**When to use**: General-purpose, any imbalance level

**Command**:
```bash
python scripts/17-normalize_kmer_counts.py \
    --input results/processed_cenhapmers/kmer_matrix.tsv \
    --output results/processed_cenhapmers/kmer_matrix_normalized.tsv \
    --marker-stats results/marker_analysis/k31_balance_marker_counts.tsv \
    --method global_mean
```

**Pros**:
- Simple and robust
- Works across all imbalance levels
- Easy to interpret

**Cons**:
- May not fully correct extreme imbalances

### Strategy 4: Regenerate Cenhapmer Database (ULTIMATE FIX)

If imbalance is >30%, consider regenerating the cenhapmer database with adjusted parameters:

**For over-represented parent** (e.g., Col-0 has too many markers):
```bash
# Use stricter filtering (higher k-mer count threshold)
# In cenhapmer generation pipeline:
--min-count 30  # Instead of default 25
```

**For under-represented parent** (e.g., Ler-0 has too few markers):
```bash
# Use more lenient filtering
--min-count 20  # Lower threshold to include more markers
```

**Goal**: Achieve <10% imbalance between parents

---

## Integration with Pipeline

### Manual Integration (Current)

1. **Run marker balance analysis** (one-time):
   ```bash
   python scripts/16-analyze_marker_balance.py \
       --cenhapmer-dir /path/to/cenhapmers_k31 \
       --output-prefix results/marker_analysis/k31_balance
   ```

2. **Check results** and decide if normalization is needed

3. **If imbalance >10%, normalize k-mer matrix**:
   ```bash
   # After process_cenhapmer_counts step completes
   python scripts/17-normalize_kmer_counts.py \
       --input results/processed_cenhapmers/threshold50/kmer_matrix.tsv \
       --output results/processed_cenhapmers/threshold50/kmer_matrix_normalized.tsv \
       --marker-stats results/marker_analysis/k31_balance_marker_counts.tsv \
       --method weighted
   ```

4. **Use normalized matrix** in downstream analysis:
   ```bash
   # Manually replace the input for segmentation
   cp results/processed_cenhapmers/threshold50/kmer_matrix_normalized.tsv \
      results/processed_cenhapmers/threshold50/kmer_matrix.tsv

   # Re-run from SEGMENT_READS with -resume
   nextflow run main.nf -resume
   ```

### Future Automated Integration (Planned)

The normalization could be integrated as an optional pipeline step:

```groovy
// In nextflow.config
params.normalize_markers = false  // Set to true to enable
params.normalization_method = 'weighted'

// In workflow
if (params.normalize_markers) {
    normalized_matrix = NORMALIZE_MARKER_COUNTS(
        kmer_matrix,
        marker_stats,
        params.normalization_method
    )
} else {
    normalized_matrix = kmer_matrix
}
```

---

## Confidence Scoring Integration

The confidence scoring system already accounts for marker balance through the **Parent Discrimination** component:

```python
# From 14-calculate_confidence_scores.py
def parent_discrimination_score(paf_data):
    """
    Measure clarity of parent assignment based on AS-NM difference.

    Higher score = clear winner (one parent much better than other)
    Lower score = ambiguous (similar AS-NM for both parents)
    """
    col_score = calculate_as_nm_score(col_paf)
    ler_score = calculate_as_nm_score(ler_paf)

    diff_pct = abs(col_score - ler_score) / max(col_score, ler_score)

    if diff_pct > 0.20:  # >20% difference
        return 1.0  # High confidence, clear winner
    elif diff_pct > 0.10:  # 10-20% difference
        return 0.7
    else:
        return 0.3  # <10% difference, ambiguous
```

**How marker balance affects this**:
- If Col-0 has 2× more markers, reads may preferentially map to Col-0
- The AS-NM difference will be smaller (both parents match well)
- Parent discrimination score will be lower
- Overall confidence score decreases → flagged for review

**This is good!** It automatically penalizes reads in imbalanced regions.

---

## Troubleshooting

### Issue 1: "No KMC databases found"

**Cause**: Cenhapmer directory structure doesn't match expected format

**Solution**:
```bash
# Check your directory structure
ls -R /path/to/cenhapmers_k31/

# Expected structure:
# Col-0_CA1.kmc_pre, Col-0_CA1.kmc_suf
# Col-0_CC1.kmc_pre, Col-0_CC1.kmc_suf
# Ler-0_LA1.kmc_pre, Ler-0_LA1.kmc_suf
# ...

# If files are in subdirectories, adjust path
python scripts/16-analyze_marker_balance.py \
    --cenhapmer-dir /path/to/cenhapmers_k31/databases \
    --output-prefix results/marker_analysis/k31_balance
```

### Issue 2: "kmc_tools not found"

**Cause**: KMC tools not in PATH

**Solution**:
```bash
# Add KMC to PATH
export PATH=/path/to/kmc/bin:$PATH

# Or use conda environment
conda activate charla_nf
```

### Issue 3: Normalization doesn't improve results

**Cause**: Imbalance may not be the primary issue

**Check**:
1. **Sequencing quality**: Low-quality reads cause poor segmentation
2. **K-mer size**: k=31 may be better than k=21 for divergent regions
3. **Segmentation thresholds**: Adjust `min_segment_length` and `min_token_count`
4. **Reference genomes**: Ensure high-quality assemblies

### Issue 4: Extreme imbalance (>50%)

**Cause**: Cenhapmer generation parameters need adjustment

**Solution**: **Regenerate cenhapmer database** with balanced filtering:
```bash
# See cenhapmer_generation/ directory for pipeline
# Adjust --min-count parameter based on marker analysis results
```

---

## Best Practices

1. **Always run marker balance analysis first**
   - Know your baseline before deciding on normalization

2. **Start with moderate normalization**
   - Try `weighted` method first
   - Only use `probability` if imbalance >20%

3. **Validate normalization effectiveness**
   - Re-run confidence scoring after normalization
   - Check if crossover calls improve
   - Compare high/medium/low confidence ratios

4. **Document your normalization**
   - Record method and parameters used
   - Include in publication methods section

5. **Consider biological interpretation**
   - Some imbalance is expected (e.g., CEN vs ARMS)
   - Don't over-normalize biological variation

---

## References

- **Cenhapmer generation pipeline**: `cenhapmer_generation/README.md`
- **Confidence scoring documentation**: `docs/CONFIDENCE_SCORING.md`
- **Segmentation parameters**: `nextflow.config` (lines 20-25)

---

## Quick Start Example

```bash
# 1. Analyze marker balance
python scripts/16-analyze_marker_balance.py \
    --cenhapmer-dir data/cenhapmers_k31 \
    --output-prefix results/marker_analysis/k31_balance

# 2. Check imbalance percentage in summary file
cat results/marker_analysis/k31_balance_balance_summary.txt

# 3. If imbalance >10%, normalize
python scripts/17-normalize_kmer_counts.py \
    --input results/processed_cenhapmers/threshold50/kmer_matrix.tsv \
    --output results/processed_cenhapmers/threshold50/kmer_matrix_normalized.tsv \
    --marker-stats results/marker_analysis/k31_balance_marker_counts.tsv \
    --method weighted

# 4. Replace original matrix and re-run
cp results/processed_cenhapmers/threshold50/kmer_matrix_normalized.tsv \
   results/processed_cenhapmers/threshold50/kmer_matrix.tsv

nextflow run main.nf -resume
```

---

**Questions or issues?** Open a GitHub issue: https://github.com/jacgonisa/charla_nf/issues
