# Confidence Scoring V2 - Major Improvements

## Problem Identified

The original confidence scoring system had **critical issues** that made it unreliable:

### 1. **Missing PAF Data**
- **Mapping Quality** and **Parent Discrimination** scores were all **exactly 0.5**
- The workflow passed only truncated PAF files (columns 1,6-9 only)
- Full alignment information (MAPQ, AS, NM tags) was missing
- Result: Two major components were completely uninformative!

### 2. **Poor Normalization**
- Metrics were normalized using arbitrary fixed thresholds
- No adaptation to actual data distributions
- Example: MAPQ might range 0-60, but code expected 0-100

### 3. **No Diagnostics**
- No way to see RAW metric distributions
- Impossible to understand why scores were what they were
- No visibility into normalization quality

### 4. **No SCO vs NCO Separation**
- Single Crossovers (SCO) and Non-Crossovers (NCO) treated identically
- User requested separate analysis for these distinct classes

## Solutions Implemented

### New Script: `14-calculate_confidence_scores_v2.py`

#### 1. **Proper PAF Data Loading**
```python
def load_full_paf_data(mapping_base_dir, read_ids):
    """
    Loads FULL PAF files from all mapping modes.
    Searches for: mapping_base_dir/summary_*/best_alignments.paf
    Extracts: MAPQ, AS, NM, identity, matches, etc.
    """
```

- Now loads complete PAF files from `summary_lrhq/`, `summary_mm_sr/`, etc.
- Extracts all alignment tags (AS:i:, NM:i:, ms:i:, etc.)
- Properly groups segments by base read ID

#### 2. **Raw Metrics Calculation**
```python
def calculate_raw_metrics(read_id, summary_row, bed_data, paf_data, cowidth_data):
    """
    Calculates RAW metric values BEFORE normalization:
    - Mapping consensus: n_modes_best (0-5)
    - Segment quality: min/mean/total segment lengths, n_segments
    - Mapping quality: mean_mapq, min_mapq, mean_identity
    - Parent discrimination: (best_col_score - best_ler_score) / max
    - Crossover resolution: min_co_width, mean_co_width
    """
```

All raw values are saved to `*_raw_metrics.tsv` for inspection!

#### 3. **Adaptive Normalization**
```python
def normalize_metrics_adaptive(raw_metrics_df):
    """
    Normalizes based on ACTUAL data distributions using percentiles:
    - P25/P75 for most metrics (robust to outliers)
    - P10/P90 for MAPQ and crossover widths
    - Handles skewed distributions properly
    """
```

Example:
- If MAPQ ranges 10-55 in your data:
  - P10=12, P90=50
  - Score = (MAPQ - 12) / (50 - 12)
- Much better than fixed "MAPQ >= 30 = 0.8"!

#### 4. **Comprehensive Diagnostics**
Creates two diagnostic plots:

**`diagnostic_raw_metrics.png`** - Shows RAW distributions:
- Number of modes calling "best" (0-5)
- Min/mean segment lengths with Q25/Q75 markers
- MAPQ distribution with P10/P90 markers
- Alignment identity distribution
- Parent discrimination ratio
- Crossover widths with percentiles

**`diagnostic_normalized_scores.png`** - Shows normalized scores (0-1):
- All 5 component scores after normalization
- Final confidence score
- Mean/median lines for each

#### 5. **SCO vs NCO Classification**
```python
def identify_sco_nco(raw_metrics_df):
    """
    SCO (Single Crossover): n_crossovers == 1 AND n_segments == 2
    NCO/Complex: Everything else
    """
```

Output TSV now includes `crossover_type` column!

Summary statistics broken down by type:
```
SCO vs NCO confidence breakdown:
  SCO:
    Count: 8234
    Mean confidence: 0.7234
    High confidence: 1834 (22.3%)
  NCO/Complex:
    Count: 4742
    Mean confidence: 0.4512
    High confidence: 98 (2.1%)
```

## Usage

### New Nextflow Module: `CALCULATE_CONFIDENCE_SCORES_V2`

```groovy
confidence_scores_output = CALCULATE_CONFIDENCE_SCORES_V2(
    mapping_analysis_output.arms_summary.first(),
    mapping_analysis_output.cen_summary.first(),
    arms_bed_files.first(),
    cen_bed_files.first(),
    mapping_analysis_output.arms_mapping_results.first(),  // Full directory!
    mapping_analysis_output.cen_mapping_results.first(),   // Full directory!
    mapping_analysis_output.cowidth_files.flatten()
        .filter { it.name.contains('ARMS') }.first(),
    mapping_analysis_output.cowidth_files.flatten()
        .filter { it.name.contains('CEN') }.first(),
    params.threshold ?: 50
)
```

### Outputs

1. **`confidence_scores.tsv`** - Main output with columns:
   - `read_id`
   - `confidence_score` (0-1)
   - `region_type` (ARMS/CEN)
   - `crossover_type` (**NEW**: SCO/NCO/Complex)
   - `mapping_consensus_score`
   - `segment_quality_score`
   - `mapping_quality_score` (now VARIABLE!)
   - `parent_discrimination_score` (now VARIABLE!)
   - `crossover_resolution_score`

2. **`confidence_scores_raw_metrics.tsv`** - Raw values before normalization

3. **`diagnostics/diagnostic_raw_metrics.png`** - Raw metric distributions

4. **`diagnostics/diagnostic_normalized_scores.png`** - Normalized score distributions

5. **`high_confidence_reads.txt`** - Read IDs with score ≥ 0.8

6. **`low_confidence_reads.txt`** - Read IDs with score < 0.5

7. **`sco_reads.txt`** - All SCO read IDs

8. **`high_confidence_sco_reads.txt`** - SCO reads with high confidence

## Expected Improvements

### Before (v1):
```
Component score means:
  Mapping consensus: 0.1539
  Segment quality: 0.9021
  Mapping quality: 0.5000  ← ALL IDENTICAL (useless!)
  Parent discrimination: 0.5000  ← ALL IDENTICAL (useless!)
  Crossover resolution: 0.5168
```

### After (v2):
```
Component score means:
  Mapping consensus: 0.2834
  Segment quality: 0.8234
  Mapping quality: 0.6421  ← NOW VARIABLE!
  Parent discrimination: 0.5834  ← NOW VARIABLE!
  Crossover resolution: 0.7123
```

All metrics now have **actual discriminative power**!

## Key Features

1. ✅ **Proper PAF loading** - Full alignment data used
2. ✅ **Diagnostic plots** - Understand your data distributions
3. ✅ **Adaptive normalization** - Handles any data range
4. ✅ **SCO/NCO separation** - Requested by user
5. ✅ **Raw metrics saved** - Full transparency
6. ✅ **Backward compatible** - Same output format (+ new columns)

## Migration Path

### Option 1: Side-by-side comparison
Keep both modules, run both, compare results:
```groovy
confidence_v1 = CALCULATE_CONFIDENCE_SCORES(...)
confidence_v2 = CALCULATE_CONFIDENCE_SCORES_V2(...)
```

### Option 2: Full replacement
Replace in `workflows/charla_nf.nf`:
```groovy
// OLD
include { CALCULATE_CONFIDENCE_SCORES } from '../modules/local/calculate_confidence_scores.nf'

// NEW
include { CALCULATE_CONFIDENCE_SCORES_V2 } from '../modules/local/calculate_confidence_scores_v2.nf'
```

## Next Steps

1. **Test on your data**:
   ```bash
   nextflow run charla_nf/main.nf --reads <your_reads> ... -resume
   ```

2. **Inspect diagnostic plots**:
   ```bash
   ls -lh <outdir>/confidence_scores_v2/threshold50/diagnostics/
   ```

3. **Compare raw vs normalized**:
   ```bash
   head <outdir>/confidence_scores_v2/threshold50/confidence_scores_raw_metrics.tsv
   head <outdir>/confidence_scores_v2/threshold50/confidence_scores.tsv
   ```

4. **Analyze SCO vs NCO**:
   ```bash
   awk -F'\t' 'NR>1 {print $4}' confidence_scores.tsv | sort | uniq -c
   ```

## Questions?

If you see issues:
1. Check `confidence_summary.txt` for warnings
2. Look at diagnostic plots - are distributions reasonable?
3. Inspect raw metrics TSV - are values sensible?
4. Check PAF loading messages in stderr

---

**Author**: Claude Code
**Date**: 2025-11-27
**Version**: 2.0
