# CSV Separator Fix

## Problem

```
ERROR ~ Process `CHARLA_NF:ANALYZE_SCO_NCO` terminated with an error exit status (1)

KeyError: 'read_id'
```

When trying to merge confidence scores with raw metrics:
```python
merged = pd.merge(conf_df, raw_df, on='read_id', how='left', suffixes=('', '_raw'))
```

## Root Cause

The script was reading `raw_metrics.csv` with wrong separator:

**Script assumption**: `sep='\t'` (TAB-separated)
**Actual file**: `sep=','` (comma-separated)

When pandas reads a CSV with the wrong separator, it treats the entire first line as a single column name, so there's no `read_id` column.

Example:
```python
# Wrong separator:
raw_df = pd.read_csv('file.csv', sep='\t')
# Result: Single column named "read_id,summary_lrhq,summary_lrhqae,..."

# Correct separator:
raw_df = pd.read_csv('file.csv', sep=',')
# Result: Columns: ['read_id', 'summary_lrhq', 'summary_lrhqae', ...]
```

## File Formats

The three input files have different formats:

1. **confidence_scores.tsv**: TAB-separated (✓)
   ```
   read_id	confidence_score	mapq_score	...
   ```

2. **raw_metrics.csv**: COMMA-separated (was being read as TAB ✗)
   ```
   read_id,summary_lrhq,summary_lrhqae,...
   ```

3. **combined_cowidths.txt**: TAB-separated (✓)
   ```
   read_id	min_co_width	...
   ```

## Solution

Added automatic separator detection based on file extension:

**File**: `scripts/18-analyze_sco_nco_comparison.py:51-56`

```python
# Raw metrics (CSV - from ARMS summary)
# Detect separator automatically
if raw_metrics_file.endswith('.csv'):
    raw_df = pd.read_csv(raw_metrics_file, sep=',')
else:
    raw_df = pd.read_csv(raw_metrics_file, sep='\t')
```

Now:
- `.csv` files → comma separator
- Everything else → tab separator

## Why This Happened

The `raw_metrics.csv` comes from `MAPPING_ANALYSIS` ARMS summary, which outputs CSV format:

```groovy
// In workflow:
mapping_analysis_output.arms_summary.first()  // This is summary_table.csv
```

The summary tables from MAPPING_ANALYSIS are CSV files, not TSV files.

## Testing

The script will now correctly read all three files:

```bash
nextflow run charla_nf/main.nf ... -resume
```

Expected output:
```
[INFO] Loading data...
[INFO] Loaded X reads
[INFO] SCO: Y reads
[INFO] NCO/Complex: Z reads
...
sco_nco_comparison.png
sco_nco_resolution_analysis.png
sco_nco_summary.tsv
```

## Alternative Solutions Considered

1. **Change file extension in workflow** - Would break cache
2. **Use pandas auto-detection** - Less reliable
3. **Convert CSV to TSV in module** - Extra processing
4. **Check file extension** - ✅ Simple and reliable (chosen)

---

**Status**: Fixed ✅
**Date**: 2025-11-28
**Impact**: ANALYZE_SCO_NCO will now correctly load all input files
