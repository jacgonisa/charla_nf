# Cache Cleanup After Bug Fix

## Problem

After fixing the classification bug in `scripts/20-calculate_confidence_v3.py`, the pipeline was still using OLD cached results from before the fix.

## Why `-resume` Didn't Work

Nextflow's `-resume` feature is smart about detecting changes:
- ✅ It WILL re-run if **input files** change
- ✅ It WILL re-run if **script content in the module** changes
- ❌ It will NOT re-run if only **external Python scripts** change (like `scripts/*.py`)

Since we modified `scripts/20-calculate_confidence_v3.py` (external to the Nextflow module), `-resume` didn't detect the change and used cached results.

## Solution

Manually delete the work directories for affected processes:

```bash
# 1. Delete CALCULATE_CONFIDENCE_V3 cache (generates confidence_scores.tsv)
rm -rf /mnt/ssd-8tb/Col_Ler_F1_pollen_q10_1kb_2023/work/77/5ffa5700eb6a9ba1014e547d596c0e

# 2. Delete ANALYZE_SCO_NCO cache (depends on confidence_scores.tsv)
rm -rf /mnt/ssd-8tb/Col_Ler_F1_pollen_q10_1kb_2023/work/d4/35c3d462e5bb4705a0835a00b666ce
```

Then run with `-resume` as usual:
```bash
nextflow run charla_nf/main.nf ... -resume
```

## How to Find Work Directories

If you need to clean other processes, you can find their work directories from the error message:

```
Work dir:
  /mnt/ssd-8tb/Col_Ler_F1_pollen_q10_1kb_2023/work/d4/35c3d462e5bb4705a0835a00b666ce
```

Or search for specific output files:
```bash
# Find where confidence_scores.tsv was generated
find /mnt/ssd-8tb/Col_Ler_F1_pollen_q10_1kb_2023/work -name "confidence_scores.tsv" -type f
```

## Alternative: Clean ALL Downstream Processes

If you want to be thorough and re-run everything from CALCULATE_CONFIDENCE_V3 onwards:

```bash
# Find all processes that depend on confidence scores
find /mnt/ssd-8tb/Col_Ler_F1_pollen_q10_1kb_2023/work -name "confidence_scores.tsv" -type l -printf "%h\n"

# Or use nextflow clean to remove failed/completed tasks
nextflow clean -f -k
```

## What Changed

The bug fix modified how segments are parsed:

**Before (WRONG)**:
```python
"325LA4"[0] = "3" → "X" → Complex
```

**After (CORRECT)**:
```python
"325LA4" → scan for 'L' → "L" → Parental (or SCO/NCO depending on pattern)
```

## Expected Results After Fix

Once the pipeline re-runs with the fixed script:

```
[INFO] Classification AFTER mapping (validated by mappers):
  Parental: ~61,751 (46.5%)
  SCO: ~64,877 (48.9%)      ← Was 0 before!
  NCO: ~5,088 (3.8%)
  DCO: ~661 (0.5%)
  Complex: ~286 (0.2%)      ← Was 2,868 (100%) before!
```

## Status

✅ Caches cleaned
✅ Pipeline ready to re-run with fixed classification
✅ Next run with `-resume` will regenerate confidence scores correctly

---

**Date**: 2025-11-28
**Affected**: CALCULATE_CONFIDENCE_V3, ANALYZE_SCO_NCO
