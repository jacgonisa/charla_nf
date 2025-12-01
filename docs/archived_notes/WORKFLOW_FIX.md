# Workflow Variable Name Fix

## Problem

```
ERROR ~ No such variable: curate_profiles_output
 -- Check script 'charla_nf/workflows/charla_nf.nf' at line: 445
```

## Root Cause

Used wrong variable name in workflow. The process output is called `curate_hybrid_profiles_output`, not `curate_profiles_output`.

## Fix

**File**: `workflows/charla_nf.nf:445`

**Before**:
```groovy
curate_profiles_output.curated_profiles.first(),  // WRONG variable name
```

**After**:
```groovy
curate_hybrid_profiles_output.cut_table.first(),  // CORRECT
```

## Output Mapping

From `CURATE_HYBRID_PROFILES` process:

```groovy
curate_hybrid_profiles_output = CURATE_HYBRID_PROFILES(...)

// Available outputs:
curate_hybrid_profiles_output.final_table        // *_final_with_combo.tsv
curate_hybrid_profiles_output.ultracurated_table // *_ultracurated.tsv
curate_hybrid_profiles_output.cut_table          // *_ultracurated_cut.tsv ← WE USE THIS
```

We use `cut_table` because it contains:
- `readname`
- `threshold`
- `ultra_curated_profile` (segments like "325LA4,319LA3")
- `hybrid_type` (SCO, NCO, etc.)
- `ectopic` (yes/no/None)

This is exactly what the confidence scoring script needs!

---

**Status**: Fixed ✅
**Date**: 2025-11-28
**Ready to**: Resume pipeline
