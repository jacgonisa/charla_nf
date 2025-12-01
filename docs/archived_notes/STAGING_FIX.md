# File Name Collision Fix

## Problem

```
ERROR ~ Process `CHARLA_NF:VISUALIZE_READ_STRUCTURES` input file name collision
-- There are multiple input files for each of the following file names:
   best_alignments.paf, summary_table.csv, inversions.paf, etc.
```

## Root Cause

Both ARMS and CEN mapping directories contain files with identical names:
- `best_alignments.paf`
- `summary_table.csv`
- `inversions.paf`
- etc.

When Nextflow stages these directories as inputs, it tries to put all files in the same work directory, causing collisions.

## Solution

Use `stageAs` to put ARMS and CEN files in separate subdirectories:

### Before (Caused Collision):
```groovy
input:
path arms_mapping_dir    // Files staged directly
path cen_mapping_dir     // Files staged directly → COLLISION!
```

Both directories' files end up in same location → duplicate names → error.

### After (No Collision):
```groovy
input:
path arms_mapping_dir, stageAs: 'arms_mapping/**'   // Files in arms_mapping/
path cen_mapping_dir, stageAs: 'cen_mapping/**'     // Files in cen_mapping/
```

Now files are in separate directories:
```
work/XX/YYYY/
├── arms_mapping/
│   └── summary_mm_sr/best_alignments.paf
│   └── summary_mm_ont/best_alignments.paf
│   └── ...
└── cen_mapping/
    └── summary_mm_sr/best_alignments.paf
    └── summary_mm_ont/best_alignments.paf
    └── ...
```

### Script Updated:
```bash
# Before:
find ${arms_mapping_dir} -name "best_alignments.paf" ...

# After:
find arms_mapping -name "best_alignments.paf" ...
find cen_mapping -name "best_alignments.paf" ...
```

## Files Modified

- `modules/local/visualize_read_structures.nf` (lines 13-14, 39-50)

## Ready to Resume

```bash
nextflow run charla_nf/main.nf ... -resume
```

The process will now stage files correctly without collisions.

---

**Status:** Fixed ✅
**Date:** 2025-11-27
