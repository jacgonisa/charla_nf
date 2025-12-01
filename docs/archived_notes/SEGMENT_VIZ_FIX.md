# Segment Visualization Fix

## Problem

```
ERROR ~ Process `CHARLA_NF:VISUALIZE_SEGMENTS` terminated with an error exit status (2)

21-visualize_read_segments.py: error: unrecognized arguments:
  CANDIDATE_reads_nonectopic_CEN.bed segments_Col_Ler_F1_pollen_q10_1kb_thr50_filtered.bed
```

## Root Cause

The workflow was passing ALL BED files from `segment_reads_output.bed_files`:
- `CANDIDATE_reads_nonectopic_ARMS.bed`
- `CANDIDATE_reads_nonectopic_CEN.bed`
- `segments_Col_Ler_F1_pollen_q10_1kb_thr50_filtered.bed`

But the script `21-visualize_read_segments.py` expects only ONE BED file (the filtered segments file).

## Solution

Filter the channel to select only the main segments BED file:

**File**: `workflows/charla_nf.nf:466-469`

```groovy
// Filter for the main segments BED file (not ARMS/CEN specific ones)
segments_bed = segment_reads_output.bed_files
    .flatten()
    .filter { it.name.contains('segments_') && it.name.contains('_filtered.bed') }
    .first()

segment_visualization_output = VISUALIZE_SEGMENTS(
    segments_bed,  // Now passes single file
    confidence_v3_output.confidence_scores.first(),
    params.max_structure_reads ?: 200
)
```

## What the Filter Does

1. `.flatten()` - Convert channel of files to individual file items
2. `.filter { ... }` - Select only files matching:
   - Contains `segments_` in filename
   - Contains `_filtered.bed` in filename
3. `.first()` - Take the first (only) matching file

This selects: `segments_Col_Ler_F1_pollen_q10_1kb_thr50_filtered.bed`

And excludes:
- `CANDIDATE_reads_nonectopic_ARMS.bed` (no "segments_" in name)
- `CANDIDATE_reads_nonectopic_CEN.bed` (no "segments_" in name)

## Ready to Resume

```bash
nextflow run charla_nf/main.nf \
    --reads reads.fq.gz \
    ... other params \
    -resume
```

The VISUALIZE_SEGMENTS process will now receive the correct single BED file.

---

**Status**: Fixed ✅
**Date**: 2025-11-28
