# PAF File Loading Fix

## Problem Identified

The `VISUALIZE_READ_STRUCTURES` process was failing with:
```
[WARNING] No PAF records loaded from paf_files
[ERROR] No PAF data loaded!
```

## Root Cause

The workflow was passing **summary PAF files** (5 columns) instead of **full PAF files** (12+ columns with MAPQ data):

### Summary PAF (5 columns - NO MAPQ):
```
7ffd019e-57a9-4faf-bb55-a5a0be8ccc70_1_CA5	Chr5_ARMS_Col_2	14992485	4988551	4989992
```

### Full PAF (12+ columns - WITH MAPQ):
```
7ffd019e-57a9-4faf-bb55-a5a0be8ccc70_1_CA5	1447	0	1447	+	Chr5_ARMS_Col_2	14992485	4988551	4989992	1438	1447	60	NM:i:9	ms:i:2843	AS:i:2840	...
                                                                                                                                                        ^^
                                                                                                                                                      MAPQ=60!
```

The summary PAF files in `plotting_crossover/${sample_id}/*.paf` are simplified files without:
- MAPQ (mapping quality)
- CIGAR strings
- Alignment identity
- NM, AS tags

These are needed for:
1. Read structure visualization
2. Confidence scoring V3
3. Parental proportion analysis

## Solution

### 1. Added New Output to MAPPING_ANALYSIS Module

**File:** `modules/local/mapping_analysis.nf` (line 28)

```groovy
output:
    // ... existing outputs ...
    path "threshold${threshold}/mapping_full_read/${params.sample_id}/nonectopic/*/summary_*/best_alignments.paf", emit: full_paf_files
```

This glob pattern captures all `best_alignments.paf` files from each mapper:
- `summary_mm_sr/best_alignments.paf` (minimap2 short-read mode)
- `summary_mm_ont/best_alignments.paf` (minimap2 ONT mode)
- `summary_lrhq/best_alignments.paf` (LRA high-quality)
- `summary_lrhqae/best_alignments.paf` (LRA high-quality all-vs-each)
- `summary_wm/best_alignments.paf` (Winnowmap)

Each file contains full PAF format with MAPQ, identity, and all alignment tags.

### 2. Updated Workflow to Use Full PAF Files

**File:** `workflows/charla_nf.nf` (line 470)

**Before:**
```groovy
read_structure_output = VISUALIZE_READ_STRUCTURES(
    mapping_analysis_output.paf_files.collect(),   // Summary PAFs (5 columns)
    ...
)
```

**After:**
```groovy
read_structure_output = VISUALIZE_READ_STRUCTURES(
    mapping_analysis_output.full_paf_files.collect(),   // Full PAFs (12+ columns)
    ...
)
```

### 3. Added Missing Parameters

**File:** `nextflow.config` (lines 96-107)

```groovy
// Non-hybrid analysis resources
nonhybrid_cpus = 8
nonhybrid_memory = '16 GB'
nonhybrid_time = '4h'

analysis_cpus = 4
analysis_memory = '8 GB'
analysis_time = '2h'

// General parameters
threads = 8                 // Default number of threads
min_quality = 10            // Minimum base quality for quality masking
```

## What's Different Now

### Old (Summary PAF):
- Only 5 fields per line
- No MAPQ information
- No alignment details
- No identity calculation possible

### New (Full PAF):
- 12+ fields per line
- MAPQ in column 12 (e.g., 60 = 1 in 1M error rate)
- NM tag = edit distance
- AS tag = alignment score
- Full CIGAR string
- Can calculate identity: `(matches / alignment_length)`

## Files Collected

The workflow will now collect ~10 full PAF files:
- ARMS_segments: 5 mappers × 1 best_alignments.paf = 5 files
- CEN_segments: 5 mappers × 1 best_alignments.paf = 5 files

Example paths:
```
threshold50/mapping_full_read/Col_Ler_F1_pollen_q10_1kb/nonectopic/
├── ARMS_segments_alntoCol_alntoLer/
│   ├── summary_mm_sr/best_alignments.paf       ← Full PAF
│   ├── summary_mm_ont/best_alignments.paf      ← Full PAF
│   ├── summary_lrhq/best_alignments.paf        ← Full PAF
│   ├── summary_lrhqae/best_alignments.paf      ← Full PAF
│   └── summary_wm/best_alignments.paf          ← Full PAF
└── CEN_segments_alntoCol_alntoLer/
    ├── summary_mm_sr/best_alignments.paf       ← Full PAF
    ├── summary_mm_ont/best_alignments.paf      ← Full PAF
    ├── summary_lrhq/best_alignments.paf        ← Full PAF
    ├── summary_lrhqae/best_alignments.paf      ← Full PAF
    └── summary_wm/best_alignments.paf          ← Full PAF
```

## What Will Work Now

With full PAF files, the pipeline can now:

1. **Extract MAPQ** - Mapping quality per alignment
2. **Calculate identity** - From CIGAR and NM tag
3. **Parental proportions** - Col vs Ler balance per read
4. **Recombination rate** - Events per Mb
5. **Confidence V3 scoring** - Using real MAPQ and identity metrics

## Expected Behavior

When you rerun with `-resume`, the `VISUALIZE_READ_STRUCTURES` process should:

```
=== Analyzing all crossover reads ===
[INFO] Loading PAF files...
[INFO] Found 10 PAF files
[INFO] Loading: best_alignments.paf (ARMS, mm_sr)
[INFO] Loading: best_alignments.paf (ARMS, mm_ont)
...
[INFO] Loaded 1,131 PAF records
[INFO] Mean MAPQ: 58.3
[INFO] Mean Identity: 0.9987
[SUCCESS] Recombination rate: 208.28 events/Mb
```

Instead of:
```
[WARNING] No PAF records loaded from paf_files
[ERROR] No PAF data loaded!
```

## Files Modified

1. ✅ `modules/local/mapping_analysis.nf` - Added `full_paf_files` output
2. ✅ `workflows/charla_nf.nf` - Changed to use `full_paf_files` instead of `paf_files`
3. ✅ `nextflow.config` - Added missing parameters (`nonhybrid_cpus`, `threads`, `min_quality`, etc.)

## Ready to Rerun

```bash
nextflow run charla_nf/main.nf ... -resume
```

The MAPPING_ANALYSIS process is cached, so the full PAF files already exist. The workflow will just use the new `full_paf_files` output channel instead of the summary `paf_files` channel.
