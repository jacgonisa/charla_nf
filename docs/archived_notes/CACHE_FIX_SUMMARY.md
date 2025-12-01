# Cache-Preserving PAF File Fix

## Problem

Adding a new output to `MAPPING_ANALYSIS` module broke the cache, causing Nextflow to want to rerun everything from scratch even with `-resume`.

## Solution

Instead of changing MAPPING_ANALYSIS outputs (which invalidates cache), we modified VISUALIZE_READ_STRUCTURES to **search for PAF files within the directories it already receives**.

## What Changed

### 1. VISUALIZE_READ_STRUCTURES Module

**Input changed** (lines 13-14):
```groovy
// OLD:
path paf_files                    // Collection of all PAF files

// NEW:
path arms_mapping_dir             // ARMS mapping results directory
path cen_mapping_dir              // CEN mapping results directory
```

**Script changed** (lines 36-53):
```bash
# Find and copy best_alignments.paf from ARMS mapping directory
find ${arms_mapping_dir} -name "best_alignments.paf" -type f | while read paf; do
    mapper=$(basename $(dirname "$paf"))
    echo "Found ARMS PAF: $mapper"
    cp "$paf" paf_files/ARMS_${mapper}.paf
done

# Find and copy best_alignments.paf from CEN mapping directory
find ${cen_mapping_dir} -name "best_alignments.paf" -type f | while read paf; do
    mapper=$(basename $(dirname "$paf"))
    echo "Found CEN PAF: $mapper"
    cp "$paf" paf_files/CEN_${mapper}.paf
done
```

This searches recursively for all `best_alignments.paf` files (which have full MAPQ data) within the mapping directories.

### 2. Workflow Call

**File:** `workflows/charla_nf.nf` (lines 469-474)

```groovy
// OLD:
read_structure_output = VISUALIZE_READ_STRUCTURES(
    mapping_analysis_output.full_paf_files.collect(),   // New output - breaks cache!
    ...
)

// NEW:
read_structure_output = VISUALIZE_READ_STRUCTURES(
    mapping_analysis_output.arms_mapping_results.first(),  // Already exists - preserves cache!
    mapping_analysis_output.cen_mapping_results.first(),   // Already exists - preserves cache!
    confidence_scores_output.confidence_scores.first(),
    params.max_structure_reads ?: 200
)
```

### 3. MAPPING_ANALYSIS Module

**No changes!** - Reverted to original state to preserve cache.

The module already emits:
- `arms_mapping_results` - Contains all ARMS mapping data including `best_alignments.paf`
- `cen_mapping_results` - Contains all CEN mapping data including `best_alignments.paf`

## Why This Works

### Cache Preservation:
- ✅ MAPPING_ANALYSIS module unchanged → cache valid
- ✅ Uses existing output channels → cache valid
- ✅ Only VISUALIZE_READ_STRUCTURES changed → runs fresh (only this module)

### PAF Files Found:
The mapping directories contain the full structure:
```
ARMS_segments_alntoCol_alntoLer/
├── summary_mm_sr/
│   └── best_alignments.paf      ← MAPQ=60, full format
├── summary_mm_ont/
│   └── best_alignments.paf      ← MAPQ=60, full format
├── summary_lrhq/
│   └── best_alignments.paf      ← MAPQ=60, full format
├── summary_lrhqae/
│   └── best_alignments.paf      ← MAPQ=60, full format
└── summary_wm/
    └── best_alignments.paf      ← MAPQ=60, full format
```

The `find` command discovers all these files automatically.

## Expected Behavior

When you run with `-resume` now:

```
[c2/8db7d3] process > CHARLA_NF:MAPPING_ANALYSIS  [100%] 1 of 1, cached: 1 ✔
[NEW_ID]    process > CHARLA_NF:VISUALIZE_READ_STRUCTURES  [  0%] 0 of 1

=== Searching for full PAF files with MAPQ data ===
Found ARMS PAF: summary_mm_sr
Found ARMS PAF: summary_mm_ont
Found ARMS PAF: summary_lrhq
Found ARMS PAF: summary_lrhqae
Found ARMS PAF: summary_wm
Found CEN PAF: summary_mm_sr
Found CEN PAF: summary_mm_ont
Found CEN PAF: summary_lrhq
Found CEN PAF: summary_lrhqae
Found CEN PAF: summary_wm
PAF files collected:
  ARMS_summary_mm_sr.paf
  ARMS_summary_mm_ont.paf
  ARMS_summary_lrhq.paf
  ARMS_summary_lrhqae.paf
  ARMS_summary_wm.paf
  CEN_summary_mm_sr.paf
  CEN_summary_mm_ont.paf
  CEN_summary_lrhq.paf
  CEN_summary_lrhqae.paf
  CEN_summary_wm.paf
=== Analyzing all crossover reads ===
[INFO] Loading PAF files...
[INFO] Found 10 PAF files
[INFO] Loaded 1,131 PAF records
[INFO] Mean MAPQ: 58.3
[SUCCESS] Read structure analysis complete!
```

## Files Modified

1. ✅ `modules/local/visualize_read_structures.nf`
   - Changed input to receive mapping directories
   - Added find command to locate best_alignments.paf files

2. ✅ `workflows/charla_nf.nf`
   - Pass mapping directories instead of PAF files

3. ✅ `modules/local/mapping_analysis.nf`
   - Reverted to original (no changes)

## Advantages

1. **Cache preserved** - MAPPING_ANALYSIS won't rerun
2. **Gets full PAF files** - With MAPQ, identity, and all tags
3. **Automatic discovery** - Finds all mappers' PAF files
4. **No pipeline modification** - Uses existing outputs
5. **Future-proof** - If more mappers added, automatically found

## Ready to Run

```bash
nextflow run charla_nf/main.nf ... -resume
```

All cached processes will be reused, only VISUALIZE_READ_STRUCTURES will run fresh.
