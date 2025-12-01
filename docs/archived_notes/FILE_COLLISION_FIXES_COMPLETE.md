# File Name Collision Fixes - Complete

## Problem

Multiple Nextflow processes were failing with input file name collision errors:

```
ERROR ~ Process input file name collision --
There are multiple input files for each of the following file names:
best_alignments.paf, summary_table.csv, etc.
```

## Root Cause

When processes receive multiple path inputs from different directories that contain files with identical names, Nextflow stages all files into the same work directory, causing collisions.

## Solution Pattern

Use `stageAs` parameter to rename files during staging, ensuring unique names:

```groovy
input:
path file1, stageAs: 'unique_name1.ext'
path file2, stageAs: 'unique_name2.ext'
```

Then update the script to reference the staged names.

---

## Fixed Modules

### 1. VISUALIZE_READ_STRUCTURES ✅

**File**: `modules/local/visualize_read_structures.nf`

**Issue**: ARMS and CEN mapping directories both contain:
- `best_alignments.paf`
- `summary_table.csv`
- `inversions.paf`
- etc.

**Fix Applied** (lines 13-14):
```groovy
input:
path arms_mapping_dir, stageAs: 'arms_mapping/**'
path cen_mapping_dir, stageAs: 'cen_mapping/**'
path confidence_scores
val max_reads
```

**Script Updated** (lines 38-62):
```bash
# Find best_alignments.paf from ARMS (stageAs flattens to numbered symlinks)
for paf in arms_mapping/*; do
    if [[ -L "$paf" ]]; then
        target=$(readlink "$paf")
        if [[ "$target" == *"best_alignments.paf" ]]; then
            mapper=$(echo "$target" | grep -oP 'summary_[^/]+' | head -1)
            cp "$paf" paf_files/ARMS_${mapper}.paf
        fi
    fi
done

# Same for CEN
for paf in cen_mapping/*; do
    ...
done
```

**Result**: Successfully loaded 19,478 alignments, generated all plots.

---

### 2. CALCULATE_CONFIDENCE_V3 ✅

**File**: `modules/local/calculate_confidence_v3.nf`

**Issue**: Both ARMS and CEN summary files named:
- `summary_table.csv`

**Fix Applied** (lines 13-17):
```groovy
input:
path parental_proportions
path arms_summary, stageAs: 'arms_summary.csv'
path cen_summary, stageAs: 'cen_summary.csv'
path arms_bed, stageAs: 'arms.bed'
path cen_bed, stageAs: 'cen.bed'
val threshold
```

**Script Updated** (lines 37-38):
```bash
python ${projectDir}/scripts/20-calculate_confidence_v3.py \
    --parental-proportions ${parental_proportions} \
    --summary-arms arms_summary.csv \
    --summary-cen cen_summary.csv \
    --output ${threshold_dir}/confidence_scores.tsv \
    --diagnostic-dir ${threshold_dir}/diagnostics
```

---

### 3. ANALYZE_SCO_NCO ✅

**File**: `modules/local/analyze_sco_nco.nf`

**Issue**: Potential collision with summary and cowidth files

**Fix Applied** (lines 12-16):
```groovy
input:
path confidence_scores                          // Confidence scores TSV from V3
path raw_metrics, stageAs: 'raw_metrics.csv'    // Raw metrics (ARMS summary)
path arms_cowidth, stageAs: 'arms.cowidths'     // ARMS cowidth file
path cen_cowidth, stageAs: 'cen.cowidths'       // CEN cowidth file
```

**Script Updated** (lines 25-33):
```bash
# Combine cowidth files
cat arms.cowidths > combined_cowidths.txt
tail -n +2 cen.cowidths >> combined_cowidths.txt

# Run SCO vs NCO analysis
python ${projectDir}/scripts/18-analyze_sco_nco_comparison.py \
    --confidence ${confidence_scores} \
    --raw-metrics raw_metrics.csv \
    --cowidth combined_cowidths.txt \
    --output-prefix sco_nco
```

---

## Technical Details

### How `stageAs` Works

1. **Without `stageAs`**:
   ```
   work/XX/YYYY/
   ├── summary_table.csv    (from ARMS)
   └── summary_table.csv    (from CEN) → COLLISION!
   ```

2. **With `stageAs`**:
   ```
   work/XX/YYYY/
   ├── arms_summary.csv    (from ARMS)
   └── cen_summary.csv     (from CEN) → No collision!
   ```

### `stageAs: '**'` with Directories

When using `stageAs: 'dirname/**'` with a directory input:
- Nextflow creates numbered symlinks (1, 2, 3, etc.)
- Each symlink points to original file location
- Use `readlink` to inspect original path
- Extract information from original path (e.g., mapper name)

Example from VISUALIZE_READ_STRUCTURES:
```bash
# After staging with stageAs: 'arms_mapping/**'
ls arms_mapping/
# Output: 1, 2, 3, 4, 5, ...

readlink arms_mapping/1
# Output: /path/to/original/summary_mm_sr/best_alignments.paf

# Extract mapper name from original path
mapper=$(echo "$target" | grep -oP 'summary_[^/]+')
```

---

## Best Practices

1. **Defensive Staging**: Even if files have different names, use `stageAs` for clarity
2. **Descriptive Names**: Use meaningful staged names that indicate source/purpose
3. **Update Scripts**: Always update script section to use staged names
4. **Document Changes**: Add comments explaining the staging strategy

---

## Files Modified

### Nextflow Modules:
- `modules/local/visualize_read_structures.nf`
- `modules/local/calculate_confidence_v3.nf`
- `modules/local/analyze_sco_nco.nf`

### No Changes Needed:
- Workflow files (channels pass correctly)
- Python scripts (receive correct filenames)

---

## Testing

To verify fixes work with cache:
```bash
nextflow run charla_nf/main.nf \
    --reads reads.fq.gz \
    --sample_id sample \
    --cenhapmer_db_dir cenhapmers/k31 \
    --reference_genomes_dir genomes/ \
    --outdir results/ \
    --col_bed_file col.bed \
    --ler_bed_file ler.bed \
    -profile singularity \
    -resume
```

All three modules should now:
1. Stage files without collisions
2. Execute successfully
3. Generate expected outputs

---

## Related Issues

- **STAGING_FIX.md**: Initial fix for VISUALIZE_READ_STRUCTURES
- **CLEANUP_SUMMARY.md**: Removal of V2, optional indel analysis
- **INDEL_MAPPING_OPTIMIZATION.md**: BAM file batch reading

---

**Status**: All collision fixes complete ✅
**Date**: 2025-11-28
**Ready for**: Pipeline execution with `-resume`
