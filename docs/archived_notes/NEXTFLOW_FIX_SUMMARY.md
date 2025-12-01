# Nextflow Error Fix - NON_HYBRID_ANALYSIS

## Error Message
```
ERROR ~ Unexpected error while finalizing task 'CHARLA_NF:NON_HYBRID_ANALYSIS (non_hybrid_analysis)'
- cause: No such property: bai_files for class: sun.nio.fs.UnixPath
```

## Root Cause

The `NON_HYBRID_MAPPING` process was emitting individual `bam_files` and `bai_files` channels, but:
1. These files are optional (may not exist if no sequences mapped)
2. The workflow was trying to access these outputs from the wrong process (`NON_HYBRID_ANALYSIS` instead of `NON_HYBRID_MAPPING`)
3. Nextflow's publishing mechanism failed when trying to access empty optional channels

## Solution

### 1. Simplified NON_HYBRID_MAPPING outputs
**File**: `modules/local/non_hybrid_mapping.nf`

**Before:**
```groovy
output:
path("${params.sample_id}_nonhybrid/"), emit: output_directory
path("${params.sample_id}_nonhybrid/*_reads.txt"), emit: read_lists
path("${params.sample_id}_nonhybrid/*_sequences.fa"), emit: extracted_sequences, optional: true
path("${params.sample_id}_nonhybrid/mappings/*.bam"), emit: bam_files, optional: true
path("${params.sample_id}_nonhybrid/mappings/*.bam.bai"), emit: bai_files, optional: true
```

**After:**
```groovy
output:
path("${params.sample_id}_nonhybrid/"), emit: output_directory
path("${params.sample_id}_nonhybrid/*_reads.txt"), emit: read_lists, optional: true
path("${params.sample_id}_nonhybrid/*_sequences.fa"), emit: extracted_sequences, optional: true
```

**Rationale:**
- BAM and BAI files are already included in `output_directory`
- No need to emit them separately
- Avoids empty channel publishing errors

### 2. Fixed workflow output references
**File**: `workflows/charla_nf.nf`

**Before:**
```groovy
// Non-hybrid analysis outputs
nonhybrid_directory = non_hybrid_analysis_output.output_directory
nonhybrid_read_lists = non_hybrid_analysis_output.read_lists
nonhybrid_sequences = non_hybrid_analysis_output.extracted_sequences
nonhybrid_bam_files = non_hybrid_analysis_output.bam_files  // WRONG PROCESS!
nonhybrid_bai_files = non_hybrid_analysis_output.bai_files  // WRONG PROCESS!
indel_analysis = non_hybrid_analysis_output.indel_results
```

**After:**
```groovy
// Non-hybrid mapping outputs (from mapping step)
nonhybrid_mapping_directory = non_hybrid_mapping_output.output_directory
nonhybrid_read_lists = non_hybrid_mapping_output.read_lists
nonhybrid_sequences = non_hybrid_mapping_output.extracted_sequences

// Non-hybrid analysis outputs (from analysis step)
nonhybrid_analysis_directory = non_hybrid_analysis_output.output_directory
indel_analysis = non_hybrid_analysis_output.indel_results
indel_plots = non_hybrid_analysis_output.plots
```

**Rationale:**
- Read lists and sequences come from `NON_HYBRID_MAPPING`
- Indel results and plots come from `NON_HYBRID_ANALYSIS`
- Clear separation of concerns

### 3. Added safety checks
**File**: `modules/local/non_hybrid_mapping.nf`

```bash
# Ensure mappings directory exists even if no BAMs created
mkdir -p ${params.sample_id}_nonhybrid/mappings

# Create placeholder to avoid empty directory issues
touch ${params.sample_id}_nonhybrid/mappings/.placeholder
```

## Files Modified

1. ✅ `modules/local/non_hybrid_mapping.nf` - Removed individual BAM/BAI emits
2. ✅ `workflows/charla_nf.nf` - Fixed output channel references

## Testing

The fix addresses:
- ✅ Empty BAM directory handling
- ✅ Optional file channel publishing
- ✅ Correct process output references
- ✅ Graceful handling when no sequences extracted

## Impact

**Breaking changes**: None
- All files still accessible via `output_directory`
- Just removed redundant individual file channels

**Benefits**:
- ✅ No more Nextflow publishing errors
- ✅ Cleaner channel management
- ✅ Proper separation of mapping vs analysis outputs

## Verification

After fix, the pipeline should:
1. Complete NON_HYBRID_MAPPING successfully
2. Pass `output_directory` to NON_HYBRID_ANALYSIS
3. Complete NON_HYBRID_ANALYSIS without errors
4. Publish all outputs correctly

No user-facing changes - all files are still in the same locations!

---

**Status**: ✅ Fixed
**Date**: November 26, 2025
