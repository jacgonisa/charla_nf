# Critical: Wrong Classification Source

## Problem Discovered

The confidence scoring and SCO/NCO analysis are using **WRONG classification**!

### Current (WRONG) Approach:
- Uses PAF mapping segments (`n_segments`, `n_switches` from alignment)
- Classifies based on mapping artifacts
- **Result**: Only 80 SCO, 2,788 NCO/Complex

### Correct Approach:
- Should use `hybrid_type` from **curated hybrid profiles**
- This is the biological classification from k-mer profiles
- **Result**: 124,209 SCO, 8,300 NCO, 136 GC-CO

## Source of Truth Files

### 1. Curated Hybrid Profiles (CLASSIFICATION)
**File**: `curated_hybrid_Col_Ler_F1_pollen_q10_1kb_ultracurated_cut.tsv`

**Columns**:
- `readname`: Read ID
- `threshold`: Threshold value (use 50)
- `ultra_curated_profile`: Actual segments (e.g., "325LA4,319LA3")
- `hybrid`: yes|N or no|N
- `hybrid_type`: **SCO, NCO, GC-CO, Unknown, None**
- `ectopic`: Ectopic recombination info

**Classification**:
- **SCO**: 124,209 reads (single crossover)
- **NCO**: 8,300 reads (non-crossover/gene conversion)
- **GC-CO**: 136 reads (gene conversion-crossover)
- **Unknown**: 26 reads
- **None**: 6,311,015 reads (parental, no crossover)

### 2. Segmentation BED (VISUALIZATION & GAPS)
**File**: `segments_*_filtered.bed`

**Purpose**:
- Visualize read structures
- **Show unassigned regions (gaps)** between segments
- Segment coordinates for plotting

**Format**:
```
read_id  start  end  segment_type
```

**Segment types**:
- **CA1-5**: Col ARMS chr1-5
- **LA1-5**: Ler ARMS chr1-5
- **CC1-5**: Col CEN chr1-5
- **LC1-5**: Ler CEN chr1-5

### 3. Parental Proportions (QUALITY METRICS)
**File**: `output_parental_proportions.tsv`

**Purpose**:
- Mapping quality (MAPQ)
- Alignment identity
- Col/Ler proportions
- Used for **confidence scoring quality metrics**

**NOT for classification** (n_segments/n_switches are from PAF mapping, not biology)

## What Each File Should Be Used For

| File | Classification | Visualization | Quality Metrics |
|------|---------------|---------------|-----------------|
| Curated Profiles | ✅ **YES** | ❌ No | ❌ No |
| Segmentation BED | ❌ No | ✅ **YES** | ❌ No |
| Parental Proportions | ❌ No | ❌ No | ✅ **YES** |

## Required Fixes

### 1. Confidence Scoring V3 (`scripts/20-calculate_confidence_v3.py`)
**Current**: Calculates `crossover_type` from PAF n_segments/n_switches
**Should**: Load `hybrid_type` from curated profiles file

**Changes needed**:
- Add input parameter: `--curated-profiles`
- Load curated profiles TSV
- Match reads by ID
- Use `hybrid_type` column directly
- Remove `identify_crossover_type()` function

### 2. SCO/NCO Analysis (`scripts/18-analyze_sco_nco_comparison.py`)
**Current**: Uses `crossover_type` from confidence scores (which is wrong)
**Should**: Still use `crossover_type` but AFTER confidence scoring is fixed

**Changes needed**:
- None (will automatically use correct classification once #1 is fixed)

### 3. Workflow Integration
**Current**: Confidence scoring doesn't receive curated profiles
**Should**: Pass curated profiles file to confidence scoring

**Changes needed in `workflows/charla_nf.nf`**:
```groovy
confidence_v3_output = CALCULATE_CONFIDENCE_V3(
    read_structure_output.proportions_table.first(),
    mapping_analysis_output.arms_summary.first(),
    mapping_analysis_output.cen_summary.first(),
    arms_bed_files.first(),
    cen_bed_files.first(),
    curate_profiles_output.curated_profiles.first(),  // ADD THIS
    params.threshold ?: 50
)
```

## Quick Fix for Current Run

To fix the double-prefix issue NOW:

**File**: `scripts/18-analyze_sco_nco_comparison.py:261`
**Change**: `f'{output_prefix}_sco_nco_comparison.png'`
**To**: `f'{output_prefix}_comparison.png'`

This will generate `sco_nco_comparison.png` instead of `sco_nco_sco_nco_comparison.png`.

## Summary

The pipeline has been using **mapping-based classification** (wrong) instead of **k-mer-based biological classification** (correct).

This explains why:
- Only 80 SCO found (should be 124,209)
- Most reads classified as "NCO/Complex" (wrong lumping)
- Numbers don't match biological expectation

**The curated profiles file is the SOURCE OF TRUTH for classification!**

---

**Priority**: HIGH - This affects all downstream analysis
**Impact**: Classification completely wrong
**Fix Required**: Update confidence scoring to use curated profiles
