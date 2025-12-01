# Classification Concordance Analysis Added

## Overview

Updated confidence scoring to use **TRUE classification** from curated k-mer profiles while also calculating **mapping-based classification** for concordance comparison.

## Changes Made

### 1. Script: `scripts/20-calculate_confidence_v3.py`

#### New Input Parameters:
```bash
--curated-profiles    # Curated hybrid profiles TSV (SOURCE OF TRUTH)
--threshold          # Threshold value to filter profiles
```

#### New Functions:
```python
normalize_hybrid_type(hybrid_type)
    # Convert curated profile types to standard categories
    # SCO → SCO
    # NCO, GC-CO → NCO
    # None → Parental
    # Unknown → Unknown

identify_crossover_type_mapping(n_segments, n_switches)
    # Renamed from identify_crossover_type
    # Calculates mapping-based classification for comparison
```

#### Classification Logic:
```python
# Get TRUE classification from k-mer profiles
crossover_type_true = normalize_hybrid_type(row['hybrid_type'])

# Calculate mapping-based classification
crossover_type_mapping = identify_crossover_type_mapping(n_segments, n_switches)

# Check concordance
concordant = (crossover_type_true == crossover_type_mapping) or \
            (crossover_type_true in ['NCO', 'Unknown'] and crossover_type_mapping == 'NCO/Complex')

# Use TRUE classification for scoring
crossover_type = crossover_type_true
```

#### New Output Columns:
```
crossover_type                # TRUE type from k-mer profiles (used for analysis)
crossover_type_mapping        # Mapping-based type (for comparison)
classification_concordant     # Boolean: do they agree?
```

### 2. Module: `modules/local/calculate_confidence_v3.nf`

#### Updated Input:
```groovy
input:
path parental_proportions
path curated_profiles        // NEW: SOURCE OF TRUTH
path arms_summary, stageAs: 'arms_summary.csv'
path cen_summary, stageAs: 'cen_summary.csv'
path arms_bed, stageAs: 'arms.bed'
path cen_bed, stageAs: 'cen.bed'
val threshold
```

#### Updated Script:
```bash
python ${projectDir}/scripts/20-calculate_confidence_v3.py \
    --parental-proportions ${parental_proportions} \
    --curated-profiles ${curated_profiles} \      # NEW
    --threshold ${threshold} \                     # NEW
    --summary-arms arms_summary.csv \
    --summary-cen cen_summary.csv \
    --output ${threshold_dir}/confidence_scores.tsv \
    --diagnostic-dir ${threshold_dir}/diagnostics
```

### 3. Workflow: `workflows/charla_nf.nf`

#### Updated Call:
```groovy
confidence_v3_output = CALCULATE_CONFIDENCE_V3(
    read_structure_output.proportions_table.first(),
    curate_profiles_output.curated_profiles.first(),  // NEW: SOURCE OF TRUTH
    mapping_analysis_output.arms_summary.first(),
    mapping_analysis_output.cen_summary.first(),
    arms_bed_files.first(),
    cen_bed_files.first(),
    params.threshold ?: 50
)
```

## Classification Sources

### K-mer Profiles (TRUE Classification):
- **Source**: `curated_hybrid_*ultracurated_cut.tsv`
- **Column**: `hybrid_type`
- **Values**: SCO, NCO, GC-CO, None, Unknown
- **Based on**: K-mer composition analysis
- **Reliability**: HIGH (biological ground truth)

### Mapping-Based (Comparison):
- **Source**: PAF alignments
- **Calculation**: Based on `n_segments` and `n_switches`
- **Logic**: 2 segments + 1 switch = SCO, else = NCO/Complex
- **Based on**: Alignment patterns
- **Reliability**: LOWER (artifacts from mapping)

## Expected Classifications

Based on your data at threshold 50:

### TRUE (K-mer profiles):
- **SCO**: ~124,000 reads
- **NCO**: ~8,300 reads
- **GC-CO**: ~136 reads
- **Parental**: ~6,311,000 reads

### Mapping-based (old method):
- **SCO**: ~80 reads (only n_segments=2, n_switches=1)
- **NCO/Complex**: ~2,788 reads (everything else)

### Concordance Analysis:
The output will now show:
1. How many reads are correctly classified by mapping
2. Which reads are misclassified
3. Whether mapping quality correlates with concordance

## Output Example

```tsv
read_id                 confidence  crossover_type  crossover_type_mapping  classification_concordant
00137367-5598-4f5a...   0.94       SCO              NCO/Complex             False
001824b6-6569-4ed7...   0.98       SCO              NCO/Complex             False
04cc321a-7d29-47a6...   0.85       SCO              SCO                     True
```

## Analysis Possibilities

With this data you can now:

1. **Check concordance rate**:
   ```bash
   awk -F'\t' 'NR>1 {if($11=="True") conc++; total++} END {print conc/total*100"%"}' confidence_scores.tsv
   ```

2. **Find discordant high-confidence reads**:
   ```bash
   awk -F'\t' 'NR>1 && $2>=0.8 && $11=="False"' confidence_scores.tsv
   ```

3. **Compare by MAPQ**:
   - Do low MAPQ reads have more discordance?
   - Are certain mappers more concordant?

4. **Identify misclassification patterns**:
   - SCO called as NCO/Complex by mapping
   - NCO called as SCO by mapping

## Biological Insight

The **curated k-mer profiles are the SOURCE OF TRUTH** because:
- K-mers represent actual biological sequence
- Not affected by mapping artifacts
- Capture true recombination patterns
- Include ARMS + CEN regions comprehensively

The **mapping-based classification** is useful for:
- Quality control (checking mapping reliability)
- Understanding which reads have complex mapping patterns
- Identifying potential mapping issues

## Usage

No changes needed to run commands:
```bash
nextflow run charla_nf/main.nf ... -resume
```

The output `confidence_scores.tsv` will now have:
- TRUE classification (from k-mer profiles)
- Mapping classification (for comparison)
- Concordance flag

## Compatibility

- ✅ Backward compatible (same output filename)
- ✅ SCO/NCO analysis uses TRUE classification
- ✅ Additional columns for comparison
- ✅ All downstream analyses use correct classification

---

**Status**: Complete ✅
**Date**: 2025-11-28
**Impact**: Now uses correct biological classification!
