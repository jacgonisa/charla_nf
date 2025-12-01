# CHARLA Pipeline Development Notes

This document contains important notes about recent development work, bug fixes, and improvements to the CHARLA pipeline.

For detailed historical session notes, see [`docs/archived_notes/`](docs/archived_notes/).

---

## Table of Contents

1. [Recent Major Updates](#recent-major-updates)
2. [Critical Bug Fixes](#critical-bug-fixes)
3. [Analysis Improvements](#analysis-improvements)
4. [Known Issues & Troubleshooting](#known-issues--troubleshooting)

---

## Recent Major Updates

### Confidence Scoring V3 (Redesigned)

**File**: `scripts/20-calculate_confidence_v3.py`

Completely redesigned confidence scoring system with weighted components:

| Component | Weight | What it measures |
|-----------|--------|------------------|
| MAPQ Score | 30% | Mapping quality (Phred-scaled) |
| Identity Score | 25% | Alignment identity % |
| Consensus Score | 20% | Agreement across 5 mappers |
| Balance Score | 15% | Parental balance (50/50 for SCO) |
| Coverage Score | 10% | Fraction of read assigned (not grey) |

**Key Features**:
- Uses **curated hybrid profiles** as source of truth for classification
- Outputs **three classifications** for concordance analysis:
  - `crossover_type`: AFTER mapping (validated by mappers) ← **USE THIS**
  - `crossover_type_before`: BEFORE mapping (k-mer only)
  - `crossover_type_mapping`: PAF-based (for comparison)
- Generates diagnostic plots automatically

### SCO vs NCO Comparison Analysis

**Files**:
- `scripts/18-analyze_sco_nco_comparison.py`
- `modules/local/analyze_sco_nco.nf`

Comprehensive comparison of Single Crossovers (SCO) vs Non-Crossovers (NCO):

**Outputs**:
- `sco_nco_comparison.png`: 9-panel comparison (counts, confidence, segments, regions)
- `sco_nco_resolution_analysis.png`: 4-panel resolution analysis (width distribution, trends)
- `sco_nco_summary.tsv`: Statistical summary table

**Key Metrics**:
- Confidence score distributions
- Crossover width/resolution
- Segment length analysis
- Mapping consensus
- ARMS vs CEN distribution

### BED Segment Analysis

**Update**: SCO/NCO analysis now extracts segment information directly from BED files.

**Why**: The PAF `n_segments` counts alignment blocks (including gaps/splits), not parental switches. BED segments give true classification segments.

**New Metrics**:
- `n_bed_segments`: Total BED entries per read
- `n_classification_segments`: Collapsed parental segments (2=SCO, 3=NCO, 4=DCO, 5+=Complex)
- `segment_pattern`: Parental pattern string ("CL", "CLC", etc.)

---

## Critical Bug Fixes

### 1. Classification Extraction Bug ⚠️ **CRITICAL**

**Problem**: Classification function was extracting the wrong character from segment strings.

**Segment format**: `<number><parental><region><chr>` (e.g., "325LA4" = 325bp + Ler + ARMS + Chr4)

**Bug**: Took first character `seg[0]` which is a digit, not the parental type!
```python
# WRONG: "325LA4"[0] = "3" → "X" → Complex
parental_pattern = ''.join([seg.strip()[0] if seg.strip()[0] in ['C', 'L'] ...
```

**Fix**: Scan for first 'C' or 'L' character:
```python
# CORRECT: "325LA4" → scan → find "L" → "L"
for char in seg:
    if char in ['C', 'L']:
        parental_types.append(char)
        break
```

**Impact**: This bug caused ALL reads to be classified as "Complex". After fix, correct SCO/NCO/DCO classifications work.

**Test**: `test_classification.py` validates all classification logic.

### 2. Classification Logic Update

**User clarification**: "3 segments is NCO ok? complex is more than 3"

**Updated logic**:
- 1 segment (C or L) → **Parental**
- 2 segments (CL or LC) → **SCO** (Single Crossover)
- 3 segments (CLC or LCL) → **NCO** (Non-Crossover / Gene Conversion)
- 4 segments (CLCL or LCLC) → **DCO** (Double Crossover)
- 5+ segments → **Complex**

### 3. Numerical Stability in Plotting

**Problem**: `log10(0)` errors and SVD convergence failures when plotting crossover widths.

**Fix**:
- Filter out zero/negative widths before log scale
- Add try/except for trend line fitting
- Require minimum 10 data points for meaningful plots

```python
widths = widths[widths > 0]  # Remove zeros
if len(widths) > 10:
    # Safe to plot
```

### 4. File Staging Collisions

**Problem**: Multiple input files with same names caused Nextflow staging collisions.

**Fix**: Use `stageAs` parameter to rename during staging:
```groovy
path arms_summary, stageAs: 'arms_summary.csv'
path cen_summary, stageAs: 'cen_summary.csv'
```

**Affected modules**:
- VISUALIZE_READ_STRUCTURES
- CALCULATE_CONFIDENCE_V3
- ANALYZE_SCO_NCO

---

## Analysis Improvements

### Marker Balance Analysis

**File**: `scripts/16-analyze_marker_balance.py`

Analyzes k-mer marker density across chromosomes and regions:
- Normalizes by genomic region size
- Compares ARMS vs CEN regions
- Helps choose optimal k-mer size (k=21, k=31, k=41)

**Key insight**: Higher k produces MORE unique markers but less error-tolerant.

### Indel Analysis Enhancements

**Files**:
- `scripts/12.2-indel_analysis.py`
- `scripts/13.2-map_indels_to_genome.py`
- `scripts/15-analyze_indel_monomer_boundaries.py`

**Improvements**:
- Optimized batch processing (10x speedup)
- Added monomer boundary analysis
- Enhanced filtering and classification
- Better visualization

### Read Structure Visualization

**File**: `scripts/19-visualize_read_structures.py`

Visualizes read structures from PAF alignments:
- Shows Col/Ler segments as colored bars
- Includes MAPQ quality scores
- Highlights crossover regions
- Displays unassigned regions (gaps)

### Segment Visualization from BED

**File**: `scripts/21-visualize_read_segments.py`

Visualizes segments directly from BED files:
- Shows actual segmented regions from mapping
- Highlights unassigned regions (gaps between segments)
- More accurate than PAF-based visualization

---

## Known Issues & Troubleshooting

### Issue: Pipeline uses cached results after script changes

**Symptom**: Modified Python scripts but pipeline still uses old results with `-resume`.

**Cause**: Nextflow doesn't detect changes in external scripts (only module files).

**Solution**: Manually delete work directories for affected processes:
```bash
# Find and delete specific process cache
rm -rf /mnt/ssd-8tb/your_run/work/77/5ffa5700...

# Or use nextflow clean
nextflow clean -f -k
```

### Issue: scipy ImportError

**Symptom**: `ModuleNotFoundError: No module named 'scipy'`

**Cause**: Container doesn't have scipy installed.

**Solution**: Made scipy optional in analysis scripts. Statistical tests will be skipped but plots still generate.

### Issue: Zero SCO reads detected

**Symptom**: Analysis shows 0 SCO, all reads classified as Complex.

**Cause**: Classification bug (see Critical Bug Fixes #1 above).

**Solution**: Update to fixed version of `scripts/20-calculate_confidence_v3.py`.

### Issue: Confusing segment counts in plots

**Symptom**: Plot says "SCO should be 2" but shows wide distribution.

**Cause**: Using PAF `n_segments` (alignment blocks) instead of classification segments.

**Solution**: Updated to extract segments from BED file (see BED Segment Analysis above).

---

## Cache Management Tips

When making script changes, you may need to invalidate caches:

```bash
# Delete specific process caches
find /path/to/work -name "confidence_scores.tsv" -printf "%h\n" | xargs rm -rf

# Clean all failed/completed tasks
nextflow clean -f -k

# Clean specific process by name (from .nextflow.log)
rm -rf /path/to/work/XX/XXXXXXXX...
```

---

## Performance Notes

### Batch Processing

Several scripts use batch processing for large datasets:
- Indel analysis: 10,000 reads per batch
- BAM reading: Batch iteration instead of loading all
- Result: 100-3000x speedup

### Parallelization

Pipeline supports parallel execution:
- Multiple mappers run in parallel
- Independent analyses run concurrently
- Configurable CPU/memory per process

---

## Testing

### Classification Tests

**File**: `test_classification.py`

Tests all classification logic with real segment patterns:
- Parental (single segment)
- SCO (2 segments)
- NCO (3 segments)
- DCO (4 segments)
- Complex (5+ segments)

Run with: `python test_classification.py`

---

## For More Details

See archived session notes in [`docs/archived_notes/`](docs/archived_notes/) for:
- Detailed debugging sessions
- Step-by-step fix descriptions
- Historical development notes
- Performance optimization details

---

**Last Updated**: 2025-11-28
**Pipeline Version**: v1.0.1
**Nextflow Version**: 23.10.0+
