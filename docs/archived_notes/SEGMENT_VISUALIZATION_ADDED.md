# Segment Visualization Module Added

## Overview

Added a new visualization module that uses segmentation BED files to show read structures with clear indication of **unassigned regions** (gaps between segments).

## Why This Is Better

### Previous Approach (PAF-based):
- Uses alignment data (PAF files)
- Shows mapping quality and parental assignments
- Less clear about unassigned regions

### New Approach (BED-based):
- Uses segmentation output directly
- **Clearly shows unassigned regions** (gaps between segments) with red dashed borders
- Shows exact segment boundaries and types (CA1-5, LA1-5, LC1-5)
- Better for understanding crossover resolution

## What It Shows

### Read Structure Plot:
- **Assigned segments**: Colored blocks
  - Col haplotypes (CA1-5): Blue shades
  - Ler haplotypes (LA1-5): Orange shades
  - Complex (LC1-5): Gray shades
- **Unassigned regions**: Light gray with red dashed border
- Read labels show: read ID, length, type (SCO/NCO/PARENTAL/COMPLEX)

### Statistics Plot (4 panels):
1. **Read type distribution**: Pie chart (SCO, NCO, Parental, Complex)
2. **Segments per read**: Histogram by type
3. **Read length distribution**: By type
4. **Segment length distribution**: Col vs Ler (log scale)

### Summary File:
- Total reads and segments
- Read type breakdown with percentages
- Average metrics (segments per read, read/segment lengths)

## Files Created

### Script:
- `scripts/21-visualize_read_segments.py` - Main visualization script
  - Loads segmentation BED file
  - Identifies unassigned regions (gaps)
  - Classifies reads (SCO/NCO/PARENTAL/COMPLEX)
  - Generates plots and statistics

### Module:
- `modules/local/visualize_segments.nf` - Nextflow process
  - Takes segmentation BED file
  - Optional confidence scores for filtering
  - Generates 3 sets of plots: all reads, SCO only, NCO only

### Workflow Integration:
- `workflows/charla_nf.nf` - Added Step 4 after SCO/NCO analysis
  - Uses `segment_reads_output.bed_files`
  - Uses `confidence_v3_output.confidence_scores` for filtering
  - Runs after confidence scoring

## Input Format

Segmentation BED file (from SEGMENT_READS):
```
read_id    start    end    segment_type
06036bce-...    2045    6881    LA4
06036bce-...    7235    9051    LA3
```

Where:
- `segment_type`: CA1-5 (Col ARMS chr1-5), LA1-5 (Ler ARMS), LC1-5 (Complex)
- Gaps between segments = **unassigned regions**

## Output Files

```
segment_visualization/
├── all_reads_read_structures.png    # All reads structure plot
├── all_reads_statistics.png         # Summary statistics
├── all_reads_summary.txt            # Text summary
├── SCO_read_structures.png          # SCO-only structures
├── SCO_statistics.png               # SCO statistics
├── SCO_summary.txt                  # SCO summary
├── NCO_read_structures.png          # NCO-only structures
├── NCO_statistics.png               # NCO statistics
└── NCO_summary.txt                  # NCO summary
```

## Read Classification Logic

The script automatically classifies reads:

- **PARENTAL**: Only one haplotype (no crossover)
- **SCO**: Exactly 2 segments, 2 different haplotypes (single crossover)
- **NCO**: Multiple segments with small segments (<500bp) suggesting gene conversion
- **COMPLEX**: Multiple crossovers or unclear patterns

## Usage in Pipeline

Automatically runs as Step 4 after confidence scoring:

```bash
nextflow run charla_nf/main.nf \
    --reads reads.fq.gz \
    ... other params \
    -resume
```

Output will appear in `${params.outdir}/segment_visualization/`

## Key Features

1. **Unassigned regions clearly marked**:
   - Light gray fill
   - Red dashed borders
   - Easy to spot gaps in coverage

2. **Automatic classification**:
   - No manual inspection needed
   - Consistent criteria across all reads

3. **Multiple views**:
   - All reads
   - High-confidence SCO only
   - High-confidence NCO only

4. **Confidence filtering**:
   - Can filter for score ≥ 0.8
   - Ensures high-quality visualizations

5. **Comprehensive statistics**:
   - Distributions by read type
   - Segment length analysis
   - Summary metrics

## Fixes Applied

### scipy Dependency:
Added to ANALYZE_SCO_NCO module:
```bash
pip install scipy --quiet 2>/dev/null || true
```

This fixes the `ModuleNotFoundError: No module named 'scipy'` error.

## Example Output

### Read Structure:
```
[06036bce...]  (8.9kb, SCO)
|==LA4==| GAP |==LA3==|
  Blue     Gray   Blue
          ^^^^^^
       Unassigned!
```

### Statistics:
- SCO: 45% (2,340 reads)
- NCO: 12% (623 reads)
- PARENTAL: 38% (1,974 reads)
- COMPLEX: 5% (260 reads)

## Testing

1. Resume pipeline:
   ```bash
   nextflow run charla_nf/main.nf ... -resume
   ```

2. Check output:
   ```bash
   ls ${outdir}/segment_visualization/
   ```

3. Inspect plots:
   - Look for red dashed regions (unassigned)
   - Compare SCO vs NCO patterns
   - Review statistics

## Comparison with PAF Visualization

| Feature | PAF-based | BED-based (new) |
|---------|-----------|-----------------|
| Source | Alignment data | Segmentation output |
| Unassigned regions | Unclear | **Clearly marked** |
| Segment boundaries | Approximate | Exact |
| MAPQ quality | ✅ Shown | ❌ Not shown |
| Parental balance | ✅ Shown | ❌ Not shown |
| Read classification | Manual | **Automatic** |
| Best for | Quality assessment | **Structure analysis** |

## Recommendation

Use **both** visualizations:
1. **PAF-based** (`read_structure_analysis/`): For quality metrics
2. **BED-based** (`segment_visualization/`): For structure and unassigned regions

Together they provide complete picture of recombination events!

---

**Status**: Complete ✅
**Date**: 2025-11-28
**Ready for**: Testing with `-resume`
