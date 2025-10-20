# Example Marker Density Analysis Results

This directory contains example outputs from marker density analysis runs to help you understand what to expect from your own analyses.

## Arabidopsis Col-0 Example

The `arabidopsis_Col-0/` directory contains real marker density analysis results from *Arabidopsis thaliana* Col-0 genome.

### Files

#### `Col-0_k21_10kb_marker_counts_chromosomes_row.png`
Genome-wide marker coverage across all 5 Arabidopsis chromosomes using k=21 and 10kb windows. This plot shows:
- **One subplot per chromosome** (Chr1-Chr5)
- **X-axis**: Genomic position (bp)
- **Y-axis**: Marker count per 10kb window
- Marker density along each chromosome showing peaks and valleys

**What to look for:**
- Coverage across entire chromosomes (no large gaps)
- Distribution patterns in both centromeric (CEN) and chromosome arm (ARMS) regions
- Peaks indicate regions with high marker density
- Valleys/dips may indicate centromeres or repetitive regions

#### `k21_ARMS_CEN_overlaid_histogram.png`, `k31_ARMS_CEN_overlaid_histogram.png`, `k41_ARMS_CEN_overlaid_histogram.png`
Overlaid histograms comparing marker density between ARMS (chromosome arms) and CEN (centromeric) regions for different k-mer sizes:
- **X-axis**: Marker count per window
- **Y-axis**: Frequency (log scale)
- **Blue**: ARMS regions
- **Red**: CEN regions

**Key observations from Arabidopsis data:**

**ARMS regions (blue - chromosome arms):**
- k=21: 25th percentile = 21 markers/kb
- k=31: 25th percentile = 41 markers/kb
- k=41: 25th percentile = 41 markers/kb
- **Observation**: Marker density remains relatively **constant** across k-mer sizes

**CEN regions (red - centromeres):**
- k=21: 25th percentile = 16 markers/kb
- k=31: 25th percentile = ~40 markers/kb
- k=41: 25th percentile = 84 markers/kb
- **Observation**: Marker density **increases dramatically** with larger k-mer sizes

**Critical biological insight**:

Centromeric regions are **highly repetitive and difficult to resolve**. The dramatic increase in CEN marker density with larger k-mer sizes shows that:

1. **Higher k provides better centromere resolution** - More unique k-mers can be found in repetitive centromeric DNA when using longer sequences
2. **ARMS regions are already well-resolved** - Chromosome arms have sufficient sequence diversity that even k=21 provides good marker coverage
3. **Tradeoff consideration**: While k=41 gives 5× better centromere resolution than k=21, it's also much less tolerant of sequencing errors

**Recommendation**: If your analysis focuses on **centromeric crossovers**, consider using k=31 or k=41 for better resolution (if your read quality allows it). For **chromosome arm crossovers**, k=21 is often sufficient and more error-tolerant.

## Interpreting Your Own Results

When you run marker density analysis on your own genomes, compare your results to these examples:

1. **Genome-wide coverage plot**: Should show relatively even coverage across chromosomes, with some expected variation in centromeric regions

2. **Overlaid histogram distributions**:
   - Should follow similar pattern: k=21 < k=31 < k=41 in terms of marker density
   - ARMS and CEN regions should show clear distributions (may overlap)
   - If you see dramatically different patterns, check your input data quality

## Choosing K-mer Size

Based on your sequencing platform and data quality:

- **Noisy long reads (ONT, PacBio CLR)**: Use **k=21**
  - Fewer markers, but more error-tolerant
  - A single error won't destroy all markers in a region

- **High-quality long reads (PacBio HiFi)**: Use **k=31**
  - Good balance of markers and error tolerance
  - Recommended for most use cases

- **Ultra high-quality data**: Use **k=41**
  - Maximum markers and specificity
  - Only if your error rate is very low (<0.1%)

## Citation

These example images were generated from publicly available *Arabidopsis thaliana* Col-0 genome assembly.
