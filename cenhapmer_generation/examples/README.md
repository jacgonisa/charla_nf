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
- **k=21**: Lower overall marker density
  - ARMS and CEN distributions are more similar
  - Fewer unique markers overall

- **k=31**: Moderate marker density
  - Clear separation between ARMS and CEN distributions
  - More markers than k=21

- **k=41**: Highest marker density
  - Broadest distributions
  - Most unique markers
  - Best discrimination between ARMS and CEN

**Important insight**: These plots clearly demonstrate that **higher k-mer sizes produce MORE unique markers** in both ARMS and CEN regions. However, remember that higher k-mer sizes are also **less resilient to sequencing errors** - a single error breaks the entire k-mer.

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
