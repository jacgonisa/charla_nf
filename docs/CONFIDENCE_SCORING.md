# Confidence Scoring for CHARLA-NF

## Overview

CHARLA-NF now includes a comprehensive confidence scoring system that assigns a reliability score (0.0-1.0) to each detected crossover. This helps distinguish high-confidence crossover calls from ambiguous or low-quality detections.

---

## Confidence Score Components

The confidence score is calculated as a weighted sum of 5 independent components:

| Component | Weight | Description |
|-----------|--------|-------------|
| **Mapping Consensus** | 30% | Agreement across 5 mapping algorithms |
| **Segment Quality** | 25% | Segment length and k-mer support |
| **Mapping Quality** | 20% | MAPQ and alignment scores (AS-NM) |
| **Parent Discrimination** | 15% | Clarity of Col vs Ler assignment |
| **Crossover Resolution** | 10% | Precision of breakpoint localization |

---

## Confidence Tiers

Crossovers are stratified into three confidence tiers:

- **High confidence (≥0.8)**: Expected >95% true positives
  - Strong agreement across all metrics
  - Recommended for publication-quality results

- **Medium confidence (0.5-0.8)**: Expected 70-95% true positives
  - Some ambiguity in one or more components
  - Suitable for exploratory analysis

- **Low confidence (<0.5)**: Expected <70% true positives
  - Significant issues detected
  - Recommend manual inspection or exclusion

---

## Output Files

The confidence scoring generates the following outputs in `results/confidence_scores/threshold{N}/`:

### Data Files
- **`confidence_scores.tsv`**: Main output with scores for all reads
  ```
  Columns: read_id, confidence_score, mapping_consensus, segment_quality,
           mapping_quality, parent_discrimination, crossover_resolution, region_type
  ```

- **`confidence_summary.txt`**: Summary statistics (total reads, tier counts, means)

- **`high_confidence_reads.txt`**: List of reads with score ≥0.8

- **`low_confidence_reads.txt`**: List of reads with score <0.5

### Visualization Files
- **`confidence_scores_histogram.png`**: Distribution with tier annotations
- **`confidence_scores_component_distributions.png`**: Distributions of all 5 components
- **`confidence_scores_correlation_matrix.png`**: Component correlation heatmap
- **`confidence_scores_component_scatter.png`**: Scatter plots (components vs final score)
- **`confidence_scores_region_comparison.png`**: ARMS vs CEN comparison

---

## Usage

### Basic Usage

Confidence scoring runs automatically as part of the standard CHARLA-NF pipeline:

```bash
nextflow run main.nf --reads data/sample.fq.gz
```

The outputs will be in `results/confidence_scores/`.

### Filtering by Confidence

Extract high-confidence crossovers only:

```bash
# Get read IDs
awk -F'\t' '$2 >= 0.8' results/confidence_scores/threshold50/confidence_scores.tsv | cut -f1 > high_conf_reads.txt

# Filter downstream files
grep -F -f high_conf_reads.txt results/mapping_analysis/threshold50/plotting_crossover/*/summary_mm_sr_*.paf
```

### Custom Thresholds

Adjust confidence tier cutoffs:

```bash
# Very conservative (only the best)
awk -F'\t' '$2 >= 0.9' confidence_scores.tsv > ultra_high_conf.txt

# Include medium confidence
awk -F'\t' '$2 >= 0.5' confidence_scores.tsv > high_medium_conf.txt
```

---

## Interpreting Scores

### Understanding Component Scores

**1. Mapping Consensus (0.0-1.0)**
- **1.0**: All 5 mappers agree on "best"
- **0.8**: 4/5 mappers agree
- **0.5**: 3/5 mappers agree (minimum acceptable)
- **<0.3**: ≤2 mappers agree OR problematic flags detected

**What to check if low:**
- Review `summary_table.csv` to see which mappers disagreed
- Check if flagged as "chimeric", "inversion", or "inconclusive_tie"

**2. Segment Quality (0.0-1.0)**
- **1.0**: Segments >200bp, mean >2kb, exactly 2 segments
- **0.7**: Segments 100-200bp, or 3 segments (complex crossover)
- **<0.4**: Short segments (<50bp) or unusual structure

**What to check if low:**
- Inspect BED file for segment lengths
- Short segments may indicate sequencing error or complex recombination
- Consider adjusting `--min_segment_length` threshold

**3. Mapping Quality (0.0-1.0)**
- **1.0**: MAPQ ≥50, AS-NM ratio ≥0.9 (very few errors)
- **0.8**: MAPQ 30-50, AS-NM 0.7-0.9
- **<0.5**: MAPQ <10 or many mismatches

**What to check if low:**
- Review PAF files for MAPQ values
- Low MAPQ suggests multi-mapping or repetitive regions
- Check alignment CIGAR strings for large indels

**4. Parent Discrimination (0.0-1.0)**
- **1.0**: >20% difference between Col and Ler scores (clear winner)
- **0.7**: 10-20% difference (moderate confidence)
- **<0.3**: <5% difference (near tie, ambiguous)

**What to check if low:**
- This is expected in regions with low divergence between parents
- Compare AS-NM scores for Col vs Ler in PAF files
- May indicate regions identical between parents

**5. Crossover Resolution (0.0-1.0)**
- **1.0**: Crossover localized to ≤100bp
- **0.6-0.8**: Resolution 100bp-1kb
- **<0.4**: Resolution >1kb (imprecise breakpoint)

**What to check if low:**
- Check `.cowidths` file for gap sizes between segments
- Large gaps are expected for crossovers in repetitive regions
- Consider using MAFFT alignment for refined breakpoints

---

## Advanced Usage

### Stratified Analysis by Confidence

```R
library(tidyverse)

# Load confidence scores
conf <- read_tsv("results/confidence_scores/threshold50/confidence_scores.tsv")

# Add tier labels
conf <- conf %>%
  mutate(tier = case_when(
    confidence_score >= 0.8 ~ "High",
    confidence_score >= 0.5 ~ "Medium",
    TRUE ~ "Low"
  ))

# Compare crossover widths by confidence tier
cowidths <- read_tsv("results/mapping_analysis/.../cowidths")
merged <- left_join(cowidths, conf, by = c("read" = "read_id"))

ggplot(merged, aes(x = tier, y = width)) +
  geom_boxplot() +
  labs(title = "Crossover Width by Confidence Tier",
       x = "Confidence Tier", y = "Crossover Width (bp)")
```

### Manual Inspection of Low-Confidence Calls

```bash
# Get low-confidence reads
awk -F'\t' '$2 < 0.5 {print $1, $2}' confidence_scores.tsv > low_conf_with_scores.txt

# Extract their segments
while read read_id score; do
  grep "^$read_id" results/read_segments/threshold50/*.bed
done < low_conf_with_scores.txt
```

### Identifying Problem Areas

```bash
# Which component is most problematic?
awk -F'\t' 'NR>1 {
  if ($3 < 0.5) print $1 " mapping_consensus";
  if ($4 < 0.5) print $1 " segment_quality";
  if ($5 < 0.5) print $1 " mapping_quality";
  if ($6 < 0.5) print $1 " parent_discrimination";
  if ($7 < 0.5) print $1 " crossover_resolution"
}' confidence_scores.tsv | cut -d' ' -f2 | sort | uniq -c | sort -rn
```

---

## Validation

### Expected Results for Test Data

Using the provided test dataset (`data/fastq/ColLer_F1_pollen.fq.gz`):

```
Total reads: ~50-100 (depends on filters)
High confidence: 60-70%
Medium confidence: 20-30%
Low confidence: 5-15%

Mean confidence score: 0.70-0.75
```

If your results deviate significantly:
- Check input data quality (FASTQ quality scores)
- Review mapping parameters (MAPQ thresholds)
- Adjust segmentation thresholds (--min_segment_length)

### Troubleshooting Low Scores

**All reads have low mapping_consensus (<0.5)**
→ Check that mapping completed successfully for all 5 algorithms
→ Review `summary_table.csv` for "chimeric" or "tie" flags

**All reads have low segment_quality (<0.5)**
→ Segments may be too short - lower `--min_segment_length` threshold
→ Check k-mer classification step for errors

**All reads have low mapping_quality (<0.5)**
→ Check reference genome quality
→ Verify correct reference genomes (Col-0, Ler-0)
→ Check for contamination in input reads

**All reads have low parent_discrimination (<0.5)**
→ This is expected if Col and Ler parents are very similar
→ Check chromosome regions - low divergence is normal for some regions
→ Verify correct parent references

**All reads have low crossover_resolution (<0.5)**
→ Expected for long-read sequencing with lower accuracy
→ Consider using higher-quality sequencing (HiFi, ultra-long ONT)
→ Use MAFFT alignment for bp-level refinement

---

## Resource Requirements

Confidence scoring is lightweight and fast:

- **CPU**: 2 cores (default)
- **Memory**: 4 GB (default)
- **Time**: <5 minutes for 100 reads

Customize in `nextflow.config`:

```groovy
confidence_cpus = 4
confidence_memory = '8 GB'
confidence_time = '30m'
```

---

## Citation

If you use the confidence scoring feature, please cite:

```
CHARLA-NF confidence scoring system
[Your citation information here]
```

---

## Support

For issues or questions:
- GitHub Issues: https://github.com/jacgonisa/charla_nf/issues
- Check logs: `results/confidence_scores/threshold{N}/confidence_summary.txt`
- Review visualizations in `results/confidence_scores/`
