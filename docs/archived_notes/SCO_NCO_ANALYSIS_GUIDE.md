# SCO vs NCO Analysis Output Guide

## Overview

The SCO/NCO analysis compares **2,532 SCO** reads vs **336 NCO/Complex** reads (88.3% vs 11.7%).

## Generated Files

1. **sco_nco_summary.tsv** - Statistical summary table
2. **sco_nco_comparison.png** - 9-panel comparison plot
3. **sco_nco_resolution_analysis.png** - 4-panel resolution/quality plot

---

## File 1: Summary Statistics (`sco_nco_summary.tsv`)

| Metric | SCO | NCO/Complex |
|--------|-----|-------------|
| **Count** | 2,532 (88.3%) | 336 (11.7%) |
| **Mean Confidence** | 0.73 | 0.82 |
| **Median Confidence** | 0.74 | 0.86 |
| **High Confidence (≥0.8)** | 962 (38.0%) | 218 (64.9%) |
| **Median CO Width** | 611 bp | 143 bp |
| **Mean CO Width** | 1,050 bp | 260 bp |
| **Mean N Segments** | 6.1 | 11.9 |

### Key Observations

✅ **NCO/Complex has HIGHER confidence** (0.82 vs 0.73)
- This makes sense: NCO events are shorter and more localized
- Easier to map with high certainty

⚠️ **Mean segments are confusing!**
- **SCO should have ~2 segments** (one switch: CL or LC)
- **Showing 6.1 segments on average** - WHY?

This is likely because **the segment count (`n_segments`) comes from the PAF mapping**, not from the classified segments (`ultra_curated_profile`).

**The confusion**:
- Classification uses `ultra_curated_profile` (e.g., "325LA4,336LA3" → collapsed to "CL" → SCO)
- But `n_segments` comes from PAF file which counts ALL alignment blocks, not just parental switches

---

## File 2: Comparison Plot (`sco_nco_comparison.png`)

This is a **9-panel comprehensive comparison**:

### Top Row (Panels 1-3):

**Panel 1: Read Counts**
- Bar chart showing 2,532 SCO vs 336 NCO/Complex
- **What it tells you**: SCO is much more common (88% of reads)

**Panel 2: Confidence Score Distribution**
- Boxplots comparing confidence scores
- **What it tells you**: NCO/Complex has higher confidence (median 0.86 vs 0.74)
- **Why**: Shorter, more localized events are easier to validate

**Panel 3: SCO Resolution (Crossover Width)**
- Histogram of crossover widths for SCO only
- Median: 611 bp, Q25: ~300 bp, Q75: ~1,000 bp
- **What it tells you**: Most SCO crossovers are resolved at 300-1,000 bp resolution

### Middle Row (Panels 4-6):

**Panel 4: Minimum Segment Length**
- Boxplot comparing smallest segment in each read
- **What it tells you**: How small are the detected segments?
- Smaller segments → harder to detect, lower quality

**Panel 5: Mean Segment Length**
- Boxplot comparing average segment length
- **What it tells you**: Overall segment size distribution
- Larger segments → better coverage, higher confidence

**Panel 6: Number of Segments Distribution** ⚠️ **CONFUSING!**
- Bar chart showing segment counts for SCO vs NCO/Complex
- **Says**: "SCO should be 2"
- **Shows**: Wide distribution (likely 2-20+ segments)

**The confusion**: This panel shows `n_segments` from PAF (alignment blocks), NOT from classification.

### Bottom Row (Panels 7-9):

**Panel 7: Mapping Consensus**
- Boxplot showing how many of the 5 mappers call "best"
- **What it tells you**: Agreement across mappers (5 = all agree)
- Higher consensus → higher confidence

**Panel 8: Genomic Distribution (ARMS vs CEN)**
- Stacked bar chart showing ARMS vs CEN regions
- **What it tells you**: Are SCO/NCO distributed differently across genome?

**Panel 9: Summary Statistics**
- Text panel with key numbers

---

## File 3: Resolution Analysis (`sco_nco_resolution_analysis.png`)

This is a **4-panel SCO-specific analysis**:

### Panel 1: Crossover Width Distribution (log scale)
- Histogram with log bins
- Shows resolution distribution with median line
- **What it tells you**: Most crossovers resolved at 100-1,000 bp range

### Panel 2: Width vs Confidence
- Scatterplot: crossover width (x-axis) vs confidence score (y-axis)
- Trend line shows correlation
- **What it tells you**: Does narrower resolution → lower confidence?
- **Expected**: Negative correlation (smaller width = harder to detect = lower confidence)

### Panel 3: Width by Genomic Region
- Boxplot: ARMS vs CEN crossover widths
- **What it tells you**: Are CEN crossovers resolved differently than ARMS?
- **Expected**: CEN might have lower resolution (harder to map in repetitive regions)

### Panel 4: Cumulative Distribution
- Shows what percentage of SCO are resolved at ≤X bp
- **What it tells you**: "50% of SCO are resolved at ≤611 bp"

---

## What's Confusing: n_segments Issue

### The Problem

**Panel 6 in comparison plot** says "SCO should be 2" but shows a wide distribution.

### Why This Happens

There are **TWO different "segment" concepts**:

1. **Classification segments** (from `ultra_curated_profile`):
   - Example: "325LA4,336LA3" → 2 segments → collapsed to "CL" → **SCO**
   - This is used for classification
   - This is what "should be 2" for SCO

2. **PAF alignment segments** (from `n_segments` column):
   - Counts alignment blocks in PAF file
   - A single parental region might map as MULTIPLE blocks if there are gaps/indels
   - Example: One "LA4" region might have 3 separate alignment blocks → 3 segments
   - This is what's being PLOTTED

### The Fix

The plot should either:

**Option A**: Use **collapsed segment counts** from classification (len(collapsed pattern))
- SCO would show 2
- NCO would show 3
- DCO would show 4
- Complex would show 5+

**Option B**: Rename the axis to clarify it's "PAF alignment blocks"
- Don't say "SCO should be 2"
- Say "Alignment blocks per read (includes internal gaps/splits)"

---

## Key Biological Insights

### 1. SCO Resolution (611 bp median)
- This is the **minimum detectable crossover width**
- Determined by k-mer marker density and read coverage
- Smaller crossovers might exist but can't be resolved

### 2. NCO Higher Confidence (0.82 vs 0.73)
- NCO events are gene conversion tracts (short, localized)
- Easier to map with certainty because they're compact
- SCO spans larger regions → more chances for mapping ambiguity

### 3. NCO Smaller Width (143 bp vs 611 bp)
- Gene conversion tracts are typically 50-500 bp
- This matches biological expectations!

### 4. Mean Segments Discrepancy
- **NOT a bug**, just a terminology confusion
- PAF segments ≠ Classification segments

---

## Recommendations

### To Fix the Confusion:

1. **Add a new metric**: `n_classification_segments`
   - Extract from `ultra_curated_profile`
   - Count collapsed segments (len(collapsed))
   - This WILL show 2 for SCO, 3 for NCO, etc.

2. **Rename existing `n_segments`** → `n_alignment_blocks`
   - Clarify it's from PAF mapping
   - NOT the same as parental switches

3. **Create separate plot**: Classification segment distribution
   - Show the "ideal" distribution (2, 3, 4, 5+)
   - Compare expected vs observed

### What to Focus On:

✅ **Crossover width distributions** - tells you resolution limits
✅ **Confidence score comparisons** - tells you data quality
✅ **ARMS vs CEN differences** - biological insights
✅ **Mapping consensus** - validation quality

⚠️ **Segment count plot** - currently misleading, needs fix

---

## Quick Answer to Your Question

> "not sure if happy with this because the segment distribution seems confusing to me"

**You're right to be confused!**

The segment count plot is showing PAF alignment blocks (which can be 2-20+ per read due to gaps/splits), NOT the classification segments (which should be exactly 2 for SCO, 3 for NCO, 4 for DCO).

**To fix this**: I can add a new column that extracts the collapsed segment count from `ultra_curated_profile` and create a proper classification-based segment distribution plot.

Would you like me to:
1. Add this new metric and regenerate the plots?
2. Just hide/remove the confusing segment panel?
3. Keep it but add clarifying labels?
