# Post-Mapping Analysis Module

This module contains scripts and tools for analyzing crossover reads **after** alignment to reference genomes. This is the "AFTER CHARLA" step that provides detailed single-read level analysis, statistics, and visualization.

## Overview

While the main CHARLA pipeline (in `charla_nf.nf`) identifies crossover reads using k-mer analysis, this module provides:
- Single-read level alignment coordinates
- Comprehensive alignment output matching collaborator formats
- Statistical analysis of crossover events
- Visualization (karyoplots, crossover maps)

## Directory Structure

```
postmapping_analysis/
├── scripts/          # R scripts for analysis and visualization
├── modules/          # Nextflow modules (optional integration)
└── README.md         # This file
```

## Scripts

### Core Analysis

#### `19-create_comprehensive_alignment_output.R`
**Purpose**: Creates single-read level comprehensive alignment output

**Key Features**:
- Combines k-mer analysis (pre-mapping) with PAF alignments (post-mapping)
- **Orders segments by reference genome position** (critical for biological interpretation)
- Handles minus-strand alignments correctly
- Filters cross-parent alignments
- Outputs in collaborator-requested format

**Input**:
- Curated confident read BED files (k-mer based segments)
- PAF alignment files for both parents
- Sample name

**Output**: TSV file with columns:
```
read_id | Mo17_read_start | Mo17_read_end | Mo17_align_length | Mo17_chr | Mo17_ref_start | Mo17_ref_end | Mo17_strand |
B73_read_start | B73_read_end | B73_read_length | B73_chr | B73_ref_start | B73_ref_end | B73_strand | recomb_type
```

**Usage**:
```bash
Rscript postmapping_analysis/scripts/19-create_comprehensive_alignment_output.R \
  --curated_dir /path/to/curated_confident_reads/ \
  --paf_b73 /path/to/segments_alnTo_B73.paf \
  --paf_mo17 /path/to/segments_alnTo_Mo17.paf \
  --sample_name sample_id \
  --output_dir /output/path/
```

**Key Technical Detail**: Segments are ordered by reference position (not read position). For example, if Mo17 aligns to chr5:156M and B73 aligns to chr5:183M, Mo17 appears first in the output regardless of where they appear on the read.

---

### Statistics & Curation

#### `14-crossover_statistics.R`
Calculates comprehensive statistics on crossover events:
- Event type distributions (SCO, DCO, GC-CO, NCO)
- Segment size distributions
- Chromosome-level statistics
- Crossover rates

#### `15-curate_confident_reads.R`
Curates and filters confident crossover reads:
- Removes ectopic recombination
- Filters by alignment quality
- Identifies high-confidence events
- Creates BED files with both pre-mapping and post-mapping coordinates

---

### Visualization

#### `15-karyoplot_crossovers.R`
Creates karyotype plots showing crossover positions on chromosomes using karyoploteR.

**Requirements**:
- karyoploteR package
- Reference genome BSgenome object

#### `16-circos_and_karyoplot_crossovers.R`
Advanced visualization combining circos-style plots with karyoplots.

#### `17-single_chr_crossover_plot.R`
Detailed single-chromosome crossover visualization showing:
- Individual crossover read positions
- Density plots
- Segment boundaries

#### `18-plot_confident_reads_comprehensive.R`
Comprehensive plotting of confident reads:
- Read structure visualization
- Segment alignment visualization
- Quality metrics plots

---

## Nextflow Modules

The `modules/` directory contains Nextflow wrappers for the scripts above. These are **optional** and provided for convenience. You can use the R scripts directly or integrate the modules into custom workflows.

### Available Modules

- `create_comprehensive_alignment_output.nf` - Wrapper for script 19
- `crossover_stats.nf` - Statistics calculation
- `curate_confident_reads.nf` - Read curation
- `karyoplot_crossovers.nf` - Karyotype visualization
- `plot_confident_reads.nf` - Comprehensive plotting

---

## Integration with Main Pipeline

This module is designed to work **after** the main CHARLA pipeline completes.

### Minimal Integration

The essential output from CHARLA for this module:
1. Curated confident read BED files (from CHARLA's curation step)
2. PAF alignment files (from CHARLA's mapping step)
3. Read classification files (SCO, NCO, etc.)

### Recommended Workflow

```
CHARLA Pipeline → Post-Mapping Analysis
     ↓                      ↓
  [K-mer Analysis]    [Script 19: Comprehensive Output]
  [Segmentation]      [Statistics & Plots]
  [Alignment]         [Visualization]
  [Curation]
```

---

## Key Findings from Analysis

From testing on maize pollen data (9,277 confident reads):

- **75.3%** of reads have non-overlapping segments
- **97.3%** of reads have both parents on the same chromosome
- **18.8%** of reads show inversions (segments align to opposite strands)
  - This indicates structural variation between parents
  - Legitimate biological phenomenon, not an error
- **1.3%** of reads have missing alignment data (NAs)
  - Due to very short segments (median 63-323 bp) that fail to align
  - Segments that fail are 5-28× smaller than successful alignments

---

## Coordinate Systems Explained

### Pre-Mapping (K-mer Based)
- Coordinates show **where on the original read** each segment is located
- Based on k-mer composition analysis
- Example: B73 segment at positions 852-12,362 on a 35kb read

### Post-Mapping (PAF Alignment)
- Coordinates show **where on the reference genome** each segment aligns
- Includes strand information (+/-)
- Example: B73 segment aligns to chr5:183,001,790-183,013,297 (-)

### Comprehensive Output (Script 19)
- Combines both coordinate systems
- **Segments ordered by reference position** (biologically meaningful)
- Includes both read positions (pre-mapping) and genome positions (post-mapping)

---

## Dependencies

### R Packages Required
```r
# Core
library(dplyr)
library(tidyr)
library(ggplot2)

# Visualization
library(karyoploteR)  # For karyoplots
library(BSgenome)     # For genome data

# I/O
library(optparse)     # Command-line parsing
```

### External Tools
- None - scripts work with CHARLA pipeline outputs

---

## Future Enhancements

Potential additions to this module:
- Interactive visualization (plotly, shiny)
- Machine learning-based crossover classification
- Comparative analysis across samples
- Integration with genetic map construction
- Crossover interference analysis
- Gene-level annotation of crossover positions

---

## Contact & Contributions

This module is part of the CHARLA pipeline. For questions about post-mapping analysis specifically, refer to the main CHARLA documentation or repository issues.

---

## License

Same as main CHARLA pipeline.
