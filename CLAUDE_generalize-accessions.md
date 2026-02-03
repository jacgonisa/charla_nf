# CLAUDE Session Summary - feature/generalize-accessions Branch

**Branch**: `feature/generalize-accessions`
**Last Updated**: 2025-02-02
**Session Focus**: Post-mapping analysis module creation and integration

---

## 🎯 What Was Accomplished

### Major Achievement: Post-Mapping Analysis Module

Created a **self-contained module** (like `cenhapmer_generation/`) for "AFTER CHARLA" analysis:

```
postmapping_analysis/
├── scripts/               # R scripts (14-19)
├── modules/               # Nextflow modules (some integrated)
├── README.md              # Main documentation
├── INTEGRATION_STATUS.md  # Which modules are integrated vs standalone
├── EXAMPLES.md            # Real data examples and results
├── COMPARISON.md          # Before/after comparison
└── TODO.md                # Future integration plans
```

---

## 🔑 Critical Fix: Comprehensive Alignment Output

### Problem Identified
Collaborators reported that segment order was **reversed** between CHARLA's output and their analysis tools.

### Root Cause
CHARLA was ordering segments by **read position** (where they appear on the sequenced read), but collaborators needed **reference position** ordering (where they appear on the genome).

### Example of the Issue
```
Read: 4201b2b4-2fd6-4c98-9a7b-120082805b05 (35,460 bp)

CHARLA (read position order):
  B73 segment:  852-12,362    (first on read)
  Mo17 segment: 17,430-35,460 (second on read)

But on reference genome (chr5):
  Mo17 aligns to: 156,805,217 (156M) ← comes FIRST
  B73 aligns to:  183,001,790 (183M) ← comes SECOND

Collaborators' tools showed: Mo17 first, B73 second (correct biological order)
```

### Solution
**Script 19** (`postmapping_analysis/scripts/19-create_comprehensive_alignment_output.R`) now:
- Orders segments by **reference genome position** (not read position)
- Combines k-mer BED (pre-mapping) + PAF (post-mapping) coordinates
- Creates single-row comprehensive output matching collaborator format
- Handles minus-strand alignments correctly

**Key code section** (lines 174-191):
```r
# Order by reference position
comprehensive_with_order <- comprehensive %>%
    group_by(read_id) %>%
    arrange(read_id, ref_start) %>%  # ← KEY: order by ref position
    mutate(
        ref_order = row_number(),
        first_parent = segment_parent[ref_order == 1],
        second_parent = segment_parent[ref_order == 2]
    )
```

---

## 📊 Output Format (Collaborator-Compatible)

**Header**:
```
read_id | Mo17_read_start | Mo17_read_end | Mo17_align_length | Mo17_chr |
Mo17_ref_start | Mo17_ref_end | Mo17_strand | B73_read_start | B73_read_end |
B73_read_length | B73_chr | B73_ref_start | B73_ref_end | B73_strand | recomb_type
```

**Key features**:
- Single row per read (easy to analyze)
- Both read positions AND genome positions
- Segments ordered by reference position
- Explicit strand information
- Complete read length for both parents

---

## 🔬 Analysis Results (Maize Data - 9,277 Reads)

### Quality Metrics
- **75.3%** have non-overlapping segments (clean boundaries)
- **97.3%** have both parents on same chromosome
- **1.3%** have missing alignment data (NAs)

### Biological Findings
- **18.8%** (1,742 reads) show **inversions**!
  - Segments align to opposite strands (+ vs -)
  - Indicates **structural variation** between B73 and Mo17
  - Not an error - real biological phenomenon

### Why NAs Occur
Segments too short to align reliably:
- Failed segments: Mean 151-743 bp, Median 63-323 bp
- Successful segments: Mean 4,004-4,299 bp, Median 1,463-1,754 bp
- **Failed segments are 5-28× smaller** than successful ones

---

## 🔄 Integration Status

### ✅ INTEGRATED (runs automatically in pipeline)

1. **CROSSOVER_STATS** (script 14)
   - Module: `postmapping_analysis/modules/crossover_stats.nf`
   - Calculates event type distributions, segment sizes

2. **CURATE_CONFIDENT_READS** (script 15)
   - Module: `postmapping_analysis/modules/curate_confident_reads.nf`
   - Filters and curates high-confidence crossover reads
   - Creates both pre-mapping and post-mapping BED files

3. **CREATE_COMPREHENSIVE_ALIGNMENT_OUTPUT** (script 19) ⭐ **MOST CRITICAL**
   - Module: `postmapping_analysis/modules/create_comprehensive_alignment_output.nf`
   - Creates single-row comprehensive TSV
   - Orders by reference position
   - **This is the main output collaborators need!**

**Workflow includes** (lines 35-40 in `workflows/charla_nf.nf`):
```groovy
// Post-mapping analysis modules (INTEGRATED into pipeline)
include { CROSSOVER_STATS as CROSSOVER_STATS_PROFILES } from '../postmapping_analysis/modules/crossover_stats.nf'
include { CURATE_CONFIDENT_READS } from '../postmapping_analysis/modules/curate_confident_reads.nf'
include { CREATE_COMPREHENSIVE_ALIGNMENT_OUTPUT } from '../postmapping_analysis/modules/create_comprehensive_alignment_output.nf'
```

### ⚠️ NOT INTEGRATED (manual use only)

4. **PLOT_CONFIDENT_READS** (script 18)
   - Comprehensive read visualization
   - Commented out in workflow (lines 666-684)

5. **KARYOPLOT_CROSSOVERS** (script 15-karyoplot)
   - Karyotype plots showing crossover positions
   - Requires karyoploteR R package
   - Commented out in workflow (lines 726-742)

6. **Single chromosome plots** (script 17)
   - No Nextflow module yet

7. **Circos plots** (script 16)
   - No Nextflow module yet

**To enable visualization**: Uncomment lines in `workflows/charla_nf.nf`

---

## 📁 File Locations

### Scripts (moved to postmapping_analysis/)
- `postmapping_analysis/scripts/14-crossover_statistics.R`
- `postmapping_analysis/scripts/15-curate_confident_reads.R`
- `postmapping_analysis/scripts/15-karyoplot_crossovers.R`
- `postmapping_analysis/scripts/16-circos_and_karyoplot_crossovers.R`
- `postmapping_analysis/scripts/17-single_chr_crossover_plot.R`
- `postmapping_analysis/scripts/18-plot_confident_reads_comprehensive.R`
- `postmapping_analysis/scripts/19-create_comprehensive_alignment_output.R` ⭐

### Modules (moved to postmapping_analysis/)
- `postmapping_analysis/modules/crossover_stats.nf`
- `postmapping_analysis/modules/curate_confident_reads.nf`
- `postmapping_analysis/modules/create_comprehensive_alignment_output.nf`
- `postmapping_analysis/modules/plot_confident_reads.nf`
- `postmapping_analysis/modules/karyoplot_crossovers.nf`

### Documentation
- `postmapping_analysis/README.md` - Main documentation
- `postmapping_analysis/INTEGRATION_STATUS.md` - Integration details ⭐
- `postmapping_analysis/EXAMPLES.md` - Real data examples
- `postmapping_analysis/COMPARISON.md` - Before/after comparison
- `postmapping_analysis/TODO.md` - Future plans

---

## 🔧 Technical Details

### Coordinate Systems

**Pre-mapping (k-mer based)**:
- From BED files: `*_premapping_readcoords.bed`
- Shows where segments are on the original read
- Example: B73 at 852-12,362, Mo17 at 17,430-35,460

**Post-mapping (alignment based)**:
- From PAF files: alignment results
- Shows where segments align on reference genome
- Example: B73 at chr5:183M, Mo17 at chr5:156M

**Comprehensive output (script 19)**:
- Combines both coordinate systems
- Orders by reference position (156M < 183M → Mo17 first)
- Single row per read

### Data Flow

```
CHARLA Pipeline
    ↓
K-mer Analysis → BED files (segments on read)
    ↓
Segmentation → Extract segments
    ↓
Alignment → PAF files (segments on genome)
    ↓
[POST-MAPPING ANALYSIS MODULE]
    ↓
Script 15: Curate confident reads
    ↓
Script 19: Create comprehensive output ⭐
    ├→ Read k-mer BED (pre-mapping coords)
    ├→ Read PAF files (post-mapping coords)
    ├→ Merge on read_id + parent
    ├→ Order by reference position
    └→ Output single-row TSV
    ↓
Script 14: Calculate statistics
    ↓
Scripts 16-18: Visualization (optional)
```

---

## 🐛 Known Issues & Fixes

### Issue 1: Overlapping Read Coordinates
**Problem**: Mo17 and B73 both showed coordinates starting at 0
**Cause**: Using PAF query coordinates (segment-relative) instead of BED (full-read)
**Fix**: Script 19 now uses k-mer BED for read positions

### Issue 2: Cross-Parent Alignments
**Problem**: Mo17 segments aligning to B73 genome appearing in results
**Cause**: Not filtering by segment parent
**Fix**: Added filter (line 81): `paf_all %>% filter(segment_parent == alignment_parent)`

### Issue 3: Segment Order Reversal
**Problem**: Segments appeared in wrong order (B73 first instead of Mo17)
**Cause**: Ordering by read position instead of reference position
**Fix**: Changed to order by `ref_start` (reference position)

### Issue 4: Workflow Broken After Module Move
**Problem**: Workflow trying to include modules from old location
**Cause**: Moved modules to `postmapping_analysis/` but didn't update includes
**Fix**: Updated include paths to `../postmapping_analysis/modules/`

---

## 🎨 Cenhapmer Generation Improvements

Also improved the cenhapmer generation module:
- Enhanced marker density analysis (script 05)
- Added marker comparison plotting (script 07)
- Various script improvements (01-06)

---

## 📝 .gitignore Updates

Added exclusions for:
```
# Build artifacts
target/
rust/*/target/

# Backup files
*.bak
*.broken

# Generated data
cenhapmer_generation/genomes_split/

# Private communications (kept local)
docs/EMAIL_TO_COLLABORATORS.md
docs/COLLABORATOR_FEEDBACK_RESPONSE.md
docs/QUESTION_FOR_COLLABORATORS.md
docs/COORDINATE_HANDLING_EXPLANATION.md
docs/PREMAPPING_VS_POSTMAPPING.md
docs/KARYOPLOTER_SETUP.md
```

---

## 🚀 How to Use

### Running the Pipeline (with integrated modules)

```bash
nextflow run main.nf \
  --input samples.csv \
  --reference_genomes_dir genomes/ \
  --parent1_name B73 \
  --parent2_name Mo17 \
  --outdir results/
```

**Automatic outputs**:
- `results/curated_confident_reads/*_comprehensive_alignments.tsv` ⭐
- `results/curated_confident_reads/*_confident_SCO_premapping_readcoords.bed`
- `results/curated_confident_reads/*_confident_SCO_postmapping_refcoords.bed`
- `results/curated_confident_reads/*_crossover_statistics.tsv`

### Running Visualization Manually

```bash
# Karyoplots
Rscript postmapping_analysis/scripts/15-karyoplot_crossovers.R \
  --bed_file crossovers.bed \
  --parent1_fai B73.fa.fai \
  --parent2_fai Mo17.fa.fai \
  --sample_id maize_pollen \
  --output_dir plots/

# Comprehensive plots
Rscript postmapping_analysis/scripts/18-plot_confident_reads_comprehensive.R \
  --confident_bed_sco SCO.bed \
  --sample_id maize_pollen \
  --output_dir plots/
```

---

## 📦 Dependencies

### Integrated Modules (required)
```r
dplyr (>= 1.0.0)
tidyr (>= 1.0.0)
optparse
```

### Visualization (optional - for manual scripts)
```r
ggplot2
karyoploteR
BSgenome
regioneR
RColorBrewer
viridis
```

---

## 🔮 Future Work

### Immediate (Next Session)
1. Test full pipeline run on real data
2. Verify comprehensive output format matches collaborators' expectations
3. Test visualization scripts manually

### Short-term
1. Enable visualization modules in pipeline (uncomment)
2. Add parameter flag: `--run_visualization true/false`
3. Create Docker container with R dependencies

### Long-term
1. Multi-sample comparison analysis
2. Interactive visualization (Shiny app)
3. Gene annotation integration
4. Machine learning classification
5. Crossover interference analysis

**See `postmapping_analysis/TODO.md` for detailed plan**

---

## 🗂️ Git Commits Summary

### Commit 1: `5fd70b1` - Initial module creation
```
Add post-mapping analysis module and fix comprehensive alignment output

- Created postmapping_analysis/ folder structure
- Moved scripts 14-19 and modules
- Fixed comprehensive alignment (ordered by ref position)
- Updated .gitignore
```

### Commit 2: `73f24ac` - Documentation
```
Add comprehensive documentation and examples for postmapping analysis

- EXAMPLES.md (real data results)
- TODO.md (integration plans)
- COMPARISON.md (before/after comparison)
```

### Commit 3: `00e4c87` - Integration fix
```
Fix workflow integration and clarify module integration status

- Updated workflow includes to postmapping_analysis/modules/
- Created INTEGRATION_STATUS.md
- Updated README with Quick Start
- Commented out visualization modules (not yet integrated)
```

---

## 🎓 Key Learnings

### 1. Coordinate Systems Matter
- K-mer analysis uses full read coordinates
- PAF uses segment-relative coordinates
- Need to combine both for complete picture
- **Reference position ordering is biologically meaningful**

### 2. Minus-Strand Handling
- Both parents can align to minus strand
- Coordinates stay in forward read orientation
- Strand information is critical for detecting inversions
- 18.8% of reads show inversions (structural variation)

### 3. Module Organization
- Self-contained modules (like cenhapmer_generation/) work well
- Clear separation: BEFORE (cenhapmer) vs AFTER (postmapping)
- Some integration is good, full integration isn't always necessary
- Documentation is critical for clarity

### 4. Collaborator Compatibility
- Single-row format much easier to work with
- Reference-based ordering matches external tools
- Explicit strand information prevents confusion
- Complete coordinates enable downstream analysis

---

## 🆘 Troubleshooting

### Pipeline fails with "module not found"
**Check**: Are postmapping modules in `postmapping_analysis/modules/`?
**Fix**: Verify include paths in `workflows/charla_nf.nf` (lines 35-40)

### Comprehensive output has NAs
**Cause**: Segments too short to align (< 100bp typically)
**Expected**: ~1-2% of reads will have NAs
**Not a bug**: These are legitimate short segments

### Segments in wrong order
**Check**: Is script 19 ordering by `ref_start`?
**Verify**: Lines 176-183 in script 19 should order by reference position

### Visualization modules don't run
**Expected**: These are commented out in workflow
**Solution**: Use scripts manually or uncomment lines 666-684 and 726-742

---

## 📞 Quick Reference

**Main output file**: `*_comprehensive_alignments.tsv`
**Integration status**: See `postmapping_analysis/INTEGRATION_STATUS.md`
**Examples**: See `postmapping_analysis/EXAMPLES.md`
**Comparison**: See `postmapping_analysis/COMPARISON.md`
**Future plans**: See `postmapping_analysis/TODO.md`

**GitHub**: https://github.com/jacgonisa/charla_nf
**Branch**: feature/generalize-accessions
**Last commit**: `00e4c87`

---

## 💡 For Future Sessions

### If starting fresh
1. Read this file first
2. Check `postmapping_analysis/INTEGRATION_STATUS.md`
3. Review the three commits (`5fd70b1`, `73f24ac`, `00e4c87`)
4. Test script 19 on example data

### If adding features
1. Check `postmapping_analysis/TODO.md` for plans
2. Decide if new script goes in postmapping_analysis/
3. Create Nextflow module if integrating
4. Update INTEGRATION_STATUS.md

### If debugging
1. Check if issue is pre-mapping or post-mapping
2. Verify coordinate systems (read vs reference)
3. Check for NAs (short segments)
4. Verify segment ordering (by ref_start)

---

*This document should be enough to pick up where we left off!*
*Created: 2025-02-02*
