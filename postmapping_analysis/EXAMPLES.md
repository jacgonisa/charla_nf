# Post-Mapping Analysis Examples

This document shows what the post-mapping analysis module produces and how it improves upon standard CHARLA output.

---

## What Does Standard CHARLA Output Provide?

CHARLA's main pipeline produces:
- **K-mer based segmentation** (BED files)
- **Classification** (SCO, DCO, NCO, GC-CO)
- **Alignment files** (PAF format)

**Example BED output** (pre-mapping coordinates):
```
read_id                                  start   end     segment  type
4201b2b4-2fd6-4c98-9a7b-120082805b05    852     12362   BH5      SCO
4201b2b4-2fd6-4c98-9a7b-120082805b05    17430   35460   MH5      SCO
```

This tells us:
- ✅ Which regions contain B73/Mo17 k-mers
- ✅ Where they are on the read
- ❌ **Missing**: Where these align on the reference genome
- ❌ **Missing**: Single-row per read format for easy analysis
- ❌ **Missing**: Biological ordering (reference position)

---

## What Does Post-Mapping Analysis Add?

### 1. Comprehensive Single-Read Output (Script 19)

**Output format** (one row per read):
```
read_id: 4201b2b4-2fd6-4c98-9a7b-120082805b05

Mo17_read_start    17430
Mo17_read_end      35460
Mo17_align_length  35460
Mo17_chr           chr5
Mo17_ref_start     156805217
Mo17_ref_end       156823253
Mo17_strand        -

B73_read_start     852
B73_read_end       12362
B73_read_length    35460
B73_chr            chr5
B73_ref_start      183001790
B73_ref_end        183013297
B73_strand         -

recomb_type        SCO
```

**Key improvements:**
- ✅ **Single row per read** - easy to analyze in spreadsheets/databases
- ✅ **Reference genome positions** - know exactly where crossovers occur
- ✅ **Biologically ordered** - Mo17 appears first because it comes first on chr5 (156M < 183M)
- ✅ **Strand information** - both segments align to minus strand
- ✅ **Complete coordinates** - both read positions AND genome positions

---

## Example Analysis Results

### From Testing on Maize Data (9,277 Confident Reads)

#### Segment Overlap Analysis
```
Non-overlapping segments:  6,988 reads (75.3%)
Overlapping segments:      2,289 reads (24.7%)
```

**What this means**: Most crossover reads (75%) have clean segment boundaries with no overlap, indicating high-quality recombination events.

---

#### Chromosome Concordance
```
Both parents on same chr:  8,913 reads (97.3%)
Parents on different chr:    364 reads ( 2.7%)
```

**What this means**: Nearly all crossovers occur within a single chromosome (expected). The 2.7% on different chromosomes may represent:
- Trans-chromosomal events (rare)
- Chimeric reads
- Assembly artifacts

---

#### Structural Variation Detection
```
Same strand orientation:   7,535 reads (81.2%)
Opposite strand (inversion): 1,742 reads (18.8%)
```

**Example read with inversion**:
```
Read: 00a7624f-d2db-4c2e-ac4a-efc6bbf9a500

Mo17: chr7:108,783,157-108,803,248 (-)  ← minus strand
B73:  chr7:15,807,312-15,815,685   (+)  ← plus strand
```

**What this means**: 18.8% of crossover reads show inversions where one parent's segment aligns to the opposite strand. This reveals **structural variation** between B73 and Mo17 genomes - regions that are inverted between parents.

---

#### Missing Alignment Data
```
Complete alignments:       9,157 reads (98.7%)
Missing data (NAs):          120 reads ( 1.3%)
```

**Why NAs occur**: Segments too short to align reliably

| Segment Status | Mean Size | Median Size |
|----------------|-----------|-------------|
| **Failed** (NA) | B73: 743 bp, Mo17: 151 bp | B73: 323 bp, Mo17: 63 bp |
| **Successful** | B73: 4,004 bp, Mo17: 4,299 bp | B73: 1,463 bp, Mo17: 1,754 bp |

**Failed segments are 5-28× smaller** than successful ones!

---

## Visualization Examples

### 1. Karyoplot (Script 15)

Shows crossover positions across all chromosomes:

```
Chr1  |====•====•==========•===|
Chr2  |=======•=====•=====•====|
Chr3  |==•=========•==•========|
...
Chr10 |=====•======•=====•=====|

• = Crossover event
```

**Use case**:
- Visualize genome-wide crossover distribution
- Identify hotspots or cold spots
- Compare between samples

---

### 2. Crossover Statistics (Script 14)

**Event Type Distribution:**
```
Event Type    Count    Percentage
SCO           6,842    73.8%
NCO           1,456    15.7%
GC-CO           654     7.0%
DCO             325     3.5%
```

**Segment Size Distribution:**
```
           Min      25%     Median    75%      Max      Mean
B73      50 bp    856 bp   1,463 bp  4,142 bp  85 kb   4,004 bp
Mo17     50 bp    932 bp   1,754 bp  4,891 bp  92 kb   4,299 bp
```

---

### 3. Single Chromosome Detailed View (Script 17)

Shows individual crossover reads on a single chromosome:

```
Chr5 Position (Mb):  0----50----100---150---200---250

Read 1  Mo17: ========•
        B73:           •========

Read 2  Mo17:     ====•
        B73:          •=======

Read 3  Mo17:               •=========
        B73:  =============•

• = Crossover point
```

**Use case**:
- Detailed examination of specific chromosomes
- Identify clustered crossovers
- Analyze interference patterns

---

## Comparison Table: Standard vs Post-Mapping Output

| Feature | Standard CHARLA | Post-Mapping Module |
|---------|----------------|---------------------|
| **K-mer segmentation** | ✅ Yes | ✅ Yes |
| **Read classification** | ✅ Yes (SCO/NCO/etc) | ✅ Yes |
| **Read coordinates** | ✅ Yes (BED) | ✅ Yes (enhanced) |
| **Reference positions** | ❌ No | ✅ **Yes** |
| **Single-row format** | ❌ No (multiple rows) | ✅ **Yes** |
| **Biological ordering** | ❌ No | ✅ **Yes** (by ref position) |
| **Strand information** | ⚠️ Implicit in PAF | ✅ **Explicit per parent** |
| **Statistics** | ⚠️ Basic | ✅ **Comprehensive** |
| **Visualization** | ❌ No | ✅ **Yes** (karyoplots, etc) |
| **Inversion detection** | ❌ No | ✅ **Yes** |
| **Quality metrics** | ⚠️ Basic | ✅ **Detailed** |
| **Collaborator format** | ❌ No | ✅ **Yes** |

---

## Real Data Example

### Input (from CHARLA pipeline):

**BED file** (k-mer segments):
```
4201b2b4    852    12362   BH5   SCO
4201b2b4   17430   35460   MH5   SCO
```

**PAF files** (alignments - separate files for each parent):
```
# B73 alignments
4201b2b4_1_BH5  11510  0  11510  -  chr5  226353449  183001790  183013297  ...

# Mo17 alignments
4201b2b4_2_MH5  18030  0  18030  -  chr5  226353449  156805217  156823253  ...
```

### Output (from Post-Mapping Analysis):

**Comprehensive alignment TSV** (single row):
```
read_id  Mo17_read_start  Mo17_read_end  Mo17_align_length  Mo17_chr  Mo17_ref_start  Mo17_ref_end  Mo17_strand  B73_read_start  B73_read_end  B73_read_length  B73_chr  B73_ref_start  B73_ref_end  B73_strand  recomb_type
4201b2b4  17430           35460          35460              chr5      156805217       156823253     -            852             12362         35460            chr5     183001790      183013297    -           SCO
```

**Key transformation:**
1. **Merged** k-mer BED and PAF data
2. **Ordered** by reference position (Mo17 @ 156M comes before B73 @ 183M)
3. **Single row** for easy downstream analysis
4. **Added** full read length for each parent
5. **Explicit** strand information

---

## Use Cases

### 1. Genetic Mapping
Use reference positions to:
- Build genetic maps
- Calculate recombination rates
- Identify recombination hotspots
- Compare to known genetic markers

### 2. Structural Variation Analysis
Use strand information to:
- Detect inversions between parents (18.8% of reads!)
- Identify rearrangements
- Study genome evolution

### 3. Quality Control
Use segment metrics to:
- Filter low-quality crossovers
- Identify chimeric reads
- Assess alignment quality
- Remove ectopic recombination

### 4. Comparative Genomics
Use complete coordinates to:
- Compare crossover patterns between samples
- Analyze parent-of-origin effects
- Study crossover interference
- Identify chromosome-specific patterns

### 5. Collaboration & Data Sharing
Single-row format allows:
- Easy import into spreadsheets
- Database storage
- Sharing with collaborators (matching their format!)
- Statistical analysis in R/Python

---

## Future Integration

**TODO**: Full integration into main Nextflow pipeline is planned. Current status:

- ✅ Scripts are functional and tested
- ✅ Nextflow modules created
- ⚠️ Integration into main workflow: **PENDING**
- ⚠️ Optional/required flag: **TO BE DECIDED**

**Target**: Make post-mapping analysis an optional pipeline step that can be enabled with a flag (e.g., `--postmapping_analysis true`)

---

## Summary

The post-mapping analysis module transforms CHARLA's k-mer based output into:

1. **Single-read level comprehensive data** with both read and genome coordinates
2. **Biologically meaningful ordering** based on reference positions
3. **Detailed statistics** on crossover patterns and quality
4. **Publication-ready visualizations** (karyoplots, chromosome maps)
5. **Structural variation detection** (inversions between parents)
6. **Collaborator-compatible format** matching external analysis tools

This bridges the gap between CHARLA's k-mer analysis (BEFORE) and downstream biological interpretation (AFTER).
