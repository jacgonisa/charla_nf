# CHARLA Standard Output vs Post-Mapping Analysis

## Visual Comparison

```
┌─────────────────────────────────────────────────────────────────────────┐
│                     CHARLA STANDARD PIPELINE                            │
└─────────────────────────────────────────────────────────────────────────┘

Input: Long reads (ONT/PacBio)
   ↓
K-mer Analysis → Segmentation → Alignment → Classification
   ↓                ↓              ↓            ↓
Output Files:
├── BED (segments)
│   read_id    start   end    segment  type
│   read123    1000    5000   BH2      SCO  ← Multiple rows per read
│   read123    6000    9000   MH2      SCO  ← Hard to analyze
│
└── PAF (alignments - separate for each parent)
    read123_1_BH2  ...  chr1  10000000  10004000  ...  ← Segment coordinates
    read123_2_MH2  ...  chr1  15000000  15003000  ...  ← Not full read

❌ Issues:
• Multiple rows per read (hard to analyze)
• No reference positions in BED
• Segments not ordered biologically
• No single comprehensive view


┌─────────────────────────────────────────────────────────────────────────┐
│              POST-MAPPING ANALYSIS MODULE (NEW!)                         │
└─────────────────────────────────────────────────────────────────────────┘

Takes CHARLA output + Adds comprehensive analysis
   ↓
Script 19: Comprehensive Alignment
   ↓
Output: Single-row per read with BOTH read and genome coordinates

read_id  | Mo17_read | Mo17_genome      | B73_read | B73_genome       | strand
─────────┼───────────┼──────────────────┼──────────┼──────────────────┼────────
read123  | 6000-9000 | chr1:15M-15.003M | 1000-5000| chr1:10M-10.004M | -/-

✅ Improvements:
• One row per read (easy analysis!)
• Both read AND genome coordinates
• Ordered by reference position (biologically meaningful)
• Strand information explicit
• Ready for downstream analysis


┌─────────────────────────────────────────────────────────────────────────┐
│                    ADDITIONAL ANALYSES                                   │
└─────────────────────────────────────────────────────────────────────────┘

Statistics (Script 14)        Karyoplots (Script 15)
┌─────────────────┐          ┌────────────────────┐
│ Event Type  %   │          │ Chr1  |•••  •  •| │
│ SCO        73.8 │          │ Chr2  |• •••   •| │
│ NCO        15.7 │          │ Chr3  |••  ••  •| │
│ GC-CO       7.0 │          │  ...              │
│ DCO         3.5 │          │ Chr10 |•  ••  ••| │
└─────────────────┘          └────────────────────┘

Quality Metrics              Single Chr View (Script 17)
┌──────────────────┐         ┌───────────────────────┐
│ Non-overlapping  │         │ Chr5 (200Mb)          │
│   segments: 75%  │         │ ├─Mo17────•          │
│                  │         │ │         └─B73────   │
│ Same chr: 97%    │         │ ├──B73───•           │
│                  │         │ │        └─Mo17────   │
│ Inversions: 19%  │         └───────────────────────┘
└──────────────────┘
```

---

## Side-by-Side Example

### Standard CHARLA Output

**File 1: BED (k-mer segments)**
```
read_id                               start   end     segment  type
4201b2b4-2fd6-4c98-9a7b-120082805b05  852     12362   BH5      SCO
4201b2b4-2fd6-4c98-9a7b-120082805b05  17430   35460   MH5      SCO
```
*Problem*: Two rows, no genome positions

**File 2: PAF (B73 alignments)**
```
4201b2b4..._1_BH5  11510  0  11510  -  chr5  ...  183001790  183013297  ...
```
*Problem*: Segment coordinates (0-11510), not full read coordinates

**File 3: PAF (Mo17 alignments)**
```
4201b2b4..._2_MH5  18030  0  18030  -  chr5  ...  156805217  156823253  ...
```
*Problem*: Separate file, hard to correlate

---

### Post-Mapping Analysis Output

**Single comprehensive file:**
```
read_id: 4201b2b4-2fd6-4c98-9a7b-120082805b05

┌─────────────┬───────────┬─────────────┬────────────────────┬────────┐
│ Parent      │ Read Pos  │ Align Length│ Genome Position    │ Strand │
├─────────────┼───────────┼─────────────┼────────────────────┼────────┤
│ Mo17 (1st)  │ 17430-    │ 35460 bp    │ chr5:156,805,217-  │   -    │
│             │  35460    │             │      156,823,253   │        │
├─────────────┼───────────┼─────────────┼────────────────────┼────────┤
│ B73 (2nd)   │ 852-      │ 35460 bp    │ chr5:183,001,790-  │   -    │
│             │  12362    │             │      183,013,297   │        │
└─────────────┴───────────┴─────────────┴────────────────────┴────────┘
Type: SCO
```

**Benefits**:
1. ✅ Single row = easy to analyze
2. ✅ Both read and genome coordinates
3. ✅ Ordered by genome position (Mo17 @ 156M comes before B73 @ 183M)
4. ✅ Complete information in one place

---

## What Information Is Added?

| Information Type | Standard CHARLA | Post-Mapping | Improvement |
|-----------------|-----------------|--------------|-------------|
| **K-mer segments** | ✅ Yes (BED) | ✅ Yes | Same |
| **Read coordinates** | ✅ Yes | ✅ Yes | Same |
| **Genome coordinates** | ⚠️ In PAF only | ✅ **In main output** | **NEW** |
| **Single-row format** | ❌ No | ✅ **Yes** | **NEW** |
| **Biological ordering** | ❌ No | ✅ **Yes** | **NEW** |
| **Strand per parent** | ⚠️ Implicit | ✅ **Explicit** | **NEW** |
| **Read length** | ⚠️ Scattered | ✅ **Explicit** | **NEW** |
| **Event statistics** | ⚠️ Basic counts | ✅ **Comprehensive** | **NEW** |
| **Visualizations** | ❌ No | ✅ **Yes** | **NEW** |
| **Quality metrics** | ⚠️ Basic | ✅ **Detailed** | **NEW** |

---

## Data Flow Diagram

```
┌──────────────┐
│ CHARLA Input │  Long Reads (FASTQ)
└──────┬───────┘
       │
       ▼
┌──────────────────────────────────────────────────┐
│          CHARLA MAIN PIPELINE                    │
├──────────────────────────────────────────────────┤
│  • K-mer counting                                │
│  • Segmentation (find B73/Mo17 regions)          │
│  • Segment extraction                            │
│  • Alignment to reference genomes                │
│  • Classification (SCO/NCO/etc)                  │
└──────┬───────────────────────────────────────────┘
       │
       ├─────────────────┬─────────────────┐
       ▼                 ▼                 ▼
┌─────────────┐  ┌──────────────┐  ┌──────────────┐
│ BED files   │  │ PAF files    │  │ Read IDs     │
│ (segments)  │  │ (alignments) │  │ (classified) │
└─────────────┘  └──────────────┘  └──────────────┘
       │                 │                 │
       └─────────────────┴─────────────────┘
                         │
                         ▼
       ┌─────────────────────────────────────┐
       │   POST-MAPPING ANALYSIS MODULE      │
       │                                     │
       │  Script 19: Comprehensive Output   │◄──── YOU ARE HERE
       │    • Merge BED + PAF               │
       │    • Order by genome position      │
       │    • Single-row format             │
       │                                    │
       │  Scripts 14-18: Analysis & Viz     │
       │    • Statistics                    │
       │    • Karyoplots                    │
       │    • Quality metrics               │
       └──────────┬─────────────────────────┘
                  │
                  ▼
       ┌───────────────────────────┐
       │  PUBLICATION-READY OUTPUT │
       ├───────────────────────────┤
       │  • Single-read TSV        │
       │  • Statistics tables      │
       │  • Karyoplots (PDF/PNG)   │
       │  • Quality reports        │
       └───────────────────────────┘
```

---

## Real Numbers: What Was Achieved

From maize pollen analysis (9,277 confident crossover reads):

### Data Completeness
- **98.7%** of reads have complete alignment data
- **1.3%** have missing data (segments too short to align)

### Biological Insights
- **75.3%** have clean, non-overlapping segments
- **97.3%** have both parents on same chromosome
- **18.8%** show **inversions** (structural variation between B73/Mo17!)

### Quality Metrics
- Mean segment size: **~4 kb** (successful alignments)
- Failed segments: **5-28× smaller** (too short to align)
- Both strands represented: minus strand dominant in this dataset

---

## Summary

| Aspect | Before (Standard) | After (Post-Mapping) |
|--------|------------------|---------------------|
| **Data format** | Multi-row BED + separate PAFs | Single-row comprehensive TSV |
| **Coordinates** | Read only | Read + Genome |
| **Ordering** | By read position | By genome position |
| **Analysis ready** | Requires merging | Immediate |
| **Collaboration** | Custom format | Matches collaborator tools |
| **Visualization** | External tools needed | Built-in plots |
| **Statistics** | Manual calculation | Automated reports |
| **Quality control** | Limited metrics | Comprehensive QC |

**Bottom line**: Post-mapping analysis transforms CHARLA's raw output into analysis-ready, publication-quality results with biological interpretation.
