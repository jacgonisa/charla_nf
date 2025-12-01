# Quick Start: Insertion Origin Analysis

## What This Does

Tests if **homologous recombination creates sequence expansions** by checking whether insertions contain k-mers from different genomic regions.

**The Question**: Can HR copy sequences from Chromosome 2 and insert them into Chromosome 1?

**The Method**: Profile insertion sequences against all cenhapmer databases to find their origin.

---

## Quick Usage

```bash
# After your non-hybrid mapping completes, run:
python scripts/14-analyze_insertion_origins.py \
    --bam-dir results/non_hybrid_mapping/sample_nonhybrid/mappings \
    --cenhapmer-dir /home/jg2070/Desktop/PhD/crossover/03-cenhapmers/k31 \
    --kmer-size 31 \
    --output-prefix results/insertion_origins/analysis
```

---

## What You'll Get

### 1. Ectopic Events Table (`*_ectopic_events.tsv`)
Shows insertions where origin ≠ mapped location:

```
read_id                 size  mapped_chr  origin_chr  is_inter_chr
read_12345_insertion_1  342   Chr1        Chr3        TRUE
read_67890_insertion_2  178   Chr2        Chr2        FALSE
```

**Inter-chr = TRUE** → Evidence of ectopic recombination!

### 2. Origin Heatmap (`*_origin_heatmap.png`)
- **Diagonal** (blue boxes): Normal local expansions
- **Off-diagonal**: Ectopic recombination events!

### 3. Summary Figure (`*_ectopic_summary.png`)
Four panels showing:
- Counts by event type
- Size distribution
- Chromosome flow matrix
- Region transitions (CEN ↔ ARMS)

---

## Interpretation

### If You See Off-Diagonal Signals:
```
🎯 SUCCESS! You found evidence of ectopic recombination!

Example: Insertion in Chr1_ARMS contains Chr3_CEN k-mers
→ Homologous recombination copied sequence from Chr3 to Chr1
→ Direct molecular evidence of HR-mediated duplication
```

### If Everything is Diagonal:
```
✅ Normal pattern - most insertions are local expansions
📊 Ectopic recombination may be rare in your sample
🔬 Try different k-mer sizes or look at larger insertions
```

### If You See CEN → ARMS:
```
🌟 Centromeric satellites expanding into chromosome arms!
🌟 Evidence of satellite repeat mobility
🌟 Look for 178bp multiples (Arabidopsis satellite)
```

---

## Expected Runtime

- **Per BAM file**: 2-10 minutes
- **Per insertion**: 1-2 seconds
- **Total**: Depends on number of insertions found

Progress is shown in real-time.

---

## Requirements

- ✅ BAM files from non-hybrid mapping
- ✅ Cenhapmer databases (same k-mer size)
- ✅ KMC tools in PATH
- ✅ Python packages: pysam, pandas, matplotlib, seaborn

---

## Quick Checks

### Check Your BAM Files
```bash
ls results/non_hybrid_mapping/*/mappings/*.bam
# Should show: CA1_*.bam, CC1_*.bam, LA1_*.bam, etc.
```

### Check Your Cenhapmers
```bash
ls /path/to/cenhapmers/k31/*.kmc_pre
# Should show many databases (CA1, CC1, etc.)
```

### Test Run (Quick)
```bash
# Test with just one BAM file
python scripts/14-analyze_insertion_origins.py \
    --bam-dir results/non_hybrid_mapping/*/mappings \
    --cenhapmer-dir /path/to/cenhapmers/k31 \
    --output-prefix test/quick_test \
    --min-size 100
# Higher min-size = fewer insertions = faster test
```

---

## Understanding the Results

### Ectopic Event Count
```
ECTOPIC RECOMBINATION EVENTS DETECTED: 15

This means: 15 insertions have k-mer signatures from different locations!
```

### Event Types
```
Inter-chromosomal: 8    → Different chromosome (Chr1 → Chr2)
Inter-regional: 7       → Different region (ARMS → CEN or vice versa)
```

### Biological Significance
- **High ectopic rate**: Active recombination, structural instability
- **Low ectopic rate**: Stable genome, normal recombination
- **CEN → ARMS**: Satellite repeat mobility
- **ARMS → CEN**: Gene capture into centromere (rare!)

---

## Troubleshooting

### "No insertions found"
→ Your reads don't have large insertions (>50bp)
→ Try `--min-size 30` to lower threshold

### "No ectopic events"
→ This is actually common! Most insertions are local
→ Try different k-mer sizes: k21, k31, k41

### "KMC error"
→ Check KMC is installed: `which kmc_tools`
→ Check cenhapmer databases exist and are readable

---

## Next Steps

1. **Run the analysis** (see command above)
2. **Check ectopic events table** for interesting cases
3. **View heatmap** for overall patterns
4. **Look for 178bp patterns** in Arabidopsis (satellite repeat)
5. **Correlate with crossover calls** if available

---

## Advanced: Compare K-mer Sizes

```bash
# Try different k-mer sizes
for k in 21 31 41; do
    python scripts/14-analyze_insertion_origins.py \
        --bam-dir results/non_hybrid_mapping/*/mappings \
        --cenhapmer-dir /path/to/cenhapmers/k${k} \
        --kmer-size $k \
        --output-prefix results/insertion_origins/k${k}_analysis
done
```

Different k-mer sizes may reveal different patterns:
- **k=21**: More sensitive, may have false positives
- **k=31**: Balanced sensitivity and specificity
- **k=41**: Very specific, fewer false positives

---

## Publications

This analysis provides empirical evidence for:
- Homologous recombination-mediated duplications
- Ectopic recombination mechanisms
- Centromeric satellite mobility
- Structural variation generation

**Use in your paper:**
> "We profiled insertion sequences against parental cenhapmer databases
> to detect ectopic recombination events, revealing X insertions (Y%)
> with k-mer signatures from non-homologous genomic regions, providing
> direct molecular evidence of recombination-mediated sequence expansions."

---

## Questions?

- **Full documentation**: `scripts/README_INSERTION_ORIGINS.md`
- **Script source**: `scripts/14-analyze_insertion_origins.py`
- **Pipeline integration**: Ask about adding to Nextflow workflow

---

**This is novel analysis!** 🚀 You're testing a fundamental question about how recombination creates structural variation. Good luck!
