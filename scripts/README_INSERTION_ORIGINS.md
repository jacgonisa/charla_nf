# Insertion Origin Analysis - Detecting Ectopic Recombination

**Script**: `14-analyze_insertion_origins.py`
**Purpose**: Test the hypothesis that homologous recombination creates sequence expansions

---

## The Biological Question

**Can homologous recombination create expansions of sequences from different genomic regions?**

This script provides empirical evidence by:
1. Extracting large insertion sequences (>50bp) from alignments
2. Profiling each insertion against ALL cenhapmer databases
3. Detecting when an insertion's k-mer signature doesn't match its mapped location
4. Revealing ectopic recombination events

---

## The Key Insight

If you find an insertion in **Chromosome 1 ARMS** that contains k-mers specific to **Chromosome 2 CEN**, this suggests:

✨ **Homologous recombination copied sequence from Chr2 CEN and inserted it into Chr1 ARMS!**

This would be direct molecular evidence of:
- Recombination-mediated duplications
- Ectopic (non-homologous) recombination
- Centromeric satellite repeat mobility
- Mechanisms of structural variation generation

---

## How It Works

### Step 1: Extract Insertions
```
For each BAM file:
  - Find all insertions >50bp (from CIGAR strings)
  - Extract the inserted sequence
  - Record the mapped location (chromosome, region)
```

### Step 2: Profile Against Cenhapmers
```
For each insertion:
  - Generate k-mers from insertion sequence
  - Compare against ALL cenhapmer databases
    (CA1, CA2, ..., CC1, CC2, ..., LA1, ..., LC5)
  - Count matching k-mers for each haplotype
  - Identify the best-matching origin
```

### Step 3: Detect Ectopic Events
```
If insertion's origin ≠ mapped location:
  → ECTOPIC RECOMBINATION EVENT!

Examples:
  - Mapped to Chr1_ARMS, origin Chr2_CEN → Inter-chromosomal
  - Mapped to Chr1_ARMS, origin Chr1_CEN → Inter-regional expansion
  - Mapped to Chr1_ARMS, origin Chr1_ARMS → Local expansion (control)
```

---

## Usage

### Basic Command
```bash
python scripts/14-analyze_insertion_origins.py \
    --bam-dir results/non_hybrid_mapping/sample_nonhybrid/mappings \
    --cenhapmer-dir /path/to/cenhapmers/k31 \
    --kmer-size 31 \
    --output-prefix results/insertion_origins/analysis
```

### Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--bam-dir` | Directory with BAM files from non-hybrid mapping | Required |
| `--cenhapmer-dir` | Directory with cenhapmer KMC databases (k21/, k31/, etc.) | Required |
| `--kmer-size` | K-mer size matching cenhapmer database | 31 |
| `--min-size` | Minimum insertion size to analyze (bp) | 50 |
| `--output-prefix` | Prefix for output files | Required |

---

## Outputs

### 1. All Insertions Table
**File**: `*_all_insertions.tsv`

Columns:
- `haplotype`: Where insertion was found (CA1, CC2, etc.)
- `read_id`: Read identifier
- `ref_chr`, `ref_pos`: Reference position
- `size`: Insertion size (bp)
- `origin`: Best-matching cenhapmer (by k-mer content)
- `origin_kmers`: Number of matching k-mers
- `is_ectopic`: TRUE if origin ≠ mapped location

### 2. Ectopic Events Table
**File**: `*_ectopic_events.tsv`

Columns:
- `read_id`, `size`
- `mapped_background`, `mapped_chr`, `mapped_region`: Where found
- `origin_background`, `origin_chr`, `origin_region`: Where it came from
- `origin_kmers`: Supporting k-mer evidence
- `is_inter_chr`: Different chromosome?
- `is_inter_region`: Different region (ARMS/CEN)?

### 3. Origin Heatmap
**Files**: `*_origin_heatmap.png/pdf`

Shows flow of insertions from origin → mapped location:
- **Diagonal** (blue boxes): Same location → Local expansions
- **Off-diagonal**: Ectopic recombination!
  - Different chromosomes → Inter-chromosomal
  - Different regions → Inter-regional (e.g., CEN → ARMS)

### 4. Ectopic Summary Figure
**Files**: `*_ectopic_summary.png/pdf`

4-panel figure showing:
- **Panel A**: Counts by event type (inter-chr, inter-region, etc.)
- **Panel B**: Size distribution of ectopic insertions
- **Panel C**: Chromosome origin matrix (which Chr → which Chr)
- **Panel D**: Region transitions (CEN → ARMS, ARMS → CEN, etc.)

---

## Interpretation Guide

### Case 1: Diagonal Dominance
```
Origin Heatmap shows:
  - Strong diagonal (same location)
  - Weak off-diagonal

Interpretation:
  ✅ Most insertions are local expansions
  ✅ Low ectopic recombination rate
  ✅ Normal recombination patterns
```

### Case 2: Off-Diagonal Signals
```
Origin Heatmap shows:
  - Signals off the diagonal
  - E.g., Chr2_CEN → Chr1_ARMS

Interpretation:
  🔥 Ectopic recombination detected!
  🔥 Sequence from Chr2 CEN copied into Chr1 ARMS
  🔥 Evidence of non-homologous recombination
```

### Case 3: CEN → ARMS Enrichment
```
Region transitions show:
  - CEN → ARMS: High count
  - ARMS → CEN: Low count

Interpretation:
  🌟 Centromeric sequences expanding into arms
  🌟 Satellite repeat mobility
  🌟 Potential mechanism: BIR (break-induced replication)
```

---

## Biological Significance

### Evidence of Recombination Mechanisms

#### 1. **Homologous Recombination (HR)**
If insertions contain sequences from **homologous chromosomes**:
- Col-0 Chr1 → Col-0 Chr1: Normal HR
- Col-0 Chr1 → Ler-0 Chr1: Allelic HR

#### 2. **Ectopic Recombination**
If insertions contain sequences from **non-homologous regions**:
- Chr1 → Chr2: Inter-chromosomal ectopic
- ARMS → CEN: Satellite repeat capture
- Evidence of promiscuous recombination

#### 3. **Break-Induced Replication (BIR)**
If large insertions (>500bp) from distant regions:
- Template switching during repair
- Can copy kilobases from donor
- Mechanism of complex rearrangements

#### 4. **Satellite Repeat Mobility**
If CEN sequences inserted into ARMS:
- 178bp satellite repeats expanding
- Centromeric instability
- Mechanism of karyotype evolution

---

## Expected Patterns in F1 Hybrids

### Homologous Chromosomes
```
Col-0 Chr1 ←→ Ler-0 Chr1
  - High recombination
  - Normal crossovers
  - Allelic exchanges
```

### Non-Homologous Chromosomes
```
Col-0 Chr1 ←→ Ler-0 Chr2
  - Rare recombination
  - Ectopic events
  - Structural rearrangements
```

### Centromeric Regions
```
CEN → ARMS expansions
  - Satellite repeat mobility
  - 178bp motif enrichment
  - Centromere instability
```

---

## Statistical Analysis Tips

### Null Hypothesis
- Insertions should match their mapped location
- Origin = Mapped location
- Diagonal dominance expected

### Alternative Hypothesis
- Ectopic recombination creates off-diagonal signals
- Insertions from different chromosomes/regions
- Evidence of sequence duplication/mobility

### Statistical Tests
```R
# Chi-square test for random expectation
observed_matrix <- read_heatmap_data()
expected <- uniform_distribution()
chisq.test(observed, p=expected)

# Fisher's exact test for inter-chr vs intra-chr
contingency_table <- create_2x2_table()
fisher.test(contingency_table)
```

---

## Integration with Pipeline

### Option 1: Run Standalone
```bash
# After non-hybrid mapping completes
python scripts/14-analyze_insertion_origins.py \
    --bam-dir results/non_hybrid_mapping/*/mappings \
    --cenhapmer-dir 03-cenhapmers/k31 \
    --output-prefix results/insertion_origins/k31
```

### Option 2: Add to Nextflow Pipeline
```groovy
process INSERTION_ORIGIN_ANALYSIS {
    input:
    path(mapping_dir)
    path(cenhapmer_dir)
    val(kmer_size)

    output:
    path("*_all_insertions.tsv")
    path("*_ectopic_events.tsv")
    path("*_origin_heatmap.png")
    path("*_ectopic_summary.png")

    script:
    """
    python ${projectDir}/scripts/14-analyze_insertion_origins.py \\
        --bam-dir ${mapping_dir} \\
        --cenhapmer-dir ${cenhapmer_dir} \\
        --kmer-size ${kmer_size} \\
        --output-prefix insertion_origins
    """
}
```

---

## Performance Considerations

### Speed
- **Per insertion**: ~1-2 seconds (KMC operations)
- **100 insertions**: ~2-5 minutes
- **Bottleneck**: KMC database queries

### Memory
- **KMC operations**: Low (~100MB per operation)
- **Temp files**: Cleaned automatically
- **Total**: <1GB RAM typically

### Optimization
```bash
# Process in parallel if many samples
for bam in *.bam; do
    python scripts/14-analyze_insertion_origins.py \
        --bam-dir $(dirname $bam) \
        --cenhapmer-dir cenhapmers/k31 \
        --output-prefix ${bam%.bam}_origins &
done
wait
```

---

## Troubleshooting

### "No insertions found"
- Check BAM files have alignments
- Lower `--min-size` threshold
- Verify CIGAR strings contain 'I' operations

### "No cenhapmer databases found"
- Check `--cenhapmer-dir` path
- Verify `.kmc_pre` and `.kmc_suf` files exist
- Ensure k-mer size matches database

### "KMC timeout"
- Some insertions may be very long
- Script skips problematic insertions
- Check logs for warnings

### "All insertions match mapped location"
- This is actually expected for many samples!
- Ectopic events may be rare
- Try different k-mer sizes (k21 vs k31)

---

## Example Results

### Sample Output
```
ECTOPIC RECOMBINATION EVENTS DETECTED: 15

Event Types:
  Inter-chromosomal: 8
  Inter-regional: 7

Examples:
  1. Chr1_ARMS insertion with Chr3_CEN signature (342bp)
     → Centromeric repeat expanded into Chr1

  2. Chr2_CEN insertion with Chr5_ARMS signature (156bp)
     → Gene sequence captured into centromere

  3. Chr4_ARMS insertion with Chr4_CEN signature (534bp)
     → Local satellite repeat expansion (3× 178bp)
```

### Biological Interpretation
```
✅ Evidence of homologous recombination-mediated expansions
✅ Centromeric satellites are mobile (CEN → ARMS)
✅ Ectopic recombination occurs between non-homologous chromosomes
✅ 178bp satellite multiples enriched (recombination hotspot)
```

---

## Citation & References

This analysis implements concepts from:

1. **Ectopic recombination**: Lichten & Haber (1989) Annu Rev Genet
2. **Break-induced replication**: Anand et al. (2013) Cell
3. **Satellite repeat dynamics**: Malik & Henikoff (2009) Nat Rev Genet

If you use this analysis, please cite CHARLA and mention this insertion origin profiling method.

---

## Future Enhancements

Potential additions:
- [ ] Sequence alignment of ectopic regions
- [ ] Motif enrichment in insertions
- [ ] Breakpoint microhomology analysis
- [ ] Integration with crossover calls
- [ ] Temporal dynamics in recombinant populations

---

**Questions?** The script is fully documented with inline comments. Check `scripts/14-analyze_insertion_origins.py` for implementation details.

**This is cutting-edge analysis!** 🚀 You're empirically testing whether HR creates sequence expansions - very cool science!
