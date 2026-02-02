# Marker Balance Analysis Suite - Complete Guide

Two powerful scripts for analyzing and visualizing cenhapmer marker databases with publication-ready figures.

---

## 📊 Script 1: Single K-mer Analysis
**Script**: `16-analyze_marker_balance.py`
**Purpose**: In-depth analysis of a single k-mer size

### Features
- Comprehensive marker count and density analysis
- Parent balance assessment (Col-0 vs Ler-0)
- Region-specific analysis (ARMS vs CEN)
- Chromosome-specific breakdown
- Statistical annotations (t-tests, p-values)
- Publication-ready figures

### Quick Start
```bash
# Basic analysis
python scripts/16-analyze_marker_balance.py \
    --cenhapmer-dir /path/to/cenhapmers/k31 \
    --output-prefix marker_analysis/k31_balance

# Full publication-ready analysis
python scripts/16-analyze_marker_balance.py \
    --cenhapmer-dir /path/to/cenhapmers/k31 \
    --output-prefix marker_analysis/k31_balance \
    --col-bed index/Col-0.bed \
    --ler-bed index/Ler-0.bed \
    --publication \
    --format all
```

### Outputs
- `*_marker_counts.tsv` - Raw data
- `*_balance_summary.txt` - Statistics
- `*_marker_distribution.png` - Overview plots
- `*_density_analysis_publication.[png/pdf/svg]` - Main figure (4-panel)
- `*_genome_landscape.[png/pdf/svg]` - Genome-wide map
- `*_density_heatmap.[png/pdf/svg]` - Statistical heatmap
- `*_normalization_recommendations.txt` - Guidance

---

## 🔬 Script 2: K-mer Size Comparison
**Script**: `17-compare_kmer_sizes.py`
**Purpose**: Compare marker databases across multiple k-mer sizes

### Features
- Auto-detection of all k-mer sizes
- Side-by-side comparison of 5+ k-mer sizes
- Trend analysis (markers vs k-mer size)
- Balance comparison across k-mer sizes
- Centromeric resolution analysis
- K-mer selection recommendations

### Quick Start
```bash
# Compare all k-mer sizes found in directory
python scripts/17-compare_kmer_sizes.py \
    --cenhapmer-dir /path/to/cenhapmers \
    --output-prefix marker_analysis/kmer_comparison \
    --col-bed index/Col-0.bed \
    --ler-bed index/Ler-0.bed \
    --format all

# Compare specific k-mer sizes only
python scripts/17-compare_kmer_sizes.py \
    --cenhapmer-dir /path/to/cenhapmers \
    --kmer-sizes 21 31 41 \
    --output-prefix marker_analysis/kmer_comparison \
    --col-bed index/Col-0.bed \
    --ler-bed index/Ler-0.bed
```

### Outputs
- `*_all_data.tsv` - Combined raw data
- `*_kmer_summary.txt` - Text summary
- `*_kmer_comparison.[png/pdf/svg]` - Main figure (6-panel)
- `*_kmer_detailed.[png/pdf/svg]` - Detailed analysis

---

## 📋 Command-Line Options

### Common Options (Both Scripts)

| Option | Description | Required |
|--------|-------------|----------|
| `--output-prefix` | Prefix for output files | Yes |
| `--col-bed` | Col-0 BED file (for density) | No* |
| `--ler-bed` | Ler-0 BED file (for density) | No* |
| `--format` | Output format: `png`, `pdf`, `svg`, `all` | No (default: png) |

\* Required for density calculations and publication figures

### Script 1 Specific

| Option | Description |
|--------|-------------|
| `--cenhapmer-dir` | Directory with KMC databases (e.g., `k31/`) |
| `--publication` | Generate publication-ready figures |

### Script 2 Specific

| Option | Description |
|--------|-------------|
| `--cenhapmer-dir` | Base directory with k-mer subdirectories |
| `--kmer-sizes` | Specific k-mer sizes to compare (optional) |

---

## 🎨 Figure Formats

### PNG (Raster)
- **Resolution**: 300 DPI
- **Best for**: Presentations, PowerPoint, web
- **File size**: Medium (200-600 KB)

### PDF (Vector)
- **Resolution**: Infinite (scalable)
- **Best for**: Journal submission, printing
- **File size**: Small (20-50 KB)

### SVG (Vector)
- **Resolution**: Infinite (scalable)
- **Best for**: Further editing in Illustrator/Inkscape
- **File size**: Medium (50-200 KB)

---

## 📊 Example Workflow

### Step 1: Compare All K-mer Sizes
```bash
python scripts/17-compare_kmer_sizes.py \
    --cenhapmer-dir /home/user/cenhapmers \
    --output-prefix results/kmer_comparison \
    --col-bed index/Col-0.bed \
    --ler-bed index/Ler-0.bed \
    --format pdf
```

**Review outputs:**
- Open `kmer_comparison_kmer_comparison.pdf`
- Check which k-mer size has best balance
- Note total marker counts

### Step 2: Detailed Analysis of Chosen K-mer
```bash
# Based on comparison, analyze k=31 in detail
python scripts/16-analyze_marker_balance.py \
    --cenhapmer-dir /home/user/cenhapmers/k31 \
    --output-prefix results/k31_detailed \
    --col-bed index/Col-0.bed \
    --ler-bed index/Ler-0.bed \
    --publication \
    --format all
```

**Review outputs:**
- Main figure: `k31_detailed_density_analysis_publication.pdf`
- Genome map: `k31_detailed_genome_landscape.pdf`
- Heatmap: `k31_detailed_density_heatmap.pdf`
- Stats: `k31_detailed_balance_summary.txt`

### Step 3: Use in Publication
- Include k-mer comparison in Methods or Supplementary
- Use detailed analysis as main marker quality figure
- Cite both analyses in manuscript

---

## 🔍 Interpreting Results

### Parent Balance
- **< 10% imbalance**: ✅ Excellent, no action needed
- **10-20% imbalance**: ⚡ Moderate, consider normalization
- **> 20% imbalance**: ⚠️ High, normalization recommended

### Marker Density
- **Normalized by region size** (markers per Mb)
- **Higher values** = more markers available for crossover detection
- **Compare ARMS vs CEN**: Expect 3-5× difference (CEN has fewer markers)

### Chromosome Variation (CV)
- **< 15%**: ✅ Good consistency across chromosomes
- **15-30%**: ⚡ Moderate variation
- **> 30%**: ⚠️ High variation, check for issues

### K-mer Size Selection
- **Low k (21-25)**: Fewer markers, more error-tolerant → ONT/noisy reads
- **Medium k (31)**: Balanced → Default recommendation
- **High k (35-41)**: Most markers, less error-tolerant → PacBio HiFi/Illumina

---

## 🎯 Use Cases

### Use Case 1: First-Time Analysis
**Goal**: Understand marker database quality

```bash
# Quick check with k=31 (common default)
python scripts/16-analyze_marker_balance.py \
    --cenhapmer-dir cenhapmers/k31 \
    --output-prefix qc/k31_check \
    --col-bed index/Col-0.bed \
    --ler-bed index/Ler-0.bed
```

**Check**: Is parent imbalance < 10%? → Proceed with analysis

---

### Use Case 2: Optimizing K-mer Size
**Goal**: Choose best k-mer size for your data

```bash
# Compare all available k-mer sizes
python scripts/17-compare_kmer_sizes.py \
    --cenhapmer-dir cenhapmers \
    --output-prefix analysis/kmer_optimization \
    --col-bed index/Col-0.bed \
    --ler-bed index/Ler-0.bed \
    --format pdf
```

**Decision criteria**:
- ONT R9/R10 (Q10-15) → Choose k=21 or k=25
- ONT R10.4+ (Q15-20) → Choose k=31
- PacBio HiFi (Q20+) → Choose k=35 or k=41
- Unsure → Choose k=31 (safe default)

---

### Use Case 3: Publication Figures
**Goal**: Generate publication-ready figures

```bash
# Generate all formats for chosen k-mer size
python scripts/16-analyze_marker_balance.py \
    --cenhapmer-dir cenhapmers/k31 \
    --output-prefix publication/marker_analysis \
    --col-bed index/Col-0.bed \
    --ler-bed index/Ler-0.bed \
    --publication \
    --format all

# Also generate k-mer comparison for supplementary
python scripts/17-compare_kmer_sizes.py \
    --cenhapmer-dir cenhapmers \
    --output-prefix publication/kmer_comparison \
    --col-bed index/Col-0.bed \
    --ler-bed index/Ler-0.bed \
    --format all
```

**Use in manuscript**:
- **Main figure**: `marker_analysis_density_analysis_publication.pdf`
- **Supplementary**: `kmer_comparison_kmer_comparison.pdf`

---

### Use Case 4: Troubleshooting Issues
**Goal**: Diagnose why crossover calls seem biased

```bash
# Check marker balance
python scripts/16-analyze_marker_balance.py \
    --cenhapmer-dir cenhapmers/k31 \
    --output-prefix troubleshoot/balance_check \
    --col-bed index/Col-0.bed \
    --ler-bed index/Ler-0.bed

# Review normalization recommendations
cat troubleshoot/balance_check_normalization_recommendations.txt
```

**Common issues**:
- High parent imbalance → May bias crossover calls toward one parent
- High chromosome variation → Some chromosomes poorly covered
- Low CEN markers → Centromeric crossovers less confident

---

## 🎨 Color Scheme (Colorblind-Safe)

Both scripts use a consistent, publication-ready color palette:

- **Col-0**: Blue (#2C7BB6)
- **Ler-0**: Red (#D7191C)
- **ARMS**: Teal (#66C2A5)
- **CEN**: Orange (#FC8D62)

These colors are:
✅ Distinguishable in grayscale
✅ Colorblind-friendly (deuteranopia/protanopia)
✅ Professional for publication

---

## 📦 Dependencies

### Required Python Packages
```bash
pip install pandas numpy matplotlib seaborn scipy
```

### Required External Tools
- **KMC** (v3.2.1+): For k-mer database operations
  - `kmc_tools` must be in PATH

### Check Installation
```bash
# Check KMC
kmc_tools --version

# Check Python packages
python -c "import pandas, numpy, matplotlib, seaborn, scipy; print('All packages OK')"
```

---

## 🐛 Troubleshooting

### "No KMC databases found"
**Solution**: Check that directory contains `.kmc_pre` and `.kmc_suf` files

```bash
ls -la /path/to/cenhapmers/k31/*.kmc_pre
```

### "Marker density not calculated"
**Solution**: Provide both BED files with `--col-bed` and `--ler-bed`

### Pandas FutureWarning messages
**Status**: Harmless warnings, already fixed in latest version

### Missing fonts in PDF
**Solution**: Install Microsoft Core Fonts
```bash
sudo apt install ttf-mscorefonts-installer
```

### Low-resolution figures
**Solution**: Use `--format pdf` or `--format svg` for vector graphics

---

## 📚 Output File Reference

### Single K-mer Analysis (Script 1)

| File | Type | Description | Size |
|------|------|-------------|------|
| `*_marker_counts.tsv` | Data | Raw counts with density | ~3-5 KB |
| `*_balance_summary.txt` | Report | Statistical summary | ~1 KB |
| `*_marker_distribution.png` | Figure | 4-panel overview | ~350 KB |
| `*_per_database.png` | Figure | Per-database breakdown | ~280 KB |
| `*_density_analysis_publication.[fmt]` | Figure | Main publication figure | PNG:370KB, PDF:31KB |
| `*_genome_landscape.[fmt]` | Figure | Genome-wide map | PNG:171KB, PDF:19KB |
| `*_density_heatmap.[fmt]` | Figure | Statistical heatmap | PNG:234KB, PDF:24KB |
| `*_normalization_recommendations.txt` | Report | Recommendations | ~1 KB |

### K-mer Comparison (Script 2)

| File | Type | Description | Size |
|------|------|-------------|------|
| `*_all_data.tsv` | Data | Combined k-mer data | ~6 KB |
| `*_kmer_summary.txt` | Report | Text summary | ~3 KB |
| `*_kmer_comparison.[fmt]` | Figure | 6-panel comparison | PNG:629KB, PDF:30KB |
| `*_kmer_detailed.[fmt]` | Figure | Detailed analysis | PNG:985KB, PDF:46KB |

---

## 🚀 Performance Tips

### Fast Analysis (Skip Density)
```bash
# Omit BED files for faster analysis (no density calculation)
python scripts/16-analyze_marker_balance.py \
    --cenhapmer-dir cenhapmers/k31 \
    --output-prefix quick/k31
```

### Multiple K-mer Sizes in Parallel
```bash
# Run analyses in parallel (if you have multiple cores)
for k in 21 25 31 35 41; do
    python scripts/16-analyze_marker_balance.py \
        --cenhapmer-dir cenhapmers/k${k} \
        --output-prefix results/k${k}_balance \
        --col-bed index/Col-0.bed \
        --ler-bed index/Ler-0.bed &
done
wait
```

### Generate Only PDFs (Smallest Files)
```bash
python scripts/17-compare_kmer_sizes.py \
    --cenhapmer-dir cenhapmers \
    --output-prefix results/comparison \
    --col-bed index/Col-0.bed \
    --ler-bed index/Ler-0.bed \
    --format pdf
```

---

## 📖 Citing

If you use these visualization tools in your publication, please cite:

1. **CHARLA**: [Your CHARLA publication]
2. **KMC**: Kokot, M., Długosz, M., & Deorowicz, S. (2017). KMC 3: counting and manipulating k-mer statistics. Bioinformatics, 33(17), 2759-2761.

---

## 💡 Tips for Publication

### Main Text Figure
Use the 4-panel density analysis as your main marker quality figure:
```
Figure X: Cenhapmer marker database characterization
(see marker_analysis_density_analysis_publication.pdf)
```

### Supplementary Material
Include:
1. K-mer size comparison (if you compared multiple)
2. Genome-wide landscape
3. Statistical heatmap
4. Raw count tables (TSV files)

### Methods Section Text (Template)
```
Marker database quality was assessed using custom visualization scripts
that calculated marker counts, density (markers per Mb), and parent balance.
All k-mer sizes tested (k=21, 25, 31, 35, 41) showed excellent parent
balance (< 1% imbalance) and consistent chromosome coverage (CV < 15%).
K-mer size k=31 was selected for crossover analysis based on balanced
performance between marker density (176K markers/Mb) and error tolerance.
```

---

## 🔗 Related Documentation

- **CHARLA Manual**: See main repository README
- **KMC Documentation**: https://github.com/refresh-bio/KMC
- **Cenhapmer Generation**: See `cenhapmer_generation/README.md`

---

## 📞 Support

For issues or questions:
1. Check this README and individual script help:
   ```bash
   python scripts/16-analyze_marker_balance.py --help
   python scripts/17-compare_kmer_sizes.py --help
   ```
2. Review output summary files (SUMMARY.md, KMER_COMPARISON_SUMMARY.md)
3. Check troubleshooting section above

---

**Last Updated**: November 26, 2025
**Version**: 2.0 (Publication-ready)
**Author**: Enhanced by Claude Code
