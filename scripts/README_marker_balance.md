# Marker Balance Analysis - Publication-Ready Figures

## Overview

The `16-analyze_marker_balance.py` script has been enhanced to generate publication-ready figures for marker database analysis. This script analyzes cenhapmer (marker) databases to identify potential biases that could affect crossover detection.

## Key Enhancements

### 1. **Publication-Ready Plot Aesthetics**
- Professional fonts (Arial/DejaVu Sans) with TrueType encoding for PDF compatibility
- Consistent color palette optimized for print and colorblind accessibility
- Clean layouts with removed chartjunk (top/right spines)
- Proper font sizes and weights for all plot elements
- High-resolution outputs (300 DPI for PNG)

### 2. **Statistical Annotations**
- Automated t-tests comparing marker densities between parents
- Significance bars with p-value annotations (*, **, ***, ns)
- Effect size calculations

### 3. **Comprehensive Visualizations**

#### **Figure 1: 4-Panel Density Analysis** (`*_density_analysis_publication.*`)
- **Panel A**: Overall marker density by parent (with significance testing)
- **Panel B**: Marker density by region type (ARMS vs CEN)
- **Panel C**: Chromosome-specific marker density comparison
- **Panel D**: Violin plots showing density distributions

#### **Figure 2: Genome-Wide Density Landscape** (`*_genome_landscape.*`)
- Horizontal bar chart showing density across the entire genome
- Color-coded by region type (ARMS/CEN)
- Chromosome boundaries clearly marked
- Density values embedded in each region

#### **Figure 3: Statistical Heatmap** (`*_density_heatmap.*`)
- Parent × Region × Chromosome breakdown
- Annotated with actual density values
- Easy identification of imbalances

### 4. **Multiple Export Formats**
- **PNG**: High-resolution (300 DPI) for presentations and web
- **PDF**: Vector format with embedded TrueType fonts for publications
- **SVG**: Scalable vector graphics for further editing in Illustrator/Inkscape
- **All**: Generate all three formats simultaneously

## Usage

### Basic Analysis
```bash
python scripts/16-analyze_marker_balance.py \
    --cenhapmer-dir 03-cenhapmers/k31/ \
    --output-prefix marker_analysis/balance
```

**Outputs:**
- `balance_marker_counts.tsv` - Raw marker counts
- `balance_balance_summary.txt` - Statistical summary
- `balance_marker_distribution.png` - Basic overview plots
- `balance_normalization_recommendations.txt` - Guidance for addressing imbalances

### Full Publication-Ready Analysis
```bash
python scripts/16-analyze_marker_balance.py \
    --cenhapmer-dir 03-cenhapmers/k31/ \
    --output-prefix marker_analysis/balance \
    --col-bed index/Col-0_renamed.bed \
    --ler-bed index/Ler-0_renamed.bed \
    --publication \
    --format all
```

**Additional Outputs:**
- `balance_density_analysis_publication.[png/pdf/svg]` - 4-panel comprehensive figure
- `balance_genome_landscape.[png/pdf/svg]` - Genome-wide density map
- `balance_density_heatmap.[png/pdf/svg]` - Statistical heatmap

### Command-Line Options

| Option | Description | Required |
|--------|-------------|----------|
| `--cenhapmer-dir` | Directory containing KMC cenhapmer databases | Yes |
| `--output-prefix` | Prefix for all output files | Yes |
| `--col-bed` | Col-0 BED file with genomic regions (for density calculation) | No* |
| `--ler-bed` | Ler-0 BED file with genomic regions (for density calculation) | No* |
| `--publication` | Generate publication-ready figures | No |
| `--format` | Output format: `png`, `pdf`, `svg`, or `all` (default: `png`) | No |

\* **Note:** BED files are required for marker density calculations and publication-ready figures.

## Example Workflow

1. **Run the analysis with publication figures:**
   ```bash
   python scripts/16-analyze_marker_balance.py \
       --cenhapmer-dir 03-cenhapmers/k31/ \
       --output-prefix marker_analysis/k31_balance \
       --col-bed index/Col-0_renamed.bed \
       --ler-bed index/Ler-0_renamed.bed \
       --publication \
       --format pdf
   ```

2. **Review the outputs:**
   - Check `k31_balance_balance_summary.txt` for statistical metrics
   - Open `k31_balance_density_analysis_publication.pdf` for main figure
   - Review `k31_balance_normalization_recommendations.txt` for guidance

3. **For presentations (high-res PNG):**
   ```bash
   python scripts/16-analyze_marker_balance.py \
       --cenhapmer-dir 03-cenhapmers/k31/ \
       --output-prefix marker_analysis/k31_balance \
       --col-bed index/Col-0_renamed.bed \
       --ler-bed index/Ler-0_renamed.bed \
       --publication \
       --format png
   ```

4. **For further editing (SVG):**
   ```bash
   python scripts/16-analyze_marker_balance.py \
       --cenhapmer-dir 03-cenhapmers/k31/ \
       --output-prefix marker_analysis/k31_balance \
       --col-bed index/Col-0_renamed.bed \
       --ler-bed index/Ler-0_renamed.bed \
       --publication \
       --format svg
   ```

## Understanding the Outputs

### Marker Density
- **Markers/Mb**: Number of unique k-mers per megabase of genomic sequence
- **Normalized by region size**: Accounts for differences in chromosome arm vs centromere sizes
- **Key metric**: More meaningful than raw counts for comparing regions

### Statistical Tests
- **t-test**: Compares marker density distributions between parents
- **Significance levels**:
  - `***` p < 0.001 (highly significant)
  - `**` p < 0.01 (very significant)
  - `*` p < 0.05 (significant)
  - `ns` p ≥ 0.05 (not significant)

### Interpretation
- **Well-balanced**: Density ratio < 1.1 (< 10% difference)
- **Moderate imbalance**: Density ratio 1.1-1.2 (10-20% difference)
- **High imbalance**: Density ratio > 1.2 (> 20% difference)

## Color Scheme

The script uses a carefully selected color palette:
- **Col-0**: Blue (#2C7BB6)
- **Ler-0**: Red (#D7191C)
- **ARMS**: Teal (#66C2A5)
- **CEN**: Orange (#FC8D62)

These colors are:
- Distinguishable in black and white printing
- Colorblind-friendly (protanopia/deuteranopia accessible)
- Professional and suitable for publication

## Dependencies

```bash
# Python packages required
pip install pandas numpy matplotlib seaborn scipy
```

External tools (must be in PATH):
- `kmc_tools` (from KMC suite v3.2.1+)

## Tips for Publication

1. **Use PDF or SVG for journals** - Vector formats scale perfectly and meet journal requirements
2. **Use PNG for presentations** - Easier to embed in PowerPoint/Keynote
3. **Check journal requirements** - Some journals prefer specific formats or DPI settings
4. **Figure legends** - Add detailed legends in your manuscript explaining each panel
5. **Accessibility** - The chosen colors work well for colorblind readers

## Troubleshooting

### "No KMC databases found"
- Ensure `--cenhapmer-dir` points to the correct directory
- Check that `.kmc_pre` and `.kmc_suf` files exist
- Verify file permissions

### "Marker density not calculated"
- Provide both `--col-bed` and `--ler-bed` files
- Verify BED files exist and are formatted correctly
- Check BED file chromosome names match KMC database names

### Missing fonts in PDF
- The script uses TrueType fonts (type 42) for maximum compatibility
- If Arial is unavailable, it falls back to DejaVu Sans or Liberation Sans
- Install Microsoft Core Fonts if needed: `sudo apt install ttf-mscorefonts-installer`

### Low-resolution figures
- For PNG: Script automatically uses 300 DPI
- For higher resolution: Use PDF or SVG formats instead
- Modify `dpi=300` in the code if you need even higher resolution

## Citation

If you use these visualizations in your publication, please cite CHARLA and mention the marker balance analysis tool.

---

**Author**: Enhanced by Claude Code
**Date**: November 2025
**Version**: 2.0 (Publication-ready)
