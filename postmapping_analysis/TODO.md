# Post-Mapping Analysis Module - TODO

## Integration Status

This module is currently **standalone** but designed for integration into the main CHARLA Nextflow pipeline.

---

## High Priority

### 🔴 CRITICAL: Nextflow Pipeline Integration

**Goal**: Integrate all post-mapping analysis scripts into the main `charla_nf.nf` workflow

**Current Status**: ⚠️ **NOT INTEGRATED**
- Scripts exist and are tested ✅
- Nextflow modules created ✅
- Main workflow integration: ❌ **PENDING**

**Integration Plan**:

1. **Add optional parameter** to `nextflow.config`:
   ```groovy
   params {
       // ... existing parameters ...

       // Post-mapping analysis
       run_postmapping_analysis = false  // Optional, default: false
       create_comprehensive_output = true  // Always create for single-read analysis
   }
   ```

2. **Integrate into main workflow** (`workflows/charla_nf.nf`):
   ```groovy
   // After curation and alignment steps
   if (params.create_comprehensive_output) {
       CREATE_COMPREHENSIVE_ALIGNMENT_OUTPUT(
           curated_confident_reads,
           paf_b73,
           paf_mo17,
           sample_name
       )
   }

   if (params.run_postmapping_analysis) {
       CROSSOVER_STATS(comprehensive_alignment)
       KARYOPLOT_CROSSOVERS(comprehensive_alignment, reference_genome)
       PLOT_CONFIDENT_READS(comprehensive_alignment)
   }
   ```

3. **Update dependencies**:
   - Ensure PAF files are properly passed from mapping step
   - Ensure curated BED files are available
   - Add reference genome as input for karyoplots

4. **Testing**:
   - Test with `--run_postmapping_analysis false` (default)
   - Test with `--run_postmapping_analysis true`
   - Verify outputs are published correctly

**Estimated effort**: 2-4 hours

---

## Medium Priority

### 🟡 Documentation & Examples

- [ ] Add example output files to `examples/` directory
- [ ] Create visual diagrams showing data flow
- [ ] Add example plots/figures to README
- [ ] Document required R packages and versions
- [ ] Add troubleshooting guide

**Estimated effort**: 1-2 hours

---

### 🟡 Containerization

- [ ] Create Docker container with all R dependencies
- [ ] Test scripts in container environment
- [ ] Add Singularity support
- [ ] Document container usage

**Estimated effort**: 3-4 hours

---

### 🟡 Performance Optimization

- [ ] Parallelize statistics calculations
- [ ] Optimize memory usage for large datasets
- [ ] Add progress bars/logging
- [ ] Profile script performance

**Estimated effort**: 2-3 hours

---

## Low Priority (Future Enhancements)

### 🟢 Additional Analysis Features

- [ ] **Crossover interference analysis**
  - Calculate interference distances
  - Statistical testing for interference
  - Visualization of interference patterns

- [ ] **Gene annotation integration**
  - Annotate crossover positions with genes
  - Calculate crossover rates per gene/intergenic region
  - Identify crossovers in specific functional regions

- [ ] **Multi-sample comparison**
  - Compare crossover patterns across samples
  - Statistical tests for differences
  - Heatmaps/clustering visualizations

- [ ] **Interactive visualization**
  - Shiny app for exploring results
  - Plotly-based interactive plots
  - IGV integration for genome browsing

- [ ] **Machine learning classification**
  - Train models to classify crossover quality
  - Predict crossover type from features
  - Identify complex rearrangements

**Estimated effort**: 1-2 weeks total

---

### 🟢 Output Format Options

- [ ] Add JSON output option
- [ ] Add BED format for genome browsers
- [ ] Add GFF3 format for annotation
- [ ] Add BigWig for signal tracks
- [ ] Add VCF-like format for crossover "variants"

**Estimated effort**: 4-6 hours

---

### 🟢 Quality Control Enhancements

- [ ] Automated outlier detection
- [ ] Quality score calculation per read
- [ ] Confidence intervals for statistics
- [ ] Batch effect detection/correction

**Estimated effort**: 1 week

---

## Technical Debt

### 🔧 Code Improvements

- [ ] Add unit tests for R scripts
- [ ] Add integration tests with test data
- [ ] Refactor duplicated code into shared functions
- [ ] Add error handling and validation
- [ ] Add command-line help messages
- [ ] Improve logging/verbose output

**Estimated effort**: 1-2 days

---

### 🔧 Module Cleanup

- [ ] Standardize Nextflow module structure
- [ ] Add module tests
- [ ] Document module inputs/outputs
- [ ] Add module parameter validation

**Estimated effort**: 4-6 hours

---

## Dependencies to Document

### R Packages Required

**Core dependencies**:
```r
dplyr (>= 1.0.0)
tidyr (>= 1.0.0)
ggplot2 (>= 3.3.0)
optparse
```

**Visualization dependencies**:
```r
karyoploteR
BSgenome
regioneR
```

**Optional dependencies**:
```r
plotly        # Interactive plots
shiny         # Web interface
RColorBrewer  # Color schemes
viridis       # Color schemes
```

### External Tools

Currently: **None** (pure R implementation)

Future possibilities:
- BEDTools (for interval operations)
- SAMtools (for BAM file handling)
- IGV (for visualization)

---

## Testing Checklist

Before marking integration as complete:

- [ ] Test on small dataset (< 100 reads)
- [ ] Test on medium dataset (1,000-10,000 reads)
- [ ] Test on large dataset (> 50,000 reads)
- [ ] Test with different species/genome sizes
- [ ] Test with missing data (NAs)
- [ ] Test with edge cases (single chromosome, inversions only, etc.)
- [ ] Test memory usage profiling
- [ ] Test on HPC cluster environment
- [ ] Test with different reference genomes
- [ ] Test error handling (malformed inputs, etc.)

---

## Known Issues

### Current Limitations

1. **Karyoplots require BSgenome objects**
   - Not available for all species
   - May need to be generated manually
   - Documentation needed for custom genomes

2. **Large datasets may require significant memory**
   - Script 19 loads entire dataset into memory
   - May need chunking for very large samples

3. **No automated parameter tuning**
   - Users need to manually adjust thresholds
   - Could benefit from automatic optimization

4. **Limited error recovery**
   - Scripts fail completely on errors
   - Could benefit from graceful degradation

---

## Community Requests

Track user-requested features here:

- [ ] *No requests yet - module just released*

---

## Version History

### v0.1.0 (Current - 2025-02-02)
- Initial release
- Core scripts 14-19 implemented
- Nextflow modules created
- Standalone functionality working
- ⚠️ **NOT YET INTEGRATED** into main pipeline

### Future Releases

**v0.2.0** (Target: TBD)
- [ ] Full Nextflow integration
- [ ] Docker container
- [ ] Complete documentation
- [ ] Example datasets

**v1.0.0** (Target: TBD)
- [ ] Production-ready integration
- [ ] Comprehensive testing
- [ ] Publication-quality documentation
- [ ] Performance optimizations

---

## Contact

For questions about integration or to contribute:
- Open an issue on GitHub
- See main CHARLA documentation
- Contact repository maintainers

---

## Priority Summary

**Immediate (Next Session)**:
1. ✅ Create this TODO file
2. ✅ Add EXAMPLES.md with comparisons
3. 🔴 Integrate Script 19 into main workflow (CRITICAL)

**Short-term (Next Week)**:
1. 🟡 Full Nextflow integration
2. 🟡 Add example outputs
3. 🟡 Documentation improvements

**Long-term (Next Month)**:
1. 🟢 Additional analysis features
2. 🟢 Interactive visualization
3. 🟢 Multi-sample support

---

*Last updated: 2025-02-02*
