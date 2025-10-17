# CHARLA Test Suite Implementation Summary

## 🎉 Overview

A comprehensive testing infrastructure has been successfully implemented for the CHARLA Nextflow pipeline using nf-test. This document summarizes what has been created and how to use it.

## ✅ What Has Been Implemented

### 1. Test Data (`tests/data/`)

**Realistic Long-Read Test Data:**
- ✅ **test_sample_longreads.fastq** - 10 ONT/PacBio-style reads (1-10kb each, ~128KB total)
  - Extracted from real simulated chimeric reads
  - Includes realistic quality scores for masking tests
  - Average read length: ~6kb (appropriate for long-read sequencing)
  - Contains variable quality profiles

**Legacy Short-Read Test Data (for backwards compatibility):**
- test_sample.fastq - 5 short reads (~250bp each)
- test_sample.fasta - 5 sequences with complex headers

### 2. Module Tests (`modules/local/tests/`)

**Preprocessing Modules:**
- ✅ `fastq_to_fasta.nf.test` - 4 comprehensive tests
  - FASTQ→FASTA without masking
  - FASTQ→FASTA with quality-based masking
  - Compressed input handling
  - Stub run testing

- ✅ `simplify_headers.nf.test` - 5 tests
  - Header simplification with complex metadata
  - Already-simple headers
  - Sequence count preservation
  - Invalid input error handling
  - Stub run testing

**K-mer Profiling Modules:**
- ✅ `generate_readmers_kmc.nf.test` - 5 tests
  - K-mer database generation (k=21)
  - K-mer database generation (k=31)
  - Different sample ID handling
  - Invalid FASTA error handling
  - Long read processing

- ✅ `get_counts_kmc.nf.test` - 2 tests
  - Per-read k-mer counting
  - Different k-mer sizes (k=21, k=31)

- ✅ `generate_histogram_kmc.nf.test` - 2 tests
  - K-mer histogram generation
  - Different k-mer sizes

**Total: 18 unit tests across 5 critical modules**

### 3. Test Configuration Files

- ✅ `nf-test.config` - Updated with test data paths and profiles
- ✅ `tests/nextflow.config` - Test-specific Nextflow configuration
- ✅ `tests/.nftignore` - Files to ignore in snapshot testing

### 4. Continuous Integration

- ✅ `.github/workflows/ci.yml` - Comprehensive CI/CD workflow
  - Module tests on push/PR/schedule
  - Pipeline integration tests
  - Code quality checks (pre-commit, linting)
  - Documentation validation
  - Container build verification
  - Multi-version Nextflow testing (24.10.5, latest-stable)
  - Test result artifacts

- ✅ `.github/markdown-link-check-config.json` - Link validation config

### 5. Documentation

- ✅ `tests/README.md` - Comprehensive testing guide
  - How to run tests
  - Test categories explained
  - Writing new tests tutorial
  - Best practices
  - Troubleshooting guide

- ✅ `TESTING.md` - This summary document
  - Quick start guide
  - Test coverage overview
  - Development workflow
  - CI/CD integration

## 📊 Test Coverage

```
Module                          Tests    Status
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
FASTQ_TO_FASTA                    4      ✅ Complete
SIMPLIFY_HEADERS                  5      ✅ Complete
generate_readmers_kmc             5      ✅ Complete
get_counts_kmc                    2      ✅ Complete
generate_histogram_kmc            2      ✅ Complete
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
TOTAL                            18      ✅ Complete
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Modules Pending Tests (can be added as needed):
- get_counts_cenhapmer_kmc
- process_cenhapmer_counts
- COMBINE_HYBRID_PROFILES
- CURATE_HYBRID_PROFILES
- SEGMENT_READS
- MAPPING_ANALYSIS
- PLOT_CROSSOVER_MAP
- NON_HYBRID_READS_ANALYSIS
```

## 🚀 Quick Start

### Run All Tests

```bash
cd /home/jg2070/Desktop/PhD/crossover/charla_nf

# Run all module tests
nf-test test tests/modules/local/*.nf.test

# Run specific module test
nf-test test tests/modules/local/fastq_to_fasta.nf.test

# Verbose mode
nf-test test --verbose

# Update snapshots after intentional changes
nf-test test --update-snapshot
```

### Test with Different Profiles

```bash
# Docker profile (default in tests)
nf-test test --profile docker

# Singularity profile
nf-test test --profile singularity

# Local profile (no containers)
nf-test test --profile local
```

## 📝 Test Data Details

### Primary Test Data: Long Reads (ONT/PacBio-style)

**File:** `tests/data/test_sample_longreads.fastq`

- **Format:** FASTQ (ONT/PacBio long reads)
- **Reads:** 10 reads
- **Size:** 128 KB
- **Average length:** ~6 kb (1-10 kb range)
- **Quality scores:** Realistic ONT quality profile
- **Source:** Subset from `simulated_ont_chimeric_reads_1pc_chimeric_error_1pc.fq`

**Why this matters:**
- CHARLA is designed for long-read sequencing (ONT/PacBio)
- Quality-based masking requires FASTQ format with Q-scores
- K-mer analysis benefits from reads >1kb for haplotype phasing
- Tests reflect real-world usage

### Test Scenarios Covered

1. **Quality Masking** (FASTQ only)
   - Low-quality bases (Q<14) masked to 'N'
   - Improves k-mer specificity
   - Reduces false hybrid calls

2. **Header Simplification**
   - Complex ONT headers → simple IDs
   - Compatibility with KMC and downstream tools

3. **K-mer Profiling**
   - Both k=21 and k=31 tested
   - Readmer generation
   - Cenhapmer counting (pending)

## 🔄 CI/CD Pipeline

### Automated Testing

**Triggers:**
- Every push to `master` or `dev`
- Pull requests
- Weekly schedule (Monday 2 AM UTC)
- Manual dispatch

**Test Matrix:**
- **Container systems:** Docker
- **Nextflow versions:** 24.10.5, latest-stable
- **Test types:** Unit, Integration, End-to-end

**Jobs:**
1. **Module Tests** - Unit tests for individual processes
2. **Pipeline Tests** - Full integration tests
3. **Lint Checks** - Code quality validation
4. **Documentation** - Link checking, completeness
5. **Containers** - Build verification
6. **Summary** - Aggregate results

### View CI Results

```bash
# GitHub Actions
https://github.com/jacgonisa/charla_nf/actions
```

## 📚 Key Files Created

```
charla_nf/
├── .github/
│   ├── workflows/
│   │   └── ci.yml                                    # CI/CD pipeline
│   └── markdown-link-check-config.json               # Doc validation
│
├── tests/
│   ├── README.md                                     # Testing guide
│   ├── data/
│   │   ├── test_sample_longreads.fastq               # PRIMARY test data
│   │   ├── test_sample.fastq                         # Legacy short reads
│   │   └── test_sample.fasta                         # FASTA test mode
│   │
│   └── modules/local/                                # Module tests
│       ├── fastq_to_fasta.nf.test
│       ├── simplify_headers.nf.test
│       ├── generate_readmers_kmc.nf.test
│       ├── get_counts_kmc.nf.test
│       └── generate_histogram_kmc.nf.test
│
├── nf-test.config                                    # nf-test config (updated)
└── TEST_SUITE_SUMMARY.md                             # This file
```

## 🎯 Next Steps for Complete Coverage

### Immediate Priorities

1. **Add Cenhapmer Tests**
   ```bash
   # Template: tests/modules/local/get_counts_cenhapmer_kmc.nf.test
   - Test with mini KMC databases
   - Test parallel processing
   - Test different accessions
   ```

2. **Add Hybrid Profile Tests**
   ```bash
   # Template: tests/modules/local/combine_hybrid_profiles.nf.test
   - Test hybrid read detection
   - Test with different thresholds
   - Test Rust binary integration
   ```

3. **Add Integration Tests**
   ```bash
   # File: tests/workflows/charla_nf.nf.test
   - Test complete preprocessing workflow
   - Test k-mer profiling workflow
   - Test hybrid detection workflow
   ```

### Future Enhancements

- [ ] Subworkflow tests for utils modules
- [ ] Performance benchmarking tests
- [ ] Resource usage validation tests
- [ ] Multi-sample processing tests
- [ ] Tissue-type specific tests (pollen vs sperm)

## 💡 Usage Tips

### For Developers

```bash
# Before committing changes
nf-test test tests/modules/local/your_module.nf.test

# After modifying module behavior
nf-test test --update-snapshot

# Check all tests pass
nf-test test --verbose
```

### For Contributors

1. **Adding a new module:**
   ```bash
   # 1. Create module in modules/local/
   # 2. Create test in tests/modules/local/
   # 3. Run test
   nf-test test tests/modules/local/new_module.nf.test
   # 4. Update this summary
   ```

2. **Test-Driven Development:**
   ```bash
   # 1. Write failing test first
   # 2. Implement module
   # 3. Run test until it passes
   # 4. Refactor
   # 5. Add edge cases
   ```

## 🐛 Troubleshooting

### Common Issues

**Issue: Test data not found**
```
Error: File not found: tests/data/test_sample_longreads.fastq
```
Solution: Ensure working directory is pipeline root

**Issue: KMC out of memory**
```
KMC killed (signal 9)
```
Solution: Increase memory in `tests/nextflow.config`:
```groovy
process.memory = '8 GB'
```

**Issue: Snapshot mismatch**
```
Snapshot mismatch for process.out
```
Solution: Review changes, update if intentional:
```bash
nf-test test --update-snapshot
```

### Getting Help

- **GitHub Issues:** https://github.com/jacgonisa/charla_nf/issues
- **Email:** jg2070@cam.ac.uk
- **nf-test Docs:** https://code.askimed.com/nf-test/

## 📈 Benefits of This Test Suite

1. **Reliability** - Catch bugs before production
2. **Confidence** - Refactor without fear
3. **Documentation** - Tests demonstrate expected behavior
4. **Collaboration** - Easy for contributors to validate changes
5. **CI/CD** - Automated validation on every change
6. **Quality** - Enforces best practices

## 🎓 Resources

- [nf-test Documentation](https://code.askimed.com/nf-test/)
- [Nextflow Testing Guide](https://www.nextflow.io/docs/latest/testing.html)
- [nf-core Testing Best Practices](https://nf-co.re/docs/contributing/modules#writing-tests)
- [GitHub Actions Docs](https://docs.github.com/en/actions)

## ✨ Conclusion

Your CHARLA pipeline now has a professional-grade testing infrastructure that:

- ✅ Uses realistic long-read test data (ONT/PacBio style, 1-10kb)
- ✅ Tests critical preprocessing and k-mer profiling modules
- ✅ Includes comprehensive CI/CD with GitHub Actions
- ✅ Provides clear documentation for contributors
- ✅ Follows nf-core and Nextflow best practices
- ✅ Ready for expansion to cover remaining modules

**To run tests now:**
```bash
cd /home/jg2070/Desktop/PhD/crossover/charla_nf
nf-test test tests/modules/local/*.nf.test
```

---

**Created:** 2025-10-17
**Pipeline Version:** 1.0.0
**Test Framework:** nf-test
**CI/CD:** GitHub Actions
