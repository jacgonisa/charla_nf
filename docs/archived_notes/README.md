# Archived Development Notes

This directory contains detailed session notes and debugging logs from the development of CHARLA pipeline.

## Contents

These files document specific bug fixes, features, and improvements made during development sessions:

### Classification & Confidence Scoring
- `CLASSIFICATION_BUG_FIX.md` - Critical bug in segment parental type extraction
- `CLASSIFICATION_UPDATE_V2.md` - Updated logic (3 segments = NCO, >3 = Complex)
- `CLASSIFICATION_CONCORDANCE_ADDED.md` - Added before/after/PAF classification comparison
- `CONFIDENCE_V3_COMPLETE.md` - Redesigned confidence scoring system
- `BEFORE_AFTER_MAPPING_CLASSIFICATION.md` - Using curated profiles as truth

### SCO/NCO Analysis
- `SCO_NCO_ANALYSIS_GUIDE.md` - Comprehensive guide to outputs and interpretation
- `BED_SEGMENT_ANALYSIS_UPDATE.md` - Extracting segments from BED instead of PAF
- `COWIDTH_MERGE_FIX.md` - Fixed crossover width merging logic

### Bug Fixes
- `FILE_COLLISION_FIXES_COMPLETE.md` - Fixed Nextflow staging collisions
- `CSV_SEPARATOR_FIX.md` - Auto-detect comma vs tab separators
- `SCIPY_OPTIONAL_FIX.md` - Made scipy dependency optional
- `STAGING_FIX.md` - Fixed file staging issues
- `WORKFLOW_FIX.md` - Workflow variable name fixes

### Visualization
- `READ_STRUCTURE_VISUALIZATION_ADDED.md` - PAF-based read structure plots
- `SEGMENT_VISUALIZATION_ADDED.md` - BED-based segment visualization
- `SEGMENT_VIZ_FIX.md` - Visualization improvements

### Analysis Features
- `INDEL_MAPPING_OPTIMIZATION.md` - 10x speedup for indel analysis
- `INSERTION_ORIGINS_QUICKSTART.md` - Insertion origin analysis guide
- `NON_ECTOPIC_FILTERING.md` - Filtering ectopic recombination

### Maintenance
- `CACHE_CLEANUP_INSTRUCTIONS.md` - How to invalidate Nextflow caches
- `CLEANUP_SUMMARY.md` - Code cleanup summary
- `FINAL_FIXES_SUMMARY.md` - Final session summary
- `SESSION_SUMMARY.md` - Overall session notes

### Tests
- `test_classification.py` - Unit tests for classification logic

## For Current Documentation

See the main documentation files in the repository root:

- [`DEVELOPMENT_NOTES.md`](../../DEVELOPMENT_NOTES.md) - Consolidated development notes
- [`README.md`](../../README.md) - Main pipeline documentation
- [`CHANGELOG.md`](../../CHANGELOG.md) - Version history

---

**Note**: These files are kept for historical reference. For the most up-to-date information, refer to the main documentation files listed above.
