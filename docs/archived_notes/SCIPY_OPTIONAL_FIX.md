# Scipy Optional Import Fix

## Problem

```
ERROR ~ Process `CHARLA_NF:ANALYZE_SCO_NCO` terminated with an error exit status (1)

Traceback (most recent call last):
  File ".../scripts/18-analyze_sco_nco_comparison.py", line 25, in <module>
    from scipy import stats
ModuleNotFoundError: No module named 'scipy'
```

## Root Cause

The Singularity container (`jacgonisa/charla_nf:v1.0.1`) does not include scipy, and the container filesystem is read-only, so `pip install scipy` cannot work inside the container.

## Solution

Made scipy **optional** in the script - if scipy is not available, statistical tests are skipped but plots are still generated.

### Changes to `scripts/18-analyze_sco_nco_comparison.py`

**Before (lines 19-25)**:
```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import sys
from scipy import stats
```

**After (lines 19-32)**:
```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import sys

# Try to import scipy, but make it optional
try:
    from scipy import stats
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    print("WARNING: scipy not available - statistical tests will be skipped")
```

**All statistical test calls wrapped** (5 locations):
```python
# Before:
if len(data_to_plot[0]) > 0 and len(data_to_plot[1]) > 0:
    t_stat, p_val = stats.ttest_ind(data_to_plot[0], data_to_plot[1])
    ...

# After:
if HAS_SCIPY and len(data_to_plot[0]) > 0 and len(data_to_plot[1]) > 0:
    t_stat, p_val = stats.ttest_ind(data_to_plot[0], data_to_plot[1])
    ...
```

This ensures that:
1. Script imports successfully even without scipy
2. All plots are generated normally
3. Statistical p-values are only added if scipy is available
4. User gets a warning if scipy is missing

### Removed Failed pip Install Attempts

Also removed the `pip install scipy` lines from modules since they don't work in read-only containers:

**Files modified**:
- `modules/local/analyze_sco_nco.nf:25` - Removed pip install line
- `modules/local/visualize_segments.nf:25` - Removed pip install line

## Impact

### Without scipy (current container):
- ✅ All plots generated
- ✅ All metrics calculated
- ❌ No p-values shown on plots
- ⚠️ Warning message printed

### With scipy (if container updated):
- ✅ All plots generated
- ✅ All metrics calculated
- ✅ P-values shown on plots
- ✅ No warnings

## Future Recommendation

Update the container to include scipy:

```dockerfile
# In the container Dockerfile:
RUN pip install scipy
```

Or add to conda environment:
```yaml
dependencies:
  - scipy
```

This will enable full statistical testing in the plots.

## Testing

The script will now run successfully without scipy:

```bash
nextflow run charla_nf/main.nf ... -resume
```

Expected output:
```
WARNING: scipy not available - statistical tests will be skipped
=== SCO vs NCO analysis complete ===
sco_nco_comparison.png
sco_nco_resolution_analysis.png
sco_nco_summary.tsv
```

All plots will be generated, just without p-value annotations.

---

**Status**: Fixed ✅
**Date**: 2025-11-28
**Impact**: Non-breaking - script works with or without scipy
