#!/usr/bin/env python3
"""Test classification function with real data examples."""

import re
import pandas as pd

def classify_read_type_from_segments(segments_str):
    """
    Classify read based on segment pattern string.
    """
    if not segments_str or pd.isna(segments_str):
        return 'Unknown'

    # Extract parental types from segments
    # Segment format: <number><parental><region><chr> e.g., "325LA4" = 325 bp, Ler, ARMS, Chr4
    # Need to find the first letter (C or L) in the segment
    segments = segments_str.split(',')
    parental_types = []
    for seg in segments:
        seg = seg.strip()
        # Find first letter that is C or L
        for char in seg:
            if char in ['C', 'L']:
                parental_types.append(char)
                break
        else:
            parental_types.append('X')

    parental_pattern = ''.join(parental_types)

    print(f"  DEBUG: segments_str='{segments_str}'")
    print(f"  DEBUG: segments={segments}")
    print(f"  DEBUG: parental_pattern='{parental_pattern}'")

    # Collapse consecutive same parental type
    collapsed = re.sub(r'C+', 'C', parental_pattern)
    collapsed = re.sub(r'L+', 'L', collapsed)

    print(f"  DEBUG: collapsed='{collapsed}'")

    # Classify based on number of segments (switches between parents)
    # User clarification: 3 segments = NCO, >3 segments = Complex
    num_segments = len(collapsed)

    if num_segments == 1:
        # Single parent (C or L)
        return 'Parental'
    elif num_segments == 2:
        # One switch (CL or LC)
        return 'SCO'
    elif num_segments == 3:
        # Two switches (CLC or LCL) - this is NCO (gene conversion)
        return 'NCO'
    elif num_segments == 4:
        # Three switches (CLCL or LCLC) - double crossover
        return 'DCO'
    else:
        # More than 4 segments = Complex
        return 'Complex'

# Test cases from real data
test_cases = [
    ("325LA4,319LA3", "Parental", "LL→L (ectopic, different chrs)"),
    ("50CA5,83LA5", "SCO", "CL (non-ectopic)"),
    ("1827CA3,336LA3", "SCO", "CL (non-ectopic)"),
    ("55CA1,316LC3,104CC5", "NCO", "CLC (3 segments = NCO)"),
    ("656LA4", "Parental", "L (single segment)"),
    ("112CA1", "Parental", "C (single segment)"),
    ("59CA3,118CA4", "Parental", "CC→C (ectopic, different chrs)"),
    ("62LA2,212LA1", "Parental", "LL→L (ectopic, different chrs)"),
    ("100CA1,200LA2,150CA3,180LA4", "DCO", "CLCL (4 segments = DCO)"),
    ("100CA1,200LA2,150CA3,180LA4,120CA5", "Complex", "CLCLC (5 segments = Complex)"),
]

print("Testing classification function:")
print("="*80)

all_correct = True
for segments, expected, note in test_cases:
    result = classify_read_type_from_segments(segments)
    status = "✓" if result == expected else "✗"
    if result != expected:
        all_correct = False
    print(f"{status} {segments:30s} → {result:10s} (expected: {expected:10s}) | {note}")

print("="*80)
if all_correct:
    print("✓ All tests passed!")
else:
    print("✗ Some tests failed!")
