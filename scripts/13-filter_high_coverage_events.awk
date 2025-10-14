#!/usr/bin/env awk -f
# ------------------------------------------------------------------
# filter_high_coverage_events.awk
# Collapse nearby alignments into clusters (±tolerance bp)
# and remove clusters with >max reads.
# Generates:
#   - filtered PAF file
#   - high-coverage BED file with problematic regions
# Also prints summary statistics to stderr.
# ------------------------------------------------------------------

BEGIN {
    tol = (tol > 0 ? tol : 50)     # default tolerance if not provided
    max = (max > 0 ? max : 3)      # default max reads per cluster
    out_prefix = (prefix ? prefix : "filtered")   # default prefix
    out_paf = out_prefix ".filtered.paf"
    out_bed = "high_coverage_regions_" out_prefix ".bed"
}

{
    # PAF simplified format: query ref_id ref_len ref_start ref_end
    ref = $2
    s   = $4
    e   = $5
    found = 0

    for (i = 1; i <= n; i++) {
        if (ref == refs[i] && (s >= starts[i] - tol && s <= starts[i] + tol) && (e >= ends[i] - tol && e <= ends[i] + tol)) {
            count[i]++
            lines[i] = lines[i] $0 ORS
            if (s < starts[i]) starts[i] = s
            if (e > ends[i]) ends[i] = e
            found = 1
            break
        }
    }

    if (!found) {
        n++
        refs[n] = ref
        starts[n] = s
        ends[n] = e
        count[n] = 1
        lines[n] = $0 ORS
    }
}

END {
    kept_clusters = 0
    filtered_clusters = 0
    kept_reads = 0
    filtered_reads = 0

    for (i = 1; i <= n; i++) {
        if (count[i] <= max) {
            printf "%s", lines[i] >> out_paf
            kept_clusters++
            kept_reads += count[i]
        } else {
            print refs[i], starts[i], ends[i], count[i] >> out_bed
            filtered_clusters++
            filtered_reads += count[i]
        }
    }

    print "Filtering summary for prefix: " out_prefix > "/dev/stderr"
    print "  Tolerance: " tol " bp" > "/dev/stderr"
    print "  Max reads per cluster: " max > "/dev/stderr"
    print "  Total clusters: " n > "/dev/stderr"
    print "  Kept clusters: " kept_clusters " (" kept_reads " reads)" > "/dev/stderr"
    print "  Filtered clusters: " filtered_clusters " (" filtered_reads " reads)" > "/dev/stderr"
}

