#!/usr/bin/env Rscript

# Script to create comprehensive alignment output matching collaborators' format
# Key insight: Use k-mer BED for read positions (where on the FULL read)
#              Use PAF for reference positions and strand

library(optparse)
library(dplyr)
library(tidyr)
library(readr)

# Parse command-line arguments
option_list <- list(
    make_option(c("--curated_dir"), type="character", help="Directory with curated confident reads"),
    make_option(c("--paf_b73"), type="character", help="Path to B73 PAF file"),
    make_option(c("--paf_mo17"), type="character", help="Path to Mo17 PAF file"),
    make_option(c("--sample_name"), type="character", help="Sample name"),
    make_option(c("--output_dir"), type="character", default=".", help="Output directory")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

curated_dir <- opt$curated_dir
paf_b73_file <- opt$paf_b73
paf_mo17_file <- opt$paf_mo17
sample_name <- opt$sample_name
output_dir <- opt$output_dir

cat("\n========================================\n")
cat("CREATING COMPREHENSIVE ALIGNMENT OUTPUT\n")
cat("========================================\n\n")
cat("Strategy:\n")
cat("  - K-mer BED: Read positions (where on the FULL read each segment is)\n")
cat("  - PAF files: Reference genome positions + strand\n\n")

# PAF format columns
paf_columns <- c(
    "query_name", "query_length", "query_start", "query_end", "strand",
    "target_name", "target_length", "target_start", "target_end",
    "num_matches", "alignment_length", "mapping_quality"
)

# Function to read and parse PAF file
read_paf <- function(paf_file, parent_name) {
    cat(sprintf("Reading PAF file for %s...\n", parent_name))

    if (!file.exists(paf_file)) {
        cat(sprintf("WARNING: PAF file not found: %s\n", paf_file))
        return(NULL)
    }

    paf <- read.delim(paf_file, header = FALSE, stringsAsFactors = FALSE, comment.char = "")
    paf <- paf[, 1:min(12, ncol(paf))]
    colnames(paf) <- paf_columns[1:ncol(paf)]

    paf <- paf %>%
        mutate(
            base_read_id = sub("_[0-9]+_[A-Z]+H[0-9]+$", "", query_name),
            # Extract segment parent from query name (BH* = B73, MH* = Mo17)
            segment_parent = ifelse(grepl("_[0-9]+_BH", query_name), "B73",
                                   ifelse(grepl("_[0-9]+_MH", query_name), "Mo17", NA)),
            alignment_parent = parent_name,
            chr = target_name
        )

    cat(sprintf("  Loaded %d alignments\n", nrow(paf)))

    return(paf)
}

# Read PAF files
paf_b73 <- read_paf(paf_b73_file, "B73")
paf_mo17 <- read_paf(paf_mo17_file, "Mo17")

# Combine PAF data
paf_all <- bind_rows(paf_b73, paf_mo17)

# Only keep alignments where segment_parent matches alignment_parent
paf_all <- paf_all %>%
    filter(segment_parent == alignment_parent)

cat(sprintf("Total correct alignments: %d\n", nrow(paf_all)))

# Read the curated confident BED files (k-mer based, FULL read coordinates)
cat("\nReading curated confident BED files (k-mer based)...\n")

if (file.info(curated_dir)$isdir) {
    bed_dir <- curated_dir
} else {
    bed_dir <- dirname(curated_dir)
}

bed_sco <- read.table(file.path(bed_dir, paste0(sample_name, "_confident_SCO.bed")),
                     header = FALSE, stringsAsFactors = FALSE,
                     col.names = c("read_id", "read_start", "read_end", "segment"))
bed_nco <- read.table(file.path(bed_dir, paste0(sample_name, "_confident_NCO.bed")),
                     header = FALSE, stringsAsFactors = FALSE,
                     col.names = c("read_id", "read_start", "read_end", "segment"))
bed_gcco <- read.table(file.path(bed_dir, paste0(sample_name, "_confident_GC-CO.bed")),
                      header = FALSE, stringsAsFactors = FALSE,
                      col.names = c("read_id", "read_start", "read_end", "segment"))

bed_sco$recomb_type <- "SCO"
bed_nco$recomb_type <- "NCO"
bed_gcco$recomb_type <- "GC-CO"

bed_all <- bind_rows(bed_sco, bed_nco, bed_gcco)

# Add parent assignment from segment (BH* = B73, MH* = Mo17)
bed_all <- bed_all %>%
    mutate(
        segment_parent = ifelse(grepl("^B", segment), "B73", "Mo17")
    )

cat(sprintf("Total BED segments: %d\n", nrow(bed_all)))

confident_reads <- unique(bed_all$read_id)
cat(sprintf("Confident reads: %d\n", length(confident_reads)))

# For each read segment, get best PAF alignment
# Group PAF by read + segment_parent + chromosome, find best
paf_summary <- paf_all %>%
    filter(base_read_id %in% confident_reads) %>%
    group_by(base_read_id, segment_parent, chr) %>%
    summarise(
        ref_start = min(target_start),
        ref_end = max(target_end),
        strand = names(sort(table(strand), decreasing = TRUE))[1],
        total_aln_length = sum(alignment_length),
        .groups = "drop"
    ) %>%
    # Keep chromosome with longest alignment per read+parent
    group_by(base_read_id, segment_parent) %>%
    arrange(desc(total_aln_length)) %>%
    slice(1) %>%
    ungroup()

cat(sprintf("PAF summary: %d read-parent pairs\n", nrow(paf_summary)))

# For each read + parent, get the segment boundaries from BED
# (these are FULL read coordinates)
bed_summary <- bed_all %>%
    group_by(read_id, segment_parent) %>%
    summarise(
        read_start = min(read_start),
        read_end = max(read_end),
        .groups = "drop"
    )

# Get full read length from BED (maximum position across all segments)
read_lengths <- bed_all %>%
    group_by(read_id) %>%
    summarise(
        read_length = max(read_end),
        .groups = "drop"
    )

# Combine BED (read positions) with PAF (reference positions + strand)
comprehensive <- bed_summary %>%
    left_join(paf_summary, by = c("read_id" = "base_read_id", "segment_parent")) %>%
    left_join(read_lengths, by = "read_id")

# Add recombination type
read_recomb_types <- bed_all %>%
    select(read_id, recomb_type) %>%
    distinct()

comprehensive <- comprehensive %>%
    left_join(read_recomb_types, by = "read_id")

cat(sprintf("Comprehensive table: %d read-parent pairs\n", nrow(comprehensive)))

# For each read, determine which parent comes first on the reference genome
# This is based on reference position (ref_start), not read position
comprehensive_with_order <- comprehensive %>%
    group_by(read_id) %>%
    arrange(read_id, ref_start) %>%
    mutate(
        ref_order = row_number(),
        first_parent = segment_parent[ref_order == 1],
        second_parent = if(n() > 1) segment_parent[ref_order == 2] else NA_character_
    ) %>%
    ungroup()

cat(sprintf("Determining reference-based order for %d reads...\n", n_distinct(comprehensive_with_order$read_id)))

# Now assign to "first" and "second" based on reference order
comprehensive_reordered <- comprehensive_with_order %>%
    mutate(
        position_label = if_else(ref_order == 1, "first", "second")
    )

# Pivot to get first and second columns
comprehensive_wide_temp <- comprehensive_reordered %>%
    pivot_wider(
        id_cols = c(read_id, read_length, recomb_type, first_parent, second_parent),
        names_from = position_label,
        values_from = c(read_start, read_end, chr, ref_start, ref_end, strand, segment_parent),
        names_sep = "_"
    )

# Now map first/second back to Mo17/B73 based on which parent was first
comprehensive_wide <- comprehensive_wide_temp %>%
    mutate(
        # Mo17 columns (could be first or second based on reference position)
        Mo17_read_start = if_else(segment_parent_first == "Mo17", read_start_first, read_start_second),
        Mo17_read_end = if_else(segment_parent_first == "Mo17", read_end_first, read_end_second),
        Mo17_align_length = read_length,
        Mo17_chr = if_else(segment_parent_first == "Mo17", chr_first, chr_second),
        Mo17_ref_start = if_else(segment_parent_first == "Mo17", ref_start_first, ref_start_second),
        Mo17_ref_end = if_else(segment_parent_first == "Mo17", ref_end_first, ref_end_second),
        Mo17_strand = if_else(segment_parent_first == "Mo17", strand_first, strand_second),

        # B73 columns (could be first or second based on reference position)
        B73_read_start = if_else(segment_parent_first == "B73", read_start_first, read_start_second),
        B73_read_end = if_else(segment_parent_first == "B73", read_end_first, read_end_second),
        B73_read_length = read_length,
        B73_chr = if_else(segment_parent_first == "B73", chr_first, chr_second),
        B73_ref_start = if_else(segment_parent_first == "B73", ref_start_first, ref_start_second),
        B73_ref_end = if_else(segment_parent_first == "B73", ref_end_first, ref_end_second),
        B73_strand = if_else(segment_parent_first == "B73", strand_first, strand_second)
    ) %>%
    select(
        read_id,
        Mo17_read_start,
        Mo17_read_end,
        Mo17_align_length,
        Mo17_chr,
        Mo17_ref_start,
        Mo17_ref_end,
        Mo17_strand,
        B73_read_start,
        B73_read_end,
        B73_read_length,
        B73_chr,
        B73_ref_start,
        B73_ref_end,
        B73_strand,
        recomb_type
    )

cat(sprintf("Final table: %d reads\n", nrow(comprehensive_wide)))

# Check for overlapping segments (should be rare)
comprehensive_wide <- comprehensive_wide %>%
    mutate(
        segments_overlap = case_when(
            is.na(Mo17_read_start) | is.na(B73_read_start) ~ NA,
            Mo17_read_end < B73_read_start ~ FALSE,  # Mo17 before B73
            B73_read_end < Mo17_read_start ~ FALSE,  # B73 before Mo17
            TRUE ~ TRUE  # Overlap
        ),
        same_chr = case_when(
            is.na(Mo17_chr) | is.na(B73_chr) ~ NA,
            Mo17_chr == B73_chr ~ TRUE,
            TRUE ~ FALSE
        )
    )

cat(sprintf("\nQuality checks:\n"))
cat(sprintf("  Reads with non-overlapping segments: %d (%.1f%%)\n",
           sum(!comprehensive_wide$segments_overlap, na.rm = TRUE),
           100 * sum(!comprehensive_wide$segments_overlap, na.rm = TRUE) / sum(!is.na(comprehensive_wide$segments_overlap))))
cat(sprintf("  Reads with overlapping segments: %d (%.1f%%)\n",
           sum(comprehensive_wide$segments_overlap, na.rm = TRUE),
           100 * sum(comprehensive_wide$segments_overlap, na.rm = TRUE) / sum(!is.na(comprehensive_wide$segments_overlap))))
cat(sprintf("  Reads with both parents on same chr: %d (%.1f%%)\n",
           sum(comprehensive_wide$same_chr, na.rm = TRUE),
           100 * sum(comprehensive_wide$same_chr, na.rm = TRUE) / sum(!is.na(comprehensive_wide$same_chr))))

# Remove QC columns before saving
comprehensive_wide <- comprehensive_wide %>%
    select(-segments_overlap, -same_chr)

# Save
outfile <- file.path(output_dir, paste0(sample_name, "_comprehensive_alignments.tsv"))
write.table(comprehensive_wide, outfile, sep = "\t", quote = FALSE,
           row.names = FALSE, col.names = TRUE)
cat(sprintf("\nSaved: %s\n", outfile))

# Show examples
cat("\n=== EXAMPLE READS (non-overlapping segments) ===\n")
examples <- comprehensive_wide %>%
    filter(!is.na(Mo17_read_start) & !is.na(B73_read_start)) %>%
    filter(Mo17_read_end < B73_read_start | B73_read_end < Mo17_read_start) %>%
    head(3)

print(examples[, c("read_id", "Mo17_read_start", "Mo17_read_end", "B73_read_start", "B73_read_end", "Mo17_chr", "B73_chr")])

cat("\nDONE!\n")
