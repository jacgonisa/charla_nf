#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(argparse)
})

## -------------------------------------------------------------------
## Parse command-line arguments
## -------------------------------------------------------------------
parser <- ArgumentParser(description = "Generate crossover statistics and width histograms")

parser$add_argument("--cowidths_file", required = FALSE, default = "NA", help = "Crossover widths file (.cowidths) - optional for post-mapping analysis")
parser$add_argument("--curated_profiles", required = TRUE, help = "Curated hybrid profiles TSV file")
parser$add_argument("--threshold", required = TRUE, type = "integer", help = "Threshold value to filter profiles")
parser$add_argument("--sample_name", required = TRUE, help = "Sample name for output files")
parser$add_argument("--analysis_stage", required = FALSE, default = "profiles", help = "Stage: 'profiles' (before mapping) or 'mapping' (after mapping)")
parser$add_argument("--parent1_name", required = FALSE, default = "Parent1", help = "Name of parent 1")
parser$add_argument("--parent2_name", required = FALSE, default = "Parent2", help = "Name of parent 2")
parser$add_argument("--output_dir", required = FALSE, default = ".", help = "Output directory")

args <- parser$parse_args()

# Assign to variables
cowidths_file <- args$cowidths_file
curated_profiles_file <- args$curated_profiles
threshold_val <- args$threshold
sample_name <- args$sample_name
analysis_stage <- args$analysis_stage
parent1_name <- args$parent1_name
parent2_name <- args$parent2_name
output_dir <- args$output_dir

# Determine if we have cowidths data
has_cowidths <- (cowidths_file != "NA" && file.exists(cowidths_file))

cat("CROSSOVER STATISTICS ANALYSIS\n")
cat("=============================\n")
cat("Analysis stage:", toupper(analysis_stage), "\n")
cat("Sample:", sample_name, "\n")
cat("Threshold:", threshold_val, "\n")
if (has_cowidths) {
    cat("Cowidths file:", cowidths_file, "\n")
} else {
    cat("Cowidths file: Not provided (k-mer profile analysis only)\n")
}
cat("Curated profiles:", curated_profiles_file, "\n\n")

## -------------------------------------------------------------------
## Read curated profiles and get hybrid type statistics
## -------------------------------------------------------------------
cat("Reading curated hybrid profiles...\n")
profiles <- read.table(curated_profiles_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Filter for the specified threshold
profiles_thr <- profiles %>% filter(threshold == threshold_val)

cat("Total reads at threshold", threshold_val, ":", nrow(profiles_thr), "\n")

# Count hybrid types
hybrid_type_counts <- profiles_thr %>%
    group_by(hybrid_type) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(proportion = count / sum(count))

cat("\nHybrid Type Statistics:\n")
print(hybrid_type_counts)

# Save hybrid type statistics
write.table(hybrid_type_counts,
            file = file.path(output_dir, paste0(sample_name, "_hybrid_type_stats.tsv")),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Create bar plot of hybrid types
p_hybrid_types <- ggplot(hybrid_type_counts, aes(x = hybrid_type, y = count, fill = hybrid_type)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0(count, "\n(", round(proportion*100, 1), "%)")),
              vjust = -0.5, size = 3) +
    labs(title = paste("Hybrid Types -", sample_name),
         subtitle = paste("Threshold:", threshold_val),
         x = "Hybrid Type", y = "Count") +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, paste0(sample_name, "_hybrid_types_barplot.svg")),
       plot = p_hybrid_types, width = 8, height = 6, dpi = 300)

cat("Saved hybrid type statistics and barplot\n\n")

## -------------------------------------------------------------------
## Read crossover widths and generate histograms (ONLY if cowidths provided)
## -------------------------------------------------------------------

if (!has_cowidths) {
    cat("\n=============================\n")
    cat("ANALYSIS COMPLETE (K-MER PROFILE STAGE)\n")
    cat("=============================\n")
    cat("\nOutput files:\n")
    cat(" -", paste0(sample_name, "_hybrid_type_stats.tsv\n"))
    cat(" -", paste0(sample_name, "_hybrid_types_barplot.svg\n"))
    cat("\nNote: Crossover width analysis will be performed after mapping step.\n\n")
    quit(save = "no", status = 0)
}

cat("\n")
cat("Reading crossover widths...\n")

co_widths <- read.table(cowidths_file, header = TRUE, stringsAsFactors = FALSE)
cat("Total width entries:", nrow(co_widths), "\n")

# Read the bed_of_best.tsv file to get actual segment information
# The file is in the same directory as the cowidths file
bed_file <- sub("\\.cowidths$", ".tsv", cowidths_file)
cat("Reading segment information from:", bed_file, "\n")

# Use first letter of parent names to identify segments
parent1_prefix <- substr(parent1_name, 1, 1)
parent2_prefix <- substr(parent2_name, 1, 1)

# Read BED file and check which reads have segments from BOTH parents
if (file.exists(bed_file)) {
    bed_data <- read.table(bed_file, header = FALSE, stringsAsFactors = FALSE,
                          col.names = c("read", "start", "end", "segment"))

    # For each read, check if it has segments from BOTH parents
    read_parent_check <- bed_data %>%
        group_by(read) %>%
        summarise(
            has_parent1 = any(grepl(paste0("^", parent1_prefix), segment)),
            has_parent2 = any(grepl(paste0("^", parent2_prefix), segment)),
            segments_list = list(segment),
            .groups = "drop"
        ) %>%
        mutate(
            has_both_parents = has_parent1 & has_parent2
        )

    cat("Total reads in BED file:", nrow(read_parent_check), "\n")
    cat("Reads with", parent1_name, "segments:", sum(read_parent_check$has_parent1), "\n")
    cat("Reads with", parent2_name, "segments:", sum(read_parent_check$has_parent2), "\n")
    cat("Reads with segments from BOTH parents:", sum(read_parent_check$has_both_parents), "\n")
    cat("Reads with only one parent (should be filtered):", sum(!read_parent_check$has_both_parents), "\n")

    # Add this information to cowidths
    co_widths <- co_widths %>%
        left_join(read_parent_check %>% select(read, has_both_parents), by = "read")

    # Store bed_data for later use in creating separate BED files
    bed_data_full <- bed_data
} else {
    cat("Warning: bed file not found\n")
    co_widths$has_both_parents <- NA
}

# Merge with hybrid type information from curated profiles
profiles_thr_unique <- profiles_thr %>%
    dplyr::select(readname, hybrid_type) %>%
    distinct()

# Join hybrid types with cowidths
co_widths <- co_widths %>%
    left_join(profiles_thr_unique, by = c("read" = "readname"))

# For reads without classification, mark as "Unknown"
co_widths$hybrid_type[is.na(co_widths$hybrid_type)] <- "Unknown"

cat("\n=== FILTERING FOR TRUE RECOMBINANT READS ===\n")
cat("Initial reads:", nrow(co_widths), "\n")

# Filter out zero-width events
co_widths_nonzero <- co_widths %>% filter(width > 0)
cat("After removing zero-width:", nrow(co_widths_nonzero), "\n")

# Filter for only SCO, NCO, and GC-CO (exclude None and Unknown)
co_widths_filtered <- co_widths_nonzero %>%
    filter(hybrid_type %in% c("SCO", "NCO", "GC-CO"))
cat("After filtering for SCO/NCO/GC-CO:", nrow(co_widths_filtered), "\n")

# CRITICAL FILTER: Only keep reads with segments from BOTH parents
if ("has_both_parents" %in% colnames(co_widths_filtered)) {
    co_widths_filtered <- co_widths_filtered %>%
        filter(has_both_parents == TRUE)
    cat("After filtering for reads with BOTH parent segments:", nrow(co_widths_filtered), "\n")
} else {
    cat("Warning: Could not check for both parents (BED file not available)\n")
}
cat("==========================================\n\n")

# Calculate overall statistics (only reads with both parents)
overall_stats <- data.frame(
    Category = "Overall",
    N = nrow(co_widths_filtered),
    Mean_width = mean(co_widths_filtered$width),
    Median_width = median(co_widths_filtered$width),
    Min_width = min(co_widths_filtered$width),
    Max_width = max(co_widths_filtered$width),
    stringsAsFactors = FALSE
)

# Calculate statistics by hybrid type
type_stats <- co_widths_filtered %>%
    group_by(hybrid_type) %>%
    summarise(
        N = n(),
        Mean_width = mean(width),
        Median_width = median(width),
        Min_width = min(width),
        Max_width = max(width),
        .groups = "drop"
    ) %>%
    rename(Category = hybrid_type)

# Combine statistics
width_stats <- bind_rows(overall_stats, type_stats)

cat("\nCrossover Width Statistics:\n")
print(width_stats)

# Save width statistics
write.table(width_stats,
            file = file.path(output_dir, paste0(sample_name, "_cowidth_stats.tsv")),
            sep = "\t", row.names = FALSE, quote = FALSE)

## -------------------------------------------------------------------
## Create separate BED files for each hybrid type with segment count filtering
## -------------------------------------------------------------------
cat("\n=== CREATING SEPARATE BED FILES BY HYBRID TYPE ===\n")
cat("Applying segment count filter:\n")
cat("  SCO: exactly 2 segments\n")
cat("  NCO: exactly 3 segments\n")
cat("  GC-CO: exactly 4 segments\n\n")

# Expected segment counts for each type
expected_segments <- list(
    "SCO" = 2,
    "NCO" = 3,
    "GC-CO" = 4
)

if (exists("bed_data_full") && nrow(bed_data_full) > 0) {
    # Count segments per read
    segment_counts <- bed_data_full %>%
        group_by(read) %>%
        summarise(n_segments = n(), .groups = "drop")

    # Get list of reads for each hybrid type (that passed all filters)
    reads_by_type <- co_widths_filtered %>%
        select(read, hybrid_type) %>%
        distinct() %>%
        left_join(segment_counts, by = "read")

    # Create curated output directory
    curated_dir <- file.path(output_dir, "curated_single_read_level")
    dir.create(curated_dir, showWarnings = FALSE, recursive = TRUE)

    all_confident_reads <- c()
    all_nonconfident_reads <- c()

    for (htype in c("SCO", "NCO", "GC-CO")) {
        expected_n <- expected_segments[[htype]]

        reads_for_type <- reads_by_type %>%
            filter(hybrid_type == htype)

        # Split into confident and non-confident
        confident_reads <- reads_for_type %>%
            filter(n_segments == expected_n) %>%
            pull(read)

        nonconfident_reads <- reads_for_type %>%
            filter(n_segments != expected_n) %>%
            pull(read)

        all_confident_reads <- c(all_confident_reads, confident_reads)
        all_nonconfident_reads <- c(all_nonconfident_reads, nonconfident_reads)

        cat(sprintf("%s: %d total reads\n", htype, nrow(reads_for_type)))
        cat(sprintf("  Confident (%d segments): %d reads\n", expected_n, length(confident_reads)))
        cat(sprintf("  Non-confident (other): %d reads\n", length(nonconfident_reads)))

        # Save confident reads
        if (length(confident_reads) > 0) {
            bed_confident <- bed_data_full %>%
                filter(read %in% confident_reads)

            outfile_confident <- file.path(curated_dir, paste0(sample_name, "_", htype, ".bed"))
            write.table(bed_confident, file = outfile_confident,
                       sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

            cat(sprintf("  Saved %s confident: %s (%d segments)\n",
                       htype, basename(outfile_confident), nrow(bed_confident)))
        }

        # Save non-confident reads
        if (length(nonconfident_reads) > 0) {
            bed_nonconfident <- bed_data_full %>%
                filter(read %in% nonconfident_reads)

            outfile_nonconfident <- file.path(curated_dir, paste0(sample_name, "_non-confident_", htype, ".bed"))
            write.table(bed_nonconfident, file = outfile_nonconfident,
                       sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

            cat(sprintf("  Saved %s non-confident: %s (%d segments)\n\n",
                       htype, basename(outfile_nonconfident), nrow(bed_nonconfident)))
        } else {
            cat("\n")
        }
    }

    cat(sprintf("Total confident reads: %d\n", length(all_confident_reads)))
    cat(sprintf("Total non-confident reads: %d\n", length(all_nonconfident_reads)))

    # Save list of confident read IDs
    writeLines(all_confident_reads, file.path(curated_dir, paste0(sample_name, "_confident_read_ids.txt")))

    ## -------------------------------------------------------------------
    ## Filter PAF files for confident reads only
    ## -------------------------------------------------------------------
    cat("\n=== FILTERING PAF FILES FOR CONFIDENT READS ===\n")

    # The cowidths file might be a symlink - resolve it first
    # The cowidths file is in threshold50/plotting_crossover/sample_name/
    # PAF files are in threshold50/mapping_full_read/sample_name/nonectopic/CHR_segments_alntoB73_alntoMo17/
    cowidths_real_path <- normalizePath(cowidths_file, mustWork = FALSE)
    cowidths_dir <- dirname(cowidths_real_path)
    threshold_dir <- dirname(dirname(cowidths_dir))  # Go up to threshold50
    paf_search_dir <- file.path(threshold_dir, "mapping_full_read", sample_name, "nonectopic", "CHR_segments_alntoB73_alntoMo17")

    cat(sprintf("Cowidths file: %s\n", cowidths_file))
    cat(sprintf("Resolved path: %s\n", cowidths_real_path))
    cat(sprintf("Threshold dir: %s\n", threshold_dir))
    cat(sprintf("Looking for PAF files in: %s\n", paf_search_dir))

    if (dir.exists(paf_search_dir)) {
        paf_files <- list.files(paf_search_dir, pattern = "\\.paf$", full.names = TRUE)

        if (length(paf_files) > 0) {
            cat(sprintf("Found %d PAF files\n", length(paf_files)))

            for (paf_file in paf_files) {
                # Read PAF file (tab-delimited, no header)
                # PAF format has variable number of columns (12 mandatory + optional tags)
                # First column contains segment IDs like "read_id_segnum_haplotype"
                # Need to extract read ID from segment ID
                tryCatch({
                    # Read all lines
                    all_lines <- readLines(paf_file)

                    # Filter lines where read ID (extracted from segment ID) is in confident reads
                    # Segment format: read_id_segnum_haplotype, e.g., "48097b0a-b3bd-49b6-a473-95e0d1a576a8_1_MH2"
                    # Read ID is everything before the last two underscores
                    confident_lines <- all_lines[sapply(strsplit(all_lines, "\t"), function(x) {
                        segment_id <- x[1]
                        # Extract read ID by removing last two underscore-separated fields
                        parts <- strsplit(segment_id, "_")[[1]]
                        if (length(parts) >= 3) {
                            read_id <- paste(parts[1:(length(parts)-2)], collapse="_")
                            return(read_id %in% all_confident_reads)
                        }
                        return(FALSE)
                    })]

                    # Save filtered PAF
                    paf_basename <- basename(paf_file)
                    outfile_paf <- file.path(curated_dir, paf_basename)
                    writeLines(confident_lines, outfile_paf)

                    cat(sprintf("  Filtered %s: %d -> %d alignments\n",
                               paf_basename, length(all_lines), length(confident_lines)))
                }, error = function(e) {
                    cat(sprintf("  Warning: Could not process %s: %s\n", basename(paf_file), e$message))
                })
            }
        } else {
            cat("No PAF files found in", paf_search_dir, "\n")
        }
    } else {
        cat("PAF directory not found:", paf_search_dir, "\n")
    }
    cat("=============================================\n\n")

} else {
    cat("Warning: BED data not available, cannot create per-type BED files\n")
}
cat("=============================================\n\n")

## -------------------------------------------------------------------
## Generate histograms
## -------------------------------------------------------------------
cat("\nGenerating histograms...\n")

# Overall histogram with median line (using filtered data with both parents)
p_overall <- ggplot(co_widths_filtered, aes(x = width)) +
    geom_histogram(binwidth = 500, fill = "steelblue", color = "white") +
    geom_vline(xintercept = median(co_widths_filtered$width),
               color = "red", linetype = "dashed", linewidth = 1) +
    annotate("text", x = median(co_widths_filtered$width),
             y = Inf,
             label = paste("Median =", round(median(co_widths_filtered$width)), "bp"),
             hjust = -0.1, vjust = 2, color = "red") +
    scale_x_continuous(limits = c(0, 30000), labels = scales::comma) +
    labs(title = paste("Crossover Width Distribution -", sample_name),
         subtitle = paste("Threshold:", threshold_val, "| N =", nrow(co_widths_filtered), "- Reads with both parent segments"),
         x = "Crossover Width (bp)", y = "Count") +
    theme_bw()

ggsave(file.path(output_dir, paste0(sample_name, "_cowidth_histogram_overall.svg")),
       plot = p_overall, width = 10, height = 6, dpi = 300)

# SCO vs NCO vs GC-CO comparison histogram (overlaid)
p_types <- ggplot(co_widths_filtered, aes(x = width, fill = hybrid_type)) +
    geom_histogram(binwidth = 500, position = "identity", alpha = 0.6) +
    scale_fill_manual(values = c("SCO" = "darkgreen", "NCO" = "darkorange", "GC-CO" = "purple")) +
    scale_x_continuous(limits = c(0, 30000), labels = scales::comma) +
    labs(title = paste("Crossover Width by Type -", sample_name),
         subtitle = paste("Threshold:", threshold_val, "- Reads with both parent segments"),
         x = "Crossover Width (bp)", y = "Count", fill = "Type") +
    theme_bw() +
    theme(legend.position = "top")

ggsave(file.path(output_dir, paste0(sample_name, "_cowidth_histogram_by_type.svg")),
       plot = p_types, width = 10, height = 6, dpi = 300)

# Faceted histograms with medians (one panel per type)
median_data <- co_widths_filtered %>%
    group_by(hybrid_type) %>%
    summarise(median_val = median(width), .groups = "drop")

p_facet <- ggplot(co_widths_filtered, aes(x = width)) +
    geom_histogram(aes(fill = hybrid_type), binwidth = 500, color = "white") +
    geom_vline(data = median_data, aes(xintercept = median_val),
               color = "red", linetype = "dashed", linewidth = 1) +
    geom_text(data = median_data,
              aes(x = median_val, y = Inf,
                  label = paste("Median =", round(median_val), "bp")),
              hjust = -0.1, vjust = 2, color = "red", size = 3) +
    facet_wrap(~ hybrid_type, scales = "free_y", ncol = 1) +
    scale_x_continuous(limits = c(0, 30000), labels = scales::comma) +
    scale_fill_manual(values = c("SCO" = "darkgreen", "NCO" = "darkorange", "GC-CO" = "purple")) +
    labs(title = paste("Crossover Width Distribution by Type -", sample_name),
         subtitle = paste("Threshold:", threshold_val, "- Reads with both parent segments"),
         x = "Crossover Width (bp)", y = "Count") +
    theme_bw() +
    theme(legend.position = "none")

ggsave(file.path(output_dir, paste0(sample_name, "_cowidth_histogram_faceted.svg")),
       plot = p_facet, width = 10, height = 10, dpi = 300)

cat("\nSaved all histograms\n")

## -------------------------------------------------------------------
## Summary report
## -------------------------------------------------------------------
cat("\n=============================\n")
cat("ANALYSIS COMPLETE (MAPPING STAGE)\n")
cat("=============================\n")
cat("\nOutput files:\n")
cat(" -", paste0(sample_name, "_hybrid_type_stats.tsv\n"))
cat(" -", paste0(sample_name, "_hybrid_types_barplot.svg\n"))
cat(" -", paste0(sample_name, "_cowidth_stats.tsv\n"))
cat(" -", paste0(sample_name, "_cowidth_histogram_overall.svg\n"))
cat(" -", paste0(sample_name, "_cowidth_histogram_by_type.svg\n"))
cat(" -", paste0(sample_name, "_cowidth_histogram_faceted.svg\n"))
cat("\n")
