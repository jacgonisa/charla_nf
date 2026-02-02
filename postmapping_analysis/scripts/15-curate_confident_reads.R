#!/usr/bin/env Rscript

# Script to curate confident recombinant reads with segment count filtering
# Creates separate BED files and filtered PAF files for confident reads only
# Also generates plots for the curated dataset

library(optparse)
library(dplyr)
library(ggplot2)
library(tidyr)

# Parse command-line arguments
option_list <- list(
    make_option(c("--cowidths_file"), type="character", help="Path to cowidths file"),
    make_option(c("--bed_file"), type="character", help="Path to BED file with segments"),
    make_option(c("--curated_profiles"), type="character", help="Path to curated hybrid profiles TSV"),
    make_option(c("--mapping_dir"), type="character", help="Directory containing PAF files for mapping"),
    make_option(c("--threshold"), type="integer", help="Threshold value"),
    make_option(c("--sample_name"), type="character", help="Sample name"),
    make_option(c("--parent1_name"), type="character", help="Parent 1 name"),
    make_option(c("--parent2_name"), type="character", help="Parent 2 name"),
    make_option(c("--output_dir"), type="character", default=".", help="Output directory")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Set parameters
cowidths_file <- opt$cowidths_file
bed_file <- opt$bed_file
curated_profiles <- opt$curated_profiles
mapping_dir <- opt$mapping_dir
threshold <- opt$threshold
sample_name <- opt$sample_name
parent1_name <- opt$parent1_name
parent2_name <- opt$parent2_name
output_dir <- opt$output_dir

# Create output directory
curated_dir <- file.path(output_dir, "curated_confident_reads")
dir.create(curated_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n========================================\n")
cat("CURATING CONFIDENT RECOMBINANT READS\n")
cat("========================================\n")
cat(sprintf("Sample: %s\n", sample_name))
cat(sprintf("Threshold: %d\n", threshold))
cat(sprintf("Parents: %s x %s\n", parent1_name, parent2_name))
cat(sprintf("Output directory: %s\n", curated_dir))
cat("\n")

# Read cowidths file
cat("Reading crossover width data...\n")
co_widths <- read.table(cowidths_file, header = TRUE, stringsAsFactors = FALSE)
cat(sprintf("Total width entries: %d\n", nrow(co_widths)))

# Read curated profiles to get hybrid types
cat("Reading curated hybrid profiles...\n")
curated_data <- read.table(curated_profiles, header = TRUE, stringsAsFactors = FALSE, sep = "\t")

# Add hybrid type information to co_widths
co_widths <- co_widths %>%
    left_join(curated_data %>% select(readname, hybrid_type), by = c("read" = "readname"))

# Filter for SCO, NCO, GC-CO only
co_widths_filtered <- co_widths %>%
    filter(hybrid_type %in% c("SCO", "NCO", "GC-CO"))

cat(sprintf("Recombinant reads (SCO/NCO/GC-CO): %d\n", length(unique(co_widths_filtered$read))))

# Read BED file
cat("Reading BED file...\n")
bed_data <- read.table(bed_file, header = FALSE, stringsAsFactors = FALSE,
                      col.names = c("read", "start", "end", "segment"))
cat(sprintf("Total segments in BED: %d\n", nrow(bed_data)))

# Get parent prefixes
parent1_prefix <- substr(parent1_name, 1, 1)
parent2_prefix <- substr(parent2_name, 1, 1)

# Check which reads have segments from BOTH parents
cat("\n=== FILTERING FOR TRUE HYBRID READS ===\n")
read_parent_check <- bed_data %>%
    group_by(read) %>%
    summarise(
        has_parent1 = any(grepl(paste0("^", parent1_prefix), segment)),
        has_parent2 = any(grepl(paste0("^", parent2_prefix), segment)),
        n_segments = n(),
        .groups = "drop"
    ) %>%
    mutate(has_both_parents = has_parent1 & has_parent2)

cat(sprintf("Reads with %s segments: %d\n", parent1_name, sum(read_parent_check$has_parent1)))
cat(sprintf("Reads with %s segments: %d\n", parent2_name, sum(read_parent_check$has_parent2)))
cat(sprintf("Reads with segments from BOTH parents: %d\n", sum(read_parent_check$has_both_parents)))

# Add parent check to cowidths
co_widths_filtered <- co_widths_filtered %>%
    left_join(read_parent_check %>% select(read, has_both_parents, n_segments), by = "read") %>%
    filter(has_both_parents == TRUE)

cat(sprintf("After filtering for both parents: %d reads\n", length(unique(co_widths_filtered$read))))

# Apply segment count filtering
cat("\n=== APPLYING SEGMENT COUNT FILTER ===\n")
cat("Expected segment counts:\n")
cat("  SCO: exactly 2 segments\n")
cat("  NCO: exactly 3 segments\n")
cat("  GC-CO: exactly 4 segments\n\n")

expected_segments <- list("SCO" = 2, "NCO" = 3, "GC-CO" = 4)

# Count segments per read and filter
reads_with_types <- co_widths_filtered %>%
    select(read, hybrid_type, n_segments) %>%
    distinct()

confident_reads_list <- list()
nonconfident_reads_list <- list()

for (htype in c("SCO", "NCO", "GC-CO")) {
    expected_n <- expected_segments[[htype]]

    reads_for_type <- reads_with_types %>% filter(hybrid_type == htype)

    confident <- reads_for_type %>% filter(n_segments == expected_n) %>% pull(read)
    nonconfident <- reads_for_type %>% filter(n_segments != expected_n) %>% pull(read)

    confident_reads_list[[htype]] <- confident
    nonconfident_reads_list[[htype]] <- nonconfident

    cat(sprintf("%s: %d total, %d confident, %d non-confident\n",
               htype, nrow(reads_for_type), length(confident), length(nonconfident)))
}

all_confident_reads <- unlist(confident_reads_list)
all_nonconfident_reads <- unlist(nonconfident_reads_list)

cat(sprintf("\nTotal confident reads: %d\n", length(all_confident_reads)))
cat(sprintf("Total non-confident reads: %d\n", length(all_nonconfident_reads)))

# Create BED files
cat("\n=== CREATING BED FILES ===\n")

for (htype in c("SCO", "NCO", "GC-CO")) {
    # Confident reads
    confident_segments <- bed_data %>%
        filter(read %in% confident_reads_list[[htype]])

    outfile <- file.path(curated_dir, paste0(sample_name, "_confident_", htype, ".bed"))
    write.table(confident_segments, file = outfile, sep = "\t",
               row.names = FALSE, col.names = FALSE, quote = FALSE)
    cat(sprintf("  %s confident: %s (%d segments from %d reads)\n",
               htype, basename(outfile), nrow(confident_segments), length(confident_reads_list[[htype]])))

    # Non-confident reads
    if (length(nonconfident_reads_list[[htype]]) > 0) {
        nonconfident_segments <- bed_data %>%
            filter(read %in% nonconfident_reads_list[[htype]])

        outfile_nc <- file.path(curated_dir, paste0(sample_name, "_nonconfident_", htype, ".bed"))
        write.table(nonconfident_segments, file = outfile_nc, sep = "\t",
                   row.names = FALSE, col.names = FALSE, quote = FALSE)
        cat(sprintf("  %s non-confident: %s (%d segments from %d reads)\n",
                   htype, basename(outfile_nc), nrow(nonconfident_segments), length(nonconfident_reads_list[[htype]])))
    }
}

# Save confident read IDs
readids_file <- file.path(curated_dir, paste0(sample_name, "_confident_read_ids.txt"))
writeLines(all_confident_reads, readids_file)
cat(sprintf("\nSaved confident read IDs: %s\n", basename(readids_file)))

# Filter PAF files
cat("\n=== FILTERING PAF FILES ===\n")
cat("Note: PAF filtering is currently disabled due to Nextflow file staging limitations.\n")
cat("The confident read BED files and plots are the primary outputs.\n")
cat("If you need filtered PAF files, you can run this script manually with the published directory paths.\n")

# Generate plots for confident reads only
cat("\n=== GENERATING PLOTS FOR CONFIDENT READS ===\n")

# Filter cowidths for confident reads only
co_widths_confident <- co_widths_filtered %>%
    filter(read %in% all_confident_reads)

# Statistics
cat("\nConfident Reads Statistics:\n")
confident_stats <- co_widths_confident %>%
    group_by(hybrid_type) %>%
    summarise(
        n_reads = n_distinct(read),
        mean_width = mean(width, na.rm = TRUE),
        median_width = median(width, na.rm = TRUE),
        min_width = min(width, na.rm = TRUE),
        max_width = max(width, na.rm = TRUE),
        .groups = "drop"
    )
print(confident_stats)

# Save statistics
stats_file <- file.path(curated_dir, paste0(sample_name, "_confident_stats.tsv"))
write.table(confident_stats, file = stats_file, sep = "\t",
           row.names = FALSE, quote = FALSE)

# Plot 1: Histogram of confident reads by type
plot_file_1 <- file.path(curated_dir, paste0(sample_name, "_confident_histogram_by_type.svg"))
svg(plot_file_1, width = 10, height = 6)
p1 <- ggplot(co_widths_confident, aes(x = width, fill = hybrid_type)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    facet_wrap(~hybrid_type, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c("SCO" = "#E41A1C", "NCO" = "#377EB8", "GC-CO" = "#4DAF4A")) +
    labs(title = paste0("Confident Recombinant Reads - ", sample_name),
         subtitle = sprintf("SCO: %d reads, NCO: %d reads, GC-CO: %d reads",
                          length(confident_reads_list$SCO),
                          length(confident_reads_list$NCO),
                          length(confident_reads_list$`GC-CO`)),
         x = "Crossover Width (bp)",
         y = "Count",
         fill = "Type") +
    theme_bw() +
    theme(legend.position = "none")
print(p1)
dev.off()
cat(sprintf("Saved plot: %s\n", basename(plot_file_1)))

# Plot 2: Barplot of confident read counts
plot_file_2 <- file.path(curated_dir, paste0(sample_name, "_confident_barplot.svg"))
svg(plot_file_2, width = 8, height = 6)
count_data <- data.frame(
    type = c("SCO", "NCO", "GC-CO"),
    confident = c(length(confident_reads_list$SCO),
                 length(confident_reads_list$NCO),
                 length(confident_reads_list$`GC-CO`)),
    nonconfident = c(length(nonconfident_reads_list$SCO),
                    length(nonconfident_reads_list$NCO),
                    length(nonconfident_reads_list$`GC-CO`))
)
count_data_long <- count_data %>%
    pivot_longer(cols = c(confident, nonconfident), names_to = "category", values_to = "count")

p2 <- ggplot(count_data_long, aes(x = type, y = count, fill = category)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("confident" = "#2ca02c", "nonconfident" = "#d62728"),
                     labels = c("Confident", "Non-confident")) +
    labs(title = paste0("Confident vs Non-confident Reads - ", sample_name),
         subtitle = sprintf("Total confident: %d, Total non-confident: %d",
                          length(all_confident_reads), length(all_nonconfident_reads)),
         x = "Recombination Type",
         y = "Number of Reads",
         fill = "Category") +
    theme_bw() +
    theme(legend.position = "top")
print(p2)
dev.off()
cat(sprintf("Saved plot: %s\n", basename(plot_file_2)))

# Plot 3: Boxplot of widths
plot_file_3 <- file.path(curated_dir, paste0(sample_name, "_confident_width_boxplot.svg"))
svg(plot_file_3, width = 8, height = 6)
p3 <- ggplot(co_widths_confident, aes(x = hybrid_type, y = width, fill = hybrid_type)) +
    geom_boxplot(alpha = 0.7) +
    scale_fill_manual(values = c("SCO" = "#E41A1C", "NCO" = "#377EB8", "GC-CO" = "#4DAF4A")) +
    scale_y_log10() +
    labs(title = paste0("Crossover Width Distribution - Confident Reads - ", sample_name),
         x = "Recombination Type",
         y = "Width (bp, log10 scale)",
         fill = "Type") +
    theme_bw() +
    theme(legend.position = "none")
print(p3)
dev.off()
cat(sprintf("Saved plot: %s\n", basename(plot_file_3)))

cat("\n========================================\n")
cat("CURATION COMPLETE\n")
cat("========================================\n")
cat(sprintf("Output directory: %s\n", curated_dir))
cat(sprintf("Confident reads: %d\n", length(all_confident_reads)))
cat(sprintf("  - SCO: %d\n", length(confident_reads_list$SCO)))
cat(sprintf("  - NCO: %d\n", length(confident_reads_list$NCO)))
cat(sprintf("  - GC-CO: %d\n", length(confident_reads_list$`GC-CO`)))
cat("\n")
