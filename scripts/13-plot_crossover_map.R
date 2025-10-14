#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(IRanges)
    library(stringr)
    library(ggplot2)
    library(optparse)
    library(argparse)
})

## -------------------------------------------------------------------
## Parse command-line arguments
## -------------------------------------------------------------------
# Define arguments
parser <- ArgumentParser(description = "Plot crossover coverage maps for Col and Ler references")

parser$add_argument("--col_bed", required = TRUE, help = "Col-0 BED file")
parser$add_argument("--ler_bed", required = TRUE, help = "Ler-0 BED file")
parser$add_argument("--col_paf", required = TRUE, help = "Filtered Col-0 PAF file")
parser$add_argument("--ler_paf", required = TRUE, help = "Filtered Ler-0 PAF file")
parser$add_argument("--sample_name", required = TRUE, help = "Sample name for output")

args <- parser$parse_args()

# Assign to variables for convenience
col_bed_file <- args$col_bed
ler_bed_file <- args$ler_bed
col_paf_file <- args$col_paf
ler_paf_file <- args$ler_paf
sample_name <- args$sample_name

cat("Running crossover plot script for sample:", sample_name, "\n")
cat("  Col BED:", col_bed_file, "\n")
cat("  Ler BED:", ler_bed_file, "\n")
cat("  Col PAF:", col_paf_file, "\n")
cat("  Ler PAF:", ler_paf_file, "\n\n")

## -------------------------------------------------------------------
## Define helper functions
## -------------------------------------------------------------------

process_bed <- function(bed_file, pattern_prefix) {
    cat("Processing BED file:", bed_file, "\n")
    bed <- read.table(bed_file, header = FALSE, stringsAsFactors = FALSE)
    colnames(bed) <- c("region", "start", "end")
    bed$length <- bed$end - bed$start
    bed$chrom <- sapply(strsplit(bed$region, "_"), `[`, 1)
    
    bed$order <- ifelse(grepl(paste0("ARMS_", pattern_prefix, "_1"), bed$region), 1,
                    ifelse(grepl(paste0("CEN_", pattern_prefix), bed$region), 2,
                    ifelse(grepl(paste0("ARMS_", pattern_prefix, "_2"), bed$region), 3, NA)))
    
    bed <- bed[order(bed$chrom, bed$order), ]
    
    bed_offsets <- bed %>%
        group_by(chrom) %>%
        arrange(order) %>%
        mutate(offset = cumsum(lag(length, default = 0))) %>%
        ungroup()
    
    return(bed_offsets)
}


process_paf <- function(paf_file, bed_offsets) {
    cat("Processing PAF file:", paf_file, "\n")
    
    if (!file.exists(paf_file)) {
        stop("ERROR: File not found: ", paf_file)
    }
    
    paf <- read.table(paf_file, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
    colnames(paf) <- c("read_id", "region", "ref_length", "start", "end")
    
    paf_merged <- merge(paf, bed_offsets[, c("region", "chrom", "offset")], by = "region")
    
    paf_merged <- paf_merged %>%
        mutate(
            new_start = start + offset,
            new_end = end + offset
        ) %>%
        mutate(
            num_code = as.integer(str_match(read_id, "_([123])_")[, 2])
        ) %>%
        filter(!is.na(num_code)) %>%
        mutate(pos = case_when(
            num_code == 1 ~ new_end,
            num_code %in% c(2, 3) ~ new_start
        ))
    
    return(paf_merged)
}


calculate_coverage <- function(paf_merged, chrom_lengths) {
    coverage_list <- list()
    chroms <- unique(paf_merged$chrom)

    for(ch in chroms) {
        ch_data <- paf_merged %>% filter(chrom == ch)
        total_length <- chrom_lengths %>% filter(chrom == ch) %>% pull(chrom_length)
        
        ranges <- IRanges(start = ch_data$pos, width = 1)
        cov <- coverage(ranges, width = total_length)
        cov_vec <- as.vector(cov)
        
        df <- data.frame(
            position = 1:total_length,
            coverage = cov_vec,
            chrom = ch
        )
        coverage_list[[ch]] <- df
    }
    
    return(bind_rows(coverage_list))
}

## -------------------------------------------------------------------
## Main logic
## -------------------------------------------------------------------

# Process BEDs
bed_offsets_col <- process_bed(col_bed_file, "Col")
bed_offsets_ler <- process_bed(ler_bed_file, "Ler")

chrom_lengths_col <- bed_offsets_col %>%
    group_by(chrom) %>%
    summarise(chrom_length = max(offset + length)) %>%
    ungroup()

chrom_lengths_ler <- bed_offsets_ler %>%
    group_by(chrom) %>%
    summarise(chrom_length = max(offset + length)) %>%
    ungroup()

# Process PAFs
paf_merged_col <- process_paf(col_paf_file, bed_offsets_col)
paf_merged_ler <- process_paf(ler_paf_file, bed_offsets_ler)

# Compute coverage
coverage_col <- calculate_coverage(paf_merged_col, chrom_lengths_col)
coverage_ler <- calculate_coverage(paf_merged_ler, chrom_lengths_ler)

coverage_col$reference <- "Col"
coverage_ler$reference <- "Ler"

# Bin coverage for plotting
binned_col <- coverage_col %>%
    group_by(chrom) %>%
    mutate(bin = cut(position, breaks = 1000, labels = FALSE)) %>%
    group_by(chrom, bin) %>%
    summarise(
        position = mean(position),
        coverage = mean(coverage),
        .groups = "drop"
    )

binned_ler <- coverage_ler %>%
    group_by(chrom) %>%
    mutate(bin = cut(position, breaks = 1000, labels = FALSE)) %>%
    group_by(chrom, bin) %>%
    summarise(
        position = mean(position),
        coverage = mean(coverage),
        .groups = "drop"
    )

# Centromere coordinates
cen_coords_col <- bed_offsets_col %>%
    filter(order == 2) %>%
    mutate(cen_start = offset + 1,
           cen_end = offset + length,
           reference = "Col") %>%
    dplyr::select(chrom, cen_start, cen_end, reference)

cen_coords_ler <- bed_offsets_ler %>%
    filter(order == 2) %>%
    mutate(cen_start = offset + 1,
           cen_end = offset + length,
           reference = "Ler") %>%
    dplyr::select(chrom, cen_start, cen_end, reference)

# Plot Col
p_col <- ggplot(binned_col, aes(x = position, y = coverage)) +
    geom_rect(data = cen_coords_col,
              aes(xmin = cen_start, xmax = cen_end, ymin = -Inf, ymax = Inf),
              fill = "grey", alpha = 0.3, inherit.aes = FALSE) +
    geom_line() +
    facet_wrap(~ chrom, scales = "free_x") +
    labs(title = paste("Crossover Coverage Plot - Col Reference"),
         subtitle = paste("Sample:", sample_name),
         x = "Chromosome Position", y = "Coverage") +
    theme_bw()

ggsave(paste0(sample_name, "_crossover_plot_Col", ".svg"),
       plot = p_col, width = 12, height = 8, dpi = 300)

# Plot Ler
p_ler <- ggplot(binned_ler, aes(x = position, y = coverage)) +
    geom_rect(data = cen_coords_ler,
              aes(xmin = cen_start, xmax = cen_end, ymin = -Inf, ymax = Inf),
              fill = "grey", alpha = 0.3, inherit.aes = FALSE) +
    geom_line() +
    facet_wrap(~ chrom, scales = "free_x") +
    labs(title = paste("Crossover Coverage Plot - Ler Reference"),
         subtitle = paste("Sample:", sample_name),
         x = "Chromosome Position", y = "Coverage") +
    theme_bw()

ggsave(paste0(sample_name, "_crossover_plot_Ler", ".svg"),
       plot = p_ler, width = 12, height = 8, dpi = 300)

# Write output tables
coverage_combined <- bind_rows(coverage_col, coverage_ler)
binned_combined <- bind_rows(binned_col, binned_ler)

write.table(coverage_combined,
            paste0(sample_name, "_coverage_data", ".txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(binned_combined,
            paste0(sample_name, "_binned_coverage", ".txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Crossover plot generation completed!\n")

