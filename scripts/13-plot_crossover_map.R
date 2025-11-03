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
parser <- ArgumentParser(description = "Plot crossover coverage maps for two parent references")

parser$add_argument("--col_bed", required = FALSE, default = "NA", help = "Parent 1 BED file (optional for centromere-aware mode)")
parser$add_argument("--ler_bed", required = FALSE, default = "NA", help = "Parent 2 BED file (optional for centromere-aware mode)")
parser$add_argument("--col_paf", required = TRUE, help = "Filtered Parent 1 PAF file")
parser$add_argument("--ler_paf", required = TRUE, help = "Filtered Parent 2 PAF file")
parser$add_argument("--sample_name", required = TRUE, help = "Sample name for output")
parser$add_argument("--parent1_name", required = FALSE, default = "Parent1", help = "Name of parent 1")
parser$add_argument("--parent2_name", required = FALSE, default = "Parent2", help = "Name of parent 2")

args <- parser$parse_args()

# Assign to variables for convenience
col_bed_file <- args$col_bed
ler_bed_file <- args$ler_bed
col_paf_file <- args$col_paf
ler_paf_file <- args$ler_paf
sample_name <- args$sample_name
parent1_name <- args$parent1_name
parent2_name <- args$parent2_name

# Detect mode
centromere_aware <- (col_bed_file != "NA" && ler_bed_file != "NA" &&
                     file.exists(col_bed_file) && file.exists(ler_bed_file))

if (centromere_aware) {
    cat("Running in CENTROMERE-AWARE mode\n")
    cat("  Sample:", sample_name, "\n")
    cat("  Parent 1 (", parent1_name, ") BED:", col_bed_file, "\n")
    cat("  Parent 2 (", parent2_name, ") BED:", ler_bed_file, "\n")
    cat("  Parent 1 PAF:", col_paf_file, "\n")
    cat("  Parent 2 PAF:", ler_paf_file, "\n\n")
} else {
    cat("Running in CENTROMERE-UNAWARE mode (chromosome-level only)\n")
    cat("  Sample:", sample_name, "\n")
    cat("  Parent 1 (", parent1_name, ") PAF:", col_paf_file, "\n")
    cat("  Parent 2 (", parent2_name, ") PAF:", ler_paf_file, "\n\n")
}

## -------------------------------------------------------------------
## Define helper functions
## -------------------------------------------------------------------

# Function for centromere-unaware mode: extract chromosome info from PAF
process_paf_simple <- function(paf_file) {
    cat("Processing PAF file (centromere-unaware mode):", paf_file, "\n")

    if (!file.exists(paf_file)) {
        stop("ERROR: File not found: ", paf_file)
    }

    # Check if file is empty
    lines <- readLines(paf_file, warn = FALSE)
    data_lines <- lines[!grepl("^#", lines)]

    if (length(data_lines) == 0) {
        cat("  Warning: PAF file is empty. Returning empty data frame.\n")
        return(data.frame(
            read_id = character(0),
            region = character(0),
            ref_length = integer(0),
            start = integer(0),
            end = integer(0),
            chrom = character(0),
            num_code = integer(0),
            pos = integer(0),
            stringsAsFactors = FALSE
        ))
    }

    paf <- read.table(paf_file, header = FALSE, stringsAsFactors = FALSE, sep = "\t", comment.char = "#")
    colnames(paf) <- c("read_id", "region", "ref_length", "start", "end")

    # Extract chromosome name from region (e.g., "CHR_chr1" -> "chr1")
    paf$chrom <- sapply(strsplit(paf$region, "_"), function(x) paste(x[-1], collapse="_"))

    # Extract num_code from read_id (positions 1, 2, or 3)
    paf$num_code <- as.integer(str_match(paf$read_id, "_([123])_")[, 2])
    paf <- paf %>% filter(!is.na(num_code))

    # Calculate position based on num_code
    paf$pos <- ifelse(paf$num_code == 1, paf$end, paf$start)

    return(paf)
}

# Function to get chromosome lengths from PAF (centromere-unaware mode)
get_chrom_lengths_from_paf <- function(paf_data) {
    chrom_lengths <- paf_data %>%
        group_by(chrom) %>%
        summarise(chrom_length = max(ref_length), .groups = "drop")
    return(chrom_lengths)
}

# Original function for centromere-aware mode
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

    # Check if file is empty or only contains comments
    lines <- readLines(paf_file, warn = FALSE)
    data_lines <- lines[!grepl("^#", lines)]

    if (length(data_lines) == 0) {
        cat("  Warning: PAF file is empty (no data). Returning empty data frame.\n")
        # Return empty data frame with correct column structure
        return(data.frame(
            read_id = character(0),
            region = character(0),
            ref_length = integer(0),
            start = integer(0),
            end = integer(0),
            chrom = character(0),
            offset = integer(0),
            new_start = integer(0),
            new_end = integer(0),
            num_code = integer(0),
            pos = integer(0),
            stringsAsFactors = FALSE
        ))
    }

    paf <- read.table(paf_file, header = FALSE, stringsAsFactors = FALSE, sep = "\t", comment.char = "#")
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
    # Handle empty PAF data
    if (nrow(paf_merged) == 0) {
        cat("  Warning: No PAF data to calculate coverage. Returning zero coverage.\n")
        # Create zero coverage for all chromosomes
        coverage_list <- list()
        for (i in 1:nrow(chrom_lengths)) {
            ch <- chrom_lengths$chrom[i]
            total_length <- chrom_lengths$chrom_length[i]
            df <- data.frame(
                position = 1:total_length,
                coverage = 0,
                chrom = ch
            )
            coverage_list[[ch]] <- df
        }
        return(bind_rows(coverage_list))
    }

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

if (centromere_aware) {
    ## CENTROMERE-AWARE MODE: Use BED files
    cat("Processing BED files and PAF files with centromere annotations...\n")

    # Process BEDs
    bed_offsets_col <- process_bed(col_bed_file, parent1_name)
    bed_offsets_ler <- process_bed(ler_bed_file, parent2_name)

    chrom_lengths_col <- bed_offsets_col %>%
        group_by(chrom) %>%
        summarise(chrom_length = max(offset + length)) %>%
        ungroup()

    chrom_lengths_ler <- bed_offsets_ler %>%
        group_by(chrom) %>%
        summarise(chrom_length = max(offset + length)) %>%
        ungroup()

    # Process PAFs with BED offsets
    paf_merged_col <- process_paf(col_paf_file, bed_offsets_col)
    paf_merged_ler <- process_paf(ler_paf_file, bed_offsets_ler)

    # Centromere coordinates for shading
    cen_coords_col <- bed_offsets_col %>%
        filter(order == 2) %>%
        mutate(cen_start = offset + 1,
               cen_end = offset + length,
               reference = parent1_name) %>%
        dplyr::select(chrom, cen_start, cen_end, reference)

    cen_coords_ler <- bed_offsets_ler %>%
        filter(order == 2) %>%
        mutate(cen_start = offset + 1,
               cen_end = offset + length,
               reference = parent2_name) %>%
        dplyr::select(chrom, cen_start, cen_end, reference)

} else {
    ## CENTROMERE-UNAWARE MODE: No BED files, use PAF directly
    cat("Processing PAF files in chromosome-level mode (no centromere annotations)...\n")

    # Process PAFs directly (simple mode)
    paf_merged_col <- process_paf_simple(col_paf_file)
    paf_merged_ler <- process_paf_simple(ler_paf_file)

    # Get chromosome lengths from PAF
    chrom_lengths_col <- get_chrom_lengths_from_paf(paf_merged_col)
    chrom_lengths_ler <- get_chrom_lengths_from_paf(paf_merged_ler)

    # No centromere coordinates to shade
    cen_coords_col <- NULL
    cen_coords_ler <- NULL
}

# Compute coverage
coverage_col <- calculate_coverage(paf_merged_col, chrom_lengths_col)
coverage_ler <- calculate_coverage(paf_merged_ler, chrom_lengths_ler)

coverage_col$reference <- parent1_name
coverage_ler$reference <- parent2_name

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

# Plot Parent 1
p_col <- ggplot(binned_col, aes(x = position, y = coverage)) +
    {if (centromere_aware && !is.null(cen_coords_col))
        geom_rect(data = cen_coords_col,
                  aes(xmin = cen_start, xmax = cen_end, ymin = -Inf, ymax = Inf),
                  fill = "grey", alpha = 0.3, inherit.aes = FALSE)
    } +
    geom_line() +
    facet_wrap(~ chrom, scales = "free_x") +
    labs(title = paste("Crossover Coverage Plot -", parent1_name, "Reference"),
         subtitle = paste("Sample:", sample_name),
         x = "Chromosome Position", y = "Coverage") +
    theme_bw()

ggsave(paste0(sample_name, "_crossover_plot_", parent1_name, ".svg"),
       plot = p_col, width = 12, height = 8, dpi = 300)

# Plot Parent 2
p_ler <- ggplot(binned_ler, aes(x = position, y = coverage)) +
    {if (centromere_aware && !is.null(cen_coords_ler))
        geom_rect(data = cen_coords_ler,
                  aes(xmin = cen_start, xmax = cen_end, ymin = -Inf, ymax = Inf),
                  fill = "grey", alpha = 0.3, inherit.aes = FALSE)
    } +
    geom_line() +
    facet_wrap(~ chrom, scales = "free_x") +
    labs(title = paste("Crossover Coverage Plot -", parent2_name, "Reference"),
         subtitle = paste("Sample:", sample_name),
         x = "Chromosome Position", y = "Coverage") +
    theme_bw()

ggsave(paste0(sample_name, "_crossover_plot_", parent2_name, ".svg"),
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

