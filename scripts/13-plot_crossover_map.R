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
parser$add_argument("--centromere_bed", required = FALSE, default = "NA", help = "Centromere coordinates BED file (optional, for highlighting in centromere-unaware mode)")

args <- parser$parse_args()

# Assign to variables for convenience
col_bed_file <- args$col_bed
ler_bed_file <- args$ler_bed
col_paf_file <- args$col_paf
ler_paf_file <- args$ler_paf
sample_name <- args$sample_name
parent1_name <- args$parent1_name
parent2_name <- args$parent2_name
centromere_bed_file <- args$centromere_bed

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

# Function to read centromere coordinates from BED file
read_centromere_bed <- function(bed_file) {
    cat("Reading centromere coordinates from:", bed_file, "\n")

    if (!file.exists(bed_file)) {
        cat("  Warning: Centromere BED file not found. Skipping centromere highlighting.\n")
        return(NULL)
    }

    # Read BED file (chromosome, start, end)
    cen_bed <- read.table(bed_file, header = FALSE, stringsAsFactors = FALSE)

    # Handle different BED formats
    if (ncol(cen_bed) >= 3) {
        colnames(cen_bed)[1:3] <- c("chrom", "cen_start", "cen_end")
        cen_bed <- cen_bed[, 1:3]
    } else {
        cat("  Warning: Invalid BED format. Expected at least 3 columns. Skipping centromere highlighting.\n")
        return(NULL)
    }

    # Add chr prefix if not present
    if (!grepl("^chr", cen_bed$chrom[1], ignore.case = TRUE)) {
        cen_bed$chrom <- paste0("chr", cen_bed$chrom)
    }

    cat("  Loaded centromere coordinates for", nrow(cen_bed), "chromosomes\n")
    return(cen_bed)
}

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

    # The region column IS the chromosome name (e.g., "chr9", "chr4")
    paf$chrom <- paf$region

    # Extract num_code from read_id (positions 1, 2, or 3)
    paf$num_code <- as.integer(str_match(paf$read_id, "_([123])_")[, 2])
    paf <- paf %>% filter(!is.na(num_code))

    # Calculate position based on num_code
    paf$pos <- ifelse(paf$num_code == 1, paf$end, paf$start)

    # Convert chromosome to factor with natural sorting for proper facet order
    # Extract numeric part for sorting (e.g., "chr9" -> 9, "chr19" -> 19)
    chrom_numbers <- as.integer(gsub("chr|Chr", "", paf$chrom))
    paf$chrom <- factor(paf$chrom, levels = unique(paf$chrom[order(chrom_numbers)]))

    return(paf)
}

# Function to get chromosome lengths from PAF (centromere-unaware mode)
get_chrom_lengths_from_paf <- function(paf_data) {
    if (nrow(paf_data) == 0) {
        cat("  Warning: Empty PAF data, returning empty chromosome lengths.\n")
        return(data.frame(
            chrom = character(0),
            chrom_length = integer(0),
            stringsAsFactors = FALSE
        ))
    }

    chrom_lengths <- paf_data %>%
        group_by(chrom) %>%
        summarise(chrom_length = max(ref_length), .groups = "drop")

    # Validate chromosome lengths
    invalid_chroms <- chrom_lengths %>%
        filter(is.na(chrom_length) | chrom_length <= 0)

    if (nrow(invalid_chroms) > 0) {
        cat("  Warning: Invalid chromosome lengths detected:\n")
        print(invalid_chroms)
        chrom_lengths <- chrom_lengths %>%
            filter(!is.na(chrom_length) & chrom_length > 0)
    }

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


calculate_coverage <- function(paf_merged, chrom_lengths, bin_size = 100000) {
    # Handle empty chromosome lengths
    if (nrow(chrom_lengths) == 0) {
        cat("  Warning: No chromosome length data. Returning empty coverage.\n")
        return(data.frame(
            position = numeric(0),
            coverage = numeric(0),
            chrom = character(0),
            stringsAsFactors = FALSE
        ))
    }

    # Handle empty PAF data - return binned zero coverage
    if (nrow(paf_merged) == 0) {
        cat("  Warning: No PAF data to calculate coverage. Returning binned zero coverage.\n")

        # Sort chromosome lengths numerically
        chrom_lengths$chrom <- as.character(chrom_lengths$chrom)
        chrom_nums <- as.integer(gsub("chr|Chr|chromosome|Chromosome", "", chrom_lengths$chrom))
        chrom_lengths <- chrom_lengths[order(chrom_nums), ]

        coverage_list <- list()
        for (i in 1:nrow(chrom_lengths)) {
            ch <- as.character(chrom_lengths$chrom[i])
            total_length <- as.numeric(chrom_lengths$chrom_length[i])

            if (is.na(total_length) || total_length <= 0) {
                cat("  Warning: Invalid chromosome length for", ch, ":", total_length, ". Skipping.\n")
                next
            }

            # Create binned positions
            num_bins <- ceiling(total_length / bin_size)
            bin_starts <- seq(1, total_length, by = bin_size)
            bin_centers <- bin_starts + (bin_size / 2)
            bin_centers <- bin_centers[1:num_bins]

            df <- data.frame(
                position = bin_centers,
                coverage = 0,
                chrom = ch,
                stringsAsFactors = FALSE
            )
            coverage_list[[ch]] <- df
        }

        if (length(coverage_list) == 0) {
            return(data.frame(
                position = numeric(0),
                coverage = numeric(0),
                chrom = character(0),
                stringsAsFactors = FALSE
            ))
        }

        return(bind_rows(coverage_list))
    }

    # Process each chromosome with binned coverage
    coverage_list <- list()

    # Convert chrom to character to avoid factor issues
    chrom_lengths$chrom <- as.character(chrom_lengths$chrom)
    paf_merged$chrom <- as.character(paf_merged$chrom)

    # Get unique chromosomes and sort them numerically
    chroms <- unique(paf_merged$chrom)
    chrom_nums <- as.integer(gsub("chr|Chr|chromosome|Chromosome", "", chroms))
    chroms <- chroms[order(chrom_nums)]

    for(ch in chroms) {
        ch_data <- paf_merged %>% filter(chrom == ch)

        # Get chromosome length
        chrom_length_subset <- chrom_lengths %>% filter(chrom == ch)

        if (nrow(chrom_length_subset) == 0) {
            cat("  Warning: No length found for chromosome", ch, ". Skipping.\n")
            next
        }

        total_length <- as.numeric(chrom_length_subset$chrom_length[1])

        if (is.na(total_length) || total_length <= 0) {
            cat("  Warning: Invalid chromosome length for", ch, ":", total_length, ". Skipping.\n")
            next
        }

        # Create bins
        num_bins <- ceiling(total_length / bin_size)
        breaks <- seq(0, total_length, by = bin_size)
        if (breaks[length(breaks)] < total_length) {
            breaks <- c(breaks, total_length)
        }

        # Bin the crossover positions
        bin_indices <- cut(ch_data$pos, breaks = breaks, labels = FALSE, include.lowest = TRUE)

        # Count crossovers per bin
        bin_counts <- table(factor(bin_indices, levels = 1:num_bins))

        # Create bin centers for plotting
        bin_centers <- breaks[-length(breaks)] + (bin_size / 2)
        bin_centers <- bin_centers[1:num_bins]

        df <- data.frame(
            position = bin_centers,
            coverage = as.numeric(bin_counts),
            chrom = ch,
            stringsAsFactors = FALSE
        )
        coverage_list[[ch]] <- df
    }

    if (length(coverage_list) == 0) {
        cat("  Warning: No valid chromosomes to process. Returning empty data frame.\n")
        return(data.frame(
            position = numeric(0),
            coverage = numeric(0),
            chrom = character(0),
            stringsAsFactors = FALSE
        ))
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

    # Debug output
    cat("Chromosome lengths for", parent1_name, ":\n")
    print(chrom_lengths_col)
    cat("Chromosome lengths for", parent2_name, ":\n")
    print(chrom_lengths_ler)

    # Optional: Load centromere coordinates for highlighting
    if (centromere_bed_file != "NA" && file.exists(centromere_bed_file)) {
        cen_coords <- read_centromere_bed(centromere_bed_file)
        if (!is.null(cen_coords)) {
            # Use same centromere coords for both references in centromere-unaware mode
            cen_coords_col <- cen_coords
            cen_coords_ler <- cen_coords
            cat("  Centromere highlighting enabled\n")
        } else {
            cen_coords_col <- NULL
            cen_coords_ler <- NULL
        }
    } else {
        # No centromere coordinates to shade
        cen_coords_col <- NULL
        cen_coords_ler <- NULL
        cat("  No centromere BED file provided - plots will not have centromere highlighting\n")
    }
}

# Compute coverage (already binned for centromere-unaware mode)
coverage_col <- calculate_coverage(paf_merged_col, chrom_lengths_col)
coverage_ler <- calculate_coverage(paf_merged_ler, chrom_lengths_ler)

coverage_col$reference <- parent1_name
coverage_ler$reference <- parent2_name

# For plotting, use coverage data directly (already binned in centromere-unaware mode)
# In centromere-aware mode with smaller regions, we may want additional binning
if (centromere_aware) {
    # Additional binning for centromere-aware mode
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
} else {
    # Use data as-is (already binned at 100kb resolution)
    binned_col <- coverage_col
    binned_ler <- coverage_ler
}

# Ensure chrom is properly ordered for faceting
if (!centromere_aware) {
    # For centromere-unaware mode, order chromosomes numerically
    binned_col$chrom <- as.character(binned_col$chrom)
    binned_ler$chrom <- as.character(binned_ler$chrom)

    # Get all unique chromosomes and sort them numerically
    all_chroms <- unique(c(binned_col$chrom, binned_ler$chrom))

    # Create a data frame to sort by numeric chromosome number
    chrom_df <- data.frame(
        chrom = all_chroms,
        num = as.integer(gsub("chr|Chr|chromosome|Chromosome", "", all_chroms)),
        stringsAsFactors = FALSE
    )
    chrom_df <- chrom_df[order(chrom_df$num), ]
    chrom_order <- chrom_df$chrom

    cat("Chromosome order for plotting:", paste(chrom_order, collapse=", "), "\n")

    binned_col$chrom <- factor(binned_col$chrom, levels = chrom_order)
    binned_ler$chrom <- factor(binned_ler$chrom, levels = chrom_order)
}

# Determine optimal plot layout based on number of chromosomes
num_chroms_col <- length(unique(binned_col$chrom))
num_chroms_ler <- length(unique(binned_ler$chrom))
max_chroms <- max(num_chroms_col, num_chroms_ler)

# Calculate optimal dimensions
# For <= 5 chromosomes: 2 columns, smaller plot
# For 6-12 chromosomes: 3 columns, medium plot
# For > 12 chromosomes: 4 columns, larger plot
if (max_chroms <= 5) {
    plot_width <- 12
    plot_height <- 8
    ncol_facet <- 2
} else if (max_chroms <= 12) {
    plot_width <- 15
    plot_height <- 12
    ncol_facet <- 3
} else {
    plot_width <- 20
    plot_height <- 15
    ncol_facet <- 4
}

cat("Detected", max_chroms, "chromosomes. Using", ncol_facet, "columns,", plot_width, "x", plot_height, "plot size.\n")

# Plot Parent 1
p_col <- ggplot(binned_col, aes(x = position, y = coverage)) +
    {if (!is.null(cen_coords_col))
        geom_rect(data = cen_coords_col,
                  aes(xmin = cen_start, xmax = cen_end, ymin = -Inf, ymax = Inf),
                  fill = "grey", alpha = 0.3, inherit.aes = FALSE)
    } +
    geom_line(color = "steelblue", size = 0.8) +
    facet_wrap(~ chrom, scales = "free_x", ncol = ncol_facet) +
    labs(title = paste("Crossover Coverage Plot -", parent1_name, "Reference"),
         subtitle = paste("Sample:", sample_name),
         x = "Chromosome Position (bp)", y = "Crossover Coverage") +
    theme_bw() +
    theme(
        strip.text = element_text(face = "bold", size = 10),
        strip.background = element_rect(fill = "lightgray"),
        panel.grid.minor = element_blank()
    )

ggsave(paste0(sample_name, "_crossover_plot_", parent1_name, ".svg"),
       plot = p_col, width = plot_width, height = plot_height, dpi = 300)

cat("Saved:", paste0(sample_name, "_crossover_plot_", parent1_name, ".svg\n"))

# Plot Parent 2
p_ler <- ggplot(binned_ler, aes(x = position, y = coverage)) +
    {if (!is.null(cen_coords_ler))
        geom_rect(data = cen_coords_ler,
                  aes(xmin = cen_start, xmax = cen_end, ymin = -Inf, ymax = Inf),
                  fill = "grey", alpha = 0.3, inherit.aes = FALSE)
    } +
    geom_line(color = "coral", size = 0.8) +
    facet_wrap(~ chrom, scales = "free_x", ncol = ncol_facet) +
    labs(title = paste("Crossover Coverage Plot -", parent2_name, "Reference"),
         subtitle = paste("Sample:", sample_name),
         x = "Chromosome Position (bp)", y = "Crossover Coverage") +
    theme_bw() +
    theme(
        strip.text = element_text(face = "bold", size = 10),
        strip.background = element_rect(fill = "lightgray"),
        panel.grid.minor = element_blank()
    )

ggsave(paste0(sample_name, "_crossover_plot_", parent2_name, ".svg"),
       plot = p_ler, width = plot_width, height = plot_height, dpi = 300)

cat("Saved:", paste0(sample_name, "_crossover_plot_", parent2_name, ".svg\n"))

# Write output tables
# Add reference column before combining
coverage_col$reference <- parent1_name
coverage_ler$reference <- parent2_name

# Ensure chromosome factors are compatible before binding
coverage_col$chrom <- as.character(coverage_col$chrom)
coverage_ler$chrom <- as.character(coverage_ler$chrom)
binned_col$chrom <- as.character(binned_col$chrom)
binned_ler$chrom <- as.character(binned_ler$chrom)

tryCatch({
    coverage_combined <- bind_rows(coverage_col, coverage_ler)
    binned_combined <- bind_rows(binned_col, binned_ler)

    write.table(coverage_combined,
                paste0(sample_name, "_coverage_data", ".txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)

    write.table(binned_combined,
                paste0(sample_name, "_binned_coverage", ".txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)

    cat("Coverage data tables saved successfully.\n")
}, error = function(e) {
    cat("Warning: Could not save combined coverage tables:", e$message, "\n")
    cat("Individual plots were saved successfully.\n")
})

cat("Crossover plot generation completed!\n")
