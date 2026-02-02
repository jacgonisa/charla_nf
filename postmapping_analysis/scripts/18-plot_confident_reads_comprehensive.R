#!/usr/bin/env Rscript

# Comprehensive plotting script for confident recombinant reads
# Combines crossover coverage plots (from script 13) and circos/karyotype plots (from script 16)
# but filters to only confident reads

suppressPackageStartupMessages({
    library(IRanges)
    library(GenomicRanges)
    library(circlize)
    library(karyoploteR)
    library(dplyr)
    library(stringr)
    library(ggplot2)
    library(argparse)
})

## -------------------------------------------------------------------
## Parse command-line arguments
## -------------------------------------------------------------------
parser <- ArgumentParser(description = "Generate comprehensive plots for confident recombinant reads")

parser$add_argument("--confident_bed_sco", required = TRUE, help = "Confident SCO BED file")
parser$add_argument("--confident_bed_nco", required = TRUE, help = "Confident NCO BED file")
parser$add_argument("--confident_bed_gcco", required = TRUE, help = "Confident GC-CO BED file")
parser$add_argument("--confident_read_ids", required = TRUE, help = "Confident read IDs file")
parser$add_argument("--parent1_paf_chr", required = TRUE, help = "Parent 1 CHR PAF file")
parser$add_argument("--parent2_paf_chr", required = TRUE, help = "Parent 2 CHR PAF file")
parser$add_argument("--parent1_genome_fai", required = TRUE, help = "Parent 1 genome FASTA index")
parser$add_argument("--parent2_genome_fai", required = TRUE, help = "Parent 2 genome FASTA index")
parser$add_argument("--cowidths_file", required = TRUE, help = "Cowidths file")
parser$add_argument("--curated_profiles", required = TRUE, help = "Curated hybrid profiles TSV")
parser$add_argument("--sample_name", required = TRUE, help = "Sample name")
parser$add_argument("--parent1_name", required = FALSE, default = "Parent1", help = "Parent 1 name")
parser$add_argument("--parent2_name", required = FALSE, default = "Parent2", help = "Parent 2 name")
parser$add_argument("--threshold", required = TRUE, type = "integer", help = "Threshold value")
parser$add_argument("--centromere_bed", required = FALSE, default = "NA", help = "Centromere BED (optional)")
parser$add_argument("--output_dir", required = FALSE, default = ".", help = "Output directory")

args <- parser$parse_args()

# Assign variables
confident_bed_sco <- args$confident_bed_sco
confident_bed_nco <- args$confident_bed_nco
confident_bed_gcco <- args$confident_bed_gcco
confident_read_ids_file <- args$confident_read_ids
parent1_paf_chr <- args$parent1_paf_chr
parent2_paf_chr <- args$parent2_paf_chr
parent1_genome_fai <- args$parent1_genome_fai
parent2_genome_fai <- args$parent2_genome_fai
cowidths_file <- args$cowidths_file
curated_profiles_file <- args$curated_profiles
sample_name <- args$sample_name
parent1_name <- args$parent1_name
parent2_name <- args$parent2_name
threshold_val <- args$threshold
centromere_bed_file <- args$centromere_bed
output_dir <- args$output_dir

## -------------------------------------------------------------------
## Setup
## -------------------------------------------------------------------
cat("\n========================================\n")
cat("COMPREHENSIVE CONFIDENT READS VISUALIZATION\n")
cat("========================================\n")
cat("Sample:", sample_name, "\n")
cat("Threshold:", threshold_val, "\n")
cat("Parents:", parent1_name, "x", parent2_name, "\n")

# Create output directory
plot_dir <- file.path(output_dir, "confident_reads_plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
cat("Output directory:", plot_dir, "\n")

# Read confident read IDs
cat("\nReading confident read IDs...\n")
confident_reads <- readLines(confident_read_ids_file)
cat("Total confident reads:", length(confident_reads), "\n")

# Read confident BED files
cat("\nReading confident BED files...\n")
bed_sco <- read.table(confident_bed_sco, header = FALSE, stringsAsFactors = FALSE)
bed_nco <- read.table(confident_bed_nco, header = FALSE, stringsAsFactors = FALSE)
bed_gcco <- read.table(confident_bed_gcco, header = FALSE, stringsAsFactors = FALSE)

# BED format: read_id, start, end, segment
colnames(bed_sco) <- c("read", "start", "end", "segment")
colnames(bed_nco) <- c("read", "start", "end", "segment")
colnames(bed_gcco) <- c("read", "start", "end", "segment")

bed_sco$recomb_type <- "SCO"
bed_nco$recomb_type <- "NCO"
bed_gcco$recomb_type <- "GC-CO"

bed_all <- rbind(bed_sco, bed_nco, bed_gcco)
cat("Total confident segments:", nrow(bed_all), "\n")

## -------------------------------------------------------------------
## PART 1: CROSSOVER COVERAGE PLOTS (based on script 13)
## -------------------------------------------------------------------
cat("\n=== GENERATING CROSSOVER COVERAGE PLOTS ===\n")

# Function to process PAF and filter for confident reads (based on script 13)
process_paf_confident <- function(paf_file, confident_read_ids) {
    cat("Processing PAF file:", basename(paf_file), "\n")

    if (!file.exists(paf_file)) {
        cat("  Warning: File not found:", paf_file, "\n")
        return(data.frame(
            read_id = character(0), region = character(0), ref_length = integer(0),
            start = integer(0), end = integer(0), chrom = character(0),
            num_code = integer(0), pos = integer(0), stringsAsFactors = FALSE
        ))
    }

    # Check if file is empty
    lines <- readLines(paf_file, warn = FALSE)
    data_lines <- lines[!grepl("^#", lines)]

    if (length(data_lines) == 0) {
        cat("  Warning: PAF file is empty\n")
        return(data.frame(
            read_id = character(0), region = character(0), ref_length = integer(0),
            start = integer(0), end = integer(0), chrom = character(0),
            num_code = integer(0), pos = integer(0), stringsAsFactors = FALSE
        ))
    }

    # Read PAF
    paf <- read.table(paf_file, header = FALSE, stringsAsFactors = FALSE, sep = "\t", comment.char = "#")
    colnames(paf) <- c("read_id", "region", "ref_length", "start", "end")

    cat("  Before filtering: ", nrow(paf), " alignments\n", sep = "")

    # Filter for confident reads only
    # Extract base read ID (handle PAF format: #UUID_NUM_CODE|best:full -> extract UUID)
    paf$base_read_id <- sub("^#", "", paf$read_id)  # Remove leading #
    paf$base_read_id <- sub("_[0-9]+.*$", "", paf$base_read_id)  # Remove _NUM and everything after
    paf_confident <- paf %>% filter(base_read_id %in% confident_read_ids)

    cat("  After filtering: ", nrow(paf_confident), " alignments (confident reads only)\n", sep = "")

    if (nrow(paf_confident) == 0) {
        return(data.frame(
            read_id = character(0), region = character(0), ref_length = integer(0),
            start = integer(0), end = integer(0), chrom = character(0),
            num_code = integer(0), pos = integer(0), stringsAsFactors = FALSE
        ))
    }

    # The region column IS the chromosome name
    paf_confident$chrom <- paf_confident$region

    # Extract num_code from read_id (positions 1, 2, or 3)
    paf_confident$num_code <- as.integer(str_match(paf_confident$read_id, "_([123])_")[, 2])
    paf_confident <- paf_confident %>% filter(!is.na(num_code))

    # Calculate position based on num_code
    paf_confident$pos <- ifelse(paf_confident$num_code == 1, paf_confident$end, paf_confident$start)

    # Convert chromosome to factor with proper numeric sorting
    # Extract numeric part and create properly ordered factor levels
    paf_confident$chrom_num <- as.integer(gsub("chr|Chr", "", paf_confident$chrom))
    # Create ordered levels based on numeric values
    chrom_levels <- paste0("chr", sort(unique(paf_confident$chrom_num)))
    paf_confident$chrom <- factor(paf_confident$chrom, levels = chrom_levels)
    paf_confident$chrom_num <- NULL  # Remove temporary column

    return(paf_confident)
}

# Function to get chromosome lengths (from script 13)
get_chrom_lengths_from_paf <- function(paf_data) {
    if (nrow(paf_data) == 0) {
        return(data.frame(chrom = character(0), chrom_length = integer(0), stringsAsFactors = FALSE))
    }

    chrom_lengths <- paf_data %>%
        group_by(chrom) %>%
        summarise(chrom_length = max(ref_length), .groups = "drop")

    return(chrom_lengths)
}

# Function to calculate coverage (from script 13)
calculate_coverage <- function(paf_data, chrom_lengths, bin_size = 1e6) {
    if (nrow(paf_data) == 0 || nrow(chrom_lengths) == 0) {
        return(data.frame(position = numeric(0), coverage = numeric(0), chrom = character(0), stringsAsFactors = FALSE))
    }

    # Store original factor levels
    chrom_levels <- levels(paf_data$chrom)

    coverage_list <- list()
    for (ch in chrom_levels) {
        ch_data <- paf_data %>% filter(chrom == ch)
        if (nrow(ch_data) == 0) next

        chrom_length_subset <- chrom_lengths %>% filter(chrom == ch)
        if (nrow(chrom_length_subset) == 0) next

        total_length <- as.numeric(chrom_length_subset$chrom_length[1])
        if (is.na(total_length) || total_length <= 0) next

        num_bins <- ceiling(total_length / bin_size)
        breaks <- seq(0, total_length, by = bin_size)
        if (breaks[length(breaks)] < total_length) {
            breaks <- c(breaks, total_length)
        }

        bin_indices <- cut(ch_data$pos, breaks = breaks, labels = FALSE, include.lowest = TRUE)
        bin_counts <- table(factor(bin_indices, levels = 1:num_bins))
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
        return(data.frame(position = numeric(0), coverage = numeric(0), chrom = character(0), stringsAsFactors = FALSE))
    }

    result <- bind_rows(coverage_list)
    # Restore factor levels in correct order
    result$chrom <- factor(result$chrom, levels = chrom_levels)
    return(result)
}

# Process PAF files for confident reads
paf_p1_confident <- process_paf_confident(parent1_paf_chr, confident_reads)
paf_p2_confident <- process_paf_confident(parent2_paf_chr, confident_reads)

# Get chromosome lengths
chrom_lengths_p1 <- get_chrom_lengths_from_paf(paf_p1_confident)
chrom_lengths_p2 <- get_chrom_lengths_from_paf(paf_p2_confident)

# Calculate coverage
coverage_p1 <- calculate_coverage(paf_p1_confident, chrom_lengths_p1)
coverage_p2 <- calculate_coverage(paf_p2_confident, chrom_lengths_p2)

# Sort coverage data by chromosome factor levels to ensure correct facet order
coverage_p1 <- coverage_p1 %>% arrange(chrom, position)
coverage_p2 <- coverage_p2 %>% arrange(chrom, position)

coverage_p1$reference <- parent1_name
coverage_p2$reference <- parent2_name

# Determine plot layout
if (nrow(coverage_p1) > 0 && nrow(coverage_p2) > 0) {
    num_chroms <- max(length(unique(coverage_p1$chrom)), length(unique(coverage_p2$chrom)))

    if (num_chroms <= 5) {
        plot_width <- 12; plot_height <- 8; ncol_facet <- 2
    } else if (num_chroms <= 12) {
        plot_width <- 15; plot_height <- 12; ncol_facet <- 3
    } else {
        plot_width <- 20; plot_height <- 15; ncol_facet <- 4
    }

    # Load centromere coordinates if provided
    cen_coords_p1 <- NULL
    cen_coords_p2 <- NULL
    if (centromere_bed_file != "NA" && file.exists(centromere_bed_file)) {
        cen_bed <- read.table(centromere_bed_file, header = FALSE, stringsAsFactors = FALSE)
        if (ncol(cen_bed) >= 3) {
            colnames(cen_bed)[1:3] <- c("chrom", "cen_start", "cen_end")
            if (!grepl("^chr", cen_bed$chrom[1], ignore.case = TRUE)) {
                cen_bed$chrom <- paste0("chr", cen_bed$chrom)
            }
            cen_coords_p1 <- cen_bed
            cen_coords_p2 <- cen_bed
            cat("  Centromere highlighting enabled\n")
        }
    }

    # Plot Parent 1
    cat("Generating coverage plot for", parent1_name, "...\n")
    if (nrow(coverage_p1) > 0) {
        p_p1 <- ggplot(coverage_p1, aes(x = position, y = coverage)) +
            {if (!is.null(cen_coords_p1))
                geom_rect(data = cen_coords_p1,
                         aes(xmin = cen_start, xmax = cen_end, ymin = -Inf, ymax = Inf),
                         fill = "grey", alpha = 0.3, inherit.aes = FALSE)
            } +
            geom_line(color = "steelblue", linewidth = 0.8) +
            facet_wrap(~chrom, scales = "free_x", ncol = ncol_facet, drop = FALSE) +
            labs(title = paste("Confident Reads Crossover Coverage -", parent1_name, "Reference"),
                 subtitle = paste("Sample:", sample_name, "| Confident reads:", length(confident_reads)),
                 x = "Chromosome Position (bp)", y = "Crossover Coverage") +
            theme_bw() +
            theme(
                strip.text = element_text(face = "bold", size = 10),
                plot.title = element_text(face = "bold", size = 14),
                axis.text.x = element_text(angle = 45, hjust = 1)
            )

        outfile_p1 <- file.path(plot_dir, paste0(sample_name, "_confident_crossover_plot_", parent1_name, ".svg"))
        ggsave(outfile_p1, plot = p_p1, width = plot_width, height = plot_height, dpi = 300)
        cat("  Saved:", basename(outfile_p1), "\n")
    }

    # Plot Parent 2
    cat("Generating coverage plot for", parent2_name, "...\n")
    if (nrow(coverage_p2) > 0) {
        p_p2 <- ggplot(coverage_p2, aes(x = position, y = coverage)) +
            {if (!is.null(cen_coords_p2))
                geom_rect(data = cen_coords_p2,
                         aes(xmin = cen_start, xmax = cen_end, ymin = -Inf, ymax = Inf),
                         fill = "grey", alpha = 0.3, inherit.aes = FALSE)
            } +
            geom_line(color = "darkorange", linewidth = 0.8) +
            facet_wrap(~chrom, scales = "free_x", ncol = ncol_facet, drop = FALSE) +
            labs(title = paste("Confident Reads Crossover Coverage -", parent2_name, "Reference"),
                 subtitle = paste("Sample:", sample_name, "| Confident reads:", length(confident_reads)),
                 x = "Chromosome Position (bp)", y = "Crossover Coverage") +
            theme_bw() +
            theme(
                strip.text = element_text(face = "bold", size = 10),
                plot.title = element_text(face = "bold", size = 14),
                axis.text.x = element_text(angle = 45, hjust = 1)
            )

        outfile_p2 <- file.path(plot_dir, paste0(sample_name, "_confident_crossover_plot_", parent2_name, ".svg"))
        ggsave(outfile_p2, plot = p_p2, width = plot_width, height = plot_height, dpi = 300)
        cat("  Saved:", basename(outfile_p2), "\n")
    }
}

## -------------------------------------------------------------------
## PART 2: CIRCOS AND KARYOTYPE PLOTS (based on script 16)
## -------------------------------------------------------------------
cat("\n=== GENERATING CIRCOS AND KARYOTYPE PLOTS ===\n")

# Read genome information
cat("Reading genome information...\n")
parent1_fai <- read.table(parent1_genome_fai, header = FALSE, stringsAsFactors = FALSE)
parent2_fai <- read.table(parent2_genome_fai, header = FALSE, stringsAsFactors = FALSE)

colnames(parent1_fai) <- c("chr", "length", "offset", "linebases", "linewidth")
colnames(parent2_fai) <- c("chr", "length", "offset", "linebases", "linewidth")

# Filter for main chromosomes
parent1_chrs <- parent1_fai %>%
    filter(grepl("^chr([0-9]+|X|Y)$", chr)) %>%
    select(chr, length) %>%
    arrange(chr)

parent2_chrs <- parent2_fai %>%
    filter(grepl("^chr([0-9]+|X|Y)$", chr)) %>%
    select(chr, length) %>%
    arrange(chr)

cat(sprintf("  %s: %d chromosomes\n", parent1_name, nrow(parent1_chrs)))
cat(sprintf("  %s: %d chromosomes\n", parent2_name, nrow(parent2_chrs)))

# Process BED data
cat("Processing confident segment data...\n")
bed_all <- bed_all %>%
    mutate(
        parent = ifelse(grepl(paste0("^", substr(parent1_name, 1, 1)), segment), parent1_name, parent2_name),
        chr_num = gsub("^[A-Z]H", "", segment),
        chr = paste0("chr", chr_num)
    )

parent1_segments <- bed_all %>% filter(parent == parent1_name)
parent2_segments <- bed_all %>% filter(parent == parent2_name)

cat(sprintf("  %s segments: %d\n", parent1_name, nrow(parent1_segments)))
cat(sprintf("  %s segments: %d\n", parent2_name, nrow(parent2_segments)))

# Prepare genome data for circos
parent1_chrs_renamed <- parent1_chrs %>%
    mutate(chr_renamed = paste0(parent1_name, "_", chr), start = 0, end = length) %>%
    select(chr_renamed, start, end)

parent2_chrs_renamed <- parent2_chrs %>%
    mutate(chr_renamed = paste0(parent2_name, "_", chr), start = 0, end = length) %>%
    select(chr_renamed, start, end)

genome_df <- rbind(parent1_chrs_renamed, parent2_chrs_renamed)
colnames(genome_df) <- c("chr", "start", "end")

chr_order <- c(
    paste0(parent1_name, "_", parent1_chrs$chr),
    rev(paste0(parent2_name, "_", parent2_chrs$chr))
)
genome_df$chr <- factor(genome_df$chr, levels = chr_order)

# Create link data
cat("Creating crossover links...\n")
reads_with_links <- intersect(parent1_segments$read, parent2_segments$read)
cat(sprintf("  Reads with segments in both parents: %d\n", length(reads_with_links)))

if (length(reads_with_links) > 0) {
    links_p1 <- parent1_segments %>%
        filter(read %in% reads_with_links) %>%
        group_by(read) %>%
        slice(1) %>%
        ungroup() %>%
        mutate(chr_renamed = paste0(parent1_name, "_", chr)) %>%
        select(chr_renamed, start, end, read)

    links_p2 <- parent2_segments %>%
        filter(read %in% reads_with_links) %>%
        group_by(read) %>%
        slice(1) %>%
        ungroup() %>%
        mutate(chr_renamed = paste0(parent2_name, "_", chr)) %>%
        select(chr_renamed, start, end, read)

    link_data <- inner_join(links_p1, links_p2, by = "read", suffix = c("_p1", "_p2"))
    cat(sprintf("  Total crossover links: %d\n", nrow(link_data)))

    # For karyoplot visualization, sample links if there are too many
    # (to make individual links visible)
    link_data_karyo <- link_data
    if (nrow(link_data) > 500) {
        set.seed(42)  # For reproducibility
        link_data_karyo <- link_data[sample(nrow(link_data), 500), ]
        cat(sprintf("  Sampling %d links for karyoplot visualization (out of %d total)\n",
                   nrow(link_data_karyo), nrow(link_data)))
    }

    ## CIRCOS PLOT (exactly like script 16)
    cat("\nGenerating CIRCOS plot...\n")
    output_circos <- file.path(plot_dir, paste0(sample_name, "_confident_circos_crossovers.pdf"))
    pdf(output_circos, width = 12, height = 12)

    circos.clear()

    col_pal_p1 <- "#E41A1C"
    col_pal_p2 <- "#377EB8"

    circos.par(
        gap.after = c(rep(2, nrow(parent1_chrs)), rep(2, nrow(parent2_chrs))),
        start.degree = 90,
        track.margin = c(0.01, 0.01),
        cell.padding = c(0.02, 0, 0.02, 0)
    )

    circos.genomicInitialize(genome_df, plotType = NULL)

    circos.track(
        ylim = c(0, 1),
        track.height = 0.08,
        bg.border = NA,
        panel.fun = function(x, y) {
            sector.index <- get.cell.meta.data("sector.index")
            color <- ifelse(grepl(parent1_name, sector.index), col_pal_p1, col_pal_p2)
            circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1, col = color, border = NA)

            clean_label <- gsub(paste0("^", parent1_name, "_|^", parent2_name, "_"), "", sector.index)
            circos.text(CELL_META$xcenter, CELL_META$ycenter, clean_label,
                       facing = "clockwise", niceFacing = TRUE, col = "white", font = 2, cex = 0.7)
        }
    )

    circos.track(
        ylim = c(0, 1),
        track.height = 0.1,
        bg.border = NA,
        panel.fun = function(x, y) {
            circos.genomicAxis(h = "top", major.by = 50e6, labels.cex = 0.4)
        }
    )

    link_colors <- colorRampPalette(c("#FF6B6B80", "#4ECDC480", "#95E1D380"))(nrow(link_data))
    for (i in 1:nrow(link_data)) {
        circos.link(
            link_data$chr_renamed_p1[i], c(link_data$start_p1[i], link_data$end_p1[i]),
            link_data$chr_renamed_p2[i], c(link_data$start_p2[i], link_data$end_p2[i]),
            col = link_colors[i], border = NA
        )
    }

    title(main = paste("Confident Reads Circos Plot -", sample_name),
          sub = paste("Threshold:", threshold_val, "| Links:", nrow(link_data)), cex.main = 1.5)

    legend("bottomleft",
           legend = c(paste(parent1_name, "chromosomes"), paste(parent2_name, "chromosomes"), "Crossover links"),
           fill = c(col_pal_p1, col_pal_p2, "#984EA3"), border = NA, bty = "n", cex = 1)

    circos.clear()
    dev.off()
    cat(sprintf("  Saved: %s\n", basename(output_circos)))

    ## KARYOTYPE PLOT (exactly like script 16)
    cat("\nGenerating KARYOTYPE plot...\n")
    output_karyo <- file.path(plot_dir, paste0(sample_name, "_confident_karyoplot_crossovers.pdf"))

    genome_granges <- GRanges(
        seqnames = c(paste0(parent1_name, "_", parent1_chrs$chr), paste0(parent2_name, "_", parent2_chrs$chr)),
        ranges = IRanges(start = 1, end = c(parent1_chrs$length, parent2_chrs$length))
    )

    parent1_granges <- GRanges(
        seqnames = paste0(parent1_name, "_", parent1_segments$chr),
        ranges = IRanges(start = parent1_segments$start, end = parent1_segments$end)
    )

    parent2_granges <- GRanges(
        seqnames = paste0(parent2_name, "_", parent2_segments$chr),
        ranges = IRanges(start = parent2_segments$start, end = parent2_segments$end)
    )

    # Use sampled links for better visibility
    link_granges_p1 <- GRanges(
        seqnames = link_data_karyo$chr_renamed_p1,
        ranges = IRanges(start = link_data_karyo$start_p1, end = link_data_karyo$end_p1)
    )

    link_granges_p2 <- GRanges(
        seqnames = link_data_karyo$chr_renamed_p2,
        ranges = IRanges(start = link_data_karyo$start_p2, end = link_data_karyo$end_p2)
    )

    pdf(output_karyo, width = 16, height = 10)

    kp <- plotKaryotype(
        genome = genome_granges,
        plot.type = 2,
        main = paste("Confident Reads Karyotype Map -", sample_name, "- Threshold:", threshold_val),
        cex = 0.8
    )

    kpPlotRegions(kp, parent1_granges, r0 = 0.05, r1 = 0.45, col = "#E41A1C99", border = "#E41A1CCC", lwd = 0.5)
    kpPlotRegions(kp, parent2_granges, r0 = 0.05, r1 = 0.45, col = "#377EB899", border = "#377EB8CC", lwd = 0.5)

    # Plot links with HIGH visibility - use bright red for maximum contrast
    kpPlotLinks(kp, data = link_granges_p1, data2 = link_granges_p2,
                col = "#FF0000BB", border = "#FF0000", lwd = 2, r0 = 0.5, r1 = 1)

    legend("topright",
           legend = c(paste(parent1_name, "segments"), paste(parent2_name, "segments"),
                     sprintf("Crossover links (showing %d of %d)", nrow(link_data_karyo), nrow(link_data))),
           fill = c("#E41A1C", "#377EB8", "#FF0000"), border = NA, bty = "n", cex = 1.2)

    dev.off()
    cat(sprintf("  Saved: %s\n", basename(output_karyo)))
}

cat("\n========================================\n")
cat("VISUALIZATION COMPLETE\n")
cat("========================================\n")
cat(sprintf("Output directory: %s\n", plot_dir))
cat(sprintf("Confident reads visualized: %d\n", length(confident_reads)))
cat("\n")
