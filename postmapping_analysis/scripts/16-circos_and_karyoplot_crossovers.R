#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(circlize)
    library(karyoploteR)
    library(GenomicRanges)
    library(dplyr)
    library(argparse)
})

## -------------------------------------------------------------------
## Parse command-line arguments
## -------------------------------------------------------------------
parser <- ArgumentParser(description = "Generate circos and karyotype plots of crossovers")

parser$add_argument("--bed_file", required = TRUE, help = "BED file with segment information (bed_of_best.tsv)")
parser$add_argument("--parent1_genome", required = TRUE, help = "Parent 1 genome FASTA index (.fai)")
parser$add_argument("--parent2_genome", required = TRUE, help = "Parent 2 genome FASTA index (.fai)")
parser$add_argument("--curated_profiles", required = TRUE, help = "Curated hybrid profiles TSV file")
parser$add_argument("--threshold", required = TRUE, type = "integer", help = "Threshold value to filter profiles")
parser$add_argument("--sample_name", required = TRUE, help = "Sample name for output files")
parser$add_argument("--parent1_name", required = FALSE, default = "Parent1", help = "Name of parent 1")
parser$add_argument("--parent2_name", required = FALSE, default = "Parent2", help = "Name of parent 2")
parser$add_argument("--output_dir", required = FALSE, default = ".", help = "Output directory")

args <- parser$parse_args()

# Assign to variables
bed_file <- args$bed_file
parent1_genome_fai <- args$parent1_genome
parent2_genome_fai <- args$parent2_genome
curated_profiles_file <- args$curated_profiles
threshold_val <- args$threshold
sample_name <- args$sample_name
parent1_name <- args$parent1_name
parent2_name <- args$parent2_name
output_dir <- args$output_dir

cat("CROSSOVER VISUALIZATION (CIRCOS + KARYOPLOT)\n")
cat("=============================================\n")
cat("Sample:", sample_name, "\n")
cat("Threshold:", threshold_val, "\n\n")

## -------------------------------------------------------------------
## Read genome information from FASTA index files
## -------------------------------------------------------------------
cat("Reading genome information...\n")

parent1_fai <- read.table(parent1_genome_fai, header = FALSE, stringsAsFactors = FALSE)
parent2_fai <- read.table(parent2_genome_fai, header = FALSE, stringsAsFactors = FALSE)

colnames(parent1_fai) <- c("chr", "length", "offset", "linebases", "linewidth")
colnames(parent2_fai) <- c("chr", "length", "offset", "linebases", "linewidth")

# Filter for main chromosomes only
parent1_chrs <- parent1_fai %>%
    filter(grepl("^chr([0-9]+|X|Y)$", chr)) %>%
    select(chr, length) %>%
    arrange(chr)

parent2_chrs <- parent2_fai %>%
    filter(grepl("^chr([0-9]+|X|Y)$", chr)) %>%
    select(chr, length) %>%
    arrange(chr)

cat("Chromosomes loaded:\n")
cat("  ", parent1_name, ":", nrow(parent1_chrs), "chromosomes\n")
cat("  ", parent2_name, ":", nrow(parent2_chrs), "chromosomes\n\n")

## -------------------------------------------------------------------
## Read BED file and curated profiles
## -------------------------------------------------------------------
cat("Reading segment and profile information...\n")

bed_data <- read.table(bed_file, header = FALSE, stringsAsFactors = FALSE,
                       col.names = c("read", "start", "end", "segment"))

profiles <- read.table(curated_profiles_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
profiles_thr <- profiles %>% filter(threshold == threshold_val)

# Filter for crossover reads
crossover_reads <- profiles_thr %>%
    filter(hybrid_type %in% c("SCO", "NCO", "GC-CO")) %>%
    pull(readname) %>%
    unique()

bed_crossovers <- bed_data %>%
    filter(read %in% crossover_reads) %>%
    mutate(
        parent = ifelse(grepl(paste0("^", substr(parent1_name, 1, 1)), segment), parent1_name, parent2_name),
        chr_num = gsub("^[A-Z]H", "", segment),
        chr = paste0("chr", chr_num)
    )

cat("Crossover reads:", length(crossover_reads), "\n")
cat("Crossover segments:", nrow(bed_crossovers), "\n\n")

# Separate by parent
parent1_segments <- bed_crossovers %>% filter(parent == parent1_name)
parent2_segments <- bed_crossovers %>% filter(parent == parent2_name)

cat("  ", parent1_name, "segments:", nrow(parent1_segments), "\n")
cat("  ", parent2_name, "segments:", nrow(parent2_segments), "\n\n")

## -------------------------------------------------------------------
## Prepare genome data for circos
## -------------------------------------------------------------------
cat("Preparing genome structure for visualization...\n")

# Create renamed chromosome data for circos (to distinguish parents)
parent1_chrs_renamed <- parent1_chrs %>%
    mutate(chr_renamed = paste0(parent1_name, "_", chr),
           start = 0,
           end = length) %>%
    select(chr_renamed, start, end)

parent2_chrs_renamed <- parent2_chrs %>%
    mutate(chr_renamed = paste0(parent2_name, "_", chr),
           start = 0,
           end = length) %>%
    select(chr_renamed, start, end)

genome_df <- rbind(parent1_chrs_renamed, parent2_chrs_renamed)
colnames(genome_df) <- c("chr", "start", "end")

# Set chromosome order (parent1 first, then parent2 reversed)
chr_order <- c(
    paste0(parent1_name, "_", parent1_chrs$chr),
    rev(paste0(parent2_name, "_", parent2_chrs$chr))
)
genome_df$chr <- factor(genome_df$chr, levels = chr_order)

## -------------------------------------------------------------------
## Create link data for circos
## -------------------------------------------------------------------
cat("Creating crossover links...\n")

# Find reads with segments in both parents
reads_with_links <- intersect(parent1_segments$read, parent2_segments$read)

# Create link pairs (first segment from each parent per read)
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

# Merge to ensure same order
link_data <- inner_join(links_p1, links_p2, by = "read", suffix = c("_p1", "_p2"))

cat("Total crossover links:", nrow(link_data), "\n\n")

if (nrow(link_data) == 0) {
    cat("Warning: No links found. Exiting.\n")
    quit(save = "no", status = 0)
}

## -------------------------------------------------------------------
## PLOT 1: CIRCOS PLOT
## -------------------------------------------------------------------
cat("Generating CIRCOS plot...\n")

output_circos <- file.path(output_dir, paste0(sample_name, "_circos_crossovers.pdf"))
pdf(output_circos, width = 12, height = 12)

# Clear circos
circos.clear()

# Color palettes
col_pal_p1 <- "#E41A1C"  # Red
col_pal_p2 <- "#377EB8"  # Blue
link_col <- "#984EA380"  # Purple with transparency

# Set circos parameters
circos.par(
    gap.after = c(rep(2, nrow(parent1_chrs)), rep(2, nrow(parent2_chrs))),
    start.degree = 90,
    track.margin = c(0.01, 0.01),
    cell.padding = c(0.02, 0, 0.02, 0)
)

# Initialize circos
circos.genomicInitialize(genome_df, plotType = NULL)

# Track 1: Chromosome labels with colors
circos.track(
    ylim = c(0, 1),
    track.height = 0.08,
    bg.border = NA,
    panel.fun = function(x, y) {
        sector.index <- get.cell.meta.data("sector.index")
        color <- ifelse(grepl(parent1_name, sector.index), col_pal_p1, col_pal_p2)
        circos.rect(CELL_META$cell.xlim[1], 0,
                    CELL_META$cell.xlim[2], 1,
                    col = color, border = NA)

        # Clean label (remove parent prefix)
        clean_label <- gsub(paste0("^", parent1_name, "_|^", parent2_name, "_"), "", sector.index)
        circos.text(CELL_META$xcenter, CELL_META$ycenter,
                    clean_label,
                    facing = "clockwise", niceFacing = TRUE,
                    col = "white", font = 2, cex = 0.7)
    }
)

# Track 2: Genomic axis
circos.track(
    ylim = c(0, 1),
    track.height = 0.1,
    bg.border = NA,
    panel.fun = function(x, y) {
        circos.genomicAxis(h = "top", major.by = 50e6, labels.cex = 0.4)
    }
)

# Draw links
link_colors <- colorRampPalette(c("#FF6B6B80", "#4ECDC480", "#95E1D380"))(nrow(link_data))
for (i in 1:nrow(link_data)) {
    circos.link(
        link_data$chr_renamed_p1[i],
        c(link_data$start_p1[i], link_data$end_p1[i]),
        link_data$chr_renamed_p2[i],
        c(link_data$start_p2[i], link_data$end_p2[i]),
        col = link_colors[i],
        border = NA
    )
}

# Add title
title(main = paste("Crossover Circos Plot -", sample_name),
      sub = paste("Threshold:", threshold_val, "| Links:", nrow(link_data)),
      cex.main = 1.5)

# Legend
legend("bottomleft",
       legend = c(paste(parent1_name, "chromosomes"),
                  paste(parent2_name, "chromosomes"),
                  "Crossover links"),
       fill = c(col_pal_p1, col_pal_p2, "#984EA3"),
       border = NA, bty = "n", cex = 1)

circos.clear()
dev.off()

cat("Circos plot saved:", output_circos, "\n\n")

## -------------------------------------------------------------------
## PLOT 2: IMPROVED KARYOPLOT
## -------------------------------------------------------------------
cat("Generating KARYOTYPE plot...\n")

output_karyo <- file.path(output_dir, paste0(sample_name, "_karyoplot_crossovers.pdf"))

# Create genome GRanges with renamed chromosomes
genome_granges <- GRanges(
    seqnames = c(
        paste0(parent1_name, "_", parent1_chrs$chr),
        paste0(parent2_name, "_", parent2_chrs$chr)
    ),
    ranges = IRanges(
        start = 1,
        end = c(parent1_chrs$length, parent2_chrs$length)
    )
)

# Create segment GRanges
parent1_granges <- GRanges(
    seqnames = paste0(parent1_name, "_", parent1_segments$chr),
    ranges = IRanges(start = parent1_segments$start, end = parent1_segments$end)
)

parent2_granges <- GRanges(
    seqnames = paste0(parent2_name, "_", parent2_segments$chr),
    ranges = IRanges(start = parent2_segments$start, end = parent2_segments$end)
)

# Create link GRanges
link_granges_p1 <- GRanges(
    seqnames = link_data$chr_renamed_p1,
    ranges = IRanges(start = link_data$start_p1, end = link_data$end_p1)
)

link_granges_p2 <- GRanges(
    seqnames = link_data$chr_renamed_p2,
    ranges = IRanges(start = link_data$start_p2, end = link_data$end_p2)
)

# Order chromosomes
pdf(output_karyo, width = 16, height = 10)

kp <- plotKaryotype(
    genome = genome_granges,
    plot.type = 2,
    main = paste("Crossover Karyotype Map -", sample_name, "- Threshold:", threshold_val),
    cex = 0.8
)

# Plot segments with better visibility
kpPlotRegions(kp, parent1_granges, r0 = 0.05, r1 = 0.45,
              col = "#E41A1C99", border = "#E41A1CCC", lwd = 0.5)
kpPlotRegions(kp, parent2_granges, r0 = 0.05, r1 = 0.45,
              col = "#377EB899", border = "#377EB8CC", lwd = 0.5)

# Plot links with MUCH better visibility
kpPlotLinks(
    kp,
    data = link_granges_p1,
    data2 = link_granges_p2,
    col = "#984EA3CC",  # More opaque purple
    border = "#984EA3",
    lwd = 2,  # Thicker lines
    r0 = 0.5, r1 = 1  # Use upper half for links
)

# Legend
legend("topright",
       legend = c(paste(parent1_name, "segments"),
                  paste(parent2_name, "segments"),
                  "Crossover links"),
       fill = c("#E41A1C", "#377EB8", "#984EA3"),
       border = NA, bty = "n", cex = 1.2)

dev.off()

cat("Karyoplot saved:", output_karyo, "\n\n")

## -------------------------------------------------------------------
## Summary
## -------------------------------------------------------------------
cat("=============================================\n")
cat("VISUALIZATION COMPLETE\n")
cat("=============================================\n")
cat("Circos plot:", output_circos, "\n")
cat("Karyotype plot:", output_karyo, "\n")
cat("Total crossover links plotted:", nrow(link_data), "\n\n")
