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
parser <- ArgumentParser(description = "Generate single-chromosome crossover plots")

parser$add_argument("--bed_file", required = TRUE, help = "BED file with segment information")
parser$add_argument("--parent1_genome", required = TRUE, help = "Parent 1 genome FASTA index (.fai)")
parser$add_argument("--parent2_genome", required = TRUE, help = "Parent 2 genome FASTA index (.fai)")
parser$add_argument("--curated_profiles", required = TRUE, help = "Curated hybrid profiles TSV file")
parser$add_argument("--threshold", required = TRUE, type = "integer", help = "Threshold value")
parser$add_argument("--sample_name", required = TRUE, help = "Sample name for output files")
parser$add_argument("--parent1_name", required = FALSE, default = "Parent1", help = "Name of parent 1")
parser$add_argument("--parent2_name", required = FALSE, default = "Parent2", help = "Name of parent 2")
parser$add_argument("--chromosome", required = FALSE, default = "1", help = "Chromosome number to plot (default: 1)")
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
chr_to_plot <- args$chromosome
output_dir <- args$output_dir

cat("SINGLE CHROMOSOME CROSSOVER VISUALIZATION\n")
cat("==========================================\n")
cat("Sample:", sample_name, "\n")
cat("Chromosome:", chr_to_plot, "\n")
cat("Threshold:", threshold_val, "\n\n")

## -------------------------------------------------------------------
## Read genome information
## -------------------------------------------------------------------
parent1_fai <- read.table(parent1_genome_fai, header = FALSE, stringsAsFactors = FALSE)
parent2_fai <- read.table(parent2_genome_fai, header = FALSE, stringsAsFactors = FALSE)
colnames(parent1_fai) <- c("chr", "length", "offset", "linebases", "linewidth")
colnames(parent2_fai) <- c("chr", "length", "offset", "linebases", "linewidth")

chr_name <- paste0("chr", chr_to_plot)

# Get chromosome info for selected chromosome
p1_chr <- parent1_fai %>% filter(chr == chr_name)
p2_chr <- parent2_fai %>% filter(chr == chr_name)

if (nrow(p1_chr) == 0 || nrow(p2_chr) == 0) {
    cat("ERROR: Chromosome", chr_name, "not found in genome files!\n")
    quit(save = "no", status = 1)
}

cat("Chromosome lengths:\n")
cat("  ", parent1_name, chr_name, ":", format(p1_chr$length, big.mark = ","), "bp\n")
cat("  ", parent2_name, chr_name, ":", format(p2_chr$length, big.mark = ","), "bp\n\n")

## -------------------------------------------------------------------
## Read data and filter for selected chromosome
## -------------------------------------------------------------------
cat("Reading segments and profiles...\n")

bed_data <- read.table(bed_file, header = FALSE, stringsAsFactors = FALSE,
                       col.names = c("read", "start", "end", "segment"))

profiles <- read.table(curated_profiles_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
profiles_thr <- profiles %>% filter(threshold == threshold_val)

crossover_reads <- profiles_thr %>%
    filter(hybrid_type %in% c("SCO", "NCO", "GC-CO")) %>%
    pull(readname) %>%
    unique()

# Filter for selected chromosome only
bed_crossovers <- bed_data %>%
    filter(read %in% crossover_reads) %>%
    mutate(
        parent = ifelse(grepl("^C", segment), parent1_name, parent2_name),
        chr_num = gsub("^[CD]H", "", segment),
        chr = paste0("chr", chr_num)
    ) %>%
    filter(chr == chr_name)  # Only selected chromosome

cat("Segments on", chr_name, ":", nrow(bed_crossovers), "\n")

# Separate by parent
parent1_segments <- bed_crossovers %>% filter(parent == parent1_name)
parent2_segments <- bed_crossovers %>% filter(parent == parent2_name)

cat("  ", parent1_name, "segments:", nrow(parent1_segments), "\n")
cat("  ", parent2_name, "segments:", nrow(parent2_segments), "\n\n")

# Create links
reads_with_links <- intersect(parent1_segments$read, parent2_segments$read)
cat("Reads with crossovers on", chr_name, ":", length(reads_with_links), "\n\n")

if (length(reads_with_links) == 0) {
    cat("No crossover links found on chromosome", chr_to_plot, "\n")
    cat("Try a different chromosome or check your data.\n")
    quit(save = "no", status = 0)
}

# Create link data
links_p1 <- parent1_segments %>%
    filter(read %in% reads_with_links) %>%
    select(read, start, end)

links_p2 <- parent2_segments %>%
    filter(read %in% reads_with_links) %>%
    select(read, start, end)

link_data <- inner_join(links_p1, links_p2, by = "read", suffix = c("_p1", "_p2"))

cat("Total crossover links to plot:", nrow(link_data), "\n\n")

## -------------------------------------------------------------------
## PLOT 1: Simple Circos Plot
## -------------------------------------------------------------------
cat("Generating CIRCOS plot...\n")

output_circos <- file.path(output_dir, paste0(sample_name, "_chr", chr_to_plot, "_circos.pdf"))
pdf(output_circos, width = 10, height = 10)

circos.clear()

# Simple genome with just two sectors (parent1 and parent2 for this chromosome)
genome_df <- data.frame(
    chr = c(paste0(parent1_name, "_", chr_name), paste0(parent2_name, "_", chr_name)),
    start = c(0, 0),
    end = c(p1_chr$length, p2_chr$length)
)

circos.par(start.degree = 90, gap.after = c(20, 20))

circos.genomicInitialize(genome_df, plotType = NULL)

# Track 1: Labels
circos.track(ylim = c(0, 1), track.height = 0.1, bg.border = NA,
    panel.fun = function(x, y) {
        sector.index <- get.cell.meta.data("sector.index")
        circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1,
                    col = ifelse(grepl(parent1_name, sector.index), "#E41A1C", "#377EB8"),
                    border = NA)
        circos.text(CELL_META$xcenter, CELL_META$ycenter, sector.index,
                    col = "white", font = 2, cex = 1.2)
    }
)

# Track 2: Axis
circos.track(ylim = c(0, 1), track.height = 0.15, bg.border = NA,
    panel.fun = function(x, y) {
        circos.genomicAxis(h = "top", major.by = 10e6, labels.cex = 0.6)
    }
)

# Draw links with different colors for each link
link_colors <- rainbow(nrow(link_data), alpha = 0.7)
for (i in 1:nrow(link_data)) {
    circos.link(
        paste0(parent1_name, "_", chr_name),
        c(link_data$start_p1[i], link_data$end_p1[i]),
        paste0(parent2_name, "_", chr_name),
        c(link_data$start_p2[i], link_data$end_p2[i]),
        col = link_colors[i],
        border = "black",
        lwd = 2
    )
}

title(main = paste("Crossovers on Chromosome", chr_to_plot, "-", sample_name),
      sub = paste(nrow(link_data), "crossover events"),
      cex.main = 1.5)

legend("bottomleft",
       legend = c(paste(parent1_name, chr_name),
                  paste(parent2_name, chr_name),
                  "Crossover links"),
       fill = c("#E41A1C", "#377EB8", NA),
       border = c(NA, NA, "black"),
       lty = c(0, 0, 1), lwd = c(0, 0, 2),
       bty = "n", cex = 1)

circos.clear()
dev.off()

cat("Circos saved:", output_circos, "\n\n")

## -------------------------------------------------------------------
## PLOT 2: Karyotype-style Linear Plot
## -------------------------------------------------------------------
cat("Generating LINEAR plot...\n")

output_linear <- file.path(output_dir, paste0(sample_name, "_chr", chr_to_plot, "_linear.pdf"))

# Create genome for just this chromosome
genome_granges <- GRanges(
    seqnames = c(paste0(parent1_name, "_", chr_name), paste0(parent2_name, "_", chr_name)),
    ranges = IRanges(start = c(1, 1), end = c(p1_chr$length, p2_chr$length))
)

# Create segment GRanges
parent1_granges <- GRanges(
    seqnames = paste0(parent1_name, "_", chr_name),
    ranges = IRanges(start = parent1_segments$start, end = parent1_segments$end)
)

parent2_granges <- GRanges(
    seqnames = paste0(parent2_name, "_", chr_name),
    ranges = IRanges(start = parent2_segments$start, end = parent2_segments$end)
)

# Link GRanges
link_granges_p1 <- GRanges(
    seqnames = paste0(parent1_name, "_", chr_name),
    ranges = IRanges(start = link_data$start_p1, end = link_data$end_p1)
)

link_granges_p2 <- GRanges(
    seqnames = paste0(parent2_name, "_", chr_name),
    ranges = IRanges(start = link_data$start_p2, end = link_data$end_p2)
)

pdf(output_linear, width = 14, height = 6)

kp <- plotKaryotype(
    genome = genome_granges,
    plot.type = 2,
    main = paste("Chromosome", chr_to_plot, "Crossovers -", sample_name, "(", nrow(link_data), "events)"),
    cex = 1.2
)

# Plot segments
kpPlotRegions(kp, parent1_granges, r0 = 0.05, r1 = 0.3,
              col = "#E41A1CDD", border = "#E41A1C", lwd = 1)
kpPlotRegions(kp, parent2_granges, r0 = 0.05, r1 = 0.3,
              col = "#377EB8DD", border = "#377EB8", lwd = 1)

# Plot links - THICK and VISIBLE
kpPlotLinks(kp,
    data = link_granges_p1,
    data2 = link_granges_p2,
    col = "#984EA3",
    border = "#984EA3",
    lwd = 3,  # Very thick
    r0 = 0.35, r1 = 0.95
)

legend("topright",
       legend = c(paste(parent1_name, "segments"),
                  paste(parent2_name, "segments"),
                  "Crossover links"),
       fill = c("#E41A1C", "#377EB8", "#984EA3"),
       border = NA, bty = "n", cex = 1.2)

dev.off()

cat("Linear plot saved:", output_linear, "\n\n")

## -------------------------------------------------------------------
## Summary
## -------------------------------------------------------------------
cat("==========================================\n")
cat("VISUALIZATION COMPLETE\n")
cat("==========================================\n")
cat("Chromosome:", chr_to_plot, "\n")
cat("Crossover events plotted:", nrow(link_data), "\n")
cat("Output files:\n")
cat("  Circos:", output_circos, "\n")
cat("  Linear:", output_linear, "\n\n")

# Print some example links for debugging
cat("Example crossover links (first 5):\n")
print(head(link_data, 5))
cat("\n")
