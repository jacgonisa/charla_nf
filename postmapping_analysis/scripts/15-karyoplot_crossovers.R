#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(karyoploteR)
    library(GenomicRanges)
    library(dplyr)
    library(argparse)
})

## -------------------------------------------------------------------
## Parse command-line arguments
## -------------------------------------------------------------------
parser <- ArgumentParser(description = "Generate karyotype plot of crossovers across genome")

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

cat("KARYOTYPE CROSSOVER VISUALIZATION\n")
cat("==================================\n")
cat("Sample:", sample_name, "\n")
cat("Threshold:", threshold_val, "\n")
cat("BED file:", bed_file, "\n")
cat("Parent1 genome:", parent1_genome_fai, "\n")
cat("Parent2 genome:", parent2_genome_fai, "\n\n")

## -------------------------------------------------------------------
## Read genome information from FASTA index files
## -------------------------------------------------------------------
cat("Reading genome information...\n")

# Read FASTA index files (.fai format: chr, length, offset, linebases, linewidth)
parent1_fai <- read.table(parent1_genome_fai, header = FALSE, stringsAsFactors = FALSE)
parent2_fai <- read.table(parent2_genome_fai, header = FALSE, stringsAsFactors = FALSE)

colnames(parent1_fai) <- c("chr", "length", "offset", "linebases", "linewidth")
colnames(parent2_fai) <- c("chr", "length", "offset", "linebases", "linewidth")

# Filter for main chromosomes only (chr1-19, chrX, chrY for mouse)
parent1_chrs <- parent1_fai %>%
    filter(grepl("^chr([0-9]+|X|Y)$", chr)) %>%
    select(chr, length) %>%
    mutate(chr = paste0(parent1_name, "_", chr))

parent2_chrs <- parent2_fai %>%
    filter(grepl("^chr([0-9]+|X|Y)$", chr)) %>%
    select(chr, length) %>%
    mutate(chr = paste0(parent2_name, "_", chr))

# Create genome GRanges
genome_granges <- GRanges(
    seqnames = c(parent1_chrs$chr, parent2_chrs$chr),
    ranges = IRanges(
        start = 1,
        end = c(parent1_chrs$length, parent2_chrs$length)
    )
)

cat("Genome chromosomes loaded:\n")
cat("  ", parent1_name, ":", nrow(parent1_chrs), "chromosomes\n")
cat("  ", parent2_name, ":", nrow(parent2_chrs), "chromosomes\n\n")

## -------------------------------------------------------------------
## Read BED file with segment information
## -------------------------------------------------------------------
cat("Reading segment information from BED file...\n")

bed_data <- read.table(bed_file, header = FALSE, stringsAsFactors = FALSE,
                       col.names = c("read", "start", "end", "segment"))

cat("Total segments:", nrow(bed_data), "\n")

## -------------------------------------------------------------------
## Read curated profiles and filter for crossover reads
## -------------------------------------------------------------------
cat("Reading curated hybrid profiles...\n")

profiles <- read.table(curated_profiles_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
profiles_thr <- profiles %>% filter(threshold == threshold_val)

# Filter for crossover reads (SCO, NCO, GC-CO)
crossover_reads <- profiles_thr %>%
    filter(hybrid_type %in% c("SCO", "NCO", "GC-CO")) %>%
    pull(readname) %>%
    unique()

cat("Crossover reads at threshold", threshold_val, ":", length(crossover_reads), "\n")

# Filter bed data for crossover reads only
bed_crossovers <- bed_data %>%
    filter(read %in% crossover_reads)

cat("Segments from crossover reads:", nrow(bed_crossovers), "\n\n")

## -------------------------------------------------------------------
## Extract chromosome information from segments
## -------------------------------------------------------------------
cat("Extracting chromosome information from segments...\n")

# Extract chromosome number from segment names (e.g., "CH13" -> "chr13", "DH7" -> "chr7")
bed_crossovers <- bed_crossovers %>%
    mutate(
        parent = ifelse(grepl("^C", segment), parent1_name, parent2_name),
        chr_num = as.numeric(gsub("^[CD]H", "", segment)),
        chr = paste0(parent, "_chr", chr_num)
    )

# Filter out segments that don't match valid chromosomes
valid_chrs <- c(parent1_chrs$chr, parent2_chrs$chr)
bed_crossovers <- bed_crossovers %>%
    filter(chr %in% valid_chrs)

cat("Valid segments after filtering:", nrow(bed_crossovers), "\n")

# Separate by parent
parent1_segments <- bed_crossovers %>% filter(parent == parent1_name)
parent2_segments <- bed_crossovers %>% filter(parent == parent2_name)

cat("  ", parent1_name, "segments:", nrow(parent1_segments), "\n")
cat("  ", parent2_name, "segments:", nrow(parent2_segments), "\n\n")

## -------------------------------------------------------------------
## Create GRanges objects for segments
## -------------------------------------------------------------------
cat("Creating GRanges objects for visualization...\n")

parent1_granges <- GRanges(
    seqnames = parent1_segments$chr,
    ranges = IRanges(start = parent1_segments$start, end = parent1_segments$end),
    read = parent1_segments$read
)

parent2_granges <- GRanges(
    seqnames = parent2_segments$chr,
    ranges = IRanges(start = parent2_segments$start, end = parent2_segments$end),
    read = parent2_segments$read
)

## -------------------------------------------------------------------
## Create links between parent segments for the same reads
## -------------------------------------------------------------------
cat("Creating links between parent segments...\n")

# Find reads that have segments in both parents
reads_with_links <- intersect(parent1_segments$read, parent2_segments$read)
cat("Reads with links between parents:", length(reads_with_links), "\n\n")

# For visualization, create link pairs (first segment from each parent per read)
links_p1 <- parent1_segments %>%
    filter(read %in% reads_with_links) %>%
    group_by(read) %>%
    slice(1) %>%
    ungroup()

links_p2 <- parent2_segments %>%
    filter(read %in% reads_with_links) %>%
    group_by(read) %>%
    slice(1) %>%
    ungroup()

# Merge to ensure same order
link_data <- inner_join(links_p1, links_p2, by = "read", suffix = c("_p1", "_p2"))

if (nrow(link_data) == 0) {
    cat("Warning: No links found between parents. Cannot create karyoplot.\n")
    cat("This may be because reads don't have segments from both parents.\n")
    quit(save = "no", status = 0)
}

link_granges_p1 <- GRanges(
    seqnames = link_data$chr_p1,
    ranges = IRanges(start = link_data$start_p1, end = link_data$end_p1)
)

link_granges_p2 <- GRanges(
    seqnames = link_data$chr_p2,
    ranges = IRanges(start = link_data$start_p2, end = link_data$end_p2)
)

## -------------------------------------------------------------------
## Create karyotype plot
## -------------------------------------------------------------------
cat("Generating karyotype plot...\n")

# Order chromosomes
seqnames_char <- as.character(seqnames(genome_granges))
chr_numbers <- as.numeric(gsub(".*chr(\\d+)$", "\\1", seqnames_char))
chr_numbers[grepl("chrX", seqnames_char)] <- 100  # X chromosome
chr_numbers[grepl("chrY", seqnames_char)] <- 101  # Y chromosome

p1_indices <- grep(paste0("^", parent1_name), seqnames_char)
p2_indices <- grep(paste0("^", parent2_name), seqnames_char)

p1_sorted <- p1_indices[order(chr_numbers[p1_indices])]
p2_sorted <- p2_indices[order(chr_numbers[p2_indices])]

# Interleave parent1 and parent2
new_order <- as.vector(rbind(p1_sorted, p2_sorted))
genome_granges <- genome_granges[new_order]
ordered_chrs <- unique(as.character(seqnames(genome_granges)))

# Create output file
output_file <- file.path(output_dir, paste0(sample_name, "_karyoplot_crossovers.pdf"))
pdf(output_file, width = 14, height = 10)

# Create karyotype plot
kp <- plotKaryotype(
    genome = genome_granges,
    plot.type = 2,  # Side-by-side genome comparison
    main = paste("Crossover Karyotype Map -", sample_name, "- Threshold:", threshold_val),
    chromosomes = ordered_chrs,
    cytobands = NULL,
    cex = 0.8
)

# Plot segments
kpPlotRegions(kp, parent1_granges, r0 = 0, r1 = 0.45, col = "#ff8d9299", border = NA)
kpPlotRegions(kp, parent2_granges, r0 = 0, r1 = 0.45, col = "#8d9aff99", border = NA)

# Plot links
kpPlotLinks(
    kp,
    data = link_granges_p1,
    data2 = link_granges_p2,
    col = "#fac7ff66",
    border = "#00000044",
    lwd = 0.5,
    r0 = 0.2, r1 = 0.8
)

# Add legend
legend("bottomright",
       legend = c(paste(parent1_name, "segments"),
                  paste(parent2_name, "segments"),
                  "Crossover links"),
       fill = c("#ff8d92", "#8d9aff", "#fac7ff"),
       border = NA,
       bty = "n",
       cex = 0.9)

dev.off()

cat("\n==================================\n")
cat("KARYOTYPE PLOT COMPLETE\n")
cat("==================================\n")
cat("Output file:", output_file, "\n")
cat("Total links plotted:", nrow(link_data), "\n\n")
