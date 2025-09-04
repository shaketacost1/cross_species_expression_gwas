#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})

# ---- Paths (EDIT the CTCF path to your file) ----
enh_bed <- "data/regulatory/macaque/derived/macaque_enhancers_merged.bed"
ctcf_bed <- "data/regulatory/macaque/derived/macaque_ctcf_merged.bed"  # <-- set this to your macaque CTCF BED
out_dir  <- "docs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Helpers ----
read_bed3 <- function(path) {
  dt <- fread(path, header = FALSE)
  stopifnot(ncol(dt) >= 3)
  setnames(dt, 1:3, c("chr","start","end"))
  dt[, chr := as.character(chr)]
  dt
}
to_gr <- function(dt) GRanges(seqnames = dt$chr, ranges = IRanges(dt$start + 1L, dt$end))  # BED: 0-based, half-open

# ---- Load ----
e_dt <- read_bed3(enh_bed)
c_dt <- read_bed3(ctcf_bed)

# normalize chr prefix if needed
if (!all(grepl("^chr", e_dt$chr))) e_dt[, chr := paste0("chr", chr)]
if (!all(grepl("^chr", c_dt$chr))) c_dt[, chr := paste0("chr", chr)]

e_gr <- to_gr(e_dt)
c_gr <- to_gr(c_dt)

# ---- Length histograms ----
png(file.path(out_dir, "macaque_enhancer_length_hist.png"), width = 1200, height = 800, res = 150)
hist(width(e_gr), breaks = 100, main = "Macaque H3K27ac enhancers (merged) — lengths", xlab = "bp")
dev.off()

png(file.path(out_dir, "macaque_ctcf_length_hist.png"), width = 1200, height = 800, res = 150)
hist(width(c_gr), breaks = 100, main = "Macaque CTCF peaks — lengths", xlab = "bp")
dev.off()

# ---- Overlap counts (enhancers vs CTCF) ----
hits <- findOverlaps(e_gr, c_gr, ignore.strand = TRUE)
has_ctcf <- logical(length(e_gr)); has_ctcf[queryHits(hits)] <- TRUE
tab <- table(`Enhancer has CTCF` = ifelse(has_ctcf, "Yes", "No"))

png(file.path(out_dir, "macaque_enhancer_ctcf_overlap_bar.png"), width = 1000, height = 800, res = 150)
barplot(tab, main = "Enhancers overlapping CTCF", ylab = "Count")
dev.off()

# ---- Tiny text summary ----
pct <- round(100 * (sum(has_ctcf) / length(has_ctcf)), 2)
cat("Enhancers:", length(e_gr), "\nCTCF peaks:", length(c_gr),
    "\nEnhancers with CTCF overlap:", sum(has_ctcf), sprintf("(%s%%)\n", pct),
    file = file.path(out_dir, "macaque_enhancer_ctcf_overlap_summary.txt"))
