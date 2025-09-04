#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  pkgs <- c("arrow","data.table","dplyr","ggplot2","GenomicRanges")
  for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) {
    if (p=="GenomicRanges") { if(!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager"); BiocManager::install("GenomicRanges", ask=FALSE, update=FALSE) }
    else install.packages(p)
  }
  library(arrow); library(data.table); library(dplyr); library(ggplot2); library(GenomicRanges)
})

args <- commandArgs(trailingOnly=TRUE)
species <- if (length(args)>=1) args[1] else "human"
cap <- paste0(toupper(substr(species,1,1)), substr(species,2,nchar(species)))

dir.create("docs/figs", recursive=TRUE, showWarnings=FALSE)

d_derived <- file.path("data","regulatory",species,"derived")
enh_parq  <- file.path(d_derived, paste0(species, "_enhancers_union.parquet"))
ctcf_parq <- file.path(d_derived, paste0(species, "_ctcf_merged.parquet"))
enh_bed   <- file.path(d_derived, paste0(species, "_enhancers_union.bed"))
ctcf_bed  <- file.path(d_derived, paste0(species, "_ctcf_merged.bed"))

# ---- Load enhancers ----
if (file.exists(enh_parq)) {
  enh <- read_parquet(enh_parq)
} else if (file.exists(enh_bed)) {
  enh <- fread(enh_bed, col.names=c("chr","start","end","name","score","strand"))
} else stop("No enhancer file found for ", species)

# keep only needed cols
need_cols <- intersect(c("chr","start","end","type","source","name"), names(enh))
enh <- as.data.table(enh)[, ..need_cols]
if (!"type" %in% names(enh)) enh[, type := NA_character_]

# ---- Load CTCF (merged) ----
if (file.exists(ctcf_parq)) {
  ctcf <- read_parquet(ctcf_parq)
} else if (file.exists(ctcf_bed)) {
  ctcf <- fread(ctcf_bed, col.names=c("chr","start","end"))
} else stop("No CTCF merged file found for ", species)

# ---- Basic counts ----
n_enh <- nrow(enh)
n_ctcf <- nrow(ctcf)

# ---- Overlaps ----
grE <- GRanges(enh$chr, IRanges(enh$start+1L, enh$end))
grC <- GRanges(ctcf$chr, IRanges(ctcf$start+1L, ctcf$end))
ov_idx <- findOverlaps(grE, grC, ignore.strand=TRUE)
n_ov <- length(unique(queryHits(ov_idx)))
pct_ov <- 100 * n_ov / n_enh

# ---- Figures ----
# A) Counts bar
counts_df <- tibble::tibble(layer=c("Enhancers (union)","CTCF (merged)"),
                            count=c(n_enh, n_ctcf))
p_counts <- ggplot(counts_df, aes(layer, count)) +
  geom_col() + theme_minimal() +
  labs(title=paste0(cap, " – counts"),
       x=NULL, y="Number of intervals")
out_counts <- file.path("docs","figs", paste0(species, "_counts.png"))
ggsave(out_counts, p_counts, width=6, height=4, dpi=150)

# B) Enhancer type composition (if available)
type_plot_path <- NA_character_
if (any(!is.na(enh$type))) {
  comp <- enh %>% filter(!is.na(type)) %>%
    count(type, name="n") %>% mutate(pct = 100*n/sum(n))
  p_type <- ggplot(comp, aes(reorder(type,-pct), pct)) +
    geom_col() + theme_minimal() +
    labs(title=paste0(cap, " – enhancer type composition"),
         x="type", y="% of enhancers")
  type_plot_path <- file.path("docs","figs", paste0(species, "_enhancer_type_composition.png"))
  ggsave(type_plot_path, p_type, width=6, height=4, dpi=150)
}

# C) Overlap by type (if type available)
ov_type_path <- NA_character_
if (any(!is.na(enh$type))) {
  has_ov <- unique(queryHits(ov_idx))
  enh[, overlapped := FALSE]
  if (length(has_ov)) enh[has_ov, overlapped := TRUE]
  bytype <- enh[!is.na(type), .(n=.N, n_ov=sum(overlapped)), by=type][
              , pct := 100*n_ov/pmax(n,1)]
  p_ovt <- ggplot(bytype, aes(reorder(type,-pct), pct)) +
    geom_col() + theme_minimal() +
    labs(title=paste0(cap, " – enhancers overlapping CTCF by type"),
         x="type", y="% overlapping CTCF")
  ov_type_path <- file.path("docs","figs", paste0(species, "_enhancers_overlap_by_type.png"))
  ggsave(ov_type_path, p_ovt, width=6, height=4, dpi=150)
}

# ---- Markdown report ----
md_path <- file.path("docs", paste0("Sprint1_", cap, "_Regulatory_Summary.md"))
cat(sprintf(
"## Sprint 1 — %s regulatory snapshot

**Enhancers (union):** %s  
**CTCF (merged):** %s  
**Enhancers overlapping CTCF:** %s (%.1f%%)

![Counts](figs/%s)
",
cap,
format(n_enh, big.mark=","),
format(n_ctcf, big.mark=","),
format(n_ov, big.mark=","),
pct_ov,
basename(out_counts)
), file=md_path)

if (!is.na(type_plot_path)) {
  cat(sprintf("\n\n![Enhancer type composition](figs/%s)\n", basename(type_plot_path)), file=md_path, append=TRUE)
}
if (!is.na(ov_type_path)) {
  cat(sprintf("\n\n![Enhancer–CTCF overlap by type](figs/%s)\n", basename(ov_type_path)), file=md_path, append=TRUE)
}

message("Wrote: ", md_path)
