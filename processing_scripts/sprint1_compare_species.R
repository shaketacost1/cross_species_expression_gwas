#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  pkgs <- c("arrow","dplyr","GenomicRanges","ggplot2","tibble","data.table")
  for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) {
    if (p=="GenomicRanges") { if(!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager"); BiocManager::install("GenomicRanges", ask=FALSE, update=FALSE) }
    else install.packages(p)
  }
  library(arrow); library(dplyr); library(GenomicRanges); library(ggplot2); library(tibble); library(data.table)
})
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) stop("Usage: sprint1_compare_species.R <speciesA> <speciesB>")
comp_stats <- function(sp){
  base <- file.path("data","regulatory",sp,"derived")
  enh  <- read_parquet(file.path(base, paste0(sp,"_enhancers_union.parquet")))
  ctcf <- read_parquet(file.path(base, paste0(sp,"_ctcf_merged.parquet")))
  grE  <- GRanges(enh$chr, IRanges(enh$start+1L, enh$end))
  grC  <- GRanges(ctcf$chr, IRanges(ctcf$start+1L, ctcf$end))
  ov   <- length(unique(queryHits(findOverlaps(grE, grC, ignore.strand=TRUE))))
  tibble(species=sp, enhancers=nrow(enh), ctcf=nrow(ctcf), overlap=ov, pct=100*ov/nrow(enh))
}
s1 <- comp_stats(args[1]); s2 <- comp_stats(args[2])
stats <- bind_rows(s1,s2)

cap <- function(x) paste0(toupper(substr(x,1,1)), substr(x,2,nchar(x)))
spA <- cap(args[1]); spB <- cap(args[2])

dir.create("docs/figs", showWarnings=FALSE, recursive=TRUE)
p <- ggplot(stats, aes(x=species, y=pct)) + geom_col() + ylim(0,100) +
     labs(title="% enhancers overlapping merged CTCF", x=NULL, y="% overlap") + theme_minimal()
fig <- sprintf("docs/figs/%s_%s_ctcf_overlap.png", tolower(spA), tolower(spB))
ggsave(fig, p, width=5, height=4, dpi=150)

md <- sprintf("docs/Sprint1_Compare_%s_%s.md", spA, spB)
cat(sprintf(
"## Sprint 1 — %s vs %s\n\n**Enhancers (union):** %s = %s, %s = %s\n\n**CTCF (merged):** %s = %s, %s = %s\n\n**Enhancers overlapping CTCF:** %s = %s (%.1f%%), %s = %s (%.1f%%)\n\n![](figs/%s)\n",
spA, spB,
spA, format(stats$enhancers[stats$species==args[1]], big.mark=","),
spB, format(stats$enhancers[stats$species==args[2]], big.mark=","),
spA, format(stats$ctcf[stats$species==args[1]], big.mark=","),
spB, format(stats$ctcf[stats$species==args[2]], big.mark=","),
spA, format(stats$overlap[stats$species==args[1]], big.mark=","), stats$pct[stats$species==args[1]],
spB, format(stats$overlap[stats$species==args[2]], big.mark=","), stats$pct[stats$species==args[2]],
basename(fig)), file=md)
message("Wrote: ", md, " and ", fig)
