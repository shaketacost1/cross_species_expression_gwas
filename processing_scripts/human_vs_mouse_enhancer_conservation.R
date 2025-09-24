#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  pkgs <- c("data.table","GenomicRanges","rtracklayer")
  for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) {
    if (p %in% c("GenomicRanges","rtracklayer")) {
      if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
      BiocManager::install(p, ask=FALSE, update=FALSE)
    } else install.packages(p)
  }
  library(data.table); library(GenomicRanges); library(rtracklayer)
})

# paths
h_bed <- "data/regulatory/human/derived/human_enhancers_union.bed"
m_bed <- "data/regulatory/mouse/derived/mouse_enhancers_union.bed"
chain_fp <- "resources/chains/mm10ToHg38.over.chain.gz"
stopifnot(file.exists(h_bed), file.exists(m_bed), file.exists(chain_fp))

out_dir  <- "data/regulatory/human/derived"
fig_dir  <- "docs/figs"
doc_path <- "docs/Sprint1_Human_EnhancerConservation_vs_Mouse.md"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)

# read minimal BEDs
read_bed <- function(fp){
  dt <- fread(fp, header=FALSE, sep="\t", data.table=TRUE, showProgress=FALSE)
  stopifnot(ncol(dt) >= 3)
  setnames(dt, 1:3, c("chr","start","end"))
  dt[, start := as.integer(start)][, end := as.integer(end)]
  dt
}
h_dt <- read_bed(h_bed)
m_dt <- read_bed(m_bed)

# GRanges (BED is 0-based half-open -> +1 for 1-based GRanges)
h_gr <- GRanges(seqnames=h_dt$chr, ranges=IRanges(start=h_dt$start+1L, end=h_dt$end))
m_gr <- GRanges(seqnames=m_dt$chr, ranges=IRanges(start=m_dt$start+1L, end=m_dt$end))

# lift mouse->human (mm10->hg38)
ch <- import(Sys.getenv("MM10_TO_HG38_CHAIN","resources/chains/mm10ToHg38.over.chain"), format="chain")
m_to_h <- liftOver(m_gr, ch)
m_to_h <- unlist(m_to_h)
m_to_h <- reduce(m_to_h, ignore.strand=TRUE)  # dedup/merge

# overlaps
hits <- findOverlaps(h_gr, m_to_h, ignore.strand=TRUE)
is_cons <- logical(length(h_gr)); is_cons[queryHits(hits)] <- TRUE
class <- ifelse(is_cons, "conserved_with_mouse", "human_specific_vs_mouse")

# write labeled human enhancer table
lab_dt <- data.table(chr=as.character(seqnames(h_gr)),
                     start=start(h_gr)-1L,
                     end=end(h_gr),
                     class=class)
fwrite(lab_dt, file.path(out_dir,"human_vs_mouse_enhancer_conservation.tsv.gz"), sep="\t")

# counts and fig
cnt <- lab_dt[, .N, by=class][order(-N)]
png(file.path(fig_dir,"human_vs_mouse_enhancer_conservation_counts.png"),
    width=1100, height=700, res=130)
barplot(cnt$N, names.arg=cnt$class, las=2, ylab="Number of human enhancers",
        main="Human enhancer conservation vs mouse (liftOver mm10→hg38)")
dev.off()

# brief markdown
cat(sprintf(
"# Sprint 1 — Human enhancer conservation vs mouse\n\n- **Total human enhancers:** %s\n- **Conserved with mouse:** %s\n- **Human-specific vs mouse:** %s\n\n![Counts](figs/%s)\n",
  format(nrow(lab_dt), big.mark=","),
  format(cnt[class=="conserved_with_mouse", N], big.mark=","),
  format(cnt[class=="human_specific_vs_mouse", N], big.mark=","),
  "human_vs_mouse_enhancer_conservation_counts.png"
), file=doc_path)

message("Wrote: ",
        file.path(out_dir,"human_vs_mouse_enhancer_conservation.tsv.gz"), " ; ",
        file.path(fig_dir,"human_vs_mouse_enhancer_conservation_counts.png"), " ; ",
        doc_path)
