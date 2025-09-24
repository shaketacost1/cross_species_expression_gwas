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

# ---- Config ----
h_bed <- "data/regulatory/human/derived/human_enhancers_union.bed"
m_bed <- "data/regulatory/mouse/derived/mouse_enhancers_union.bed"

h2m_chain_path <- Sys.getenv("HG38_TO_MM10_CHAIN", "resources/chains/hg38ToMm10.over.chain.gz")
m2h_chain_path <- Sys.getenv("MM10_TO_HG38_CHAIN", "resources/chains/mm10ToHg38.over.chain.gz")
rbh_min_recip_frac <- as.numeric(Sys.getenv("RBH_MIN_RECIP_FRAC", "0.5"))

out_dir <- "data/regulatory/human/derived"
fig_dir <- "docs/figs"
doc_path <- "docs/Sprint1_Human_EnhancerConservation_RBH_vs_Mouse.md"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)

stopifnot(file.exists(h_bed), file.exists(m_bed), file.exists(h2m_chain_path), file.exists(m2h_chain_path))

# ---- Helpers ----
bed_to_gr <- function(dt) {
  GRanges(seqnames = dt$chr,
          ranges   = IRanges(start = as.integer(dt$start) + 1L, end = as.integer(dt$end)))
}

# ---- Load union beds ----
h_dt <- fread(h_bed, header=FALSE, sep="\t", showProgress=FALSE)
# expect: chr start end name score strand (6 cols)
stopifnot(ncol(h_dt) >= 3)
setnames(h_dt, 1:3, c("chr","start","end"))

# we don't actually need mouse GRanges for RBH classification, but read to ensure file exists
m_dt <- fread(m_bed, header=FALSE, sep="\t", showProgress=FALSE) # sanity only
rm(m_dt)

h_gr <- bed_to_gr(h_dt)

# ---- Import chains (gz OK) ----
h2m <- import(h2m_chain_path, format="chain")
m2h <- import(m2h_chain_path, format="chain")

# ---- Liftover: human -> mouse ----
h2m_list <- liftOver(h_gr, h2m)
one_map <- elementNROWS(h2m_list) == 1L
idx_h_one <- which(one_map)
m_gr <- unlist(h2m_list[one_map], use.names=FALSE)
m_gr$hid <- idx_h_one

# ---- Back: mouse -> human ----
m2h_list <- liftOver(m_gr, m2h)
one_back <- elementNROWS(m2h_list) == 1L
back_gr  <- unlist(m2h_list[one_back], use.names=FALSE)
back_gr$hid <- m_gr$hid[one_back]

# ---- RBH filter: sufficient pairwise overlap after round-trip ----
orig <- h_gr[back_gr$hid]
ov   <- pintersect(back_gr, orig, ignore.strand=TRUE)
frac <- width(ov) / width(orig)
keep <- (!is.na(frac)) & (frac >= rbh_min_recip_frac)

hid_conserved <- unique(back_gr$hid[keep])

# ---- Label all human enhancers ----
label <- rep("human_specific_vs_mouse_RBH", length(h_gr))
label[hid_conserved] <- "conserved_with_mouse_RBH"

res <- data.table(chr = h_dt$chr, start = h_dt$start, end = h_dt$end, class = label)

# ---- Write TSV.GZ (compress via system gzip for portability) ----
out_tsv_gz <- file.path(out_dir, "human_vs_mouse_enhancer_conservation_RBH.tsv.gz")
tmp_tsv    <- sub("\\.gz$", "", out_tsv_gz)
fwrite(res, tmp_tsv, sep="\t")
system(paste("gzip -f", shQuote(tmp_tsv)))

# ---- Plot simple counts ----
cnt <- res[, .N, by=class][order(-N)]
png(file.path(fig_dir, "human_vs_mouse_enhancer_conservation_RBH_counts.png"),
    width=900, height=600, res=120)
par(mar=c(6,6,3,1))
bp <- barplot(cnt$N, names.arg=gsub("_", "\n", cnt$class), las=2, ylab="count (log10)",
              main="Human enhancer conservation vs mouse (RBH)", log="y")
text(bp, cnt$N, labels=format(cnt$N, big.mark=","), pos=3, cex=0.8)
dev.off()

# ---- Tiny markdown summary ----
cat(sprintf(
"# Sprint 1 — Human enhancer conservation vs mouse (RBH)\n\n- **Total human enhancers:** %s\n- **Conserved (RBH):** %s\n- **Human-specific vs mouse (RBH):** %s\n\n![Counts](figs/%s)\n",
  format(nrow(res), big.mark=","),
  format(cnt[class=='conserved_with_mouse_RBH', N], big.mark=","),
  format(cnt[class=='human_specific_vs_mouse_RBH', N], big.mark=","),
  "human_vs_mouse_enhancer_conservation_RBH_counts.png"
), file=doc_path)

message("Wrote: ", out_tsv_gz, " ; ",
        file.path(fig_dir,"human_vs_mouse_enhancer_conservation_RBH_counts.png"), " ; ",
        doc_path)
