#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  pkgs <- c("data.table","arrow","GenomicRanges")
  for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p, repos="https://cloud.r-project.org")
  library(data.table); library(arrow); library(GenomicRanges)
})

gwas_pq <- Sys.getenv("DEMO_GWAS","data/gwas/clean/bmi.parquet")
col_chr <- Sys.getenv("DEMO_GWAS_CHR","")
col_pos <- Sys.getenv("DEMO_GWAS_POS","")
col_p   <- Sys.getenv("DEMO_GWAS_P","p_num")
p_sig   <- as.numeric(Sys.getenv("DEMO_PSIG","5e-8"))

fig_dir <- "docs/figs"; dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)
doc_path <- "docs/Sprint1_EndGoal_Demo.md"

en_dir <- "data/regulatory/human/derived"
con_fp <- file.path(en_dir, "human_enhancers_conserved_with_mouse_RBH.bed.gz")
hs_fp  <- file.path(en_dir, "human_enhancers_human_specific_vs_mouse_RBH.bed.gz")

read_bed <- function(fp){
  if (!file.exists(fp)) stop("Missing enhancer BED: ", fp)
  dt <- fread(cmd=paste("gzip -dc", shQuote(fp)), sep="\t", header=FALSE,
              col.names=c("chr","start","end"))
  dt[, start := as.integer(start)]; dt[, end := as.integer(end)]
  dt
}
con_dt <- read_bed(con_fp); hs_dt <- read_bed(hs_fp)
gr_con <- GRanges(seqnames=con_dt$chr, ranges=IRanges(con_dt$start+1L, con_dt$end))
gr_hs  <- GRanges(seqnames=hs_dt$chr,  ranges=IRanges(hs_dt$start+1L,  hs_dt$end))

tbl <- as.data.frame(read_parquet(gwas_pq))
nms <- names(tbl)
if (col_chr=="" || col_pos=="") stop("Set DEMO_GWAS_CHR and DEMO_GWAS_POS (e.g. CHR_ID / CHR_POS).")
for (qq in c(col_chr,col_pos,col_p)) if (!(qq %in% nms))
  stop("Column not found: ", qq, "  (available: ", paste(nms, collapse=", "), ")")

chr <- as.character(tbl[[col_chr]])
chr <- ifelse(grepl("^chr", chr, ignore.case=TRUE), chr, paste0("chr", chr))
pos <- as.integer(tbl[[col_pos]])
pv  <- suppressWarnings(as.numeric(tbl[[col_p]]))
pv[!is.finite(pv)] <- NA_real_

gwas <- data.table(chr=chr, pos=pos, p=pv)
gwas <- gwas[!is.na(chr) & !is.na(pos) & !is.na(p)]
sig  <- gwas[p <= p_sig]
if (nrow(sig)==0) stop("No variants pass p <= ", p_sig, " in ", gwas_pq)

gr_sig <- GRanges(seqnames=sig$chr, ranges=IRanges(sig$pos, sig$pos))

h_con <- length(unique(queryHits(findOverlaps(gr_sig, gr_con))))
h_hs  <- length(unique(queryHits(findOverlaps(gr_sig, gr_hs))))
len_con <- sum(width(reduce(gr_con))); len_hs <- sum(width(reduce(gr_hs)))
r_con <- h_con/(len_con/1e6); r_hs <- h_hs/(len_hs/1e6)

png(file.path(fig_dir,"demo_end_goal_hits_per_Mb.png"), width=800, height=500)
vals <- c(r_con, r_hs)
bp <- barplot(vals, names.arg=c("Conserved","Human-specific"),
              ylab="GWAS hits per Mb (p \u2264 5e-8)",
              main="GWAS enrichment in enhancer classes")
text(bp, vals, labels=round(vals,3), pos=3, cex=0.9)
dev.off()

cat(sprintf(
"# End-goal demo: GWAS enrichment in enhancer classes\n
- GWAS file: `%s`
- Significant variants: %s
- Conserved span: %.1f Mb; hits: %s (%.3f / Mb)
- Human-specific span: %.1f Mb; hits: %s (%.3f / Mb)\n
![Demo](figs/%s)\n",
  gwas_pq, format(nrow(sig), big.mark=","),
  len_con/1e6, format(h_con, big.mark=","), r_con,
  len_hs/1e6,  format(h_hs,  big.mark=","), r_hs,
  "demo_end_goal_hits_per_Mb.png"
), file=doc_path)

message("Wrote: ", doc_path, " and docs/figs/demo_end_goal_hits_per_Mb.png")
