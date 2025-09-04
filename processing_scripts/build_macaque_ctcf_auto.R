#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); library(GenomicRanges) })

# ---- Config ----
input_root <- normalizePath(".", mustWork = TRUE)          # search in project root
out_dir    <- "data/regulatory/macaque/derived"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_merged <- file.path(out_dir, "macaque_ctcf_merged.bed") # what the viz expects
out_union  <- file.path(out_dir, "macaque_ctcf_union.bed")  # optional

# macaque filename prefixes to keep
macaque_prefix_whitelist <- c("rheMac","mmul","Mmul","Macaca_mulatta","macMul",
                              "macFas","Macaca_fascicularis","mfasc","macaque")

# ---- Helpers ----
stem <- function(x) sub("(?:\\.bed(?:\\.gz)?|\\.narrowPeak(?:\\.gz)?)+$", "", x, perl=TRUE, ignore.case=TRUE)
is_header_like <- function(dt){
  (tolower(dt[[1]]) %in% c("chr","chrom","chromosome")) |
    (tolower(dt[[2]]) %in% "start") | (tolower(dt[[3]]) %in% "end")
}
read_bed_any <- function(fp){
  # data.table::fread can read .gz directly when given a file path
  dt <- tryCatch(
    fread(fp, sep="\t", header=FALSE, data.table=TRUE, fill=TRUE, select=1:3, showProgress=FALSE),
    error=function(e) tryCatch(
      fread(fp, sep=" ", header=FALSE, data.table=TRUE, fill=TRUE, select=1:3, showProgress=FALSE),
      error=function(e2) NULL
    )
  )
  if (is.null(dt) || ncol(dt) < 3L) return(NULL)
  setnames(dt, 1:3, c("chr","start","end"))
  dt <- dt[!is_header_like(dt)]
  suppressWarnings({ dt[, start := as.integer(start)]; dt[, end := as.integer(end)] })
  dt <- dt[!is.na(start) & !is.na(end) & end > start & start >= 0]
  if (!nrow(dt)) return(NULL)
  dt[, chr := as.character(chr)]
  dt
}
to_gr <- function(dt){
  # BED is 0-based, half-open -> GRanges is 1-based closed
  GRanges(seqnames = dt$chr, ranges = IRanges(dt$start + 1L, dt$end))
}
write_bed3 <- function(gr, path){
  dt <- data.table(chr = as.character(seqnames(gr)),
                   start = start(gr) - 1L,  # back to 0-based
                   end   = end(gr))
  fwrite(dt, path, sep="\t", col.names=FALSE)
}
matches_any_prefix <- function(prefix, whitelist) {
  any(startsWith(prefix, whitelist)) || any(startsWith(tolower(prefix), tolower(whitelist)))
}

# ---- Find macaque CTCF files (accept BED/bed.gz/narrowPeak/narrowPeak.gz) ----
all_files <- list.files(input_root, recursive=TRUE, full.names=TRUE)
cand <- all_files[
  grepl("CTCF", basename(all_files), ignore.case=TRUE) &
    grepl("\\.(bed|bed\\.gz|narrowPeak|narrowPeak\\.gz)$", basename(all_files), ignore.case=TRUE)
]

if (!length(cand)) {
  message("No files matched '*CTCF*.(bed|narrowPeak)[.gz]'.")
  # Helpful diagnostic: show any paths containing 'CTCF'
  diag <- all_files[grepl("CTCF", basename(all_files), ignore.case=TRUE)]
  if (length(diag)) {
    message("CTCF-like names seen (but not BED/narrowPeak):")
    print(head(basename(diag), 20))
  }
  quit(save="no", status=1)
}

pref <- sub("-.*", "", basename(cand))
is_macaque <- vapply(pref, matches_any_prefix, logical(1), whitelist = macaque_prefix_whitelist)
ctcf_files <- cand[is_macaque]
if (!length(ctcf_files)) {
  stop("Found CTCF files, but none look macaque by prefix. Prefixes seen: ",
       paste(unique(pref), collapse=", "))
}

# Deduplicate by stem
stems <- stem(basename(ctcf_files))
keep  <- !duplicated(stems)
if (any(!keep)) message("Dropping duplicate stems: ", paste(basename(ctcf_files[!keep]), collapse=", "))
ctcf_files <- ctcf_files[keep]

# ---- Read, union, merge ----
parts <- list()
for (p in ctcf_files) {
  dt <- read_bed_any(p)
  if (is.null(dt)) { message("Skip (no valid rows): ", basename(p)); next }
  if (!all(grepl("^chr", dt$chr))) dt[, chr := paste0("chr", chr)]
  parts[[length(parts)+1]] <- dt
}
if (!length(parts)) stop("No valid rows in any CTCF file.")

union_dt <- rbindlist(parts, use.names=TRUE, fill=TRUE)

# write UNION (optional)
fwrite(union_dt[, .(chr,start,end)], out_union, sep="\t", col.names=FALSE)

# MERGE (reduce overlaps)
gr <- to_gr(union_dt)
grm <- reduce(gr, ignore.strand=TRUE)
write_bed3(grm, out_merged)

message("✅ Wrote: ", out_union)
message("✅ Wrote: ", out_merged, " (", length(grm), " intervals)")
