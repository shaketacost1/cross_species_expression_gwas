#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# ----------------------------
# Simple arg parser
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)
argval <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0L) return(default)
  if (i == length(args)) return(TRUE) # boolean flag
  nxt <- args[i + 1L]
  if (grepl("^--", nxt)) return(TRUE) # boolean flag
  nxt
}

# Defaults
species        <- argval("--species", "macaque")
mark_pattern   <- argval("--mark", "H3K27[Aa]c")     # regex
normalize_chr_names <- isTRUE(as.logical(argval("--chr-norm", "TRUE")))
promoters_bed  <- argval("--promoters", NULL)
blacklist_bed  <- argval("--blacklist", NULL)
max_len        <- argval("--max-len", NULL)

if (!is.null(max_len)) {
  if (!grepl("^[0-9]+$", max_len)) stop("--max-len must be an integer number of bp")
  max_len <- as.integer(max_len)
}

# ----------------------------
# Configuration (paths)
# ----------------------------
input_root <- normalizePath(".", mustWork = TRUE)
out_bed_dir    <- file.path("data", "regulatory", species, "enhancers_raw")
derived_dir    <- file.path("data", "regulatory", species, "derived")
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)

union_bed_path  <- file.path(derived_dir, sprintf("%s_enhancers_union.bed", species))
merged_bed_path <- file.path(derived_dir, sprintf("%s_enhancers_merged.bed", species))
qc_txt   <- file.path(derived_dir, sprintf("%s_enhancers_qc.txt", species))
hist_png <- file.path(derived_dir, sprintf("%s_enhancers_len_hist.png", species))
union_script   <- "processing_scripts/enhancers_build_union.R"

# ----------------------------
# Helpers
# ----------------------------
must_normpath <- function(p) normalizePath(p, mustWork = TRUE)
stem <- function(x) sub("(?:\\.bed(?:\\.gz)?)+$", "", x, perl = TRUE, ignore.case = TRUE)
is_header_like <- function(dt) {
  (tolower(dt[[1]]) %in% c("chr", "chrom", "chromosome")) |
    (tolower(dt[[2]]) %in% c("start")) |
    (tolower(dt[[3]]) %in% c("end"))
}
read_min_bed <- function(lines) {
  if (!length(lines)) return(NULL)
  lines <- lines[!grepl("^(#|track|browser)", lines)]
  if (!length(lines)) return(NULL)
  
  dt <- tryCatch(
    fread(text = lines, sep = "\t", header = FALSE, data.table = TRUE,
          fill = TRUE, showProgress = FALSE, select = 1:3),
    error = function(e) tryCatch(
      fread(text = lines, sep = " ", header = FALSE, data.table = TRUE,
            fill = TRUE, showProgress = FALSE, select = 1:3),
      error = function(e2) NULL
    )
  )
  if (is.null(dt) || ncol(dt) < 3L) return(NULL)
  setnames(dt, 1:3, c("chr", "start", "end"))
  dt <- dt[!is_header_like(dt)]
  suppressWarnings({
    dt[, start := as.integer(start)]
    dt[,   end := as.integer(end)]
  })
  dt <- dt[!is.na(start) & !is.na(end) & end > start & start >= 0]
  if (nrow(dt) == 0L) return(NULL)
  dt[, chr := as.character(chr)]
  dt
}
write_min_bed <- function(dt, out_path) {
  fwrite(dt[, .(chr, start, end)], out_path, sep = "\t", col.names = FALSE)
}
normalize_chr_prefix <- function(bed_path) {
  if (!file.exists(bed_path)) return(invisible(FALSE))
  dt <- fread(bed_path, header = FALSE)
  if (ncol(dt) < 3L) return(invisible(FALSE))
  setnames(dt, 1:3, c("chr","start","end"))
  dt[, chr := as.character(chr)]
  dt[!grepl("^chr", chr), chr := paste0("chr", chr)]
  dt[chr %in% c("chrM","chrMT"), chr := "chrMT"]
  fwrite(dt, bed_path, sep = "\t", col.names = FALSE)
  invisible(TRUE)
}
bedtools_available <- (system("bedtools --version >/dev/null 2>&1") == 0)
matches_any_prefix <- function(prefix, whitelist) {
  any(startsWith(prefix, whitelist)) ||
    any(startsWith(tolower(prefix), tolower(whitelist)))
}
macaque_prefix_whitelist <- c(
  "rheMac", "mmul", "Mmul", "Macaca_mulatta", "macMul",
  "macFas", "Macaca_fascicularis", "mfasc",
  "macaque"
)

# ----------------------------
# 1) Discover input files
# ----------------------------
ae_dir <- must_normpath(input_root)
all_files <- list.files(ae_dir, recursive = TRUE, full.names = TRUE)
peak_files <- all_files[
  grepl(mark_pattern, basename(all_files), ignore.case = TRUE) &
    grepl("replicated-peaks_macs", basename(all_files), ignore.case = TRUE)
]

if (!length(peak_files)) {
  message("No H3K27ac replicated-peak files found under: ", ae_dir)
  quit(save = "no", status = 0)
}

all_prefixes <- sub("-.*", "", basename(peak_files))
is_macaque   <- vapply(all_prefixes, matches_any_prefix, logical(1),
                       whitelist = macaque_prefix_whitelist)
mac_files    <- peak_files[is_macaque]
message("Found ", length(peak_files), " H3K27ac peak file(s); ",
        sum(is_macaque), " appear macaque by prefix.")

if (!length(mac_files)) {
  message("No macaque-prefixed H3K27ac files detected. Prefixes present: ",
          paste(unique(all_prefixes), collapse = ", "))
  quit(save = "no", status = 0)
}

# ----------------------------
# 2) De-duplicate by stem
# ----------------------------
mac_stems <- stem(basename(mac_files))
keep_idx  <- !duplicated(mac_stems)
if (any(!keep_idx)) {
  message("De-duplicating by stem; dropping duplicates: ",
          paste(basename(mac_files[!keep_idx]), collapse = ", "))
}
mac_files <- mac_files[keep_idx]

# ----------------------------
# 3) Parse -> clean -> write minimal BEDs
# ----------------------------
dir.create(out_bed_dir, recursive = TRUE, showWarnings = FALSE)
wrote <- 0L
seen_bases <- character(0)
for (p in mac_files) {
  base <- gsub("[^A-Za-z0-9_.-]", "_", stem(basename(p)))
  if (base %in% seen_bases) {
    message("⏭️  Skipping duplicate stem: ", basename(p))
    next
  }
  
  txt <- readLines(p, warn = FALSE)
  dt  <- read_min_bed(txt)
  if (is.null(dt)) {
    message("⚠️  Skip (no valid rows): ", basename(p))
    next
  }
  
  out <- file.path(out_bed_dir, paste0(base, ".bed"))
  write_min_bed(dt, out)
  
  seen_bases <- c(seen_bases, base)
  wrote <- wrote + 1L
  message("✅ Wrote: ", out, " (", nrow(dt), " rows)")
}
message("✅ Wrote ", wrote, " BED file(s) to: ", out_bed_dir)

if (wrote == 0L) {
  message("No macaque BEDs written—nothing to build. Exiting.")
  quit(save = "no", status = 0)
}

# ----------------------------
# 4) Cleanup legacy duplicates
# ----------------------------
files <- list.files(out_bed_dir, full.names = TRUE)
bases <- basename(files)
stems <- stem(bases)
to_keep_idx <- !duplicated(stems)
to_delete   <- files[!to_keep_idx]
if (length(to_delete)) {
  message("Cleaning duplicates in enhancers_raw:\n  - ",
          paste(basename(to_delete), collapse = "\n  - "))
  invisible(file.remove(to_delete))
}

# ----------------------------
# 5) Build union/merged
# ----------------------------
cmd <- sprintf("Rscript %s %s", union_script, species)
message("Running: ", cmd)
status <- system(cmd)
if (status != 0) stop("Union/merge script failed with status ", status)
message("✅ Done: enhancer UNION and MERGED layers written for species: ", species)
message("   UNION:  ", union_bed_path, " (+ tsv.gz/parquet where applicable)")
message("   MERGED: ", merged_bed_path, " (+ parquet where applicable)")

# ----------------------------
# 6) Normalize chr names
# ----------------------------
if (normalize_chr_names) {
  message("Normalizing chromosome names to UCSC-style for UNION and MERGED...")
  normalize_chr_prefix(union_bed_path)
  normalize_chr_prefix(merged_bed_path)
}

# ----------------------------
# 7) QC: counts, lengths, histogram, report
# ----------------------------
u_dt <- fread(union_bed_path, header = FALSE, select = 1:3)
setnames(u_dt, c("chr","start","end"))
m_dt <- fread(merged_bed_path, header = FALSE, select = 1:3)
setnames(m_dt, c("chr","start","end"))
stopifnot(all(u_dt$end > u_dt$start, na.rm = TRUE), all(u_dt$start >= 0, na.rm = TRUE))
stopifnot(all(m_dt$end > m_dt$start, na.rm = TRUE), all(m_dt$start >= 0, na.rm = TRUE))
u_n <- nrow(u_dt); m_n <- nrow(m_dt)
collapsed <- u_n - m_n; collapsed_pct <- round(100 * collapsed / max(1, u_n), 2)
m_dt[, len := end - start]

png(hist_png, width = 1200, height = 800, res = 150)
hist(m_dt$len, breaks = 100, main = sprintf("Merged peak lengths: %s H3K27ac", species),
     xlab = "Length (bp)")
dev.off()

capped_path <- NULL
if (!is.null(max_len)) {
  capped_path <- file.path(derived_dir, sprintf("%s_enhancers_merged.max%dbp.bed", species, max_len))
  fwrite(m_dt[len <= max_len, .(chr, start, end)], capped_path, sep = "\t", col.names = FALSE)
  message("Wrote length-capped merged BED: ", capped_path)
}

sink(qc_txt)
cat("QC report for", species, "H3K27ac\n")
cat("Union intervals:", u_n, "\n")
cat("Merged intervals:", m_n, "\n")
cat("Collapsed overlaps:", collapsed, sprintf("(%s%%)", collapsed_pct), "\n\n")
print(summary(m_dt$len))
if (!is.null(max_len)) {
  cat("\nMax length cap:", max_len, "bp\n")
  cat("Capped file:", capped_path, "\n")
}
sink()
message("QC written: ", qc_txt)
message("Histogram written: ", hist_png)

# ----------------------------
# 8) Optional promoter/distal split
# ----------------------------
if (!is.null(promoters_bed) && file.exists(promoters_bed) && bedtools_available) {
  pro_out <- file.path(derived_dir, sprintf("%s_enhancers_promoter.bed", species))
  dis_out <- file.path(derived_dir, sprintf("%s_enhancers_distal.bed", species))
  system(sprintf("bedtools intersect -u -a %s -b %s > %s", merged_bed_path, promoters_bed, pro_out))
  system(sprintf("bedtools intersect -v -a %s -b %s > %s", merged_bed_path, promoters_bed, dis_out))
  message("Promoter/distal split written:\n  - ", pro_out, "\n  - ", dis_out)
} else if (!is.null(promoters_bed)) {
  if (!file.exists(promoters_bed)) message("Promoters BED not found: ", promoters_bed)
  if (!bedtools_available) message("bedtools not found in PATH; skipping promoter/distal split.")
}

# ----------------------------
# 9) Optional blacklist filtering
# ----------------------------
if (!is.null(blacklist_bed) && file.exists(blacklist_bed) && bedtools_available) {
  bl_out <- file.path(derived_dir, sprintf("%s_enhancers_merged.blacklistFiltered.bed", species))
  system(sprintf("bedtools intersect -v -a %s -b %s > %s", merged_bed_path, blacklist_bed, bl_out))
  message("Blacklist-filtered merged BED:\n  - ", bl_out)
} else if (!is.null(blacklist_bed)) {
  if (!file.exists(blacklist_bed)) message("Blacklist BED not found: ", blacklist_bed)
  if (!bedtools_available) message("bedtools not found in PATH; skipping blacklist filtering.")
}

message("All checks passed. Done.")
