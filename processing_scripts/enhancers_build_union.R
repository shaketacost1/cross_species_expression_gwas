#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)   # fread/fwrite
  library(dplyr)        # joins, summarise
  library(arrow)        # parquet
  library(GenomicRanges)
  library(stringr)
  library(rlang)
  library(GenomeInfoDb)
  library(BiocGenerics)
})

# ---------------- CLI & Paths ----------------
args <- commandArgs(trailingOnly = TRUE)
species <- if (length(args) >= 1) args[1] else "human"

base_dir        <- file.path("data", "regulatory", species)
enh_raw_dir     <- file.path(base_dir, "enhancers_raw")   # main input dir
h3k27ac_raw_dir <- file.path(base_dir, "h3k27ac_raw")     # optional input dir
derived_dir     <- file.path(base_dir, "derived")
meta_file       <- file.path(base_dir, "enhancers_report.tsv")  # optional metadata table

dir.create(derived_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------- Helpers ----------------
read_bed_any <- function(fp) {
  # Robust reader for .bed and .bed.gz
  if (grepl("\\.gz$", fp, ignore.case = TRUE)) {
    dt <- data.table::fread(
      cmd = paste("gzip -dc", shQuote(fp)),
      sep = "\t", header = FALSE, showProgress = FALSE
    )
  } else {
    dt <- data.table::fread(
      fp, sep = "\t", header = FALSE, showProgress = FALSE
    )
  }
  if (ncol(dt) < 3) stop(paste("File has <3 columns or failed to parse:", fp))
  data.table::setnames(dt, 1:3, c("chr","start","end"))
  dt[]
}

infer_type <- function(fn_base) {
  nm <- tolower(fn_base)
  dplyr::case_when(
    grepl("els",    nm) ~ "ELS",
    grepl("ccre",   nm) ~ "cCREs",
    grepl("vista",  nm) ~ "vista",
    grepl("villar", nm) ~ "villar2015",
    grepl("h3k27ac",nm) ~ "h3k27ac",
    TRUE                ~ "enhancer_other"
  )
}

# ---------------- Gather Files ----------------
files_enh <- if (dir.exists(enh_raw_dir)) list.files(enh_raw_dir, pattern="\\.bed(\\.gz)?$", full.names = TRUE) else character(0)
files_h3k <- if (dir.exists(h3k27ac_raw_dir)) list.files(h3k27ac_raw_dir, pattern="\\.bed(\\.gz)?$", full.names = TRUE) else character(0)
bed_files <- c(files_enh, files_h3k)

if (length(bed_files) == 0) {
  stop(sprintf("No enhancer BED files found in:\n  %s\n  %s", enh_raw_dir, h3k27ac_raw_dir))
}

# ---------------- Read & Concatenate ----------------
message(sprintf("Found %d input BED file(s) for %s", length(bed_files), species))
all_list <- vector("list", length(bed_files))

for (i in seq_along(bed_files)) {
  fp    <- bed_files[i]
  base  <- sub("\\.bed(\\.gz)?$", "", basename(fp), ignore.case = TRUE)
  dt    <- read_bed_any(fp)
  dt[, `:=`(
    start  = as.integer(start),
    end    = as.integer(end),
    source = basename(fp),
    type   = infer_type(base),
    enh_id = paste0(base, "_row", .I)
  )]
  all_list[[i]] <- dt
}

enh_union <- data.table::rbindlist(all_list, use.names = TRUE, fill = TRUE)

# ---------------- Optional Metadata Join ----------------
# If you keep a table with columns like file_accession/source and biosample names
if (file.exists(meta_file)) {
  meta <- data.table::fread(meta_file, sep = "\t")
  if ("file_accession" %in% names(meta)) meta <- dplyr::rename(meta, source = file_accession)
  tissue_col <- intersect(c("biosample_term_name","tissue","cell_type","biosample_summary"), names(meta))
  if (length(tissue_col) > 0) {
    meta_small <- meta %>% dplyr::select(source, !!sym(tissue_col[1])) %>% dplyr::rename(tissue = !!sym(tissue_col[1]))
    enh_union  <- enh_union %>% dplyr::left_join(meta_small, by = "source")
  } else {
    enh_union[, tissue := NA_character_]
  }
} else {
  enh_union[, tissue := NA_character_]
}

# ---------------- Write UNION Outputs ----------------
union_prefix <- file.path(derived_dir, paste0(species, "_enhancers_union"))

# Full table (TSV + Parquet)
data.table::fwrite(enh_union, paste0(union_prefix, ".tsv.gz"), sep = "\t")
arrow::write_parquet(enh_union, paste0(union_prefix, ".parquet"))

# BED6 for broad tool compatibility (chr, start, end, name, score, strand)
bed6 <- enh_union[, .(chr, start, end, name = enh_id, score = 0L, strand = ".")]
data.table::fwrite(bed6, paste0(union_prefix, ".bed"), sep = "\t", col.names = FALSE)

# ---------------- Build MERGED (Reduced) Intervals ----------------
gr <- GenomicRanges::makeGRangesFromDataFrame(
  enh_union, keep.extra.columns = FALSE,
  seqnames.field = "chr", start.field = "start", end.field = "end",
  ignore.strand = TRUE
)
grm <- GenomicRanges::reduce(gr, ignore.strand = TRUE)

merged_dt <- data.table::data.table(
  chr   = as.character(GenomeInfoDb::seqnames(grm)),
  start = BiocGenerics::start(grm) - 1L,  # convert back to 0-based BED
  end   = BiocGenerics::end(grm)
)

merged_prefix <- file.path(derived_dir, paste0(species, "_enhancers_merged"))
data.table::fwrite(merged_dt, paste0(merged_prefix, ".bed"), sep = "\t", col.names = FALSE)
arrow::write_parquet(merged_dt, paste0(merged_prefix, ".parquet"))

# ---------------- Quick QC Tables ----------------
qc_by_source <- enh_union %>% dplyr::summarise(n_intervals = dplyr::n(), .by = source) %>% dplyr::arrange(desc(n_intervals))
qc_by_type   <- enh_union %>% dplyr::summarise(n_intervals = dplyr::n(), .by = type)   %>% dplyr::arrange(desc(n_intervals))
qc_by_tissue <- enh_union %>% dplyr::summarise(n_intervals = dplyr::n(), .by = tissue) %>% dplyr::arrange(desc(n_intervals))

data.table::fwrite(qc_by_source, file.path(derived_dir, paste0(species, "_enhancers_qc_by_source.tsv")), sep = "\t")
data.table::fwrite(qc_by_type,   file.path(derived_dir, paste0(species, "_enhancers_qc_by_type.tsv")),   sep = "\t")
data.table::fwrite(qc_by_tissue, file.path(derived_dir, paste0(species, "_enhancers_qc_by_tissue.tsv")), sep = "\t", quote = FALSE)

message("✅ Done: enhancer UNION and MERGED layers written for species: ", species)
message("   UNION:  ", union_prefix, ".[bed|tsv.gz|parquet]")
message("   MERGED: ", merged_prefix, ".[bed|parquet]")

