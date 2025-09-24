#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  pkgs <- c("data.table","arrow","dplyr","stringr","GenomicRanges")
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) {
    if (p == "GenomicRanges") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install("GenomicRanges", ask = FALSE, update = FALSE)
    } else install.packages(p)
  }
  library(data.table); library(arrow); library(dplyr); library(stringr); library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
species <- if (length(args) >= 1) args[1] else "human"

base_dir    <- file.path("data","regulatory", species)
raw_dir     <- file.path(base_dir, "ctcf_raw")
derived_dir <- file.path(base_dir, "derived")
meta_file   <- file.path(base_dir, "ctcf_report.tsv") # optional
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)

read_bed_any <- function(fp) {
  if (grepl("\\.gz$", fp, ignore.case = TRUE)) {
    dt <- data.table::fread(cmd = paste("gzip -dc", shQuote(fp)),
                            sep = "\t", header = FALSE, data.table = TRUE, showProgress = FALSE)
  } else {
    dt <- data.table::fread(fp, sep = "\t", header = FALSE, data.table = TRUE, showProgress = FALSE)
  }
  stopifnot(ncol(dt) >= 3)
  data.table::setnames(dt, names(dt)[1:3], c("chr","start","end"))
  dt[]
}

bed_files <- list.files(raw_dir, pattern = "\\.bed(\\.gz)?$", full.names = TRUE)
if (length(bed_files) == 0) stop("No BED files found in: ", raw_dir)

all_list <- vector("list", length(bed_files))
for (i in seq_along(bed_files)) {
  fp    <- bed_files[i]
  encff <- sub("\\.bed(\\.gz)?$", "", basename(fp), ignore.case = TRUE) # ENCFF...
  dt    <- read_bed_any(fp)
  dt[, source := encff]
  dt[, ctcf_id := paste0(encff, "_peak", .I)]
  all_list[[i]] <- dt
}
ctcf_union <- data.table::rbindlist(all_list, use.names = TRUE, fill = TRUE)
ctcf_union[, start := as.integer(start)]
ctcf_union[, end   := as.integer(end)]

# optional: join metadata (if present)
if (file.exists(meta_file)) {
  meta <- fread(meta_file)
  if ("file_accession" %in% names(meta)) meta <- dplyr::rename(meta, source = file_accession)
  tcol <- intersect(c("biosample_term_name","tissue","cell_type","biosample_summary"), names(meta))
  if (length(tcol)) {
    meta_small <- meta %>% dplyr::select(source, !!rlang::sym(tcol[1])) %>% dplyr::rename(tissue = !!rlang::sym(tcol[1]))
    ctcf_union <- dplyr::left_join(ctcf_union, meta_small, by = "source")
  } else {
    ctcf_union[, tissue := NA_character_]
  }
} else ctcf_union[, tissue := NA_character_]

# write union
bed_union_path <- file.path(derived_dir, paste0(species, "_ctcf_union.bed"))
bed_union_dt   <- ctcf_union[, .(chr, start, end, name = ctcf_id, score = 0L, strand = ".")]
fwrite(bed_union_dt, bed_union_path, sep = "\t", col.names = FALSE)
tsv_union_path <- file.path(derived_dir, paste0(species, "_ctcf_union.tsv.gz"))
fwrite(ctcf_union, tsv_union_path, sep = "\t")
parquet_union_path <- file.path(derived_dir, paste0(species, "_ctcf_union.parquet"))
arrow::write_parquet(ctcf_union, parquet_union_path)

# merge intervals (reduce)
gr <- GenomicRanges::GRanges(seqnames = ctcf_union$chr,
                             ranges = IRanges::IRanges(start = ctcf_union$start + 1L, end = ctcf_union$end))
gr_merged <- GenomicRanges::reduce(gr, ignore.strand = TRUE)
merged_dt <- data.table(chr   = as.character(GenomeInfoDb::seqnames(gr_merged)),
                        start = BiocGenerics::start(gr_merged) - 1L,
                        end   = BiocGenerics::end(gr_merged))
bed_merged_path     <- file.path(derived_dir, paste0(species, "_ctcf_merged.bed"))
parquet_merged_path <- file.path(derived_dir, paste0(species, "_ctcf_merged.parquet"))
fwrite(merged_dt, bed_merged_path, sep = "\t", col.names = FALSE)
arrow::write_parquet(merged_dt, parquet_merged_path)

# tiny QC tables
ctcf_union[, .N, by = source][order(-N)] |>
  fwrite(file.path(derived_dir, paste0(species, "_ctcf_qc_by_source.tsv")), sep="\t")
if ("tissue" %in% names(ctcf_union))
  ctcf_union[, .N, by = tissue][order(-N)] |>
  fwrite(file.path(derived_dir, paste0(species, "_ctcf_qc_by_tissue.tsv")), sep="\t")

message("✅ Done: ", species,
        " — union/merged written to ", file.path("data","regulatory",species,"derived"))

