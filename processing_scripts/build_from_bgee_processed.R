suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(stringr) })
species <- c("human")
in_root  <- "data/raw/bgee"; out_root <- "data/expression"; dir.create(out_root, recursive=TRUE, showWarnings=FALSE)

normalize_names <- function(nms){
  n <- tolower(nms); n <- gsub("[^a-z0-9]+","_",n); gsub("^_|_$","",n)
}
has_col <- function(df, key){ if (key %in% names(df)) key else NA_character_ }

standardize_tbl <- function(df_raw){
  colnames(df_raw) <- normalize_names(colnames(df_raw))
  df <- df_raw
  geneIdCol <- has_col(df,"gene_id"); if (is.na(geneIdCol)) {
    cand <- intersect(c("ensemblgeneid","ensembl_gene_id","bgeegeneid"), names(df)); if (length(cand)) geneIdCol <- cand[1]
  }
  libCol <- has_col(df,"library_id"); if (is.na(libCol)) {
    cand <- intersect(c("sample_id","rnaseqlibraryid","rna_seq_library_id"), names(df)); if (length(cand)) libCol <- cand[1]
  }
  tpmCol <- if ("tpm" %in% names(df)) "tpm" else {
    cand <- grep("tpm", names(df), value=TRUE); if (length(cand)) cand[1] else NA_character_
  }
  if (any(is.na(c(geneIdCol, libCol, tpmCol)))) return(tibble())

  geneNameCol <- (intersect(c("gene_name","symbol","associated_gene_name","gene_symbol"), names(df)))[1]
  anatIdCol   <- (intersect(c("anatomical_entity_id","uberon_id","anat_entity_id"), names(df)))[1]
  anatNameCol <- (intersect(c("anatomical_entity_name","tissue","organ","anat_entity_name"), names(df)))[1]
  stageIdCol  <- (intersect(c("stage_id","dev_stage_id"), names(df)))[1]
  stageCol    <- (intersect(c("stage_name","dev_stage_name","stage"), names(df)))[1]
  sexCol      <- (intersect(c("sex"), names(df)))[1]
  strainCol   <- (intersect(c("strain"), names(df)))[1]
  platformCol <- (intersect(c("platform","library_platform","libraryplatform","library_type"), names(df)))[1]

  tibble(
    geneId    = as.character(df[[geneIdCol]]),
    geneName  = if (!is.na(geneNameCol)) as.character(df[[geneNameCol]]) else NA_character_,
    libraryId = as.character(df[[libCol]]),
    tpm_raw   = as.character(df[[tpmCol]]),
    anatId    = if (!is.na(anatIdCol)) as.character(df[[anatIdCol]]) else NA_character_,
    anatName  = if (!is.na(anatNameCol)) as.character(df[[anatNameCol]]) else NA_character_,
    stageId   = if (!is.na(stageIdCol)) as.character(df[[stageIdCol]]) else NA_character_,
    stage     = if (!is.na(stageCol)) as.character(df[[stageCol]]) else NA_character_,
    sex       = if (!is.na(sexCol)) as.character(df[[sexCol]]) else NA_character_,
    strain    = if (!is.na(strainCol)) as.character(df[[strainCol]]) else NA_character_,
    platform  = if (!is.na(platformCol)) as.character(df[[platformCol]]) else NA_character_
  ) |>
    dplyr::filter(!is.na(geneId), geneId != "", !is.na(libraryId), libraryId != "")
}

process_species <- function(sp){
  message("== ", sp, " ==")
  in_dir <- file.path(in_root, sp, "RNA_SEQ")
  out_dir <- file.path(out_root, sp)
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

  # RECURSIVE search (human extracts into subdirs)
  tsvs <- list.files(in_dir, pattern="\\.(tsv|txt)(\\.gz)?$", full.names=TRUE, ignore.case=TRUE, recursive=TRUE)
  if (length(tsvs) == 0) {
    tarball <- list.files(in_dir, pattern="_RNA-Seq_read_counts_TPM\\.tar\\.gz$", full.names=TRUE, recursive=FALSE)
    if (length(tarball) == 1) {
      message("   Extracting: ", basename(tarball))
      utils::untar(tarball, exdir=in_dir)
      tsvs <- list.files(in_dir, pattern="\\.(tsv|txt)(\\.gz)?$", full.names=TRUE, ignore.case=TRUE, recursive=TRUE)
    } else stop("No TSVs and no tarball in ", in_dir)
  }

  all_rows <- list()
  for (f in tsvs) {
    df <- try(suppressMessages(readr::read_tsv(f, guess_max=100000, progress=FALSE,
                                               col_types=readr::cols(.default="c"))), silent=TRUE)
    if (inherits(df, "try-error") || !is.data.frame(df)) next
    std <- try(standardize_tbl(df), silent=TRUE)
    if (!inherits(std, "try-error") && nrow(std) > 0) all_rows[[length(all_rows) + 1]] <- std
  }
  if (length(all_rows) == 0) stop("No usable expression rows for ", sp)

  expr <- dplyr::bind_rows(all_rows) |>
          dplyr::distinct(geneId, geneName, libraryId, .keep_all=TRUE) |>
          dplyr::mutate(tpm = suppressWarnings(as.numeric(tpm_raw))) |>
          dplyr::select(-tpm_raw)

  tpm_mat <- expr |>
    dplyr::select(geneId, libraryId, tpm) |>
    dplyr::group_by(geneId, libraryId) |>
    dplyr::summarize(tpm = mean(tpm, na.rm=TRUE), .groups="drop") |>
    tidyr::pivot_wider(names_from=libraryId, values_from=tpm, values_fill=0) |>
    dplyr::arrange(geneId)

  sample_meta <- expr |>
    dplyr::select(libraryId, anatId, anatName, stageId, stage, sex, strain, platform) |>
    dplyr::distinct() |>
    dplyr::arrange(libraryId)

  readr::write_csv(tpm_mat, file.path(out_dir, "tpm_matrix.csv"))
  readr::write_csv(sample_meta, file.path(out_dir, "sample_metadata.csv"))
  message("   Wrote: ", file.path(out_dir, "tpm_matrix.csv"))
  message("   Wrote: ", file.path(out_dir, "sample_metadata.csv"))
}

for (sp in species) process_species(sp)
message("All species done.")
