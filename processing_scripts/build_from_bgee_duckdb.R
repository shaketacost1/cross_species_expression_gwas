# install once
install.packages(c("duckdb","DBI","vroom","dplyr","tidyr","stringr","cli"))

suppressPackageStartupMessages({
  library(duckdb); library(DBI); library(vroom)
  library(dplyr);  library(tidyr); library(stringr); library(cli)
})

species   <- c("human","mouse","macaque")   # add "pig","chicken" later
in_root   <- "data/raw/bgee"
out_root  <- "data/expression"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

normalize_names <- function(nms) {
  n <- tolower(nms); n <- gsub("[^a-z0-9]+","_",n); gsub("^_|_$","",n)
}
pick_col <- function(df, options) {
  cols <- names(df); low <- tolower(cols)
  for (opt in options) { j <- which(low == tolower(opt)); if (length(j)) return(cols[j[1]]) }
  ncols <- normalize_names(cols); nopts <- normalize_names(options)
  for (opt in nopts) { j <- which(ncols == opt); if (length(j)) return(cols[j[1]]) }
  for (opt in nopts) { j <- which(grepl(opt, ncols, fixed=TRUE)); if (length(j)) return(cols[j[1]]) }
  NA_character_
}
ensure_extracted <- function(in_dir){
  tsvs <- list.files(in_dir, pattern="\\.(tsv|txt)(\\.gz)?$", full.names=TRUE, recursive=TRUE, ignore.case=TRUE)
  if (length(tsvs)>0) return(invisible())
  tars <- list.files(in_dir, pattern="\\.tar(\\.gz)?$", full.names=TRUE, recursive=FALSE, ignore.case=TRUE)
  if (length(tars)==0) return(invisible())
  cli::cli_alert_info("Extracting: {basename(tars[1])}")
  utils::untar(tars[1], exdir=in_dir)
}
standardize_tbl <- function(df_raw){
  colnames(df_raw) <- normalize_names(colnames(df_raw)); df <- df_raw
  geneIdCol <- pick_col(df, c("Gene ID","EnsemblGeneId","GeneID","gene_id","ensembl_gene_id","bgeegeneid"))
  libCol    <- pick_col(df, c("Library ID","RNA-seq library ID","rnaSeqLibraryId","Sample ID","library_id","sample_id","rna_seq_library_id"))
  tpmCol    <- pick_col(df, c("TPM","tpm","TPM value","expression TPM","tpm_value"))
  if (any(is.na(c(geneIdCol,libCol,tpmCol)))) return(tibble())
  geneNameCol <- pick_col(df, c("gene_name","symbol","associated_gene_name","gene_symbol"))
  anatIdCol   <- pick_col(df, c("Anatomical entity ID","Uberon ID","anat_id","anatomical_entity_id"))
  anatNameCol <- pick_col(df, c("Anatomical entity name","Tissue","Organ","anatomical_entity_name","tissue","organ"))
  stageIdCol  <- pick_col(df, c("Stage ID","Developmental stage ID","stage_id","dev_stage_id"))
  stageCol    <- pick_col(df, c("Stage name","Stage","Developmental stage","stage_name","dev_stage_name"))
  sexCol      <- pick_col(df, c("Sex","sex"))
  strainCol   <- pick_col(df, c("Strain","strain"))
  platCol     <- pick_col(df, c("Library type","Platform","Library platform","library_platform","platform","library type"))
  tibble(
    geneId    = as.character(df[[geneIdCol]]),
    geneName  = if (!is.na(geneNameCol)) as.character(df[[geneNameCol]]) else NA_character_,
    libraryId = as.character(df[[libCol]]),
    tpm       = suppressWarnings(as.numeric(df[[tpmCol]])),
    anatId    = if (!is.na(anatIdCol)) as.character(df[[anatIdCol]]) else NA_character_,
    anatName  = if (!is.na(anatNameCol)) as.character(df[[anatNameCol]]) else NA_character_,
    stageId   = if (!is.na(stageIdCol)) as.character(df[[stageIdCol]]) else NA_character_,
    stage     = if (!is.na(stageCol)) as.character(df[[stageCol]]) else NA_character_,
    sex       = if (!is.na(sexCol)) as.character(df[[sexCol]]) else NA_character_,
    strain    = if (!is.na(strainCol)) as.character(df[[strainCol]]) else NA_character_,
    platform  = if (!is.na(platCol)) as.character(df[[platCol]]) else NA_character_
  ) |>
    filter(!is.na(geneId), geneId!="", !is.na(libraryId), libraryId!="", !is.na(tpm))
}

process_species_duckdb <- function(sp){
  cli::cli_h1(sp)
  in_dir  <- file.path(in_root, sp, "RNA_SEQ")
  out_dir <- file.path(out_root, sp); dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  ensure_extracted(in_dir)

  tsvs <- list.files(in_dir, pattern="\\.(tsv|txt)(\\.gz)?$", full.names=TRUE, recursive=TRUE, ignore.case=TRUE)
  cli::cli_alert_info("Found {length(tsvs)} TSV files under {in_dir}")
  if (!length(tsvs)) stop("No TSV files for ", sp)

  # DuckDB file for this species
  dbfile <- file.path(out_dir, "bgee.duckdb")
  con <- dbConnect(duckdb::duckdb(dbdir=dbfile, read_only=FALSE))
  on.exit({ dbDisconnect(con, shutdown=TRUE) }, add=TRUE)

  # Prepare tables
  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS expr_long (
      geneId TEXT, libraryId TEXT, tpm DOUBLE
    );
  ")
  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS sample_meta (
      libraryId TEXT, anatId TEXT, anatName TEXT,
      stageId TEXT, stage TEXT, sex TEXT, strain TEXT, platform TEXT
    );
  ")

  # Stream each file: standardize, then append straight to DuckDB
  for (i in seq_along(tsvs)) {
    f <- tsvs[i]
    if (i %% 25 == 1 || i == length(tsvs)) cli::cli_alert("Reading {i}/{length(tsvs)}: {basename(f)}")
    df <- try(vroom::vroom(f, delim="\t", col_types = vroom::cols(.default="c"), progress=FALSE), silent=TRUE)
    if (inherits(df, "try-error")) next
    std <- try(standardize_tbl(df), silent=TRUE)
    rm(df); gc()
    if (inherits(std, "try-error") || !nrow(std)) next
    # write to DuckDB (append)
    DBI::dbWriteTable(con, "expr_long", std |> select(geneId, libraryId, tpm), append=TRUE)
    DBI::dbWriteTable(con, "sample_meta", std |> select(libraryId, anatId, anatName, stageId, stage, sex, strain, platform) |> distinct(), append=TRUE)
    rm(std); gc()
  }

  # De-dup metadata in DB
  dbExecute(con, "CREATE OR REPLACE TABLE sample_meta AS SELECT DISTINCT * FROM sample_meta;")

  # Write metadata CSV directly from DuckDB to disk (no big R dataframe)
  DBI::dbExecute(con, sprintf("
    COPY (
      SELECT * FROM sample_meta ORDER BY libraryId
    ) TO '%s' (HEADER, DELIMITER ',')
  ", file.path(out_dir, "sample_metadata.csv")))

  # TPM aggregation to one row per (geneId, libraryId)
  dbExecute(con, "
    CREATE OR REPLACE TABLE agg AS
    SELECT geneId, libraryId, AVG(tpm) AS tpm
    FROM expr_long
    GROUP BY geneId, libraryId;
  ")

  # Build a wide matrix *in DuckDB* and write it straight to CSV.
  # 1) get all libraryIds
  libs <- dbGetQuery(con, "SELECT libraryId FROM sample_meta ORDER BY libraryId;")$libraryId
  if (!length(libs)) stop("No libraries found.")

  # 2) dynamic SELECT with CASE WHEN to pivot on disk
  cols_sql <- paste(
    sprintf("COALESCE(MAX(CASE WHEN libraryId = '%s' THEN tpm END), 0) AS \"%s\"", libs, libs),
    collapse = ",\n    "
  )
  wide_sql <- sprintf("
    COPY (
      SELECT geneId,
             %s
      FROM agg
      GROUP BY geneId
      ORDER BY geneId
    ) TO '%s' (HEADER, DELIMITER ',');
  ", cols_sql, file.path(out_dir, "tpm_matrix.csv"))

  DBI::dbExecute(con, wide_sql)

  cli::cli_alert_success("Wrote: {file.path(out_dir,'tpm_matrix.csv')}")
  cli::cli_alert_success("Wrote: {file.path(out_dir,'sample_metadata.csv')}")
}

for (sp in species) process_species_duckdb(sp)


