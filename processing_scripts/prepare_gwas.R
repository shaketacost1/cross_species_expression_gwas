## --- CRAN mirror + user library (non-interactive safe) ---
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!dir.exists(Sys.getenv("R_LIBS_USER"))) {
  dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE, showWarnings = FALSE)
}
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

## Install missing packages once
need <- setdiff(
  c("vroom","dplyr","tidyr","stringr","cli","duckdb","DBI"),
  rownames(installed.packages())
)
if (length(need)) install.packages(need)

suppressPackageStartupMessages({
  library(vroom)
  library(dplyr); library(tidyr); library(stringr); library(cli)
  library(duckdb); library(DBI)
})

## --------- SETTINGS ----------
dir.create("data/gwas/raw",   recursive = TRUE, showWarnings = FALSE)
dir.create("data/gwas/clean", recursive = TRUE, showWarnings = FALSE)

## Traits to extract (edit this list any time)
trait_specs <- list(
  alcohol = c("alcohol", "alcohol consumption", "alcohol use", "alcohol dependence", "AUD"),
  inflammation = c("C-reactive protein", "CRP", "inflammation", "inflammatory", "immune"),
  bmi = c("body mass index", "BMI")
)

## P-value threshold for genome-wide significance
pval_threshold <- 5e-8

## ---- Helper: write a data.frame to Parquet via DuckDB (no Arrow needed) ----
write_parquet_duckdb <- function(df, out_path) {
  con <- DBI::dbConnect(duckdb::duckdb(dbdir = ":memory:", read_only = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  DBI::dbWriteTable(con, "t", df, overwrite = TRUE)
  DBI::dbExecute(con, sprintf(
    "COPY (SELECT * FROM t) TO '%s' (FORMAT 'parquet');", normalizePath(out_path, mustWork = FALSE)
  ))
}

## ---- Step 1: Download the current full associations TSV from GWAS Catalog ----
## Primary documented endpoint for the curated associations dump:
full_urls <- c(
  "https://www.ebi.ac.uk/gwas/api/search/downloads/full",         # main
  "https://www.ebi.ac.uk/gwas/api/search/downloads/alternative"   # fallback (schema varies slightly)
)

#!/usr/bin/env Rscript

## --- CRAN mirror + user library (non-interactive safe) -----------------------
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!dir.exists(Sys.getenv("R_LIBS_USER"))) {
  dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE, showWarnings = FALSE)
}
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

## --- Install missing packages -------------------------------------------------
need <- setdiff(
  c("vroom","dplyr","stringr","arrow","cli"),
  rownames(installed.packages())
)
if (length(need)) install.packages(need)

suppressPackageStartupMessages({
  library(vroom)
  library(dplyr)
  library(stringr)
  library(arrow)
  library(cli)
})

## --- Paths -------------------------------------------------------------------
home_dir <- Sys.getenv("HOME")
raw_dir  <- file.path(home_dir, "data/gwas/raw")
out_dir  <- file.path(home_dir, "data/gwas/clean")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gwas_url <- "https://www.ebi.ac.uk/gwas/api/search/downloads/full"
raw_tsv  <- file.path(raw_dir, "gwas_catalog_full.tsv.gz")

## --- Download (cached) -------------------------------------------------------
if (!file.exists(raw_tsv)) {
  cli::cli_alert_info("Attempting download: {gwas_url}")
  ok <- tryCatch({
    download.file(gwas_url, destfile = raw_tsv, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) FALSE)
  if (!ok || !file.exists(raw_tsv)) {
    stop("Could not download GWAS Catalog file.")
  }
  cli::cli_alert_success("Downloaded GWAS Catalog full associations to {raw_tsv}")
} else {
  cli::cli_alert_info("Using cached GWAS Catalog at {raw_tsv}")
}

## --- Read TSV (as character, we'll coerce selectively) -----------------------
cli::cli_alert_info("Reading GWAS TSV (this can take a minute)...")
gwas <- vroom::vroom(
  raw_tsv,
  delim = "\t",
  col_types = vroom::cols(.default = "c"),
  progress = FALSE
)

## --- Normalize convenient columns --------------------------------------------
# Standardize an easy "trait" alias for DISEASE/TRAIT (column name includes a slash)
if (!"DISEASE/TRAIT" %in% names(gwas)) {
  stop("Column 'DISEASE/TRAIT' not found in GWAS file.")
}
gwas <- gwas |> mutate(trait = `DISEASE/TRAIT`)

## --- Robust p-value parsing ---------------------------------------------------
# Numeric column often exists but may be NA when p-value is stored as text.
num_col  <- "P-VALUE"
text_col <- "P-VALUE (TEXT)"
if (!(num_col %in% names(gwas)))  gwas[[num_col]]  <- NA_character_
if (!(text_col %in% names(gwas))) gwas[[text_col]] <- NA_character_

# Convert obvious numeric strings
gwas$p_num <- suppressWarnings(as.numeric(gwas[[num_col]]))

# Parse scientific notation embedded in text like "3 x 10^-9", "3×10^-9", "1e-8", "p=3e-9"
parse_p_text <- function(x) {
  if (is.na(x) || !nzchar(x)) return(NA_real_)
  y <- tolower(trimws(x))
  # unify "x" or "×" into "e" format via 10^k -> 1e{k}
  y <- gsub("\\s*[×x]\\s*", " x ", y, perl = TRUE)
  y <- gsub("10\\^\\s*([-+]?\\d+)", "1e\\1", y, perl = TRUE)   # 10^-9 -> 1e-9
  y <- gsub("\\s*x\\s*", "e", y, perl = TRUE)                  # "3 x 1e-9" -> "3e1e-9" (fix next)
  y <- gsub("e1e", "e", y, fixed = TRUE)                       # collapse "3e1e-9" -> "3e-9"
  y <- gsub("\\s+", "", y)
  # try direct numeric
  val <- suppressWarnings(as.numeric(y))
  if (!is.na(val)) return(val)
  # extract 1e-9 pattern if embedded
  m <- regmatches(y, regexpr("[-+]?\\d+(?:\\.\\d+)?e[-+]?\\d+", y, perl = TRUE))
  if (length(m) == 1) return(suppressWarnings(as.numeric(m)))
  # sometimes plain decimals like 0.00000001
  m2 <- suppressWarnings(as.numeric(y))
  if (!is.na(m2)) return(m2)
  NA_real_
}

# Fill p_num from text column when missing
fill_from_text <- is.na(gwas$p_num) & !is.na(gwas[[text_col]]) & nzchar(gwas[[text_col]])
if (any(fill_from_text)) {
  gwas$p_num[fill_from_text] <- vapply(gwas[[text_col]][fill_from_text], parse_p_text, numeric(1))
}

## --- Trait regex patterns ----------------------------------------------------
# Alcohol-related traits (consumption, dependence/AUD, etc.)
pat_alcohol <- "(?i)(\\balcohol( use| consumption| intake| dependence)?\\b|\\bAUD\\b|ethanol|drinking)"

# Inflammation umbrella + common biomarkers/phenotypes used across species
pat_inflam <- paste(
  "(?i)inflamm",                                       # inflammation, inflammatory
  "\\bCRP\\b|c-?reactive protein",                      # biomarker
  "\\bTNF\\b|tumou?r necrosis factor",                  # cytokine
  "interleukin|\\bIL-?6\\b|\\bIL6\\b",                  # cytokines
  "inflammatory bowel disease|\\bIBD\\b|crohn|ulcerative colitis",
  sep = "|"
)

# BMI (aliases)
pat_bmi <- "(?i)\\bBMI\\b|body mass index"

## --- Common output columns (keep what's available) ---------------------------
preferred_cols <- c(
  "PUBMEDID","DATE","STUDY","JOURNAL","LINK",
  "DISEASE/TRAIT","trait",
  "CHR_ID","CHR_POS","SNPS","STRONGEST SNP-RISK ALLELE",
  "RISK ALLELE FREQUENCY","OR or BETA","95% CI (TEXT)",
  "P-VALUE","P-VALUE (TEXT)",
  "REPORTED GENE(S)","MAPPED_GENE","MAPPED_TRAIT","MAPPED_TRAIT_URI",
  "CONTEXT","INTERGENIC","SNP_ID_CURRENT"
)
keep_cols <- intersect(preferred_cols, names(gwas))

## --- Helper: write filtered trait set ----------------------------------------
write_trait <- function(df, out_file, label) {
  if (nrow(df) == 0) {
    cli::cli_alert_warning("No rows matched for {label}.")
    return(invisible(NULL))
  }
  res <- df |>
    filter(!is.na(p_num), p_num <= 5e-8) |>
    select(all_of(keep_cols), p_num)

  n_all <- nrow(df); n_sig <- nrow(res)
  cli::cli_rule(label)
  cli::cli_alert_info("Hits total: {n_all}")
  cli::cli_alert_info("Hits passing p <= 5e-08: {n_sig}")

  if (n_sig > 0) {
    arrow::write_parquet(res, out_file)
    cli::cli_alert_success("Saved {n_sig} rows -> {out_file}")
  } else {
    cli::cli_alert_warning("No rows at genome-wide significance for {label}.")
  }
}

## --- Filter & write three trait groups ---------------------------------------
write_trait(
  gwas |> filter(str_detect(trait, pat_alcohol)),
  file.path(out_dir, "alcohol.parquet"),
  "alcohol"
)

write_trait(
  gwas |> filter(str_detect(trait, pat_inflam)),
  file.path(out_dir, "inflammation.parquet"),
  "inflammation"
)

write_trait(
  gwas |> filter(str_detect(trait, pat_bmi)),
  file.path(out_dir, "bmi.parquet"),
  "bmi"
)

cli::cli_h2("GWAS trait extraction complete.")
raw_tsv_gz <- "data/gwas/raw/gwas_catalog_full.tsv.gz"

download_ok <- FALSE
for (u in full_urls) {
  cli::cli_alert_info("Attempting download: {u}")
  ok <- tryCatch({
    utils::download.file(u, destfile = raw_tsv_gz, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)
  if (ok && file.exists(raw_tsv_gz) && file.info(raw_tsv_gz)$size > 1e6) {  # >1MB sanity check
    download_ok <- TRUE
    cli::cli_alert_success("Downloaded GWAS Catalog full associations to {raw_tsv_gz}")
    break
  }
}
if (!download_ok) stop("Could not download GWAS Catalog associations file from known endpoints.")

## ---- Step 2: Read with vroom (fast) and normalize column names ----
normalize_names <- function(nms) {
  x <- tolower(gsub("[^A-Za-z0-9]+", "_", nms))
  gsub("^_|_$", "", x)
}

cli::cli_alert_info("Reading GWAS TSV (this can take a minute)...")
raw_df <- vroom::vroom(raw_tsv_gz, delim = "\t", altrep = TRUE, progress = FALSE)
names(raw_df) <- normalize_names(names(raw_df))

## We’ll try to harmonize commonly used columns across format variants.
## The “full” file typically includes (among many): snps, chromosome, chromosomal_region, chr_id, chr_pos,
## p_value, p_value_text, strongest_snp_risk_allele, effect_allele, other_allele, beta, or, se,
## mapped_trait, mapped_trait_uri, reported_trait
pick <- function(df, candidates, default = NA_character_) {
  have <- intersect(candidates, names(df))
  if (length(have)) have[1] else NA_character_
}
colmap <- list(
  snps      = pick(raw_df, c("snps","rsid","variant_id")),
  chr       = pick(raw_df, c("chr_id","chromosome","chromosome_name")),
  pos       = pick(raw_df, c("chr_pos","chromosome_position","position")),
  pval      = pick(raw_df, c("p_value")),
  pval_txt  = pick(raw_df, c("p_value_text","pvalue_text")),
  eff_alle  = pick(raw_df, c("effect_allele","strongest_snp_risk_allele")),
  other_al  = pick(raw_df, c("other_allele")),
  beta      = pick(raw_df, c("beta")),
  or        = pick(raw_df, c("or","odds_ratio")),
  se        = pick(raw_df, c("standard_error","se")),
  mapped    = pick(raw_df, c("mapped_trait")),
  mapped_uri= pick(raw_df, c("mapped_trait_uri","mapped_trait_uri_s")),
  reported  = pick(raw_df, c("reported_trait"))
)

## Minimal clean frame
clean_base <- raw_df %>%
  transmute(
    rsid            = if (!is.na(colmap$snps)) .data[[colmap$snps]] else NA_character_,
    chr             = if (!is.na(colmap$chr)) .data[[colmap$chr]] else NA_character_,
    pos             = suppressWarnings(as.integer(if (!is.na(colmap$pos)) .data[[colmap$pos]] else NA)),
    pval            = suppressWarnings(as.numeric(if (!is.na(colmap$pval)) .data[[colmap$pval]] else NA)),
    pval_text       = if (!is.na(colmap$pval_txt)) .data[[colmap$pval_txt]] else NA_character_,
    effect_allele   = if (!is.na(colmap$eff_alle)) .data[[colmap$eff_alle]] else NA_character_,
    other_allele    = if (!is.na(colmap$other_al)) .data[[colmap$other_al]] else NA_character_,
    beta            = suppressWarnings(as.numeric(if (!is.na(colmap$beta)) .data[[colmap$beta]] else NA)),
    odds_ratio      = suppressWarnings(as.numeric(if (!is.na(colmap$or)) .data[[colmap$or]] else NA)),
    se              = suppressWarnings(as.numeric(if (!is.na(colmap$se)) .data[[colmap$se]] else NA)),
    mapped_trait    = if (!is.na(colmap$mapped)) .data[[colmap$mapped]] else NA_character_,
    mapped_trait_uri= if (!is.na(colmap$mapped_uri)) .data[[colmap$mapped_uri]] else NA_character_,
    reported_trait  = if (!is.na(colmap$reported)) .data[[colmap$reported]] else NA_character_
  ) %>%
  ## fix common weirdness
  mutate(
    rsid = ifelse(!is.na(rsid), str_extract(rsid, "rs[0-9]+"), rsid),
    chr  = as.character(chr)
  ) %>%
  ## keep rows with at least rsid & coordinates if present
  filter(!is.na(rsid) | (!is.na(chr) & !is.na(pos)))

## ---- Step 3: Per-trait filtering + write Parquet ----
select_cols <- c("rsid","chr","pos","effect_allele","other_allele",
                 "beta","odds_ratio","se","pval","pval_text",
                 "mapped_trait","mapped_trait_uri","reported_trait")

for (trait_name in names(trait_specs)) {
  keywords <- trait_specs[[trait_name]]
  rx <- paste0(keywords, collapse = "|")
  cli::cli_h1(trait_name)

  trait_df <- clean_base %>%
    filter(
      !is.na(pval) & pval <= pval_threshold,
      (
        (!is.na(mapped_trait)  & str_detect(tolower(mapped_trait),  tolower(rx))) |
        (!is.na(reported_trait)& str_detect(tolower(reported_trait),tolower(rx)))
      )
    ) %>%
    distinct(across(all_of(select_cols)))

  cli::cli_alert_info("Hits passing p <= {format(pval_threshold, scientific = TRUE)}: {nrow(trait_df)}")

  out_parquet <- file.path("data/gwas/clean", paste0(trait_name, ".parquet"))
  out_csv     <- file.path("data/gwas/clean", paste0(trait_name, ".csv"))

  if (nrow(trait_df)) {
    # write Parquet (via DuckDB) and also a CSV for quick inspection
    write_parquet_duckdb(trait_df, out_parquet)
    readr::write_csv(trait_df, out_csv)
    cli::cli_alert_success("Wrote: {out_parquet}")
    cli::cli_alert_success("Wrote: {out_csv}")
  } else {
    cli::cli_alert_warning("No rows matched for trait: {trait_name}.")
  }
}

cli::cli_h2("GWAS trait extraction complete.")

