#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DBI); library(RSQLite); library(readr); library(stringr)
})

app_dir  <- normalizePath("~/Desktop/Project Alpha", mustWork = FALSE)
src_root <- normalizePath("~/my_project", mustWork = FALSE)
db_path  <- file.path(app_dir, "data", "regland.sqlite")
expr_out <- file.path(app_dir, "data", "expression_tpm.tsv")
dir.create(dirname(db_path), showWarnings = FALSE, recursive = TRUE)

con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
on.exit(try(DBI::dbDisconnect(con), silent = TRUE), add = TRUE)

DBI::dbExecute(con,"CREATE TABLE IF NOT EXISTS genes(gene_id INTEGER PRIMARY KEY, symbol TEXT, species_id TEXT, chrom TEXT, start INTEGER, end INTEGER);")
DBI::dbExecute(con,"CREATE TABLE IF NOT EXISTS enhancers(enh_id INTEGER PRIMARY KEY, species_id TEXT, chrom TEXT, start INTEGER, end INTEGER, tissue TEXT, score REAL, source TEXT);")
DBI::dbExecute(con,"CREATE TABLE IF NOT EXISTS enhancer_class(enh_id INTEGER, class TEXT);")
DBI::dbExecute(con,"CREATE TABLE IF NOT EXISTS gwas_snps(snp_id INTEGER PRIMARY KEY, rsid TEXT UNIQUE, chrom TEXT, pos INTEGER, trait TEXT, pval REAL, source TEXT, category TEXT);")
DBI::dbExecute(con,"CREATE TABLE IF NOT EXISTS snp_to_enhancer(snp_id INTEGER, enh_id INTEGER, UNIQUE(snp_id,enh_id));")
DBI::dbExecute(con,"CREATE TABLE IF NOT EXISTS gene_to_enhancer(gene_id INTEGER, enh_id INTEGER, UNIQUE(gene_id,enh_id));")
DBI::dbExecute(con,"CREATE TABLE IF NOT EXISTS ctcf_sites(site_id INTEGER PRIMARY KEY, species_id TEXT, chrom TEXT, start INTEGER, end INTEGER, score REAL, motif_p REAL, cons_class TEXT);")
DBI::dbExecute(con,"CREATE TABLE IF NOT EXISTS tad_domains(tad_id INTEGER PRIMARY KEY, species_id TEXT, chrom TEXT, start INTEGER, end INTEGER, source TEXT);")

invisible(try(DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_genes_species_chrom_start ON genes(species_id,chrom,start);"), TRUE))
invisible(try(DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_enh_species_chrom_bounds  ON enhancers(species_id,chrom,start,end);"), TRUE))
invisible(try(DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_enh_tissue ON enhancers(tissue);"), TRUE))
invisible(try(DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_ec_enh ON enhancer_class(enh_id);"), TRUE))
invisible(try(DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_snp_chrom_pos ON gwas_snps(chrom,pos);"), TRUE))
invisible(try(DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_snp_rsid ON gwas_snps(rsid);"), TRUE))
invisible(try(DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_ctcf_species_chrom_bounds ON ctcf_sites(species_id,chrom,start,end);"), TRUE))
invisible(try(DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_ctcf_cons_class ON ctcf_sites(cons_class);"), TRUE))
invisible(try(DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_tad_species_chrom_bounds ON tad_domains(species_id,chrom,start,end);"), TRUE))

map_species_dir <- function(path){
  nm <- basename(path)
  if (grepl("human|hsa|hg38", nm, TRUE)) return("human_hg38")
  if (grepl("mouse|mmus|mm39", nm, TRUE)) return("mouse_mm39")
  if (grepl("macaque|rhemac|rheMac", nm, TRUE)) return("macaque_rheMac10")
  if (grepl("pig|susScr", nm, TRUE)) return("susScr")
  NA_character_
}
sql_write <- function(tbl, df){
  if (!is.null(df) && nrow(df)) {
    DBI::dbBegin(con); on.exit(try(DBI::dbRollback(con), silent=TRUE), add=TRUE)
    DBI::dbWriteTable(con, tbl, df, append = TRUE)
    DBI::dbCommit(con); on.exit(NULL, add=FALSE)
  }
}

# -------- 1) CTCF (base R only) --------
ctcf_dir <- file.path(src_root, "data", "capstone data")
if (dir.exists(ctcf_dir)) {
  files <- list.files(ctcf_dir, pattern = "CTCF\\.(bed|bed\\.gz)$", full.names = TRUE)
  message("CTCF files: ", length(files))
  for (f in files) {
    sp <- if (grepl("human", basename(f), TRUE)) "human_hg38"
          else if (grepl("mouse", basename(f), TRUE)) "mouse_mm39"
          else NA_character_
    if (is.na(sp)) next

    bed <- tryCatch(readr::read_tsv(f, col_names = FALSE, show_col_types = FALSE, progress = FALSE),
                    error = function(e) NULL)
    if (is.null(bed) || !nrow(bed)) next

    n <- ncol(bed)
    cols <- c("chrom","start","end","name","score","strand")[seq_len(min(6, n))]
    names(bed)[seq_along(cols)] <- cols
    if (!"score" %in% names(bed)) bed$score <- NA_real_

    df <- data.frame(
      species_id = rep(sp, nrow(bed)),
      chrom = as.character(bed$chrom),
      start = as.integer(bed$start),
      end   = as.integer(bed$end),
      score = suppressWarnings(as.numeric(bed$score)),
      stringsAsFactors = FALSE
    )
    df <- df[!is.na(df$chrom) & !is.na(df$start) & !is.na(df$end), , drop = FALSE]
    df$motif_p    <- NA_real_
    df$cons_class <- NA_character_
    message("CTCF df cols: ", paste(names(df), collapse=","))

    tryCatch({ message("CTCF writing rows: ", nrow(df)); DBI::dbWriteTable(con, "ctcf_sites", as.data.frame(df, stringsAsFactors=FALSE), append=TRUE); message("CTCF wrote rows OK") }, error=function(e){ message("CTCF write error -> ", conditionMessage(e)); traceback(2); quit(status=1) })
  }
} else {
  message("CTCF directory missing: ", ctcf_dir)
}

# -------- 2) Enhancers (base R) --------
reg_dir <- file.path(src_root, "data", "regulatory")
if (dir.exists(reg_dir)) {
  sp_dirs <- list.dirs(reg_dir, recursive = FALSE)
  for (sd in sp_dirs) {
    sp <- map_species_dir(sd); if (is.na(sp)) next
    beds <- list.files(sd, pattern="\\.bed(\\.gz)?$", full.names = TRUE, recursive = TRUE)
    if (!length(beds)) next
    message("Enhancers: ", basename(sd), " (", sp, ") files=", length(beds))
    for (f in beds) {
      bed <- tryCatch(readr::read_tsv(f, col_names = FALSE, show_col_types = FALSE, progress = FALSE),
                      error=function(e) NULL)
      if (is.null(bed) || !nrow(bed)) next
      n <- ncol(bed)
      cols <- c("chrom","start","end","name","score","strand")[seq_len(min(6, n))]
      names(bed)[seq_along(cols)] <- cols
      if (!"score" %in% names(bed)) bed$score <- NA_real_

      df <- data.frame(
        species_id = rep(sp, nrow(bed)),
        chrom = as.character(bed$chrom),
        start = as.integer(bed$start),
        end   = as.integer(bed$end),
        tissue = NA_character_,
        score  = suppressWarnings(as.numeric(bed$score)),
        source = basename(f),
        stringsAsFactors = FALSE
      )
      df <- df[!is.na(df$chrom) & !is.na(df$start) & !is.na(df$end), , drop = FALSE]
      sql_write("enhancers", df)
    }
  }
} else {
  message("Regulatory dir missing: ", reg_dir)
}

# -------- 3) GWAS (base R-ish) --------
gwas_dirs  <- file.path(src_root, "data", "gwas", c("clean","raw"))
gwas_files <- unique(unlist(lapply(gwas_dirs, function(d) if (dir.exists(d)) list.files(d, pattern="\\.tsv(\\.gz)?$", full.names=TRUE) else character(0))))
ingest_gwas <- function(f){
  df <- tryCatch(readr::read_tsv(f, show_col_types = FALSE, progress = FALSE, guess_max = 100000), error=function(e) NULL)
  if (is.null(df) || !nrow(df)) return(NULL)
  nm <- names(df); pick <- function(x){ y <- intersect(x, nm); if (length(y)) y[1] else NA_character_ }
  rs  <- pick(c("rsid","SNP","snp","variant_id"))
  chr <- pick(c("chrom","chr","CHR","chromosome"))
  pos <- pick(c("pos","position","BP","bp"))
  pv  <- pick(c("p","pval","P","PVAL","p_value"))
  tr  <- pick(c("trait","phenotype","TRAIT"))
  if (is.na(chr) || is.na(pos)) return(NULL)
  out <- data.frame(
    rsid  = if (!is.na(rs)) as.character(df[[rs]]) else NA_character_,
    chrom = as.character(df[[chr]]),
    pos   = as.integer(df[[pos]]),
    trait = if (!is.na(tr)) as.character(df[[tr]]) else basename(f),
    pval  = suppressWarnings(as.numeric(if (!is.na(pv)) df[[pv]] else NA_real_)),
    source = basename(f),
    category = NA_character_,
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$chrom) & !is.na(out$pos), , drop = FALSE]
  if (!("rsid" %in% names(out))) out$rsid <- NA_character_
  out <- out[!duplicated(out[c("rsid","chrom","pos","trait")]), , drop = FALSE]
  out
}
if (length(gwas_files)) {
  message("GWAS files: ", length(gwas_files))
  gwas_list <- lapply(gwas_files, ingest_gwas)
  gwas_all <- do.call(rbind, Filter(function(z) !is.null(z) && nrow(z), gwas_list))
  if (!is.null(gwas_all) && nrow(gwas_all)) sql_write("gwas_snps", gwas_all) else message("No usable GWAS rows.")
} else {
  message("No GWAS files found under data/gwas/{clean,raw}")
}

# -------- 4) Expression -> expression_tpm.tsv (base R) --------
expr_root <- file.path(src_root, "data", "expression")
collect_expr <- function(){
  if (!dir.exists(expr_root)) return(data.frame())
  spp <- list.dirs(expr_root, recursive = FALSE)
  files <- unlist(lapply(spp, function(d) c(list.files(d, pattern="\\.tsv$", full.names=TRUE),
                                            list.files(d, pattern="\\.csv$", full.names=TRUE))))
  if (!length(files)) return(data.frame())
  out <- list()
  for (f in files) {
    delim <- if (grepl("\\.csv$", f)) "," else "\t"
    df <- tryCatch(readr::read_delim(f, delim=delim, show_col_types=FALSE, progress=FALSE, guess_max=50000),
                   error=function(e) NULL)
    if (is.null(df) || !nrow(df)) next
    nm <- names(df)
    symbol <- intersect(c("symbol","gene","gene_symbol","Gene","SYMBOL"), nm)
    tpm    <- intersect(c("tpm","TPM","expression","expr"), nm)
    tissue <- intersect(c("tissue","Tissue","tissue_name"), nm)
    if (!length(symbol) || !length(tpm)) next
    out[[length(out)+1]] <- data.frame(
      symbol = as.character(df[[symbol[1]]]),
      tissue = if (length(tissue)) as.character(df[[tissue[1]]]) else "All",
      tpm    = suppressWarnings(as.numeric(df[[tpm[1]]])),
      source = basename(f),
      stringsAsFactors = FALSE
    )
  }
  if (!length(out)) return(data.frame())
  do.call(rbind, out)
}
expr_df <- collect_expr()
if (nrow(expr_df)) {
  expr_df$symbol <- toupper(expr_df$symbol)
  expr_df$tissue <- ifelse(grepl("brain", expr_df$tissue, TRUE), "Brain",
                     ifelse(grepl("heart|card", expr_df$tissue, TRUE), "Heart",
                     ifelse(grepl("liver|hep", expr_df$tissue, TRUE), "Liver", "Other")))
  expr_df <- expr_df[expr_df$tissue %in% c("Brain","Heart","Liver"), , drop=FALSE]
  if (nrow(expr_df)) {
    # median by symbol+tissue
    key <- paste(expr_df$symbol, expr_df$tissue, sep="|")
    split_idx <- split(seq_len(nrow(expr_df)), key)
    med_rows <- lapply(split_idx, function(ix){
      data.frame(symbol = expr_df$symbol[ix[1]],
                 tissue = expr_df$tissue[ix[1]],
                 tpm = median(expr_df$tpm[ix], na.rm=TRUE),
                 stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, med_rows)
    readr::write_tsv(out, expr_out)
    message("Wrote ", expr_out)
  } else {
    readr::write_tsv(data.frame(symbol=character(), tissue=character(), tpm=double()), expr_out)
    message("No expression rows after filtering; wrote empty ", expr_out)
  }
} else {
  readr::write_tsv(data.frame(symbol=character(), tissue=character(), tpm=double()), expr_out)
  message("No expression rows found; wrote empty ", expr_out)
}

# -------- 5) Optional TADs (base R) --------
tad_files <- if (dir.exists(reg_dir)) list.files(reg_dir, pattern="tad.*\\.bed$", full.names=TRUE, recursive=TRUE) else character(0)
if (length(tad_files)) {
  message("TAD files: ", length(tad_files))
  for (f in tad_files){
    sp <- map_species_dir(dirname(f)); if (is.na(sp)) next
    bed <- tryCatch(readr::read_tsv(f, col_names=FALSE, show_col_types=FALSE, progress=FALSE), error=function(e) NULL)
    if (is.null(bed) || !nrow(bed)) next
    n <- ncol(bed); cols <- c("chrom","start","end","name")[seq_len(min(4, n))]
    names(bed)[seq_along(cols)] <- cols
    df <- data.frame(
      species_id = sp,
      chrom = as.character(bed$chrom),
      start = as.integer(bed$start),
      end   = as.integer(bed$end),
      source = basename(f),
      stringsAsFactors = FALSE
    )
    df <- df[!is.na(df$chrom) & !is.na(df$start) & !is.na(df$end), , drop=FALSE]
    sql_write("tad_domains", df)
  }
}

# -------- summary --------
n_tbl <- function(tbl) tryCatch(DBI::dbGetQuery(con, paste0("SELECT COUNT(*) AS n FROM ", tbl))$n, error=function(e) NA_integer_)
cat("\n=== Bootstrap summary ===\n",
    "ctcf_sites:   ", n_tbl("ctcf_sites"), "\n",
    "enhancers:    ", n_tbl("enhancers"), "\n",
    "gwas_snps:    ", n_tbl("gwas_snps"), "\n",
    "tad_domains:  ", n_tbl("tad_domains"), "\n",
    "expression_tpm.tsv rows: ", tryCatch(nrow(readr::read_tsv(expr_out, show_col_types=FALSE)), error=function(e) 0), "\n", sep="")