f <- "scripts/bootstrap_data.R"
x <- readLines(f)

# locate the CTCF block
start <- grep('^ctcf_dir <- file\\.path', x)
end   <- grep('^reg_dir <- file\\.path', x)
if (!length(start) || !length(end)) stop("Could not locate CTCF or next block marker")
end <- end[1L] - 1L

# new CTCF block (base R only)
new <- c(
  'ctcf_dir <- file.path(src_root, "data", "capstone data")',
  'if (dir.exists(ctcf_dir)) {',
  '  files <- list.files(ctcf_dir, pattern = "CTCF\\\\.(bed|bed\\\\.gz)$", full.names = TRUE)',
  '  message("CTCF files: ", length(files))',
  '  for (f in files) {',
  '    sp <- if (grepl("human", basename(f), TRUE)) "human_hg38" else if (grepl("mouse", basename(f), TRUE)) "mouse_mm39" else NA_character_',
  '    if (is.na(sp)) next',
  '    bed <- tryCatch(readr::read_tsv(f, col_names = FALSE, show_col_types = FALSE, progress = FALSE), error = function(e) NULL)',
  '    if (is.null(bed) || !nrow(bed)) next',
  '    n <- ncol(bed)',
  '    cols <- c("chrom","start","end","name","score","strand")[seq_len(min(6, n))]',
  '    names(bed)[seq_along(cols)] <- cols',
  '    if (!"score" %in% names(bed)) bed$score <- NA_real_',
  '    # build df with base R only',
  '    df <- data.frame(',
  '      species_id = rep(sp, nrow(bed)),',
  '      chrom = as.character(bed$chrom),',
  '      start = as.integer(bed$start),',
  '      end   = as.integer(bed$end),',
  '      score = suppressWarnings(as.numeric(bed$score)),',
  '      stringsAsFactors = FALSE',
  '    )',
  '    df <- df[!is.na(df$chrom) & !is.na(df$start) & !is.na(df$end), , drop = FALSE]',
  '    df$motif_p    <- NA_real_',
  '    df$cons_class <- NA_character_',
  '    message("CTCF df cols: ", paste(names(df), collapse=","))',
  '    sql_write("ctcf_sites", df)',
  '  }',
  '} else {',
  '  message("CTCF directory missing: ", ctcf_dir)',
  '}'
)

# splice in the new block
x <- c(x[seq_len(start-1L)], new, x[seq(from=end+1L, to=length(x))])
writeLines(x, f)
