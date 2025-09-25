suppressPackageStartupMessages({ library(DBI); library(RSQLite); library(readr) })

app_dir  <- normalizePath("~/Desktop/Project Alpha", mustWork = FALSE)
src_root <- normalizePath("~/my_project", mustWork = FALSE)
db_path  <- file.path(app_dir, "data", "regland.sqlite")
con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
on.exit(try(DBI::dbDisconnect(con), silent=TRUE), add=TRUE)

ctcf_dir <- file.path(src_root, "data", "capstone data")
stopifnot(dir.exists(ctcf_dir))
files <- list.files(ctcf_dir, pattern="CTCF\\.(bed|bed\\.gz)$", full.names=TRUE)
if (!length(files)) stop("No CTCF .bed files in: ", ctcf_dir)

f <- files[1]
message("Testing file: ", f)

sp <- if (grepl("human", basename(f), TRUE)) "human_hg38" else if (grepl("mouse", basename(f), TRUE)) "mouse_mm39" else NA_character_
stopifnot(!is.na(sp))

bed <- tryCatch(readr::read_tsv(f, col_names=FALSE, show_col_types=FALSE, progress=FALSE), error=function(e) NULL)
stopifnot(!is.null(bed), nrow(bed) > 0)

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
df <- df[!is.na(df$chrom) & !is.na(df$start) & !is.na(df$end), , drop=FALSE]

df$motif_p    <- NA_real_
df$cons_class <- NA_character_

message("Built df with columns: ", paste(names(df), collapse=","))
str(df, list.len=5)

# Try a direct write to a temp table
tbl <- "ctcf_sites_debug"
try(DBI::dbExecute(con, paste0("DROP TABLE IF EXISTS ", tbl)), silent=TRUE)
DBI::dbExecute(con, sprintf('CREATE TABLE %s(
  site_id INTEGER PRIMARY KEY,
  species_id TEXT, chrom TEXT, start INTEGER, end INTEGER,
  score REAL, motif_p REAL, cons_class TEXT
)', tbl))

DBI::dbWriteTable(con, tbl, df, append=TRUE)
cat("Write OK. Row count in debug table: ",
    DBI::dbGetQuery(con, sprintf("SELECT COUNT(*) n FROM %s", tbl))$n, "\n", sep="")
