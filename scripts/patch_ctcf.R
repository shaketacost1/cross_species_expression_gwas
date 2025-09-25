f <- "scripts/bootstrap_data.R"
x <- readLines(f)

# 1) Remove the brittle subsetting line in the CTCF block
rx <- '^\\s*df\\s*<-\\s*df\\s*\\[,\\s*c\\("species_id","chrom","start","end","score","motif_p","cons_class"\\)\\)\\]\\s*$'
x <- x[!grepl(rx, x)]

# 2) Insert a debug print AFTER the cons_class assignment (if not already there)
i <- grep('df\\$cons_class\\s*<-\\s*NA_character_', x)
if (length(i)) {
  ins <- 'message("CTCF df cols: ", paste(names(df), collapse=","))'
  need <- rep(TRUE, length(i))
  for (k in seq_along(i)) {
    if (i[k] < length(x) && x[i[k]+1] == ins) need[k] <- FALSE
  }
  if (any(need)) x <- append(x, ins, after=i[need])
}

writeLines(x, f)
