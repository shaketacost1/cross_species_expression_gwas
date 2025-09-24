# app.R — genome track + conservation UX + expression + CTCF/3D + hardening

library(shiny)
library(bslib)
library(DT)
library(DBI)
library(RSQLite)
library(ggplot2)
library(dplyr)
library(tidyr)
library(thematic)
library(shinyjs)
library(scales)
library(readr)
library(stringr)
library(rlang)

thematic::thematic_on()
`%||%` <- function(a, b) if (!is.null(a)) a else b
kb_to_bp <- function(kb) kb * 1000

db_path <- "data/regland.sqlite"
con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
onStop(function() try(DBI::dbDisconnect(con), silent = TRUE))

has_table <- function(tbl) DBI::dbExistsTable(con, tbl)
has_col   <- function(tbl, col) (has_table(tbl) && (col %in% DBI::dbListFields(con, tbl)))
has_tissue_col <- function() has_col("enhancers", "tissue")

# ---- schema & indexes (non-destructive) ----
ensure_schema <- function(con){
  dbExecute(con,"CREATE TABLE IF NOT EXISTS genes(gene_id INTEGER PRIMARY KEY,symbol TEXT,species_id TEXT,chrom TEXT,start INTEGER,end INTEGER);")
  dbExecute(con,"CREATE TABLE IF NOT EXISTS enhancers(enh_id INTEGER PRIMARY KEY,species_id TEXT,chrom TEXT,start INTEGER,end INTEGER,tissue TEXT,score REAL,source TEXT);")
  dbExecute(con,"CREATE TABLE IF NOT EXISTS enhancer_class(enh_id INTEGER,class TEXT);")
  dbExecute(con,"CREATE TABLE IF NOT EXISTS gwas_snps(snp_id INTEGER PRIMARY KEY,rsid TEXT UNIQUE,chrom TEXT,pos INTEGER,trait TEXT,pval REAL,source TEXT,category TEXT);")
  dbExecute(con,"CREATE TABLE IF NOT EXISTS snp_to_enhancer(snp_id INTEGER,enh_id INTEGER,UNIQUE(snp_id,enh_id));")
  dbExecute(con,"CREATE TABLE IF NOT EXISTS gene_to_enhancer(gene_id INTEGER,enh_id INTEGER,UNIQUE(gene_id,enh_id));")
  dbExecute(con,"CREATE TABLE IF NOT EXISTS ctcf_sites(site_id INTEGER PRIMARY KEY,species_id TEXT,chrom TEXT,start INTEGER,end INTEGER,score REAL,motif_p REAL,cons_class TEXT);")
  dbExecute(con,"CREATE TABLE IF NOT EXISTS tad_domains(tad_id INTEGER PRIMARY KEY,species_id TEXT,chrom TEXT,start INTEGER,end INTEGER,source TEXT);")
  invisible(TRUE)
}
ensure_schema(con)

ensure_indexes <- function(con){
  for (q in c(
    "CREATE INDEX IF NOT EXISTS idx_genes_species_chrom_start ON genes(species_id,chrom,start);",
    "CREATE INDEX IF NOT EXISTS idx_enh_species_chrom_bounds  ON enhancers(species_id,chrom,start,end);",
    "CREATE INDEX IF NOT EXISTS idx_enh_tissue             ON enhancers(tissue);",
    "CREATE INDEX IF NOT EXISTS idx_ec_enh                 ON enhancer_class(enh_id);",
    "CREATE INDEX IF NOT EXISTS idx_ge_gene                ON gene_to_enhancer(gene_id);",
    "CREATE INDEX IF NOT EXISTS idx_ge_enh                 ON gene_to_enhancer(enh_id);",
    "CREATE INDEX IF NOT EXISTS idx_se_enh                 ON snp_to_enhancer(enh_id);",
    "CREATE INDEX IF NOT EXISTS idx_snp_chrom_pos          ON gwas_snps(chrom,pos);",
    "CREATE INDEX IF NOT EXISTS idx_snp_rsid               ON gwas_snps(rsid);",
    "CREATE INDEX IF NOT EXISTS idx_ctcf_species_chrom_bounds ON ctcf_sites(species_id,chrom,start,end);",
    "CREATE INDEX IF NOT EXISTS idx_ctcf_cons_class           ON ctcf_sites(cons_class);",
    "CREATE INDEX IF NOT EXISTS idx_tad_species_chrom_bounds  ON tad_domains(species_id,chrom,start,end);"
  )) try(dbExecute(con, q), silent = TRUE)
  invisible(TRUE)
}
ensure_indexes(con)

# ---- hydrate GWAS categories (optional TSV) ----
try({
  if (file.exists("data/gwas_hg38_with_category.tsv")) {
    gc <- readr::read_tsv("data/gwas_hg38_with_category.tsv", show_col_types = FALSE, progress = FALSE)
    nm <- names(gc)
    rsid_col <- intersect(c("rsid","RSID","snp","SNP"), nm)[1]
    cat_col  <- intersect(c("category","Category","cat","CAT"), nm)[1]
    if (!is.na(rsid_col) && !is.na(cat_col)) {
      gc2 <- gc %>% transmute(rsid = .data[[rsid_col]], category = .data[[cat_col]]) %>%
        filter(!is.na(rsid), !is.na(category))
      if (nrow(gc2)) {
        dbExecute(con, "CREATE TEMP TABLE IF NOT EXISTS tmp_gwas_cat(rsid TEXT, category TEXT)")
        dbExecute(con, "DELETE FROM tmp_gwas_cat")
        dbWriteTable(con, "tmp_gwas_cat", gc2, temporary = TRUE, append = TRUE)
        dbExecute(con, "
          UPDATE gwas_snps
             SET category = (SELECT category FROM tmp_gwas_cat WHERE tmp_gwas_cat.rsid = gwas_snps.rsid)
           WHERE rsid IN (SELECT rsid FROM tmp_gwas_cat)
        ")
        dbExecute(con, "DROP TABLE tmp_gwas_cat")
      }
    }
  }
}, silent = TRUE)

# ---- helpers ----
vneed <- function(cond, msg = " ") {
  if (isTRUE(cond)) return(invisible())
  validate(need(FALSE, as.character(if (length(msg) == 0 || is.na(msg)) " " else msg)[1]))
}
sql_in <- function(col, values){
  if (length(values)==0) return(list(sql=" AND 1=0 ", params=list()))
  ph <- paste(rep("?", length(values)), collapse=",")
  list(sql = sprintf(" AND %s IN (%s) ", col, ph), params = as.list(values))
}
ucsc_db_for <- function(species_id){
  c(
    human_hg38       = "hg38",
    mouse_mm39       = "mm39",
    macaque_rheMac10 = "rheMac10",
    chicken_galGal6  = "galGal6",
    pig_susScr11     = "susScr11"
  )[species_id] %||% NA_character_
}

# ---- expression data (force Brain/Heart/Liver only) ----
expr_path <- "data/expression_tpm.tsv"
expr_tbl <- tryCatch({
  raw <- readr::read_tsv(expr_path, show_col_types = FALSE, progress = FALSE)
  if (!all(c("symbol","tissue","tpm") %in% names(raw))) {
    sym_col <- intersect(c("symbol","gene","Gene","SYMBOL"), names(raw))[1]
    stopifnot(!is.na(sym_col))
    raw <- raw %>%
      rename(symbol = all_of(sym_col)) %>%
      pivot_longer(-symbol, names_to = "tissue", values_to = "tpm")
  }
  raw %>%
    mutate(tissue = case_when(
      str_detect(tissue, regex("brain|cortex|cereb|hippo|amyg|putamen|nucleus_acc|caudate", TRUE)) ~ "Brain",
      str_detect(tissue, regex("heart|atrial|ventricle|cardiac|aorta", TRUE))                       ~ "Heart",
      str_detect(tissue, regex("liver", TRUE))                                                      ~ "Liver",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(tissue)) %>%
    group_by(symbol, tissue) %>%
    summarize(tpm = mean(as.numeric(tpm), na.rm = TRUE), .groups="drop") %>%
    mutate(
      symbol = as.character(symbol),
      tissue = factor(tissue, levels = c("Brain","Heart","Liver"))
    )
}, error = function(e) NULL)

# ---- UI ----
reg_theme <- bs_theme(version=5, bootswatch="flatly",
                      base_font = font_google("Inter"),
                      heading_font = font_google("Inter"),
                      "enable-rounded"=TRUE, "enable-shadows"=TRUE)

ui <- page_navbar(
  title = tags$span("Regulatory Landscapes"), theme = reg_theme,
  navbar_options = navbar_options(collapsible = TRUE), id = "topnav",

  nav_panel("Home", layout_columns(card(card_header("Welcome"),
                                        p("Use Explore Genes to visualize enhancer conservation, expression and GWAS overlap.")))),

  nav_panel(
    "Explore Genes",
    layout_sidebar(
      sidebar = sidebar(
        position = "left", sticky = TRUE, width = 320,
        useShinyjs(),
        textInput("gene", NULL, placeholder = "Search gene (e.g., BDNF)", value = "BDNF"),
        selectInput(
          "species","Species",
          c("Human"="human_hg38",
            "Mouse"="mouse_mm39",
            "Macaque"="macaque_rheMac10",
            "Chicken"="chicken_galGal6",
            "Pig"="pig_susScr11"
          ),
          selected="human_hg38"
        ),
        selectInput("tissue","Tissue", c("Liver","Brain","Heart","Other"), selected="Liver"),
        selectInput("preset","Preset",
          c("None"="",
            "Brain (BDNF, SCN1A, GRIN2B, DRD2, APOE)"="brain",
            "Heart (TTN, MYH6, MYH7, PLN, KCNQ1)"="heart",
            "Liver (ALB, APOB, CYP3A4, HNF4A, PCSK9)"="liver"),
          selected=""),
        uiOutput("gene_suggestions"),

        sliderInput("tss_kb","Distance to TSS", min=0,max=1000,value=100,post=" kb"),
        checkboxGroupInput("cls","Enhancer Class",
                           choices=c(conserved="conserved",gained="gained",lost="lost",unlabeled="unlabeled"),
                           selected=c("conserved","gained","lost","unlabeled")),
        sliderInput("nbins","Bins across window", min=10, max=60, value=30, step=1),

        checkboxInput("norm_rows","Normalize rows (0–100%)", FALSE),
        checkboxInput("show_counts","Show bin counts", FALSE),
        checkboxInput("mark_tss","Mark TSS", TRUE),

        tags$hr(),
        checkboxInput("track_stack",    "Stack lanes by class", TRUE),
        checkboxInput("track_show_gene","Show gene body", TRUE),
        checkboxInput("track_show_snps","Show GWAS SNPs", TRUE),

        actionButton("apply","Apply", class="btn-primary btn-lg w-100"),
        actionButton("boost_cov","Boost coverage (±250 kb)", class="btn-outline-secondary btn-sm w-100 mt-2"),

        tags$hr(),
        checkboxGroupInput("gwas_cat","GWAS categories",
                           choices=c("Alcohol","BMI","Inflammation"),
                           selected=character(0)),
        tags$small(
          "Clear: ",
          actionLink("clr_cats","Inflammation"), " ",
          actionLink("clr_bmi","BMI"), " ",
          actionLink("clr_alc","Alcohol")
        ),

        br(), tags$small(textOutput("applied_status"))
      ),
      layout_columns(
        col_widths = c(12),
        card(card_header(div(class="d-flex justify-content-between align-items-center",
                             div("Genome Tracks"),
                             downloadButton("dl_tracks","PNG", class="btn-secondary btn-sm"))),
             plotOutput("track_bar", height = 260)),
        card(card_header("Mini genome browser"), uiOutput("ucsc_widget")),
        layout_columns(
          col_widths=c(6,6),
          card(card_header(div(class="d-flex justify-content-between align-items-center",
                               div("Conservation Matrix"),
                               downloadButton("dl_heat","PNG", class="btn-secondary btn-sm"))),
               plotOutput("conservation_heat", height = 360)),
          card(card_header(div(class="d-flex justify-content-between align-items-center",
                               div("Expression (per tissue)"),
                               tagList(
                                 checkboxInput("expr_log","log10(TPM+1)", FALSE, width="auto"),
                                 checkboxInput("expr_vals","Show values", TRUE, width="auto"),
                                 downloadButton("dl_expr","CSV", class="btn-secondary btn-sm")
                               ))),
               plotOutput("expr_bar", height = 280))
        ),
        card(card_header(div(class="d-flex justify-content-between align-items-center",
                             div("GWAS hits (near gene)"),
                             downloadButton("dl_gwas","CSV", class="btn-secondary btn-sm"))),
             DTOutput("tbl_gwas"))
      )
    )
  ),

  nav_panel("Conservation Map", layout_columns(card(card_header("Across-species heatmap"),
                                                    plotOutput("conservation_overview", height=420)))),

  nav_panel(
    "CTCF & 3D",
    layout_sidebar(
      sidebar = sidebar(
        width = 320, sticky = TRUE,
        radioButtons("link_mode", "Link enhancers to gene by",
                     c("TSS window" = "tss",
                       "CTCF-bounded domain" = "ctcf"),
                     selected = "tss"),
        conditionalPanel(
          "input.link_mode == 'tss'",
          sliderInput("tss_kb_ctcf", "Window around TSS", min=10, max=1000, value=250, post=" kb")
        ),
        conditionalPanel(
          "input.link_mode == 'ctcf'",
          checkboxInput("domain_snap_tss", "Snap to domain containing the TSS", TRUE),
          helpText("If off, the Explore Genes window is used; any TADs overlapping it are considered.")
        ),
        tags$hr(),
        checkboxGroupInput("ctcf_cons_groups", "CTCF conservation",
                           choices=c("conserved","human_specific","other"),
                           selected=c("conserved","human_specific")),
        checkboxGroupInput("enh_cons_groups", "Enhancer classes",
                           choices=c("conserved","gained","lost","unlabeled"),
                           selected=c("conserved","gained","lost","unlabeled")),
        sliderInput("ctcf_dist_cap_kb", "Cap distance plot at", min=25, max=1000, value=250, post=" kb"),
        tags$hr(),
        checkboxGroupInput("assoc_outcomes", "Associate conservation with",
                           choices=c("RNA expression","GWAS hits in enhancers","CTCF strength")),
        actionButton("apply_ctcf", "Analyze", class="btn-primary w-100")
      ),
      layout_columns(
        col_widths = c(12),
        card(
          card_header(div(class="d-flex justify-content-between align-items-center",
                          div("3D domain view"),
                          downloadButton("dl_ctcf_tracks","PNG", class="btn-secondary btn-sm"))),
          plotOutput("ctcf_tracks", height = 280)
        ),
        layout_columns(
          col_widths = c(6,6),
          card(card_header("Distances to nearest CTCF"), plotOutput("ctcf_dist_plot", height=320)),
          card(card_header("Enhancers per region (by class)"), plotOutput("enh_per_domain", height=320))
        ),
        layout_columns(
          col_widths = c(6,6),
          card(card_header("Conserved enhancers vs RNA"), plotOutput("assoc_expr", height=320)),
          card(card_header("GWAS over-representation"), DTOutput("tbl_partition"))
        ),
        card(card_header("Top CTCF sites near gene (ranked by score)"), DTOutput("tbl_ctcf"))
      )
    )
  ),

  nav_panel("Expression", layout_columns(card(card_header("Gene expression by tissue/species"),
                                              plotOutput("expr_multi", height=420)))),
  nav_panel("GWAS / Heritability", layout_columns(card(card_header("Trait overlap & enrichment"),
                                                       plotOutput("gwas_enrich", height=420)))),
  nav_panel("Downloads", layout_columns(card(card_header("Results"), DTOutput("tbl_downloads")))),
  nav_menu("About", nav_item(a(href="https://doi.org/10.1038/s41559-017-0377-2", target="_blank","Berthelot et al., 2018")))
)

# ---- server ----
server <- function(input, output, session){

  # presets & suggestions
  brain_genes <- c("BDNF","SCN1A","GRIN2B","DRD2","APOE")
  heart_genes <- c("TTN","MYH6","MYH7","PLN","KCNQ1")
  liver_genes <- c("ALB","APOB","CYP3A4","HNF4A","PCSK9")

  observeEvent(input$preset, ignoreInit = TRUE, {
    if (identical(input$preset,"brain")) {
      updateSelectInput(session, "tissue", selected = "Brain")
      updateTextInput(session, "gene", value = brain_genes[[1]])
    } else if (identical(input$preset,"heart")) {
      updateSelectInput(session, "tissue", selected = "Heart")
      updateTextInput(session, "gene", value = heart_genes[[1]])
    } else if (identical(input$preset,"liver")) {
      updateSelectInput(session, "tissue", selected = "Liver")
      updateTextInput(session, "gene", value = liver_genes[[1]])
    }
    shinyjs::click("apply")
  })

  output$gene_suggestions <- renderUI({
    gg <- switch(input$preset, "brain"=brain_genes, "heart"=heart_genes, "liver"=liver_genes, NULL)
    if (is.null(gg)) return(NULL)
    tags$div(
      tags$small("Suggestions: "),
      lapply(seq_along(gg), function(i) {
        gid <- paste0("sugg_", i)
        actionLink(gid, gg[[i]], class = "me-2")
      })
    )
  })
  observe({
    gg <- switch(input$preset, "brain"=brain_genes, "heart"=heart_genes, "liver"=liver_genes, NULL)
    if (is.null(gg)) return()
    lapply(seq_along(gg), function(i){
      gid <- paste0("sugg_", i)
      observeEvent(input[[gid]], {
        updateTextInput(session, "gene", value = gg[[i]])
        shinyjs::click("apply")
      }, ignoreInit = TRUE)
    })
  })

  observeEvent(input$boost_cov, { updateSliderInput(session, "tss_kb", value = 250); shinyjs::click("apply") })
  observeEvent(input$clr_cats, { updateCheckboxGroupInput(session,"gwas_cat",selected=setdiff(input$gwas_cat,"Inflammation")) })
  observeEvent(input$clr_bmi,  { updateCheckboxGroupInput(session,"gwas_cat",selected=setdiff(input$gwas_cat,"BMI")) })
  observeEvent(input$clr_alc,  { updateCheckboxGroupInput(session,"gwas_cat",selected=setdiff(input$gwas_cat,"Alcohol")) })

  gene_typed <- debounce(reactive(input$gene), 400)
  observeEvent(gene_typed(), ignoreInit=TRUE, { if (nzchar(gene_typed())) shinyjs::click("apply") })
  observeEvent(input$species, ignoreInit=TRUE, { shinyjs::click("apply") })
  observeEvent(input$tss_kb,  ignoreInit=TRUE, { shinyjs::click("apply") })
  observeEvent(input$cls,     ignoreInit=TRUE, { shinyjs::click("apply") })
  observeEvent(input$nbins,   ignoreInit=TRUE, { shinyjs::click("apply") })
  observeEvent(input$norm_rows, ignoreInit=TRUE, { shinyjs::click("apply") })
  observeEvent(input$show_counts, ignoreInit=TRUE, { shinyjs::click("apply") })
  observeEvent(input$mark_tss, ignoreInit=TRUE, { shinyjs::click("apply") })
  observeEvent(input$track_stack, ignoreInit=TRUE, { shinyjs::click("apply") })
  observeEvent(input$track_show_gene, ignoreInit=TRUE, { shinyjs::click("apply") })
  observeEvent(input$track_show_snps, ignoreInit=TRUE, { shinyjs::click("apply") })

  region <- eventReactive(input$apply, {
    sym <- toupper(if (nzchar(input$gene)) input$gene else "BDNF")
    g <- tryCatch(DBI::dbGetQuery(con,
                                  "SELECT gene_id,symbol,chrom,start,end
                                     FROM genes
                                    WHERE UPPER(symbol)=? AND species_id=?
                                    LIMIT 1",
                                  params=list(sym,input$species)),
                  error=function(e) data.frame())
    if (!nrow(g)) return(NULL)
    tss  <- as.numeric(g$start[1])
    half <- kb_to_bp(input$tss_kb)
    data.frame(gene_id=g$gene_id[1], chrom=g$chrom[1],
               start=max(0, tss-half), end=tss+half, tss=tss)
  }, ignoreInit = TRUE)

  output$applied_status <- renderText({
    r <- isolate(region())
    if (is.null(r)) sprintf("Waiting… species=%s | ±%skb", input$species, input$tss_kb)
    else sprintf("Applied: gene=%s | species=%s | tissue=%s | ±%skb",
                 toupper(input$gene %||% "BDNF"), input$species, input$tissue, input$tss_kb)
  })

  # ---- GENOME TRACKS ----
  tracks_plot <- reactive({
    req(region())
    r <- region()
    cls_filter <- sql_in("COALESCE(ec.class,'unlabeled')", input$cls)

    base_sql <- "
      SELECT e.enh_id, e.chrom, e.start, e.end,
             COALESCE(ec.class,'unlabeled') AS class
        FROM enhancers e
        LEFT JOIN enhancer_class ec USING(enh_id)
       WHERE e.species_id = ?"
    # treat NULL tissue as 'Other'
    tissue_sql <- if (has_tissue_col()) " AND (COALESCE(e.tissue,'Other') = ? OR ? = 'Other')" else ""
    tail_sql   <- " AND e.chrom = ? AND e.start < ? AND e.end > ? "
    q <- paste0(base_sql, tissue_sql, tail_sql, cls_filter$sql)

    params <- c(list(input$species),
                if (has_tissue_col()) list(input$tissue, input$tissue) else NULL,
                list(r$chrom, r$end, r$start),
                cls_filter$params)

    enh <- tryCatch(DBI::dbGetQuery(con, q, params=params), error=function(e) data.frame())

    if (isTRUE(input$track_stack)) {
      enh$class[is.na(enh$class)] <- "unlabeled"
      enh$class <- factor(enh$class, levels=c("conserved","gained","lost","unlabeled"))
      ymap <- data.frame(
        class = factor(c("conserved","gained","lost","unlabeled"), levels=c("conserved","gained","lost","unlabeled")),
        ymin  = c(0.58, 0.46, 0.34, 0.22),
        ymax  = c(0.78, 0.66, 0.54, 0.42)
      )
      enh <- merge(enh, ymap, by="class", all.x=TRUE)
    } else {
      enh$ymin <- 0.56
      enh$ymax <- 0.80
      enh$class[is.na(enh$class)] <- "unlabeled"
      enh$class <- factor(enh$class, levels=c("conserved","gained","lost","unlabeled"))
    }

    snps <- NULL
    if (isTRUE(input$track_show_snps) &&
        has_table("gwas_snps") && has_table("snp_to_enhancer") && has_table("gene_to_enhancer")) {
      snps <- tryCatch(DBI::dbGetQuery(
        con, "
          WITH g AS (
            SELECT gene_id, chrom FROM genes
             WHERE gene_id = ? AND species_id = ?
          )
          SELECT s.rsid, s.pval, s.chrom, s.pos
            FROM g, gwas_snps s
            JOIN snp_to_enhancer se ON se.snp_id = s.snp_id
            JOIN enhancers e        ON e.enh_id = se.enh_id
            JOIN gene_to_enhancer ge ON ge.enh_id = e.enh_id
           WHERE ge.gene_id = g.gene_id
             AND s.chrom = g.chrom
             AND s.pos BETWEEN ? AND ?
        ",
        params = list(r$gene_id, input$species, r$start, r$end)
      ), error=function(e) NULL)
      if (!is.null(snps) && nrow(snps)) {
        snps$mlog10p <- ifelse(is.na(snps$pval), NA_real_, pmin(-log10(snps$pval), 30))
      }
    }

    gbody <- NULL
    if (isTRUE(input$track_show_gene)) {
      gbody <- tryCatch(DBI::dbGetQuery(
        con,"SELECT start AS gstart, end AS gend FROM genes WHERE gene_id=? LIMIT 1",
        params=list(r$gene_id)), error=function(e) NULL)
    }

    xmin <- r$start; xmax <- r$end
    pal <- c(conserved="#31c06a", gained="#ffcf33", lost="#8f9aa7", unlabeled="#4ea4ff")
    min_pad <- max(50, (xmax - xmin) * 0.0002)

    step_bp <- max(2.5e6, (xmax - xmin) / 16)
    bands <- data.frame(
      xmin = seq(floor(xmin/step_bp)*step_bp, xmax, by = 2*step_bp),
      xmax = seq(floor(xmin/step_bp)*step_bp, xmax, by = 2*step_bp) + step_bp
    )
    bands$xmax <- pmin(bands$xmax, xmax)

    if (nrow(enh)) {
      enh$xmin <- pmax(enh$start, xmin)
      enh$xmax <- pmin(enh$end,   xmax)
      w <- enh$xmax - enh$xmin
      too_thin <- w < (2*min_pad)
      if (any(too_thin)) {
        mid <- (enh$xmin[too_thin] + enh$xmax[too_thin]) / 2
        enh$xmin[too_thin] <- pmax(xmin, mid - min_pad)
        enh$xmax[too_thin] <- pmin(xmax, mid + min_pad)
      }
    }

    ggplot() +
      { if (nrow(bands))
          geom_rect(data=bands,
                    aes(xmin=xmin, xmax=xmax, ymin=0.20, ymax=0.86),
                    inherit.aes=FALSE, fill="#f5f7fa") } +
      annotate("rect", xmin=xmin, xmax=xmax, ymin=0.86, ymax=0.90, fill="#e1e6ec") +
      { if (!is.null(gbody) && nrow(gbody) > 0)
          geom_rect(data=gbody, aes(xmin=pmax(gstart,xmin), xmax=pmin(gend,xmax),
                                    ymin=0.10, ymax=0.16),
                    inherit.aes = FALSE, fill="#c7d0da") } +
      { if (nrow(enh))
          geom_rect(data=enh,
                    aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=class),
                    color="#26303a", linewidth=0.15, alpha=0.96, inherit.aes=FALSE) } +
      geom_segment(aes(x=r$tss, xend=r$tss, y=0.16, yend=0.92),
                   linetype="dashed", linewidth=0.6, color="#e8590c") +
      geom_label(data = data.frame(x = r$tss, y = 0.93, label = "TSS"),
                 aes(x=x, y=y, label=label),
                 vjust = 0, size = 3, label.size = NA, fill = NA, color = "#e8590c") +
      { if (!is.null(snps) && nrow(snps))
          geom_segment(data=snps,
                       aes(x=pos, xend=pos, y=0.90, yend=pmin(0.90 + (mlog10p/35), 0.98)),
                       color="#212529", linewidth=0.5) } +
      { if (!is.null(snps) && nrow(snps))
          geom_point(data=snps, aes(x=pos, y=pmin(0.90 + (mlog10p/35), 0.98)),
                     size=1.9, color="#212529") } +
      coord_cartesian(xlim=c(xmin,xmax), ylim=c(0.08,1.02), expand=0) +
      scale_x_continuous(labels=function(x) number(x/1e6, accuracy=0.1, suffix=" Mb")) +
      scale_fill_manual(values=pal, limits=c("conserved","gained","lost","unlabeled"),
                        drop=FALSE, name=NULL) +
      labs(x=NULL,y=NULL) +
      theme_minimal(base_size=12) +
      theme(axis.text.y=element_blank(), axis.ticks=element_blank(),
            panel.grid=element_blank(), legend.position="bottom",
            plot.margin = margin(4,8,4,8))
  })

  output$track_bar <- renderPlot({
    vneed(!is.null(region()), " ")
    print(tracks_plot())
  }, res = 96)

  # ---- UCSC ----
  output$ucsc_widget <- renderUI({
    if (is.null(region())) return(div(class="text-muted","Select a gene or click Apply."))
    r <- region(); db <- ucsc_db_for(input$species)
    loc <- sprintf("%s:%d-%d", r$chrom, r$start, r$end)
    if (is.na(db)) return(div(class="text-muted","UCSC link not available for this species."))
    href <- sprintf("https://genome.ucsc.edu/cgi-bin/hgTracks?db=%s&position=%s",
                    db, URLencode(loc, reserved = TRUE))
    tagList(
      tags$div(tags$strong("Locus: "), loc),
      a(href=href, target="_blank", class="btn btn-outline-secondary btn-sm mt-2",
        "Open in UCSC Genome Browser")
    )
  })

  # ---- heatmap (conservation matrix) ----
  heat_plot <- reactive({
    req(region())
    r <- region()
    nbins <- input$nbins %||% 30
    binw <- max(1,(r$end-r$start)/nbins)
    bins <- tibble(bin=1:nbins,
                   bin_start=r$start + (0:(nbins-1))*binw,
                   bin_end  =r$start + (1:nbins)*binw)

    cls <- sql_in("COALESCE(ec.class,'unlabeled')", input$cls)

    base_sql <- "
      SELECT e.start,e.end,COALESCE(ec.class,'unlabeled') AS class
        FROM enhancers e
        LEFT JOIN enhancer_class ec USING(enh_id)
       WHERE e.species_id=?"
    # treat NULL tissue as 'Other'
    tissue_sql <- if (has_tissue_col()) " AND (COALESCE(e.tissue,'Other') = ? OR ?='Other')" else ""
    tail_sql   <- " AND e.chrom=? AND e.start<? AND e.end>?"
    q <- paste0(base_sql, tissue_sql, tail_sql, cls$sql)
    params <- c(list(input$species),
                if (has_tissue_col()) list(input$tissue, input$tissue) else NULL,
                list(r$chrom, r$end, r$start),
                cls$params)
    df <- tryCatch(DBI::dbGetQuery(con, q, params=params), error=function(e) data.frame())

    if (!nrow(df)) {
      return(ggplot() + theme_void() + labs(title="No enhancers in window") +
               theme(plot.margin=margin(8,8,8,8)))
    }

    df2 <- bins |>
      rowwise() |>
      mutate(data = list(dplyr::filter(df, end > bin_start, start < bin_end))) |>
      tidyr::unnest(data) |>
      mutate(bin=factor(bin)) |>
      count(bin, class, name="n") |>
      tidyr::complete(bin, class=c("conserved","gained","lost","unlabeled"), fill=list(n=0)) |>
      group_by(class) |>
      mutate(n_norm = if (isTRUE(input$norm_rows)) ifelse(max(n)>0, n / max(n), 0) else n) |>
      ungroup()

    bin_centers <- (0:(nbins-1)) + 0.5
    pos_bp      <- r$start + bin_centers * binw
    tick_every  <- max(1, floor(nbins / 6))
    ticks_i     <- seq(1, nbins, by = tick_every)
    tick_labels <- scales::number(pos_bp[ticks_i]/1e6, accuracy = 0.1, suffix = " Mb")

    p <- ggplot(df2, aes(x=bin, y=class, fill=n_norm)) +
      geom_tile() +
      scale_fill_gradient(low="#e9ecef", high="#51cf66",
                          name = if (isTRUE(input$norm_rows)) "% of row max" else "Count") +
      labs(x=NULL,y=NULL) +
      scale_x_discrete(breaks=factor(ticks_i), labels=tick_labels, expand=c(0,0)) +
      theme_minimal(base_size=12) +
      theme(panel.grid=element_blank(), axis.text.x=element_text(size=9),
            plot.margin=margin(8,8,8,8))

    if (isTRUE(input$show_counts)) {
      p <- p + geom_text(aes(label=ifelse(n>0, n, "")), size=3, color="#334155")
    }
    if (isTRUE(input$mark_tss)) {
      tss_bin <- floor((r$tss - r$start)/binw) + 1
      if (is.finite(tss_bin) && tss_bin>=1 && tss_bin<=nbins) {
        p <- p + geom_vline(xintercept = tss_bin, linetype="dashed", color="#6c757d")
      }
    }
    p
  })

  output$conservation_heat <- renderPlot({
    vneed(!is.null(region()), " ")
    print(heat_plot())
  }, res = 96)

  # ---- expression (per tissue) ----
  output$expr_bar <- renderPlot({
    vneed(!is.null(region()), " ")
    vneed(!is.null(expr_tbl) && nrow(expr_tbl) > 0, "Expression file not found.")

    gene_sym <- toupper(input$gene %||% "BDNF")
    d <- expr_tbl |> dplyr::filter(toupper(symbol) == gene_sym)

    d <- dplyr::right_join(
      d,
      data.frame(tissue = factor(c("Brain","Heart","Liver"), levels = c("Brain","Heart","Liver"))),
      by = "tissue"
    ) |> dplyr::mutate(symbol = gene_sym, tpm = tpm %||% 0)

    vneed(nrow(d) > 0, sprintf("No expression rows for %s.", gene_sym))

    dd <- d |> arrange(tpm)
    ylab <- if (isTRUE(input$expr_log)) "log10(TPM+1)" else "TPM (GTEx v10 median)"
    if (isTRUE(input$expr_log)) dd <- dd |> mutate(tpm = log10(tpm + 1))

    g <- ggplot(dd, aes(x=tissue, y=tpm)) +
      geom_col() +
      { if (isTRUE(input$expr_vals)) geom_text(aes(label = scales::number(tpm, accuracy = 0.01)),
                                               vjust=-0.3, size=3.2) } +
      labs(x=NULL, y=ylab, subtitle=gene_sym) +
      theme_minimal(base_size=12) +
      theme(axis.text.x = element_text(angle = 20, hjust = 1),
            plot.margin = margin(8,8,8,8))
    print(g)
  }, res = 96)

  output$dl_expr <- downloadHandler(
    filename = function() paste0("expression_", toupper(input$gene %||% "GENE"), ".csv"),
    content = function(file){
      if (is.null(expr_tbl)) { write.csv(data.frame(), file, row.names=FALSE); return() }
      gene_sym <- toupper(input$gene %||% "BDNF")
      df <- expr_tbl |> filter(toupper(symbol) == gene_sym)
      if (!nrow(df)) write.csv(data.frame(), file, row.names=FALSE) else readr::write_csv(df, file)
    }
  )

  # ---- GWAS table ----
  output$tbl_gwas <- DT::renderDT({
    if (is.null(region())) {
      return(datatable(data.frame(rsid=character(), trait=character(), pval=numeric(), category=character()),
                       options=list(pageLength=5), rownames=FALSE))
    }
    r <- region()
    if (!(has_table("gwas_snps") && has_table("snp_to_enhancer") && has_table("gene_to_enhancer"))) {
      return(datatable(data.frame(rsid=character(), trait=character(), pval=numeric(), category=character()),
                       options=list(pageLength=5), rownames=FALSE))
    }

    qbase <- "
        SELECT DISTINCT s.rsid, s.trait, s.pval, s.category
          FROM gwas_snps s
          JOIN snp_to_enhancer se  ON se.snp_id = s.snp_id
          JOIN enhancers e         ON e.enh_id = se.enh_id
          JOIN gene_to_enhancer ge ON ge.enh_id = e.enh_id
         WHERE ge.gene_id = ? AND e.chrom = ? AND e.start < ? AND e.end > ?
         ORDER BY COALESCE(s.pval, 1e99) ASC, s.rsid ASC
      "

    dat <- tryCatch(DBI::dbGetQuery(con, qbase, params=list(r$gene_id, r$chrom, r$end, r$start)),
                    error=function(e) data.frame())

    if (!"category" %in% names(dat)) dat$category <- NA_character_

    cats <- input$gwas_cat %||% character(0)
    if (length(cats) && nrow(dat)) {
      if (any(!is.na(dat$category))) {
        dat <- dplyr::filter(dat, !is.na(category) & category %in% cats)
      } else {
        keymap <- list(
          Alcohol      = "(alcohol|alcoholism|drinking)",
          BMI          = "(\\bBMI\\b|body mass index|obesity|adiposity)",
          Inflammation = "(\\bCRP\\b|C-reactive|inflamm|\\bIL[- ]?6\\b|interleukin|\\bTNF\\b)"
        )
        rx <- paste(unlist(keymap[cats]), collapse="|")
        dat <- dplyr::filter(dat, grepl(rx, trait, ignore.case = TRUE))
      }
    }

    datatable(dat, options=list(pageLength=5), rownames=FALSE)
  })

  # ---- downloads for Explore Genes ----
  png_device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"

  output$dl_tracks <- downloadHandler(
    filename = function() paste0("tracks_", toupper(input$gene %||% "BDNF"), "_", input$species, "_", input$tss_kb, "kb.png"),
    content  = function(file) ggsave(file, plot = tracks_plot(), width=11, height=2.8, dpi=160, device=png_device)
  )
  output$dl_heat <- downloadHandler(
    filename = function() paste0("heatmap_", toupper(input$gene %||% "BDNF"), "_", input$species, "_", input$tss_kb, "kb.png"),
    content  = function(file) ggsave(file, plot = heat_plot(), width=7.5, height=3.6, dpi=160, device=png_device)
  )
  output$dl_gwas <- downloadHandler(
    filename = function() paste0("gwas_", toupper(input$gene %||% "BDNF"), "_", input$species, "_", input$tss_kb, "kb.csv"),
    content  = function(file) {
      if (is.null(region())) { write.csv(data.frame(), file, row.names=FALSE); return() }
      r <- region()
      base <- tryCatch(DBI::dbGetQuery(con, "
        SELECT DISTINCT s.rsid, s.trait, s.pval, s.category
          FROM gwas_snps s
          JOIN snp_to_enhancer se  ON se.snp_id = s.snp_id
          JOIN enhancers e         ON e.enh_id = se.enh_id
          JOIN gene_to_enhancer ge ON ge.enh_id = e.enh_id
         WHERE ge.gene_id = ? AND e.chrom = ? AND e.start < ? AND e.end > ?
         ORDER BY COALESCE(s.pval, 1e99) ASC, s.rsid ASC",
        params = list(r$gene_id, r$chrom, r$end, r$start)
      ), error=function(e) data.frame())
      if (!"category" %in% names(base)) base$category <- NA_character_

      cats <- input$gwas_cat %||% character(0)
      if (length(cats) && nrow(base)) {
        if (any(!is.na(base$category))) {
          base <- dplyr::filter(base, !is.na(category) & category %in% cats)
        } else {
          keymap <- list(
            Alcohol      = "(alcohol|alcoholism|drinking)",
            BMI          = "(\\bBMI\\b|body mass index|obesity|adiposity)",
            Inflammation = "(\\bCRP\\b|C-reactive|inflamm|\\bIL[- ]?6\\b|interleukin|\\bTNF\\b)"
          )
          rx <- paste(unlist(keymap[cats]), collapse="|")
          base <- dplyr::filter(base, grepl(rx, trait, ignore.case = TRUE))
        }
      }
      if (requireNamespace("readr", quietly=TRUE)) readr::write_csv(base, file) else write.csv(base, file, row.names=FALSE)
    }
  )

  # =========================
  # CTCF & 3D TAB
  # =========================

  observeEvent(list(input$link_mode, input$tss_kb_ctcf, input$domain_snap_tss,
                    input$ctcf_cons_groups, input$enh_cons_groups, input$ctcf_dist_cap_kb), {
    if (!is.null(region())) shinyjs::click("apply_ctcf")
  })

  domain_region <- eventReactive(input$apply_ctcf, {
    r <- region(); req(r)
    if (identical(input$link_mode, "tss")) {
      half <- kb_to_bp(input$tss_kb_ctcf %||% 250)
      return(transform(r, start = max(0, r$tss - half), end = r$tss + half))
    }
    tad <- tryCatch(DBI::dbGetQuery(con, "
             SELECT * FROM tad_domains
              WHERE species_id=? AND chrom=? AND start<? AND end>?
              ORDER BY (end-start) ASC LIMIT 1",
             params = list(input$species, r$chrom, r$tss, r$tss)), error=function(e) NULL)
    if (isTRUE(input$domain_snap_tss) && !is.null(tad) && nrow(tad)) {
      data.frame(gene_id=r$gene_id, chrom=r$chrom, start=tad$start[1], end=tad$end[1], tss=r$tss)
    } else {
      r
    }
  }, ignoreInit=TRUE)

  ctcf_data <- reactive({
    r <- domain_region(); req(r)
    if (!has_table("ctcf_sites")) return(data.frame())
    cons <- input$ctcf_cons_groups %||% character(0)
    cons_filter <- if (length(cons)) sql_in("COALESCE(cons_class,'other')", cons) else list(sql="",params=list())
    q <- paste0("
      SELECT site_id, chrom, start, end, score, COALESCE(cons_class,'other') AS cons_class
        FROM ctcf_sites
       WHERE species_id = ? AND chrom = ? AND start < ? AND end > ? ",
       cons_filter$sql)
    params <- c(list(input$species, r$chrom, r$end, r$start), cons_filter$params)
    tryCatch(dbGetQuery(con, q, params=params), error=function(e) data.frame())
  })

  enh_in_region <- reactive({
    r <- domain_region(); req(r)
    cls <- input$enh_cons_groups %||% c("conserved","gained","lost","unlabeled")
    cls_filter <- sql_in("COALESCE(ec.class,'unlabeled')", cls)
    q <- paste0("
       SELECT e.enh_id, e.chrom, e.start, e.end, COALESCE(ec.class,'unlabeled') AS class
         FROM enhancers e
         LEFT JOIN enhancer_class ec USING(enh_id)
        WHERE e.species_id=? AND e.chrom=? AND e.start<? AND e.end>? ", cls_filter$sql)
    params <- c(list(input$species, (domain_region())$chrom, (domain_region())$end, (domain_region())$start),
                cls_filter$params)
    tryCatch(dbGetQuery(con, q, params=params), error=function(e) data.frame())
  })

  gwas_in_enh <- reactive({
    r <- domain_region(); req(r)
    if (!(has_table("gwas_snps") && has_table("snp_to_enhancer"))) return(data.frame())
    tryCatch(dbGetQuery(con, "
      SELECT DISTINCT s.rsid, s.trait, s.pval, s.chrom, s.pos
        FROM gwas_snps s
        JOIN snp_to_enhancer se ON se.snp_id = s.snp_id
        JOIN enhancers e        ON e.enh_id = se.enh_id
       WHERE e.species_id=? AND e.chrom=? AND e.start<? AND e.end>? ",
      params=list(input$species, r$chrom, r$end, r$start)), error=function(e) data.frame())
  })

  ctcf_tracks_plot <- reactive({
    r <- domain_region(); req(r)
    enh <- enh_in_region(); ctcf <- ctcf_data()
    xmin <- r$start; xmax <- r$end

    pal_enh <- c(conserved="#31c06a", gained="#ffcf33", lost="#8f9aa7", unlabeled="#4ea4ff")
    pal_ctcf <- c(conserved="#1f77b4", human_specific="#d62728", other="#7f7f7f")
    pad <- max(50, (xmax-xmin)*0.0002)

    enh2 <- NULL
    if (nrow(enh)) {
      enh2 <- within(enh, {
        xmin <- pmax(start, xmin); xmax <- pmin(end, xmax)
        thin <- (pmin(end, xmax) - pmax(start, xmin)) < (2*pad)
        mid  <- (pmax(start, xmin) + pmin(end, xmax)) / 2
        xmin2 <- ifelse(thin, pmax(xmin, mid - pad), pmax(start, xmin))
        xmax2 <- ifelse(thin, pmin(xmax, mid + pad), pmin(end, xmax))
      })
    }

    ggplot() +
      annotate("rect", xmin=xmin, xmax=xmax, ymin=0.58, ymax=0.62, fill="#e1e6ec") +
      { if (!is.null(enh2) && nrow(enh2))
          geom_rect(data=enh2,
                    aes(xmin=xmin2, xmax=xmax2, ymin=0.64, ymax=0.80, fill=class),
                    color="#26303a", linewidth=0.15, alpha=0.96) } +
      { if (nrow(ctcf))
          geom_segment(data=ctcf,
            aes(x=(start+end)/2, xend=(start+end)/2, y=0.44, yend=0.58, color=cons_class),
            linewidth=0.6) } +
      geom_segment(aes(x=r$tss, xend=r$tss, y=0.82, yend=0.90),
                   linetype="dashed", linewidth=0.6, color="#e8590c") +
      geom_label(data = data.frame(x = r$tss, y = 0.91, label = "TSS"),
                 aes(x=x, y=y, label=label),
                 vjust = 0, size = 3, label.size = NA, fill = NA, color = "#e8590c") +
      coord_cartesian(xlim=c(xmin,xmax), ylim=c(0.40,0.94), expand=0) +
      scale_x_continuous(labels=function(x) scales::number(x/1e6, accuracy=0.1, suffix = " Mb")) +
      scale_fill_manual(values=pal_enh, limits=names(pal_enh), drop=FALSE, name="Enhancer") +
      scale_color_manual(values=pal_ctcf, limits=names(pal_ctcf), drop=FALSE, name="CTCF") +
      theme_minimal(base_size=12) +
      labs(x=NULL,y=NULL) +
      theme(axis.text.y=element_blank(), panel.grid=element_blank(),
            legend.position="bottom", plot.margin = margin(4,8,4,8))
  })

  output$ctcf_tracks <- renderPlot({
    req(domain_region())
    print(ctcf_tracks_plot())
  }, res=96)

  output$dl_ctcf_tracks <- downloadHandler(
    filename = function() paste0("ctcf_tracks_", toupper(input$gene %||% "GENE"), ".png"),
    content  = function(file) ggsave(file, plot = ctcf_tracks_plot(), width=11, height=2.8, dpi=160, device=png_device)
  )

  output$ctcf_dist_plot <- renderPlot({
    r <- domain_region(); req(r)
    enh <- enh_in_region(); ctcf <- ctcf_data()
    vneed(nrow(enh) && nrow(ctcf), "Need enhancers and CTCF sites in view.")

    enh$mid <- (pmax(enh$start, r$start) + pmin(enh$end, r$end)) / 2
    ctcf$mid <- (ctcf$start + ctcf$end) / 2

    df <- merge(enh, ctcf, by="chrom")
    df$absdist <- abs(df$mid.x - df$mid.y)
    nn <- df %>% group_by(enh_id) %>% slice_min(absdist, n = 1, with_ties = FALSE) %>%
      ungroup() %>% transmute(class, cons_class, dist_kb = pmin(absdist/1000, input$ctcf_dist_cap_kb))

    ggplot(nn, aes(x=class, y=dist_kb, fill=class)) +
      geom_violin(trim=TRUE, scale="width") +
      geom_boxplot(width=0.15, outlier.size=0.6) +
      facet_wrap(~cons_class, nrow=1) +
      labs(x=NULL, y="Nearest CTCF distance (kb, capped)") +
      theme_minimal(base_size=12) + theme(legend.position="none")
  }, res=96)

  output$enh_per_domain <- renderPlot({
    enh <- enh_in_region(); vneed(nrow(enh), "No enhancers in region.")
    dd <- enh %>% count(class)
    ggplot(dd, aes(x=class, y=n, fill=class)) +
      geom_col() +
      geom_text(aes(label=n), vjust=-0.2) +
      labs(x=NULL, y="# Enhancers in region") +
      theme_minimal(base_size=12) + theme(legend.position="none")
  }, res=96)

  output$assoc_expr <- renderPlot({
    vneed(!is.null(expr_tbl) && nrow(expr_tbl)>0, "Load expression_tpm.tsv")
    r <- domain_region(); req(r)
    enh <- enh_in_region()
    cons_count <- sum(enh$class == "conserved", na.rm=TRUE)

    gene_sym <- toupper(input$gene %||% "BDNF")
    d <- expr_tbl %>% filter(toupper(symbol)==gene_sym) %>%
      group_by(symbol) %>% summarize(tpm_mean = mean(tpm, na.rm=TRUE))
    vneed(nrow(d)>0, "No expression for gene.")
    df <- data.frame(conserved_enh = cons_count, tpm = d$tpm_mean[1])
    ggplot(df, aes(x=conserved_enh, y=tpm)) +
      geom_point(size=3) +
      labs(x="# conserved enhancers in region", y="Mean TPM", subtitle=gene_sym) +
      theme_minimal(base_size=12)
  }, res=96)

  # ---- FIXED: robust GWAS partition table ----
  output$tbl_partition <- DT::renderDT({
    r <- domain_region(); req(r)
    enh <- enh_in_region()
    if (!nrow(enh)) return(DT::datatable(data.frame(), options=list(dom='t')))

    has_snp <- tryCatch({
      dbGetQuery(con, "
        SELECT DISTINCT e.enh_id
          FROM enhancers e
          JOIN snp_to_enhancer se USING(enh_id)
          JOIN gwas_snps s       ON s.snp_id = se.snp_id
         WHERE e.species_id=? AND e.chrom=? AND e.start<? AND e.end>?",
        params=list(input$species, r$chrom, r$end, r$start))$enh_id
    }, error=function(e) integer(0))

    enh$has_gwas <- enh$enh_id %in% has_snp
    enh$bucket   <- ifelse(enh$class=="conserved", "conserved", "non_conserved")

    tab <- enh %>%
      count(bucket, has_gwas) %>%
      tidyr::complete(bucket, has_gwas = c(FALSE, TRUE), fill = list(n = 0)) %>%
      tidyr::pivot_wider(names_from = has_gwas, values_from = n, values_fill = 0)

    if (!"FALSE" %in% names(tab)) tab$`FALSE` <- 0
    if (!"TRUE"  %in% names(tab)) tab$`TRUE`  <- 0
    tab <- dplyr::rename(tab, no_gwas = `FALSE`, with_gwas = `TRUE`)
    tab$total <- tab$no_gwas + tab$with_gwas
    tab$prop  <- ifelse(tab$total>0, tab$with_gwas/tab$total, NA_real_)

    if (nrow(tab)==2) {
      m <- matrix(c(tab$with_gwas[tab$bucket=="conserved"], tab$no_gwas[tab$bucket=="conserved"],
                    tab$with_gwas[tab$bucket!="conserved"], tab$no_gwas[tab$bucket!="conserved"]), nrow=2, byrow=TRUE)
      ft <- tryCatch(fisher.test(m), error=function(e) NULL)
      or <- if (!is.null(ft)) unname(ft$estimate) else NA_real_
      p  <- if (!is.null(ft)) ft$p.value else NA_real_
    } else { or <- NA_real_; p <- NA_real_ }

    out <- tab %>% select(bucket, total, with_gwas, prop) %>%
      mutate(odds_ratio = ifelse(bucket=="conserved", or, NA_real_),
             p_value = ifelse(bucket=="conserved", p, NA_real_))
    DT::datatable(out, options=list(pageLength=5), rownames=FALSE)
  })

  output$tbl_ctcf <- DT::renderDT({
    r <- domain_region(); req(r)
    ctcf <- ctcf_data()
    if (!nrow(ctcf)) return(DT::datatable(data.frame(), options=list(dom='t')))
    df <- ctcf %>%
      mutate(mid = round((start+end)/2),
             width = end-start) %>%
      arrange(desc(score)) %>%
      select(cons_class, score, chrom, start, end, width, mid)
    DT::datatable(df, options=list(pageLength=6), rownames=FALSE)
  })

  # ---- placeholders / other tabs ----
  gg_blank <- function(tt) ggplot() + theme_void() + labs(title=tt) + theme(plot.margin=margin(8,8,8,8))
  output$conservation_overview <- renderPlot({ print(gg_blank("Conservation overview (todo)")) }, res=96)
  output$expr_multi            <- renderPlot({ print(gg_blank("Expression by tissue/species (todo)")) }, res=96)
  output$gwas_enrich           <- renderPlot({ print(gg_blank("GWAS enrichment (todo)")) }, res=96)
  output$tbl_downloads         <- renderDT({ datatable(data.frame(File=character(), Link=character())) })
}

shinyApp(ui, server)
