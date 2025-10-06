"""
Python Shiny port of app.R (best-effort, full-feature parity attempt)

Notes:
- Requires: shiny, pandas, plotnine, numpy, scipy, pillow
    Install with: pip install shiny pandas plotnine numpy scipy pillow
- This is a direct translation of the server/UI logic in `app.R` to Python Shiny.
- Some UI styling and download handlers may need small adjustments when running.

To run:
    pip install shiny pandas plotnine numpy scipy pillow
    shiny run --reload python_app.py

"""

from shiny import App, ui, render, reactive
import pandas as pd
import numpy as np
import sqlite3
import logging
from plotnine import (
    ggplot, aes, geom_rect, geom_segment, geom_label, geom_point, coord_cartesian,
    scale_x_continuous, scale_fill_manual, scale_fill_gradient, geom_tile, geom_text,
    labs, theme_minimal, theme, element_blank, element_text, facet_wrap, geom_violin,
    geom_boxplot, geom_col
)
from scipy.stats import fisher_exact
import os
import math
from utils import sql_utils, plot_utils, export_utils

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
# Primary DB path (app-local)
DB_PATH = os.path.join(BASE_DIR, "data", "regland.sqlite")
# Fallback to parent project data directory (Project Alpha 2/data/regland.sqlite)
if not os.path.exists(DB_PATH):
    alt = os.path.normpath(os.path.join(BASE_DIR, "..", "data", "regland.sqlite"))
    if os.path.exists(alt):
        DB_PATH = alt
conn = None
try:
    conn = sqlite3.connect(DB_PATH, check_same_thread=False)
except Exception:
    conn = None

# Simple logger for debug/diagnostics
logger = logging.getLogger("regland_app")
if not logger.handlers:
    h = logging.StreamHandler()
    h.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    logger.addHandler(h)
# Default logger level (set to INFO in normal runs). To capture SQL/debug logs set
# environment variable REGLAND_LOG_SQL=1 before starting the server.
logger.setLevel(logging.INFO)

# Toggle SQL parameter logging via env var (useful for debugging schema/query issues)
LOG_SQL_PARAMS = os.environ.get('REGLAND_LOG_SQL', '0') == '1'

def log_sql(q, params=None):
    """Log SQL query and parameters at DEBUG level when enabled."""
    try:
        if LOG_SQL_PARAMS:
            # Keep queries compact in log output
            qshort = (q.replace('\n', ' ').strip()[:1000] + ("..." if len(q) > 1000 else ""))
            logger.debug("SQL: %s | params: %s", qshort, params)
    except Exception:
        pass

# If SQL parameter logging requested, raise logger level so debug messages are emitted
if LOG_SQL_PARAMS:
    logger.setLevel(logging.DEBUG)

# Determine enhancers table name (some DBs use 'enhancers_all')
ENH_TBL = 'enhancers'
if conn is not None:
    ENH_TBL = sql_utils.detect_enh_table(conn)

# Helpers are provided by utils/sql_utils
sql_in = sql_utils.sql_in
has_table = lambda tbl: sql_utils.has_table(conn, tbl)
has_col = lambda tbl, col: sql_utils.has_col(conn, tbl, col)
species_col = lambda tbl: sql_utils.species_col(conn, tbl)
standardize_coords = sql_utils.standardize_coords
chrom_variants = sql_utils.chrom_variants

ucsc_map = {
    "human_hg38": "hg38",
    "mouse_mm39": "mm39",
    "macaque_rheMac10": "rheMac10",
    "chicken_galGal6": "galGal6",
    "pig_susScr11": "susScr11",
}

# Map UI labels to species_id used in the DB
SPECIES_MAP = {
    'Human': 'human_hg38',
    'Mouse': 'mouse_mm39',
    'Macaque': 'macaque_rheMac10',
    'Chicken': 'chicken_galGal6',
    'Pig': 'pig_susScr11',
}

# ---- Expression loading (robust) ----
expr_path = os.path.join(BASE_DIR, "data", "expression_tpm.tsv")
expr_tbl = None
try:
    raw = pd.read_csv(expr_path, sep="\t")
    # Try to find a symbol column
    sym_col = next((c for c in raw.columns if c.lower() in ("symbol","gene","gene_symbol","geneid","gene_id")), None)
    if sym_col is None:
        raise ValueError("No gene symbol column detected.")

    # If already long
    if set([sym_col, "tissue", "tpm"]).issubset(raw.columns):
        df = raw.rename(columns={sym_col: "symbol"})
    else:
        # wide -> long
        df = raw.rename(columns={sym_col: "symbol"}).melt(
            id_vars=["symbol"], var_name="tissue", value_name="tpm"
        )

    # Clean
    df["symbol"] = df["symbol"].astype(str)
    df["tissue"] = df["tissue"].astype(str)
    df["tpm"] = pd.to_numeric(df["tpm"], errors="coerce")

    # Optional: collapse noisy tissue names to 3 big buckets *when possible*
    t_low = df["tissue"].str.lower()
    # Broaden matching to capture common variants (ventricular, cardiomyocyte, frontal, etc.)
    # Use non-capturing groups (?:...) to avoid pandas UserWarning about regex match groups
    bucket = np.where(t_low.str.contains(r"heart|atrial|ventric(?:le|ular)?|cardiac|aorta|cardiomyocyte"),
             "Heart",
             np.where(t_low.str.contains(r"brain|cortex|hippo|amyg|putamen|nucleus|caudate|cerebral|frontal"),
             "Brain",
             np.where(t_low.str.contains(r"liver"), "Liver", None)))
    # Keep both the raw tissue and the bucket; the plot will prefer bucket if present
    df["tissue_bucket"] = bucket
    # Map any missing/None bucket values to explicit 'Other' to avoid dropping rows later
    try:
        df["tissue_bucket"] = df["tissue_bucket"].fillna('Other')
    except Exception:
        pass

    expr_tbl = df
except Exception:
    expr_tbl = None

## UI: restored full navbar layout (inputs/outputs mapped to server)
app_ui = ui.page_navbar(
    ui.nav_panel("Home", ui.h3("Welcome"), ui.p("Use Explore Genes to visualize enhancer conservation, expression and GWAS overlap.")),
    ui.nav_panel("Explore Genes",
                 ui.layout_sidebar(
                     ui.sidebar(
                         ui.input_text("gene", None, placeholder="Search gene (e.g., BDNF)", value="BDNF"),
                         ui.input_select("species","Species", {"Human":"human_hg38","Mouse":"mouse_mm39","Macaque":"macaque_rheMac10","Chicken":"chicken_galGal6","Pig":"pig_susScr11"}, selected="human_hg38"),
                         ui.input_select("tissue","Tissue", ["Liver","Brain","Heart","Other"], selected="Liver"),
                         ui.input_select("preset","Preset", ["","brain","heart","liver"], selected=""),
                         ui.output_ui("gene_suggestions"),
                         ui.input_slider("tss_kb","Distance to TSS", 0, 1000, 100, step=1, post=" kb"),
                         ui.input_checkbox_group("cls","Enhancer Class", {"conserved":"conserved","gained":"gained","lost":"lost","unlabeled":"unlabeled"}, selected=["conserved","gained","lost","unlabeled"]),
                         ui.input_slider("nbins","Bins across window", 10, 60, 30, step=1),
                         ui.input_checkbox("norm_rows","Normalize rows (0–100%)", False),
                         ui.input_checkbox("show_counts","Show bin counts", False),
                         ui.input_checkbox("mark_tss","Mark TSS", True),
                         ui.hr(),
                         ui.input_checkbox("track_stack",    "Stack lanes by class", True),
                        ui.input_checkbox("expr_log", "Log scale expression (log10 TPM+1)", False),
                        ui.input_checkbox("expr_vals", "Show bar labels", True),
                         ui.input_checkbox("track_show_gene","Show gene body", True),
                         ui.input_checkbox("track_show_snps","Show GWAS SNPs", True),
                         ui.input_action_button("apply","Apply"),
                         ui.input_action_button("boost_cov","Boost coverage (±250 kb)"),
                        ui.input_checkbox("auto_zoom", "Auto-zoom to data span", True),
                        ui.input_numeric("auto_pad_kb", "Zoom padding (kb)", 25, min=0, max=1000),
                         ui.hr(),
                         ui.input_checkbox_group("gwas_cat","GWAS categories", {"Alcohol":"Alcohol","BMI":"BMI","Inflammation":"Inflammation"}, selected=[]),
                         ui.br(), ui.output_text("applied_status")
                             , ui.output_text("debug_counts")
                        , ui.hr(),
                        ui.h5("Ask AI about this view"),
                        ui.input_checkbox("ai_include_gwas", "Include GWAS table context", True),
                        ui.output_ui("ai_prompt_box")
                     ),
                     # Main content area: genome tracks, mini browser, conservation + expression side-by-side, GWAS
                    ui.tags.div(
                        # Genome Tracks header + PNG download button
                        ui.tags.div(
                            ui.h4("Genome Tracks"),
                            ui.download_button("dl_tracks_png", "PNG", class_="btn btn-sm", style="float:right;"),
                        ),
                        ui.output_plot("track_bar"),
                        ui.br(),
                        ui.h5("Mini genome browser"),
                        ui.output_ui("mini_browser"),
                        ui.br(),
                       # two-column row: conservation (left) and expression (right)
                       ui.tags.div(
                           # left column (conservation)
                           ui.tags.div(
                               ui.tags.div(
                                   ui.tags.div(
                                       ui.h5("Conservation Matrix"),
                                       ui.download_button("dl_cons_png", "PNG", class_="btn btn-sm", style="float:right;"),
                                   ),
                                   ui.output_plot("conservation_heat"),
                                   class_="card-body"
                               ),
                               class_="card"
                           , style="flex:2;"),
                           # right column (expression)
                           ui.tags.div(
                               ui.tags.div(
                                   ui.tags.div(
                                       ui.h5("Expression (per tissue)"),
                                       # keep toggles in the sidebar but expose a download button here
                                       ui.download_button("dl_expr_csv", "CSV", class_="btn btn-sm", style="float:right;"),
                                   ),
                                   ui.output_plot("expr_bar"),
                                   class_="card-body"
                               ),
                               class_="card"
                           , style="flex:1;"),
                       style="display:flex; gap:16px; align-items:stretch;"
                       ),
                        ui.br(),
                        ui.tags.div(
                            ui.h5("GWAS hits (near gene)"),
                            ui.download_button("dl_gwas_csv", "CSV", class_="btn btn-sm", style="float:right;"),
                        ),
                        ui.output_ui("tbl_gwas")
                    )
                 )
    ),
    ui.nav_panel("Conservation Map", ui.output_plot("conservation_overview")),
    ui.nav_panel("CTCF & 3D",
                 ui.layout_sidebar(
                     ui.sidebar(
                         ui.input_radio_buttons("link_mode", "Link enhancers to gene by", {"tss":"TSS window", "ctcf":"CTCF-bounded domain"}, selected="tss"),
                         ui.input_slider("tss_kb_ctcf", "Window around TSS", 10, 1000, 250, post=" kb"),
                         ui.input_checkbox("domain_snap_tss", "Snap to domain containing the TSS", True),
                         ui.hr(),
                        ui.input_checkbox_group("ctcf_cons_groups", "CTCF conservation", {"conserved":"conserved","human_specific":"human_specific","other":"other"}, selected=["conserved","human_specific"]),
                        # Hidden proxy input used by automated badges to ensure server sees selection
                        ui.tags.div(ui.input_text("ctcf_cons_proxy", None, value=""), style="display:none;"),
                         ui.input_checkbox_group("enh_cons_groups", "Enhancer classes", {"conserved":"conserved","gained":"gained","lost":"lost","unlabeled":"unlabeled"}, selected=["conserved","gained","lost","unlabeled"]),
                         ui.input_slider("ctcf_dist_cap_kb", "Cap distance plot at", 25, 1000, 250, post=" kb"),
                         ui.hr(),
                         ui.input_checkbox_group("assoc_outcomes", "Associate conservation with", {"RNA expression":"RNA expression","GWAS hits in enhancers":"GWAS hits in enhancers","CTCF strength":"CTCF strength"}),
                         ui.input_action_button("apply_ctcf", "Analyze")
                     ),
                     ui.output_ui("ctcf_notice"),
                     ui.output_ui("ctcf_cons_in_region"),
                     ui.output_plot("ctcf_tracks"),
                     ui.output_plot("ctcf_dist_plot"),
                     ui.output_plot("enh_per_domain"),
                     ui.output_plot("assoc_expr"),
                     ui.output_ui("tbl_partition"),
                     ui.output_ui("tbl_ctcf")
                 )
    ),
    ui.nav_panel(
        "UCSC",
        ui.layout_sidebar(
            ui.sidebar(
                ui.p("Open the current region in the UCSC Genome Browser."),
                ui.input_text("ucsc_extra_tracks", "Custom track URL(s) (optional)", placeholder="Paste custom track hub / session URL(s)"),
                ui.input_checkbox("ucsc_show_gene", "Highlight TSS", True),
                ui.input_numeric("ucsc_pad_kb", "Padding (kb)", 25, min=0, max=2000),
                ui.input_action_button("open_ucsc_btn", "Build link")
            ),
            ui.output_ui("ucsc_panel")
        )
    ),
    ui.nav_panel("Expression", ui.output_plot("expr_multi")),
    ui.nav_panel("GWAS / Heritability", ui.output_plot("gwas_enrich")),
    ui.nav_panel("Downloads", ui.output_ui("tbl_downloads")),
    ui.nav_panel("About", ui.a("Berthelot et al., 2018", href="https://doi.org/10.1038/s41559-017-0377-2", target="_blank"))
)

# Server

def server(input, output, session):

    # presets
    brain_genes = ["BDNF","SCN1A","GRIN2B","DRD2","APOE"]
    heart_genes = ["TTN","MYH6","MYH7","PLN","KCNQ1"]
    liver_genes = ["ALB","APOB","CYP3A4","HNF4A","PCSK9"]

    @output()
    @render.ui
    def gene_suggestions():
        preset = input.preset()
        gg = None
        if preset == "brain":
            gg = brain_genes
        elif preset == "heart":
            gg = heart_genes
        elif preset == "liver":
            gg = liver_genes
        if conn is not None:
            ENH_TBL = sql_utils.detect_enh_table(conn)
        sp = input.species()
        if sp is None:
            return None
        # If input is already a species_id (contains underscore), return it
        if isinstance(sp, str) and '_' in sp:
            return sp
        return SPECIES_MAP.get(sp, sp)
    # Helper: normalize the UI's species selection to species_id used in DB
    def selected_species():
        sp = input.species()
        if sp is None:
            return None
        # If input is already a species_id (contains underscore), return it
        if isinstance(sp, str) and '_' in sp:
            return sp
        return SPECIES_MAP.get(sp, sp)

    # Helper to compute tight x-limits from available data (enhancers, snps, gene body)
    def compute_xlim_from_ctx(ctx, default_start, default_end, pad_kb=25):
        xs = []
        enh = ctx.get('enh') if ctx is not None else None
        snps = ctx.get('snps') if ctx is not None else None
        gbody = ctx.get('gbody') if ctx is not None else None
        try:
            if enh is not None and not enh.empty:
                xs.extend([int(enh['start'].min()), int(enh['end'].max())])
        except Exception:
            pass
        try:
            if snps is not None and not snps.empty and 'pos' in snps.columns:
                xs.extend([int(snps['pos'].min()), int(snps['pos'].max())])
        except Exception:
            pass
        try:
            if gbody is not None and not gbody.empty:
                xs.extend([int(gbody['gstart'].min()), int(gbody['gend'].max())])
        except Exception:
            pass
        if not xs:
            return default_start, default_end
        pad = int((pad_kb or 0) * 1000)
        xmin = max(0, min(xs) - pad)
        xmax = max(xs) + pad
        if xmax <= xmin:
            xmax = xmin + 1
        return xmin, xmax

    # chrom_variants and standardize_coords are provided by utils.sql_utils
    # (available via chrom_variants, standardize_coords aliases at module top)

    # --- UCSC / download helpers ---
    UCSC_COLORS = {
        "conserved":  (49, 192, 106),
        "gained":     (255, 207, 51),
        "lost":       (143, 154, 167),
        "unlabeled":  (78, 164, 255),
    }

    def _enhancers_in_view_for_ucsc():
        ctx = tracks_plot_df()
        if ctx is None:
            return pd.DataFrame(), None
        r = ctx.get('r', None)
        enh = ctx.get('enh', pd.DataFrame())
        if enh is None or enh.empty:
            return pd.DataFrame(), r
        cols = [c for c in ['chrom','start','end','class','enh_id'] if c in enh.columns]
        out = enh[cols].copy()
        out['class'] = out['class'].fillna('unlabeled')
        return out, r

    def _custom_track_text(name="Enhancers", visibility="dense"):
        enh, r = _enhancers_in_view_for_ucsc()
        if r is None or enh.empty:
            return ""
        header = [
            f'track name="{name}" description="{name} near { (input.gene() or "").upper() }" visibility={visibility} itemRgb="On"',
            f"browser position {r['chrom']}:{r['start']}-{r['end']}"
        ]
        rows = []
        for _, row in enh.iterrows():
            chrom = row['chrom']
            start = int(row['start'])
            end = int(row['end'])
            cl = str(row.get('class','unlabeled'))
            r_, g_, b_ = UCSC_COLORS.get(cl, UCSC_COLORS['unlabeled'])
            # include enh_id to make UCSC hover labels unique/useful
            name_ = f"{cl}|{row.get('enh_id','')}"
            rows.append(f"{chrom}\t{start}\t{end}\t{name_}\t0\t.\t{r_},{g_},{b_}")
        return "\n".join(header + rows)

    def current_locus_and_url(pad_kb=0, mark_tss=True):
        r = region()
        if r is None:
            return None, None
        db = ucsc_map.get(selected_species(), None)
        if db is None:
            return None, None
        pad = int((pad_kb or 0) * 1000)
        start = max(0, r['start'] - pad)
        end = r['end'] + pad
        locus = f"{r['chrom']}:{start}-{end}"
        url = f"https://genome.ucsc.edu/cgi-bin/hgTracks?db={db}&position={r['chrom']}:{start}-{end}"
        if mark_tss and isinstance(r.get('tss', None), int):
            url += f"&hgt.highlight={r['chrom']}:{r['tss']-1}-{r['tss']}"
        return locus, url

    # ---- Pure figure builders (safe to call from downloads) ----
    def build_track_bar_figure():
        # Delegate to utils.plot_utils for clarity and reuse
        ctx = tracks_plot_df()
        return plot_utils.build_track_bar_figure(ctx, input, selected_species, compute_xlim_from_ctx)


    def build_conservation_heat_figure():
        ctx = heat_plot_df()
        return plot_utils.build_conservation_heat_figure(ctx)

    # Mini genome browser UI (inline UCSC link)
    @output()
    @render.ui
    def mini_browser():
        locus, url = current_locus_and_url(pad_kb=0, mark_tss=True)
        if locus is None:
            return ui.tags.div("Locus: (no region)")
        return ui.tags.div(
            ui.tags.div(f"Locus: {locus}", style="margin-bottom:6px;"),
            ui.tags.a("Open in UCSC Genome Browser", href=url, target="_blank", class_="btn btn-outline-secondary btn-sm", style="width:100%;"),
        )

    @reactive.effect
    @reactive.event(input.open_ucsc_btn)
    def _noop_ucsc_build():
        # trigger recompute on button click
        pass

    @output()
    @render.ui
    def ucsc_panel():
        locus, base_url = current_locus_and_url(pad_kb=(input.ucsc_pad_kb() or 0), mark_tss=bool(input.ucsc_show_gene()))
        if base_url is None:
            return ui.tags.div("Select a gene and click Apply.")
        extra = (input.ucsc_extra_tracks() or "").strip()
        url = base_url
        if extra:
            url += f"&hubUrl={extra}"
        ct_text = _custom_track_text(name="RegLand Enhancers", visibility="dense")
        # Build a proper UI using Shiny-for-Python tags instead of raw HTML strings
        form = None
        if ct_text:
            form = ui.tags.form(
                ui.tags.input_(type="hidden", name="db", value=ucsc_map.get(selected_species(), "")),
                ui.tags.input_(type="hidden", name="position", value=locus),
                ui.tags.textarea(ct_text, name="hgt.customText", style="display:none"),
                id="ucscPostForm",
                action="https://genome.ucsc.edu/cgi-bin/hgTracks",
                method="POST",
                target="_blank",
            )

        # Compose the panel with proper tag objects and a real Shiny download button
        buttons = [
            ui.tags.a("Open in UCSC (base)", href=url, target="_blank", class_="btn btn-outline-secondary"),
        ]
        if form is not None:
            buttons.append(form)
            buttons.append(ui.tags.button("Open with enhancers (custom track)", class_="btn btn-primary", onclick="document.getElementById('ucscPostForm').submit()"))
        buttons.append(ui.download_button("download_ucsc_bed", "Download BED", class_="btn btn-outline-secondary"))

        return ui.div(
            ui.tags.div(ui.tags.b("Locus:"), f" {locus}", style="margin-bottom:8px"),
            ui.tags.div(*buttons, style="display:flex; gap:8px; flex-wrap:wrap; align-items:center"),
            ui.tags.div(ui.tags.small("Tip: Paste a public hub or session URL above to load additional tracks."), style="margin-top:6px"),
        )

    # Export helpers
    def _build_tracks_figure():
        fig = track_bar()
        return fig

    def _build_cons_heat_figure():
        fig = conservation_heat()
        return fig

    def _expr_export_df():
        if region() is None or expr_tbl is None or expr_tbl.empty:
            return pd.DataFrame()
        gene_sym = (input.gene() or "BDNF").upper()
        d = expr_tbl[expr_tbl['symbol'].str.upper() == gene_sym].copy()
        if d.empty:
            return pd.DataFrame()
        if 'tissue_bucket' in d.columns and d['tissue_bucket'].notna().any():
            dd = (d.dropna(subset=['tissue_bucket']).groupby('tissue_bucket', as_index=False)['tpm'].mean().rename(columns={'tissue_bucket':'tissue'}))
            cats = [c for c in ['Brain','Heart','Liver'] if c in dd['tissue'].unique()]
            if cats:
                dd['tissue'] = pd.Categorical(dd['tissue'], categories=cats, ordered=True)
                dd = dd.sort_values('tissue')
        else:
            dd = d.groupby('tissue', as_index=False)['tpm'].mean().sort_values('tpm', ascending=False).head(10)
        dd.insert(0, 'gene', gene_sym)
        return dd

    def _gwas_export_df():
        r = region()
        if r is None or not (has_table('gwas_snps') and has_table('snp_to_enhancer') and has_table('gene_to_enhancer')):
            return pd.DataFrame()

        q = f"""
            SELECT DISTINCT s.rsid, s.trait, s.pval, s.category
            FROM gwas_snps s
            JOIN snp_to_enhancer se ON se.snp_id = s.snp_id
            JOIN {ENH_TBL} e ON e.enh_id = se.enh_id
            JOIN gene_to_enhancer ge ON ge.enh_id = e.enh_id
            WHERE ge.gene_id = ? AND REPLACE(e.chrom,'chr','') = REPLACE(?,'chr','')
              AND e.start < ? AND e.end > ?
            ORDER BY COALESCE(s.pval, 1e99) ASC, s.rsid ASC
        """
        try:
            chrom_param = r['chrom']
            log_sql(q, params=(r['gene_id'], chrom_param, r['end'], r['start']))
            return pd.read_sql_query(q, conn, params=(r['gene_id'], chrom_param, r['end'], r['start']))
        except Exception:
            return pd.DataFrame()

    def _stamp(prefix, ext):
        g = (input.gene() or "GENE").upper()
        sp = selected_species() or "species"
        return f"{prefix}_{g}_{sp}.{ext}"

    @output()
    @render.download(filename=lambda: _stamp("tracks", "png"))
    def dl_tracks_png():
        fig = build_track_bar_figure()
        data = export_utils.fig_to_png_bytes(fig, dpi=200)
        yield data

    @output()
    @render.download(filename=lambda: _stamp("conservation", "png"))
    def dl_cons_png():
        fig = build_conservation_heat_figure()
        data = export_utils.fig_to_png_bytes(fig, dpi=200)
        yield data

    @output()
    @render.download(filename=lambda: _stamp("expression", "csv"))
    def dl_expr_csv():
        dd = _expr_export_df()
        yield export_utils.df_to_csv_bytes(dd)

    @output()
    @render.download(filename=lambda: _stamp("gwas_hits", "csv"))
    def dl_gwas_csv():
        dat = _gwas_export_df()
        if dat.empty:
            dat = pd.DataFrame(columns=["rsid","trait","pval","category"])
        yield export_utils.df_to_csv_bytes(dat)

    @output()
    @render.download(filename=lambda: f"enhancers_{(input.gene() or 'GENE').upper()}_{selected_species() or 'species'}.bed")
    def download_ucsc_bed():
        enh, r = _enhancers_in_view_for_ucsc()
        txt = export_utils.bed_text_from_enhancers(enh, r, UCSC_COLORS)
        if not txt:
            txt = "# no enhancers in view\n"
        yield txt.encode('utf-8')

    # --- Client-side-only AI prompt: assemble text for user to copy into ChatGPT ---
    import html

    def _assemble_ai_prompt(ctx, expr_df, gwas_df):
        r = ctx.get('r', {}) if ctx else {}
        g = (input.gene() or 'NA').upper()
        sp = selected_species() or 'NA'
        enh = ctx.get('enh') if ctx else None
        enh_summary = (enh['class'].value_counts().to_dict() if (enh is not None and not enh.empty) else {})
        expr_summary = (expr_df[['tissue','tpm']].to_dict(orient='records') if expr_df is not None and not expr_df.empty else [])
        gwas_summary = (gwas_df.head(25).to_dict(orient='records') if gwas_df is not None and not gwas_df.empty else [])
        return f"""You are analyzing a regulatory-locus snapshot.

Gene: {g}
Species: {sp}
Locus: {r.get('chrom')}:{r.get('start')}-{r.get('end')} (TSS={r.get('tss')})

Enhancer classes (counts): {enh_summary}
Expression summary (TPM means): {expr_summary}
GWAS hits (top 25 by p): {gwas_summary}

Please write a concise interpretation: patterns in enhancer conservation vs. expression,
notable GWAS traits, and suggested next checks.
"""

    @output()
    @render.ui
    def ai_prompt_box():
        try:
            ctx = tracks_plot_df()
            expr_df = _expr_export_df()
            gwas_df = _gwas_export_df() if input.ai_include_gwas() else None
            prompt = _assemble_ai_prompt(ctx, expr_df, gwas_df)
        except Exception:
            prompt = "(No region selected)"
        prompt_html = html.escape(prompt)
        return ui.HTML(f"""
          <div style="margin-top:6px">
            <textarea id="aiPromptText" style="width:100%;height:180px;white-space:pre-wrap;">{prompt_html}</textarea>
            <div style="display:flex;gap:8px;margin-top:6px;flex-wrap:wrap">
              <button class="btn btn-sm btn-outline-secondary"
                onclick="
                  const el = document.getElementById('aiPromptText');
                  if (navigator.clipboard && window.isSecureContext) {{
                    navigator.clipboard.writeText(el.value);
                  }} else {{
                    el.focus(); el.select(); document.execCommand('copy');
                  }}
                ">
                Copy prompt
              </button>
              <a class="btn btn-sm btn-primary" target="_blank" rel="noopener" href="https://chat.openai.com/">
                Open ChatGPT
              </a>
            </div>
            <small>Copy the text, open ChatGPT, paste, and run. (Clipboard API works on HTTPS or localhost.)</small>
          </div>
        """)

    # region reactive: triggered by apply button
    @reactive.Calc
    def region():
        # depend on apply count
        _ = input.apply()
        sym = (input.gene() or "BDNF").upper()
        species = selected_species()
        try:
            q = "SELECT gene_id, symbol, chrom, start, end FROM genes WHERE UPPER(symbol)=? AND species_id=? LIMIT 1"
            log_sql(q, params=(sym, species))
            df = pd.read_sql_query(q, conn, params=(sym, species))
            if df.empty:
                return None
            tss = int(df.loc[0, 'start'])
            half = int((input.tss_kb() or 100) * 1000)
            return {
                'gene_id': int(df.loc[0, 'gene_id']),
                'chrom': df.loc[0, 'chrom'],
                'start': max(0, tss - half),
                'end': tss + half,
                'tss': tss
            }
        except Exception:
            return None

    @output()
    @render.text
    def applied_status():
        r = region()
        if r is None:
            return f"Waiting… species={input.species()} | ±{input.tss_kb()}kb"
        return f"Applied: gene={ (input.gene() or 'BDNF').upper() } | species={input.species()} | tissue={input.tissue()} | ±{input.tss_kb()}kb"

    @output()
    @render.text
    def debug_counts():
        ctx = tracks_plot_df()
        if ctx is None:
            return "region: None"
        enh = ctx.get('enh', None)
        snps = ctx.get('snps', None)
        gbody = ctx.get('gbody', None)
        n_enh = int(len(enh)) if enh is not None else 0
        n_snps = int(len(snps)) if snps is not None else 0
        n_g = int(len(gbody)) if gbody is not None else 0
        return f"enhancers={n_enh} | snps={n_snps} | gene_body={n_g}"

    # Tracks plot
    def tracks_plot_df():
        r = region()
        if r is None:
            return None

        spc = species_col(ENH_TBL) or 'species_id'
        cls_filter_sql, cls_params = sql_in("COALESCE(ec.class,'unlabeled')", input.cls() or [])
        has_tissue = has_col(ENH_TBL, 'tissue')

        base_sql = f"""
            SELECT e.enh_id, e.chrom, e.start, e.end, COALESCE(ec.class,'unlabeled') AS class
            FROM {ENH_TBL} e LEFT JOIN enhancer_class ec USING(enh_id)
            WHERE e.{spc} = ?
        """
        tissue_sql = " AND (COALESCE(e.tissue,'Other') = ? OR ? = 'Other')" if has_tissue else ""
        tail_sql = " AND REPLACE(e.chrom,'chr','') = REPLACE(?,'chr','') AND e.start < ? AND e.end > ? "
        q = base_sql + tissue_sql + tail_sql + cls_filter_sql

        params = [selected_species()]
        if has_tissue:
            params += [input.tissue(), input.tissue()]
        chrom_param = r['chrom']
        params += [chrom_param, r['end'], r['start']] + cls_params

        try:
            log_sql(q, params=params)
            enh = pd.read_sql_query(q, conn, params=params)
            enh = standardize_coords(enh)
            logger.info(f"tracks_plot_df: enh rows={0 if enh is None else len(enh)} (chrom={r['chrom']}, species={selected_species()})")
        except Exception:
            enh = pd.DataFrame()
        # snps
        snps = pd.DataFrame()
        if input.track_show_snps() and has_table('gwas_snps') and has_table('snp_to_enhancer') and has_table('gene_to_enhancer'):
            try:
                snps_q = f"""
                    WITH g AS (SELECT gene_id, chrom FROM genes WHERE gene_id = ? AND species_id = ?)
                    SELECT s.rsid, s.pval, s.chrom, s.pos
                    FROM g, gwas_snps s
                    JOIN snp_to_enhancer se ON se.snp_id = s.snp_id
                    JOIN {ENH_TBL} e ON e.enh_id = se.enh_id
                    JOIN gene_to_enhancer ge ON ge.enh_id = e.enh_id
                    WHERE ge.gene_id = g.gene_id
                      AND REPLACE(s.chrom,'chr','') = REPLACE(g.chrom,'chr','')
                      AND s.pos BETWEEN ? AND ?
                """
                log_sql(snps_q, params=(r['gene_id'], selected_species(), r['start'], r['end']))
                snps = pd.read_sql_query(snps_q, conn, params=(r['gene_id'], selected_species(), r['start'], r['end']))
                if not snps.empty:
                    snps['mlog10p'] = snps['pval'].apply(lambda v: min(-math.log10(v), 30) if pd.notna(v) else np.nan)
            except Exception:
                snps = pd.DataFrame()
        # gene body
        gbody = pd.DataFrame()
        if input.track_show_gene():
            try:
                qg = "SELECT start AS gstart, end AS gend FROM genes WHERE gene_id=? LIMIT 1"
                log_sql(qg, params=(r['gene_id'],))
                gbody = pd.read_sql_query(qg, conn, params=(r['gene_id'],))
            except Exception:
                gbody = pd.DataFrame()

        return dict(r=r, enh=enh, snps=snps, gbody=gbody)

    @output()
    @render.plot
    def track_bar():
        return build_track_bar_figure()

    # heatmap
    def heat_plot_df():
        r = region()
        if r is None:
            return None
        ctx = tracks_plot_df()
        if ctx is None:
            return None
        xmin_req = r['start']; xmax_req = r['end']
        if getattr(input, 'auto_zoom', False):
            xmin, xmax = compute_xlim_from_ctx(ctx, xmin_req, xmax_req, input.auto_pad_kb())
        else:
            xmin, xmax = xmin_req, xmax_req

        nbins = input.nbins() or 30
        binw = max(1, (xmax - xmin) / nbins)
        bins = pd.DataFrame(
            {'bin': np.arange(1, nbins+1),
             'bin_start': (xmin + np.arange(0, nbins)*binw).astype(int),
             'bin_end':   (xmin + np.arange(1, nbins+1)*binw).astype(int)}
        )

        cls_sql, cls_params = sql_in("COALESCE(ec.class,'unlabeled')", input.cls() or [])
        has_tissue = has_col(ENH_TBL,'tissue')
        spc = species_col(ENH_TBL) or 'species_id'
        base_sql = f"SELECT e.start, e.end, COALESCE(ec.class,'unlabeled') AS class FROM {ENH_TBL} e LEFT JOIN enhancer_class ec USING(enh_id) WHERE e.{spc}=?"
        tissue_sql = " AND (COALESCE(e.tissue,'Other') = ? OR ?='Other')" if has_tissue else ""
        # Use normalized chrom comparison (handles 'chr' prefix variants)
        tail_sql = " AND REPLACE(e.chrom,'chr','') = REPLACE(?,'chr','') AND e.start<? AND e.end>?"
        q = base_sql + tissue_sql + tail_sql + cls_sql

        params = [selected_species()]
        if has_tissue:
            params += [input.tissue(), input.tissue()]
        ch_match = r['chrom']
        # pass region end/start as params matching e.start < end AND e.end > start
        params += [ch_match, xmax, xmin] + cls_params
        logger.info(f"heat_plot_df: region chrom={r['chrom']} xmin={xmin} xmax={xmax} nbins={nbins}")
        try:
            log_sql(q, params=params)
            df = pd.read_sql_query(q, conn, params=params)
            try:
                logger.info(f"heat_plot_df: enhancer query returned {0 if df is None else len(df)} rows (chrom={r['chrom']}, species={selected_species()})")
            except Exception:
                logger.debug("heat_plot_df: could not log df row count")
        except Exception:
            df = pd.DataFrame()
        if df.empty:
            return None
        rows = []
        for idx, b in bins.iterrows():
            sub = df[(df['end'] > b['bin_start']) & (df['start'] < b['bin_end'])]
            counts = sub.groupby('class').size().to_dict()
            for cl in ['conserved','gained','lost','unlabeled']:
                n = counts.get(cl,0)
                rows.append({'bin': b['bin'], 'class': cl, 'n': n})
        df2 = pd.DataFrame(rows)
        if input.norm_rows():
            df2['n_norm'] = df2.groupby('class')['n'].transform(lambda x: x / (x.max() if x.max()>0 else 1))
        else:
            df2['n_norm'] = df2['n']
        df2['bin_f'] = df2['bin'].astype(str)
        return dict(df2=df2, r=r, nbins=nbins, binw=binw)

    @output()
    @render.plot
    def conservation_heat():
        return build_conservation_heat_figure()

    @output()
    @render.plot
    def expr_bar():
        import matplotlib.pyplot as plt

        if region() is None:
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, "Select a gene or click Apply.", ha='center')
            ax.axis('off')
            return fig

        if expr_tbl is None or expr_tbl.empty:
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, "Expression file not found.", ha='center')
            ax.axis('off')
            return fig

        gene_sym = (input.gene() or 'BDNF').upper()
        d = expr_tbl[expr_tbl['symbol'].str.upper() == gene_sym].copy()

        if d.empty:
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, f"No expression rows for {gene_sym}.", ha='center')
            ax.axis('off')
            return fig

        # Prefer canonical tissue buckets if present; otherwise use raw tissues
        if 'tissue_bucket' in d.columns and d['tissue_bucket'].notna().any():
            dd = (
                d.dropna(subset=['tissue_bucket'])
                 .groupby('tissue_bucket', as_index=False)['tpm'].mean()
                 .rename(columns={'tissue_bucket':'tissue'})
            )
            cats = [c for c in ['Brain','Heart','Liver'] if c in dd['tissue'].unique()]
            if cats:
                dd['tissue'] = pd.Categorical(dd['tissue'], categories=cats, ordered=True)
                dd = dd.sort_values('tissue')
        else:
            dd = d.groupby('tissue', as_index=False)['tpm'].mean().sort_values('tpm', ascending=False).head(10)

        # Optional log scale
        do_log = bool(getattr(input, 'expr_log', False))
        if do_log:
            dd['tpm'] = np.log10(dd['tpm'].clip(lower=0) + 1)
            ylab = 'log10(TPM + 1)'
        else:
            ylab = 'TPM'

        # Build chart
        p = ggplot(dd, aes(x='tissue', y='tpm')) + geom_col()
        if bool(getattr(input, 'expr_vals', False)):
            if do_log:
                dd['label'] = dd['tpm'].round(3)
            else:
                dd['label'] = dd['tpm'].map(lambda v: float(f"{v:.3g}"))
            # Only show labels for bars with measurable height to avoid clutter
            try:
                thresh = 0.02 if not do_log else 0.01
                dd['show_label'] = dd['tpm'] > thresh
                p = p + geom_text(aes(label='label'), va='bottom', data=dd[dd['show_label']])
            except Exception:
                p = p + geom_text(aes(label='label'), va='bottom')

        p = p + labs(x=None, y=ylab, subtitle=gene_sym) + theme_minimal()

        # Always return a Matplotlib Figure
        try:
            fig = p.draw()
            import matplotlib.figure
            if isinstance(fig, matplotlib.figure.Figure):
                return fig
        except Exception:
            pass

        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, 'Error rendering expression plot', ha='center')
        ax.axis('off')
        return fig

    # GWAS table as HTML
    @output()
    @render.ui
    def tbl_gwas():
        r = region()
        if r is None:
            return ui.tags.div("No region selected")
        if not (has_table('gwas_snps') and has_table('snp_to_enhancer') and has_table('gene_to_enhancer')):
            return ui.tags.div("GWAS data not available")

        q = f"""
            SELECT DISTINCT s.rsid, s.trait, s.pval, s.category
            FROM gwas_snps s
            JOIN snp_to_enhancer se ON se.snp_id = s.snp_id
            JOIN {ENH_TBL} e ON e.enh_id = se.enh_id
            JOIN gene_to_enhancer ge ON ge.enh_id = e.enh_id
            WHERE ge.gene_id = ?
              AND REPLACE(e.chrom,'chr','') = REPLACE(?,'chr','')
              AND e.start < ? AND e.end > ?
            ORDER BY COALESCE(s.pval, 1e99) ASC, s.rsid ASC
        """
        try:
            chrom_param = r['chrom']
            log_sql(q, params=(r['gene_id'], chrom_param, r['end'], r['start']))
            dat = pd.read_sql_query(q, conn, params=(r['gene_id'], chrom_param, r['end'], r['start']))
        except Exception:
            dat = pd.DataFrame()
        if dat.empty:
            return ui.tags.div("No GWAS hits")
        cats = input.gwas_cat() or []
        if cats:
            if 'category' in dat.columns and dat['category'].notna().any():
                dat = dat[dat['category'].isin(cats)]
            else:
                keymap = {
                    'Alcohol': '(alcohol|alcoholism|drinking)',
                    'BMI': '(\\bBMI\\b|body mass index|obesity|adiposity)',
                    'Inflammation': '(\\bCRP\\b|C-reactive|inflamm|\\bIL[- ]?6\\b|interleukin|\\bTNF\\b)'
                }
                rx = '|'.join([keymap[c] for c in cats if c in keymap])
                dat = dat[dat['trait'].str.contains(rx, case=False, na=False)]
        return ui.HTML(dat.to_html(index=False, classes='table table-sm'))

    # placeholders for other tabs
    @output()
    @render.plot
    def conservation_overview():
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.text(0.5,0.5,"Conservation overview (todo)", ha='center')
        ax.axis('off')
        return fig

    @output()
    @render.plot
    def expr_multi():
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.text(0.5,0.5,"Expression by tissue/species (todo)", ha='center')
        ax.axis('off')
        return fig

    @output()
    @render.plot
    def gwas_enrich():
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.text(0.5,0.5,"GWAS enrichment (todo)", ha='center')
        ax.axis('off')
        return fig

    @output()
    @render.ui
    def tbl_downloads():
                # Build a simple downloads panel that triggers the existing download
                # buttons (they're present elsewhere in the DOM). Using onclick to
                # invoke those elements keeps the behaviour portable and avoids
                # needing to know internal Shiny download URLs.
                html = """
                <div>
                    <ul>
                        <li><a href="#" onclick="document.getElementById('dl_tracks_png')?.click(); return false;">Download tracks (PNG)</a></li>
                        <li><a href="#" onclick="document.getElementById('dl_cons_png')?.click(); return false;">Download conservation heat (PNG)</a></li>
                        <li><a href="#" onclick="document.getElementById('dl_expr_csv')?.click(); return false;">Download expression (CSV)</a></li>
                        <li><a href="#" onclick="document.getElementById('dl_gwas_csv')?.click(); return false;">Download GWAS hits (CSV)</a></li>
                        <li><a href="#" onclick="document.getElementById('download_ucsc_bed')?.click(); return false;">Download UCSC BED (BED9)</a></li>
                    </ul>
                </div>
                """
                return ui.HTML(html)

    # CTCF & 3D: simplified versions
    def domain_region():
        _ = input.apply_ctcf()
        r = region()
        if r is None:
            return None
        if input.link_mode() == 'tss':
            half = int((input.tss_kb_ctcf() or 250) * 1000)
            return {**r, 'start': max(0, r['tss'] - half), 'end': r['tss'] + half}
        else:
            # attempt to fetch tad containing tss
            try:
                spc = species_col('tad_domains') or 'species_id'
                tad_q = (
                    f"SELECT * FROM tad_domains WHERE {spc}=? AND REPLACE(chrom,'chr','') = REPLACE(?,'chr','') "
                    f"AND start<? AND end>? ORDER BY (end-start) ASC LIMIT 1"
                )
                chrom_param = r['chrom']
                log_sql(tad_q, params=(selected_species(), chrom_param, r['tss'], r['tss']))
                tad = pd.read_sql_query(tad_q, conn, params=(selected_species(), chrom_param, r['tss'], r['tss']))
                tad = standardize_coords(tad)
                logger.info(f"domain_region: tad rows={0 if tad is None else len(tad)} for chrom={r['chrom']} species={selected_species()}")
                if input.domain_snap_tss() and tad is not None and not tad.empty:
                    return {**r, 'start': int(tad.loc[0,'start']), 'end': int(tad.loc[0,'end'])}
            except Exception:
                pass
            return r

    def ctcf_data_df():
        r = domain_region()
        if r is None or not has_table('ctcf_sites'):
            return pd.DataFrame()
        # reset fallback flag for this call
        session._ctcf_fallback_used = False
        # Prefer explicit checkbox-group input; fall back to proxy text input which
        # automated badges set when checkbox bindings don't propagate in time.
        cons = input.ctcf_cons_groups() or ([] if not input.ctcf_cons_proxy() else [input.ctcf_cons_proxy()])
        cons_sql, cons_params = sql_in("COALESCE(cons_class,'other')", cons) if cons else ("", [])
        spc = species_col('ctcf_sites') or 'species_id'
        chrom_param = r['chrom']
        q = f"SELECT site_id, chrom, start, end, score, COALESCE(cons_class,'other') AS cons_class FROM ctcf_sites WHERE {spc} = ? AND REPLACE(chrom,'chr','') = REPLACE(?,'chr','') AND start < ? AND end > ? " + cons_sql
        params = [selected_species(), chrom_param, r['end'], r['start']] + cons_params
        try:
            # Primary query (respects selected conservation groups)
            log_sql(q, params=params)
            ct = pd.read_sql_query(q, conn, params=params)
            ct = standardize_coords(ct)
            # normalize conservation labels
            if 'cons_class' not in ct.columns:
                ct['cons_class'] = 'other'
            ct['cons_class'] = (
                ct['cons_class']
                  .fillna('other')
                  .replace({'human-specific': 'human_specific',
                            'humanSpecific': 'human_specific',
                            'Conserved': 'conserved'})
            )
            nrows = 0 if ct is None else len(ct)
            logger.info(f"ctcf_data_df: query returned {nrows} rows (chrom={r['chrom']}, species={selected_species()})")

            # Store present CTCF classes/counts in session for UI displays
            try:
                counts = {}
                if ct is not None and not ct.empty and 'cons_class' in ct.columns:
                    tmp = ct['cons_class'].fillna('other').astype(str).apply(lambda s: s.replace('-', '_').replace('humanSpecific', 'human_specific').replace('Conserved', 'conserved'))
                    counts = tmp.value_counts().to_dict()
                session._ctcf_present = counts
            except Exception:
                session._ctcf_present = {}

            # Smarter fallback: if primary returned zero rows AND the user selected
            # conservation groups, first check what CTCF cons_class values exist in
            # the region. Only retry without the cons_class filter when CTCF sites
            # exist but NONE of the selected groups are present. This avoids
            # overriding user intent when selected groups are actually present
            # but the primary query returned zero for another reason.
            if nrows == 0 and cons:
                try:
                    distinct_q = (
                        f"SELECT DISTINCT COALESCE(cons_class,'other') AS cons_class "
                        f"FROM ctcf_sites WHERE {spc} = ? AND REPLACE(chrom,'chr','') = REPLACE(?,'chr','') "
                        f"AND start < ? AND end > ?"
                    )
                    params_dist = (selected_species(), chrom_param, r['end'], r['start'])
                    log_sql(distinct_q, params=params_dist)
                    df_dist = pd.read_sql_query(distinct_q, conn, params=params_dist)
                    present = set()
                    if not df_dist.empty and 'cons_class' in df_dist.columns:
                        # normalize labels to match UI choices
                        present = set(df_dist['cons_class'].fillna('other').astype(str).apply(
                            lambda s: s.replace('-', '_').replace('humanSpecific', 'human_specific').replace('Conserved', 'conserved')
                        ).tolist())

                    # If there are CTCF rows but none of the user-selected groups are present,
                    # then retry without the cons_class filter to show available CTCF data.
                    if present and not (set(cons) & present):
                        logger.info("ctcf_data_df: selected cons groups absent in-region; retrying ctcf query without cons_class filter")
                        # mark that we widened the query for UI feedback
                        try:
                            session._ctcf_fallback_used = True
                        except Exception:
                            pass
                        q2 = (
                            f"SELECT site_id, chrom, start, end, score, COALESCE(cons_class,'other') AS cons_class "
                            f"FROM ctcf_sites WHERE {spc} = ? AND REPLACE(chrom,'chr','') = REPLACE(?,'chr','') "
                            f"AND start < ? AND end > ?"
                        )
                        params2 = (selected_species(), chrom_param, r['end'], r['start'])
                        log_sql(q2, params=params2)
                        ct2 = pd.read_sql_query(q2, conn, params=params2)
                        ct2 = standardize_coords(ct2)
                        if 'cons_class' not in ct2.columns:
                            ct2['cons_class'] = 'other'
                        ct2['cons_class'] = ct2['cons_class'].fillna('other')
                        # update present counts for fallback result
                        try:
                            tmp2 = ct2['cons_class'].fillna('other').astype(str).apply(lambda s: s.replace('-', '_').replace('humanSpecific', 'human_specific').replace('Conserved', 'conserved'))
                            session._ctcf_present = tmp2.value_counts().to_dict()
                        except Exception:
                            session._ctcf_present = {}
                        logger.info(f"ctcf_data_df: fallback query returned {0 if ct2 is None else len(ct2)} rows (chrom={r['chrom']}, species={selected_species()})")
                        return ct2
                    else:
                        logger.debug(f"ctcf_data_df: distinct cons_class present in region: {present}; no fallback needed")
                except Exception:
                    # if anything in this heuristic fails, silently continue and
                    # return the original (possibly empty) primary result
                    logger.debug("ctcf_data_df: fallback heuristic failed; returning primary result")

            return ct
        except Exception:
            return pd.DataFrame()

    def enh_in_region_df():
        r = domain_region()
        if r is None:
            return pd.DataFrame()
        cls = input.enh_cons_groups() or ['conserved','gained','lost','unlabeled']
        cls_sql, cls_params = sql_in("COALESCE(ec.class,'unlabeled')", cls)
        spc = species_col(ENH_TBL) or 'species_id'
        chrom_param = r['chrom']
        q = f"SELECT e.enh_id, e.chrom, e.start, e.end, COALESCE(ec.class,'unlabeled') AS class FROM {ENH_TBL} e LEFT JOIN enhancer_class ec USING(enh_id) WHERE e.{spc} = ? AND REPLACE(e.chrom,'chr','') = REPLACE(?,'chr','') AND e.start<? AND e.end>? " + cls_sql
        params = [selected_species(), chrom_param, r['end'], r['start']] + cls_params
        try:
            log_sql(q, params=params)
            eh = pd.read_sql_query(q, conn, params=params)
            eh = standardize_coords(eh)
            logger.info(f"enh_in_region_df: query returned {0 if eh is None else len(eh)} enh rows (chrom={r['chrom']}, species={selected_species()})")
            return eh
        except Exception:
            return pd.DataFrame()

    @output()
    @render.plot
    def ctcf_tracks():
        r = domain_region()
        if r is None:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(); ax.text(0.5,0.5,'No region'); ax.axis('off'); return fig
        enh = enh_in_region_df(); ctcf = ctcf_data_df()
        xmin = r['start']; xmax = r['end']
        pal_enh = {"conserved":"#31c06a","gained":"#ffcf33","lost":"#8f9aa7","unlabeled":"#4ea4ff"}
        pal_ctcf = {"conserved":"#1f77b4","human_specific":"#d62728","other":"#7f7f7f"}
        p = ggplot()
        if not enh.empty:
            enh['xmin'] = enh['start'].clip(lower=xmin)
            enh['xmax'] = enh['end'].clip(upper=xmax)
            p = p + geom_rect(enh, aes(xmin='xmin', xmax='xmax', ymin=0.64, ymax=0.80, fill='class'), color='#26303a', alpha=0.96)
        if not ctcf.empty:
            ctcf['mid'] = (ctcf['start'] + ctcf['end'])/2
            # plotnine uses 'size' for line width (not 'linewidth')
            p = p + geom_segment(ctcf, aes(x='mid', xend='mid', y=0.44, yend=0.58, color='cons_class'), size=0.6)
        # Apply manual palettes so legends appear and match expected colors, only when data exist
        from plotnine import scale_fill_manual, scale_color_manual
        if not enh.empty:
            p = p + scale_fill_manual(values=pal_enh)
        if not ctcf.empty:
            p = p + scale_color_manual(values=pal_ctcf)
        logger.info(f"ctcf_tracks: enh={0 if enh is None else len(enh)} ctcf={0 if ctcf is None else len(ctcf)} species={selected_species()} chrom={r['chrom']}")
        p = p + geom_segment(aes(x=r['tss'], xend=r['tss'], y=0.82, yend=0.90), linetype='dashed')
        p = p + coord_cartesian(xlim=(xmin,xmax), ylim=(0.40,0.94), expand=False)
        p = p + scale_x_continuous(labels=lambda x: [f"{xi/1e6:.1f} Mb" for xi in x])
        p = p + labs(x=None, y=None) + theme_minimal()
        return p.draw()

    @output()
    @render.ui
    def ctcf_notice():
        # Show a small notice if the ctcf query had to be widened (fallback used)
        try:
            used = getattr(session, '_ctcf_fallback_used', False)
            # To make the notice decision deterministic (not dependent on delayed
            # client-side input propagation), query the DB to see if any CTCF
            # sites exist in-region. If there are CTCF sites, suppress the
            # fallback warning because the quick-select badges already show them.
            present = {}
            try:
                r = domain_region()
                if r is not None and has_table('ctcf_sites'):
                    spc = species_col('ctcf_sites') or 'species_id'
                    distinct_q = (
                        f"SELECT COALESCE(cons_class,'other') AS cons_class, COUNT(*) AS n "
                        f"FROM ctcf_sites WHERE {spc} = ? AND REPLACE(chrom,'chr','') = REPLACE(?,'chr','') "
                        f"AND start < ? AND end > ? GROUP BY COALESCE(cons_class,'other')"
                    )
                    params_dist = (selected_species(), r['chrom'], r['end'], r['start'])
                    log_sql(distinct_q, params=params_dist)
                    df_dist = pd.read_sql_query(distinct_q, conn, params=params_dist)
                    if not df_dist.empty and 'cons_class' in df_dist.columns:
                        tmp = df_dist['cons_class'].fillna('other').astype(str).apply(lambda s: s.replace('-', '_').replace('humanSpecific', 'human_specific').replace('Conserved', 'conserved'))
                        present = tmp.value_counts().to_dict()
            except Exception:
                present = {}

            # If fallback was used but present classes exist, suppress the notice.
            if used:
                if present:
                    return None
                else:
                    return ui.tags.div(ui.tags.div(f"Note: no CTCF sites matched your selected conservation groups; showing available CTCF sites in region.", class_='alert alert-warning', role='alert'))
        except Exception:
            pass
        return None

    @output()
    @render.ui
    def ctcf_cons_in_region():
        """Show distinct CTCF conservation classes present in the current domain and
        render quick-select badges so users can set `ctcf_cons_groups` to those values.
        """
        r = domain_region()
        if r is None or not has_table('ctcf_sites'):
            return None
        spc = species_col('ctcf_sites') or 'species_id'
        try:
            # Query counts per cons_class so we can display counts on badges
            counts_q = (
                f"SELECT COALESCE(cons_class,'other') AS cons_class, COUNT(*) AS n "
                f"FROM ctcf_sites WHERE {spc} = ? AND REPLACE(chrom,'chr','') = REPLACE(?,'chr','') "
                f"AND start < ? AND end > ? GROUP BY COALESCE(cons_class,'other')"
            )
            params = (selected_species(), r['chrom'], r['end'], r['start'])
            log_sql(counts_q, params=params)
            df_dist = pd.read_sql_query(counts_q, conn, params=params)
            if df_dist.empty or 'cons_class' not in df_dist.columns:
                return ui.tags.div(ui.tags.small("No CTCF conservation classes in region."))
            # Normalize labels and prepare (value, count) pairs
            df_dist['cons_class'] = df_dist['cons_class'].fillna('other').astype(str).apply(lambda s: s.replace('-', '_').replace('humanSpecific', 'human_specific').replace('Conserved', 'conserved'))
            vals = [(row['cons_class'], int(row['n'])) for _, row in df_dist.iterrows()]
            # Render badges/buttons that set the checkbox group to the selected value(s)
            buttons = []
            for v, n in vals:
                # Also set ctcf_cons_proxy (hidden text input) so automated flows can
                # propagate a simple value that the server will read if checkbox
                # binding hasn't updated yet.
                js = (
                    f"Shiny.setInputValue('ctcf_cons_groups', ['{v}'], {{priority:'event'}}); "
                    f"Shiny.setInputValue('ctcf_cons_proxy', '{v}', {{priority:'event'}}); "
                    f"Shiny.setInputValue('apply_ctcf', Math.random(), {{priority:'event'}});"
                )
                label = f"{v} ({n})"
                buttons.append(
                    ui.tags.button(label, class_='btn btn-outline-secondary btn-sm', style='margin-right:6px;margin-bottom:4px', onclick=js)
                )
            # Also offer a button to select all present values
            all_js_vals = ",".join([f"'{v}'" for v, _ in vals])
            # use first value as proxy placeholder when selecting all
            proxy_val = vals[0][0] if vals else ''
            all_js = (
                f"Shiny.setInputValue('ctcf_cons_groups', [{all_js_vals}], {{priority:'event'}}); "
                f"Shiny.setInputValue('ctcf_cons_proxy', '{proxy_val}', {{priority:'event'}}); "
                f"Shiny.setInputValue('apply_ctcf', Math.random(), {{priority:'event'}});"
            )
            buttons.append(ui.tags.button('Select all', class_='btn btn-primary btn-sm', style='margin-left:8px', onclick=all_js))
            return ui.tags.div(ui.tags.small('CTCF classes in region:'), ui.tags.div(*buttons, style='margin-top:6px'))
        except Exception:
            return None

    @output()
    @render.plot
    def ctcf_dist_plot():
        r = domain_region()
        enh = enh_in_region_df(); ctcf = ctcf_data_df()
        if enh.empty or ctcf.empty:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(); ax.text(0.5,0.5,'Need enhancers and CTCF sites in view.'); ax.axis('off'); return fig
        enh['mid'] = (enh['start'].clip(lower=r['start']) + enh['end'].clip(upper=r['end'])) / 2
        ctcf['mid'] = (ctcf['start'] + ctcf['end']) / 2
        df = enh.merge(ctcf, on='chrom', suffixes=('.x','.y'))
        df['absdist'] = (df['mid.x'] - df['mid.y']).abs()
        # guard against NaNs and empty groups
        df = df.dropna(subset=['absdist'])
        if df.empty or df['absdist'].isna().all():
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(); ax.text(0.5,0.5,'Need enhancers and CTCF sites in view.'); ax.axis('off'); return fig

        # safer nearest-CTCF per enhancer: take row with min absdist per enh_id
        try:
            nn_idx = df.groupby('enh_id')['absdist'].idxmin()
            nn = df.loc[nn_idx.dropna()].copy()
        except Exception:
            nn = pd.DataFrame()
        if nn.empty:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(); ax.text(0.5,0.5,'Need enhancers and CTCF sites in view.'); ax.axis('off'); return fig
        nn = nn.assign(dist_kb=lambda d: d['absdist']/1000)
        nn['dist_kb'] = nn['dist_kb'].clip(upper=input.ctcf_dist_cap_kb())
        p = ggplot(nn, aes(x='class', y='dist_kb', fill='class')) + geom_violin(trim=True)
        p = p + geom_boxplot(width=0.15)
        p = p + facet_wrap('~cons_class', nrow=1) + labs(x=None, y='Nearest CTCF distance (kb, capped)') + theme_minimal()
        return p.draw()

    @output()
    @render.plot
    def enh_per_domain():
        enh = enh_in_region_df()
        if enh.empty:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(); ax.text(0.5,0.5,'No enhancers'); ax.axis('off'); return fig
        dd = enh.groupby('class').size().reset_index(name='n')
        pal_enh = {"conserved":"#31c06a","gained":"#ffcf33","lost":"#8f9aa7","unlabeled":"#4ea4ff"}
        p = (
            ggplot(dd, aes(x='class', y='n', fill='class'))
            + geom_col()
            + geom_text(aes(label='n'), va='bottom')
            + scale_fill_manual(values=pal_enh)
            + theme_minimal()
        )
        return p.draw()

    @output()
    @render.plot
    def assoc_expr():
        if expr_tbl is None or expr_tbl.empty:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(); ax.text(0.5,0.5,'Load expression_tpm.tsv'); ax.axis('off'); return fig
        r = domain_region(); enh = enh_in_region_df()
        cons_count = int((enh['class'] == 'conserved').sum()) if not enh.empty else 0
        gene_sym = (input.gene() or 'BDNF').upper()
        d = expr_tbl[expr_tbl['symbol'].str.upper() == gene_sym]
        if d.empty:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(); ax.text(0.5,0.5,'No expression for gene.'); ax.axis('off'); return fig
        tpm_mean = d['tpm'].mean()
        dfp = pd.DataFrame({'conserved_enh':[cons_count], 'tpm':[tpm_mean]})
        p = ggplot(dfp, aes(x='conserved_enh', y='tpm')) + geom_point(size=3) + labs(x='# conserved enhancers in region', y='Mean TPM', subtitle=gene_sym) + theme_minimal()
        return p.draw()

    @output()
    @render.ui
    def tbl_partition():
        r = domain_region()
        enh = enh_in_region_df()
        if enh.empty:
            return ui.HTML('<em>No enhancers in region</em>')
        try:
            spc = species_col(ENH_TBL) or 'species_id'
            has_snp_q = f"SELECT DISTINCT e.enh_id FROM {ENH_TBL} e JOIN snp_to_enhancer se USING(enh_id) JOIN gwas_snps s ON s.snp_id = se.snp_id WHERE e.{spc} = ? AND (e.chrom = ? OR e.chr = ?) AND e.start < ? AND e.end > ?"
            ch1, ch2 = chrom_variants(r['chrom'])
            log_sql(has_snp_q, params=(selected_species(), ch1, ch2, r['end'], r['start']))
            has_snp = pd.read_sql_query(has_snp_q, conn, params=(selected_species(), ch1, ch2, r['end'], r['start']))['enh_id'].tolist()
        except Exception:
            has_snp = []
        enh['has_gwas'] = enh['enh_id'].isin(has_snp)
        enh['bucket'] = enh['class'].apply(lambda c: 'conserved' if c=='conserved' else 'non_conserved')
        tab = enh.groupby(['bucket','has_gwas']).size().unstack(fill_value=0).reset_index()
        if True not in tab.columns: tab[True] = 0
        if False not in tab.columns: tab[False] = 0
        tab['total'] = tab[False] + tab[True]
        tab['prop'] = tab[True] / tab['total'].replace({0: np.nan})
        # compute Fisher exact test for conserved vs non_conserved
        or_val = np.nan
        p_val = np.nan
        try:
            conserved_row = tab[tab['bucket']=='conserved']
            noncons_row = tab[tab['bucket']!='conserved']
            if not conserved_row.empty and not noncons_row.empty:
                a = int(conserved_row[True].iloc[0]) if True in conserved_row.columns else 0
                b = int(conserved_row[False].iloc[0]) if False in conserved_row.columns else 0
                c = int(noncons_row[True].iloc[0]) if True in noncons_row.columns else 0
                d = int(noncons_row[False].iloc[0]) if False in noncons_row.columns else 0
                table = [[a, b], [c, d]]
                or_res, p_res = fisher_exact(table)
                or_val = float(or_res)
                p_val = float(p_res)
        except Exception:
            or_val = np.nan; p_val = np.nan

        tab_out = tab.copy()
        tab_out['odds_ratio'] = np.nan
        tab_out['p_value'] = np.nan
        mask = tab_out['bucket']=='conserved'
        if mask.any():
            tab_out.loc[mask, 'odds_ratio'] = or_val
            tab_out.loc[mask, 'p_value'] = p_val

        return ui.HTML(tab_out.to_html(index=False))

    @output()
    @render.ui
    def tbl_ctcf():
        r = domain_region(); ctcf = ctcf_data_df()
        if ctcf.empty:
            return ui.HTML('<em>No CTCF sites</em>')
        df = ctcf.copy(); df['mid'] = ((df['start'] + df['end'])/2).round().astype(int); df['width'] = df['end'] - df['start']
        df = df.sort_values('score', ascending=False)[['cons_class','score','chrom','start','end','width','mid']]
        return ui.HTML(df.to_html(index=False))


app = App(app_ui, server)
