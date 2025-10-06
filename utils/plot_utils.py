"""Plotting helpers for Project Alpha.

Functions accept dataframes/context dictionaries and return matplotlib.Figure objects.
These are small wrappers around the plotnine-based builders used in the app.
"""
import math
import numpy as np
import pandas as pd
from plotnine import ggplot, aes, geom_rect, geom_segment, geom_point, coord_cartesian, scale_x_continuous, scale_fill_manual, labs, theme_minimal, geom_tile, scale_fill_gradient, geom_text, geom_col, geom_violin, geom_boxplot, facet_wrap
import matplotlib.pyplot as plt
from PIL import Image


UCSC_COLORS = {
    "conserved":  (49, 192, 106),
    "gained":     (255, 207, 51),
    "lost":       (143, 154, 167),
    "unlabeled":  (78, 164, 255),
}


def build_track_bar_figure(ctx, input, selected_species, compute_xlim_from_ctx):
    if ctx is None:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "Select a gene or click Apply.", ha='center')
        ax.axis('off')
        return fig

    r = ctx['r']; enh = ctx['enh']; snps = ctx['snps']; gbody = ctx['gbody']
    xmin_req = r['start']; xmax_req = r['end']
    if getattr(input, 'auto_zoom', False):
        xmin, xmax = compute_xlim_from_ctx(ctx, xmin_req, xmax_req, input.auto_pad_kb())
    else:
        xmin, xmax = xmin_req, xmax_req

    pal = {"conserved":"#31c06a","gained":"#ffcf33","lost":"#8f9aa7","unlabeled":"#4ea4ff"}

    if enh is not None and not enh.empty:
        enh = enh.copy()
        enh['xmin'] = enh['start'].clip(lower=xmin)
        enh['xmax'] = enh['end'].clip(upper=xmax)
        enh['class'] = enh['class'].fillna('unlabeled')
        if input.track_stack():
            ymap = {"conserved":(0.58,0.78),"gained":(0.46,0.66),"lost":(0.34,0.54),"unlabeled":(0.22,0.42)}
            enh['ymin'] = enh['class'].map(lambda c: ymap.get(c,(0.56,0.80))[0])
            enh['ymax'] = enh['class'].map(lambda c: ymap.get(c,(0.56,0.80))[1])
        else:
            enh['ymin'] = 0.56; enh['ymax']=0.80

    p = ggplot()

    if gbody is not None and not gbody.empty:
        gb = gbody.copy()
        gb['xmin'] = gb['gstart'].clip(lower=xmin); gb['xmax'] = gb['gend'].clip(upper=xmax)
        p = p + geom_rect(gb, aes(xmin='xmin', xmax='xmax', ymin=0.10, ymax=0.16), fill='#c7d0da')

    if enh is not None and not enh.empty:
        p = p + geom_rect(enh, aes(xmin='xmin', xmax='xmax', ymin='ymin', ymax='ymax', fill='class'), color='#26303a', alpha=0.96)

    p = p + geom_segment(aes(x=r['tss'], xend=r['tss'], y=0.16, yend=0.92), linetype='dashed')

    if snps is not None and not snps.empty:
        snps_plot = snps.copy()
        snps_plot['mlog10p'] = snps_plot['pval'].apply(lambda v: min(-math.log10(v), 30) if pd.notna(v) else np.nan)
        snps_plot['yend'] = snps_plot['mlog10p'].apply(lambda v: min(0.90 + (v/35.0), 0.98) if pd.notna(v) else 0.90)
        p = p + geom_segment(snps_plot, aes(x='pos', xend='pos', y=0.90, yend='yend'), color='#212529')
        p = p + geom_point(snps_plot, aes(x='pos', y='yend'), size=1.9, color='#212529')

    p = p + coord_cartesian(xlim=(xmin,xmax), ylim=(0.08,1.02), expand=False)
    p = p + scale_x_continuous(labels=lambda x: [f"{xi/1e6:.1f} Mb" for xi in x])

    if enh is not None and not enh.empty and 'class' in enh.columns and enh['class'].notna().any():
        present = [c for c in ['conserved','gained','lost','unlabeled'] if c in enh['class'].unique()]
        vals = {k: pal.get(k, '#4ea4ff') for k in present}
        if vals:
            p = p + scale_fill_manual(name='Enhancer class', values=vals)

    p = p + labs(x=None, y=None) + theme_minimal()

    try:
        fig = p.draw()
        import matplotlib.figure
        if isinstance(fig, matplotlib.figure.Figure):
            return fig
    except Exception:
        pass

    buf = None
    try:
        import io
        buf = io.BytesIO()
        p.save(buf, format='png', dpi=150)
        buf.seek(0)
        img = Image.open(buf).convert('RGBA')
        fig, ax = plt.subplots(figsize=(6.4, 4.0), dpi=150)
        ax.imshow(img); ax.axis('off')
        return fig
    finally:
        if buf is not None:
            try: buf.close()
            except Exception: pass


def build_conservation_heat_figure(ctx):
    import matplotlib.pyplot as plt
    if ctx is None:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No enhancers in window", ha='center'); ax.axis('off')
        return fig

    df2 = ctx['df2']
    p = ggplot(df2, aes(x='bin_f', y='class', fill='n_norm')) + geom_tile()
    p = p + scale_fill_gradient(low="#e9ecef", high="#51cf66") + labs(x=None, y=None) + theme_minimal()
    try:
        fig = p.draw()
        return fig
    except Exception:
        import io
        buf = io.BytesIO()
        p.save(buf, format='png', dpi=150)
        buf.seek(0)
        img = Image.open(buf).convert('RGBA')
        fig, ax = plt.subplots(figsize=(6.4, 4.0), dpi=150)
        ax.imshow(img); ax.axis('off')
        return fig
