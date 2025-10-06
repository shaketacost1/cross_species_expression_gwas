"""Export helpers for Project Alpha (figures, BED text, CSV bytes)"""
import io
import pandas as pd


def fig_to_png_bytes(fig, dpi=200):
    buf = io.BytesIO()
    try:
        fig.savefig(buf, format='png', dpi=dpi, bbox_inches='tight')
        buf.seek(0)
        return buf.read()
    finally:
        try:
            buf.close()
        except Exception:
            pass


def df_to_csv_bytes(df):
    buf = io.StringIO()
    try:
        df.to_csv(buf, index=False)
        buf.seek(0)
        return buf.read().encode('utf-8')
    finally:
        try:
            buf.close()
        except Exception:
            pass


def bed_text_from_enhancers(enh_df, r, colors_map):
    if enh_df is None or enh_df.empty or r is None:
        return ""
    header = [
        f'track name="Enhancers" description="Enhancers near region" visibility=dense itemRgb="On"',
        f"browser position {r['chrom']}:{r['start']}-{r['end']}"
    ]
    rows = []
    for _, row in enh_df.iterrows():
        chrom = row.get('chrom')
        start = int(row.get('start',0))
        end = int(row.get('end',0))
        cl = str(row.get('class','unlabeled'))
        r_,g_,b_ = colors_map.get(cl, colors_map['unlabeled'])
        name_ = f"{cl}|{row.get('enh_id','')}"
        rows.append(f"{chrom}\t{start}\t{end}\t{name_}\t0\t.\t{r_},{g_},{b_}")
    return "\n".join(header + rows)
