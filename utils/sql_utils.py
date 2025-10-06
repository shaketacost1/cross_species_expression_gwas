"""SQL helper utilities for Project Alpha

This module provides small helpers to normalize SQL usage across database schemas and
to log SQL when requested via environment variable REGLAND_LOG_SQL.
"""
import os
import logging
import pandas as pd

logger = logging.getLogger("regland_app")

LOG_SQL_PARAMS = os.environ.get('REGLAND_LOG_SQL', '0') == '1'
if LOG_SQL_PARAMS:
    logger.setLevel(logging.DEBUG)

_SPECIES_COL_CACHE = {}


def log_sql(q, params=None):
    try:
        if LOG_SQL_PARAMS:
            qshort = (q.replace('\n', ' ').strip()[:1000] + ("..." if len(q) > 1000 else ""))
            logger.debug("SQL: %s | params: %s", qshort, params)
    except Exception:
        pass


def sql_in(col, values):
    if not values:
        return " AND 1=0 ", []
    ph = ",".join(["?"] * len(values))
    return f" AND {col} IN ({ph}) ", list(values)


def has_table(conn, tbl):
    try:
        cur = conn.execute("SELECT name FROM sqlite_master WHERE type='table' AND name=?", (tbl,))
        return cur.fetchone() is not None
    except Exception:
        return False


def has_col(conn, tbl, col):
    try:
        cur = conn.execute(f"PRAGMA table_info({tbl})")
        cols = [r[1] for r in cur.fetchall()]
        return col in cols
    except Exception:
        return False


def species_col(conn, tbl):
    if tbl in _SPECIES_COL_CACHE:
        return _SPECIES_COL_CACHE[tbl]
    try:
        if has_col(conn, tbl, 'species_id'):
            logger.info(f"species_col: using 'species_id' for table {tbl}")
            _SPECIES_COL_CACHE[tbl] = 'species_id'
            return 'species_id'
        if has_col(conn, tbl, 'species'):
            logger.info(f"species_col: using 'species' for table {tbl}")
            _SPECIES_COL_CACHE[tbl] = 'species'
            return 'species'
    except Exception:
        pass
    logger.info(f"species_col: no species column found for table {tbl}")
    _SPECIES_COL_CACHE[tbl] = None
    return None


def standardize_coords(df):
    if df is None or df.empty:
        return df
    cols = {c.lower(): c for c in df.columns}
    chr_col = cols.get('chrom') or cols.get('chr') or cols.get('chromosome')
    start_col = cols.get('start') or cols.get('chromstart') or cols.get('txstart')
    end_col = cols.get('end') or cols.get('chromend') or cols.get('txend')
    ren = {}
    if chr_col and chr_col != 'chrom': ren[chr_col] = 'chrom'
    if start_col and start_col != 'start': ren[start_col] = 'start'
    if end_col and end_col != 'end': ren[end_col] = 'end'
    if ren:
        try:
            df = df.rename(columns=ren)
            logger.info(f"standardize_coords: renamed columns {ren}")
        except Exception:
            pass
    return df


def chrom_variants(chrom):
    if chrom is None:
        return (None, None)
    c = str(chrom)
    if c.lower().startswith('chr'):
        return (c, c[3:])
    else:
        return (c, 'chr' + c)


def detect_enh_table(conn):
    """Return enhancers table name present in DB and log PRAGMA for key tables.
    Defaults to 'enhancers' when detection fails.
    """
    try:
        cur = conn.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='enhancers'")
        if cur.fetchone() is None:
            cur2 = conn.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='enhancers_all'")
            if cur2.fetchone() is not None:
                tbl = 'enhancers_all'
            else:
                tbl = 'enhancers'
        else:
            tbl = 'enhancers'
    except Exception:
        tbl = 'enhancers'
    # log PRAGMA table_info for debugging
    try:
        def _log_pragma(t):
            try:
                cur = conn.execute(f"PRAGMA table_info({t})")
                cols = [r[1] for r in cur.fetchall()]
                logger.info(f"PRAGMA {t}: {cols}")
            except Exception as e:
                logger.info(f"PRAGMA {t} failed: {e}")
        _log_pragma(tbl)
        _log_pragma('ctcf_sites')
        _log_pragma('tad_domains')
    except Exception:
        pass
    return tbl
