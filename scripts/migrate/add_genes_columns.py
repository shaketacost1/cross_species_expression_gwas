#!/usr/bin/env python3
import sqlite3

DB = "ipa.db"
missing_cols = {
    "strand": "TEXT",
    "gene_name": "TEXT",
    "biotype": "TEXT",
    "source": "TEXT"
}

con = sqlite3.connect(DB)
cur = con.cursor()
cur.execute("PRAGMA table_info(genes);")
present = {row[1] for row in cur.fetchall()}

to_add = [(c, t) for c, t in missing_cols.items() if c not in present]
if not to_add:
    print("✅ No migration needed — all columns already exist.")
else:
    for col, typ in to_add:
        cur.execute(f"ALTER TABLE genes ADD COLUMN {col} {typ};")
        print(f"✅ Added column {col} {typ}")
    con.commit()
    print("✅ Migration complete — genes table updated.")
con.close()
