#!/usr/bin/env python3
import sqlite3, argparse, csv, sys
from pathlib import Path

def upsert_enhancers(cur, species_id, bed_path, source, method):
    with open(bed_path) as f:
        r = csv.reader(f, delimiter='\t')
        for row in r:
            if len(row) < 3: 
                continue
            chrom, start, end = row[0], int(row[1]), int(row[2])
            # Optional columns: name, score, strand… we only try to parse a numeric score if available
            score = None
            if len(row) > 4:
                try: score = float(row[4])
                except: pass
            cur.execute("""INSERT INTO enhancers(species_id,chrom,start,end,source,score,method)
                           VALUES (?,?,?,?,?,?,?)""",
                        (species_id, chrom, start, end, source, score, method))

def upsert_ctcf(cur, species_id, bed_path, method):
    with open(bed_path) as f:
        r = csv.reader(f, delimiter='\t')
        for row in r:
            if len(row) < 3: 
                continue
            chrom, start, end = row[0], int(row[1]), int(row[2])
            score = None
            if len(row) > 4:
                try: score = float(row[4])
                except: pass
            cur.execute("""INSERT INTO ctcf_sites(species_id,chrom,start,end,score,method)
                           VALUES (?,?,?,?,?,?)""", (species_id, chrom, start, end, score, method))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True)
    ap.add_argument("--species", required=True, choices=["human","mouse","pig","chicken","dog"])
    ap.add_argument("--type", required=True, choices=["enhancer","ctcf"])
    ap.add_argument("--bed", required=True)
    ap.add_argument("--source", default="")
    ap.add_argument("--method", default="")
    args = ap.parse_args()

    conn = sqlite3.connect(args.db)
    cur = conn.cursor()
    cur.execute("PRAGMA foreign_keys = ON;")
    cur.execute("SELECT species_id FROM species WHERE common_name=?", (args.species,))
    row = cur.fetchone()
    if not row:
        sys.exit(f"Unknown species {args.species}")
    species_id = row[0]

    if args.type == "enhancer":
        upsert_enhancers(cur, species_id, args.bed, args.source, args.method)
    else:
        upsert_ctcf(cur, species_id, args.bed, args.method)

    conn.commit(); conn.close()
    print("OK")

if __name__ == "__main__":
    main()
