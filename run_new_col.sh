#!/usr/bin/env python3
import sqlite3
import argparse
import gzip
import sys

def parse_attrs(attr_str):
    """Parse the attributes column from a GTF line."""
    attrs = {}
    for entry in attr_str.strip().split(";"):
        entry = entry.strip()
        if not entry:
            continue
        if " " not in entry:
            continue
        key, val = entry.split(" ", 1)
        attrs[key] = val.strip('"')
    gene_id = attrs.get("gene_id")
    gene_name = attrs.get("gene_name")
    biotype = (
        attrs.get("gene_biotype") or
        attrs.get("gene_type") or
        attrs.get("biotype")
    )
    return gene_id, gene_name, biotype


def open_any(path):
    """Open plain or gzipped file transparently."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def main():
    ap = argparse.ArgumentParser(
        description="Load genes + biotypes from GTF into SQLite database."
    )
    ap.add_argument("--db", required=True, help="SQLite database file (ipa.db)")
    ap.add_argument("--species", required=True,
                    choices=["human", "mouse", "pig", "chicken", "dog"],
                    help="Species common name (must exist in species table)")
    ap.add_argument("--gtf", required=True, help="Input GTF (can be .gz)")
    ap.add_argument("--source", required=True, help="Annotation source, e.g., GENCODE or Ensembl")
    args = ap.parse_args()

    # Connect to database
    conn = sqlite3.connect(args.db)
    cur = conn.cursor()

    # Get species ID
    cur.execute("SELECT species_id FROM species WHERE common_name=?", (args.species,))
    row = cur.fetchone()
    if not row:
        sys.exit(f"Unknown species: {args.species}")
    species_id = row[0]

    n = 0
    seen = set()

    # Open and parse GTF
    with open_any(args.gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            chrom, src, feature, start, end, score, strand, frame, attrs = parts
            if feature != "gene":
                continue

            gid, gname, biotype = parse_attrs(attrs)
            if not gid:
                continue

            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue

            # Avoid duplicates
            if gid in seen:
                continue
            seen.add(gid)

            # Insert with explicit type coercion
            cur.execute("""
                INSERT OR REPLACE INTO genes
                  (gene_stable_id, species_id, chrom, start, end, strand, gene_name, biotype, source)
                VALUES (?,?,?,?,?,?,?,?,?)
            """, (
                str(gid),
                int(species_id),
                str(chrom),
                int(start_i),
                int(end_i),
                str(strand),
                str(gname) if gname else None,
                str(biotype) if biotype else None,
                str(args.source)
            ))

            n += 1
            if n % 100000 == 0:
                conn.commit()
                print(f"... processed {n:,} genes", flush=True)

    conn.commit()
    conn.close()
    print(f"✅ OK - inserted {n:,} genes for {args.species}")


if __name__ == "__main__":
    main()

