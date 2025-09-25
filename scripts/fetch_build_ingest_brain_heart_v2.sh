#!/usr/bin/env bash
set -euo pipefail

DB_PATH="data/regland.sqlite"
SPECIES_ID="human_hg38"
PROMOTER_BP=2000
LINK_KB=250
WORKDIR="data/tissue_build"
BRAIN_H3="$WORKDIR/brain/h3k27ac"
BRAIN_DN="$WORKDIR/brain/dnase"
HEART_H3="$WORKDIR/heart/h3k27ac"
HEART_DN="$WORKDIR/heart/dnase"
GENCODE_GTF="data/gencode.v46.annotation.gtf.gz"
BLACKLIST_GZ="data/hg38-blacklist.v2.bed.gz"
SOURCE_LABEL="ENCODE_build"

mkdir -p "$BRAIN_H3" "$BRAIN_DN" "$HEART_H3" "$HEART_DN" "$WORKDIR" data

# Ensure tools
for t in curl jq bedtools sqlite3 gzip awk sort; do command -v "$t" >/dev/null || { echo "Missing $t"; exit 1; }; done
HAS_BB2B="$(command -v bigBedToBed || true)"

# Get GENCODE & blacklist
[ -f "$GENCODE_GTF" ] || curl -L -o "$GENCODE_GTF" "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz"
[ -f "$BLACKLIST_GZ" ] || curl -L -o "$BLACKLIST_GZ" "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz"

# Helper: download ENCODE File records (BED + bigBed) using the File endpoint with accession-based @@download URLs
download_encode_files () {
  organ="$1"; assay="$2"; outdir="$3"
  mkdir -p "$outdir"
  base="https://www.encodeproject.org/search/?type=File&status=released&assembly=GRCh38&output_type=peaks&format=json&limit=all&biosample_ontology.organ_slims=${organ}"
  if [ "$assay" = "H3K27ac" ]; then
    q_bed="${base}&file_format=bed&assay_title=Histone%20ChIP-seq&target.label=H3K27ac"
    q_bb="${base}&file_format=bigBed&assay_title=Histone%20ChIP-seq&target.label=H3K27ac"
  else
    q_bed="${base}&file_format=bed&assay_title=DNase-seq"
    q_bb="${base}&file_format=bigBed&assay_title=DNase-seq"
  fi

  curl -sL "$q_bed" \
  | jq -r '.["@graph"][] | select(.file_format=="bed" and (.output_type|test("peaks";"i"))) | "https://www.encodeproject.org/files/" + .accession + "/@@download/" + .accession + ".bed.gz"' \
  | while read -r url; do [ -n "$url" ] || continue; fn="$outdir/$(basename "$url")"; [ -f "$fn" ] || curl -sL -o "$fn" "$url" || true; done

  curl -sL "$q_bb" \
  | jq -r '.["@graph"][] | select(.file_format=="bigBed" and (.output_type|test("peaks";"i"))) | "https://www.encodeproject.org/files/" + .accession + "/@@download/" + .accession + ".bigBed"' \
  | while read -r url; do
      [ -n "$url" ] || continue
      fn="$outdir/$(basename "$url")"
      [ -f "$fn" ] || curl -sL -o "$fn" "$url" || true
      if [ -n "$HAS_BB2B" ] && [ -s "$fn" ]; then
        base="${fn%.bigBed}"; tmp="${base}.bed"
        bigBedToBed "$fn" "$tmp" || true
        [ -s "$tmp" ] && gzip -f "$tmp"
      fi
    done

  if [ -z "$HAS_BB2B" ] && ls "$outdir"/*.bigBed 1>/dev/null 2>&1; then
    echo "warning: bigBed files present in $outdir but bigBedToBed not found; skipping conversion"
  fi
}

download_encode_files brain H3K27ac "$BRAIN_H3"
download_encode_files brain DNase   "$BRAIN_DN"
download_encode_files heart H3K27ac "$HEART_H3"
download_encode_files heart DNase   "$HEART_DN"

# Union helpers: merge many peak files into one BED
union_dir () {
  indir="$1"; out="$2"
  if ls "$indir"/* 1>/dev/null 2>&1; then
    ( gzip -cd "$indir"/*.gz 2>/dev/null || cat "$indir"/*.bed 2>/dev/null || true ) \
    | awk 'BEGIN{OFS="\t"} $1~/^chr[0-9XYM]+$/ {print $1,$2,$3}' \
    | sort -k1,1 -k2,2n | bedtools merge -i - > "$out"
  else
    : > "$out"
  fi
}

union_dir "$BRAIN_H3"  "$WORKDIR/brain_h3.bed"
union_dir "$BRAIN_DN"  "$WORKDIR/brain_dnase.bed"
union_dir "$HEART_H3"  "$WORKDIR/heart_h3.bed"
union_dir "$HEART_DN"  "$WORKDIR/heart_dnase.bed"

# Fallbacks if one assay missing
[ -s "$WORKDIR/brain_h3.bed"    ] || cp "$WORKDIR/brain_dnase.bed" "$WORKDIR/brain_h3.bed"
[ -s "$WORKDIR/brain_dnase.bed" ] || cp "$WORKDIR/brain_h3.bed"    "$WORKDIR/brain_dnase.bed"
[ -s "$WORKDIR/heart_h3.bed"    ] || cp "$WORKDIR/heart_dnase.bed" "$WORKDIR/heart_h3.bed"
[ -s "$WORKDIR/heart_dnase.bed" ] || cp "$WORKDIR/heart_h3.bed"    "$WORKDIR/heart_dnase.bed"

# Promoters (±PROMOTER_BP around TSS) and blacklist
prom="$WORKDIR/promoters.bed"
gzip -cd "$GENCODE_GTF" \
| awk -v PBP="$PROMOTER_BP" 'BEGIN{OFS="\t"} $0!~/^#/ && $3=="gene" && $1~/^chr[0-9XYM]+$/ {tss0=($7=="+"?$4-1:$5-1); s=tss0-PBP; if(s<0)s=0; e=tss0+PBP; print $1,s,e}' \
| sort -k1,1 -k2,2n | bedtools merge -i - > "$prom"

bl="$WORKDIR/blacklist.bed"
gzip -cd "$BLACKLIST_GZ" | awk 'BEGIN{OFS="\t"} $1~/^chr/ {print $1,$2,$3}' | sort -k1,1 -k2,2n > "$bl"

# Final enhancers
bedtools intersect -a "$WORKDIR/brain_h3.bed" -b "$WORKDIR/brain_dnase.bed" -u \
| bedtools subtract -a - -b "$prom" \
| bedtools subtract -a - -b "$bl" \
| sort -k1,1 -k2,2n > data/brain_enhancers.bed

bedtools intersect -a "$WORKDIR/heart_h3.bed" -b "$WORKDIR/heart_dnase.bed" -u \
| bedtools subtract -a - -b "$prom" \
| bedtools subtract -a - -b "$bl" \
| sort -k1,1 -k2,2n > data/heart_enhancers.bed

gzip -f data/brain_enhancers.bed
gzip -f data/heart_enhancers.bed

# Ingest into SQLite
sqlite3 "$DB_PATH" <<SQL
PRAGMA journal_mode=WAL;
CREATE TABLE IF NOT EXISTS genes(gene_id INTEGER PRIMARY KEY,symbol TEXT,species_id TEXT,chrom TEXT,start INTEGER,end INTEGER);
CREATE TABLE IF NOT EXISTS enhancers(enh_id INTEGER PRIMARY KEY,species_id TEXT,chrom TEXT,start INTEGER,end INTEGER,tissue TEXT,score REAL,source TEXT);
CREATE TABLE IF NOT EXISTS enhancer_class(enh_id,class TEXT);
CREATE TABLE IF NOT EXISTS gwas_snps(snp_id INTEGER PRIMARY KEY,rsid TEXT UNIQUE,chrom TEXT,pos INTEGER,trait TEXT,pval REAL,source TEXT);
CREATE TABLE IF NOT EXISTS snp_to_enhancer(snp_id INTEGER,enh_id INTEGER,UNIQUE(snp_id,enh_id));
CREATE TABLE IF NOT EXISTS gene_to_enhancer(gene_id INTEGER,enh_id INTEGER,UNIQUE(gene_id,enh_id));
CREATE INDEX IF NOT EXISTS idx_enh_loc    ON enhancers(species_id,chrom,start,end);
CREATE INDEX IF NOT EXISTS idx_enh_tissue ON enhancers(tissue);
CREATE UNIQUE INDEX IF NOT EXISTS idx_enh_unique ON enhancers(species_id,chrom,start,end,tissue);
CREATE INDEX IF NOT EXISTS idx_gene_loc ON genes(species_id,chrom,start);
CREATE INDEX IF NOT EXISTS idx_snp_loc  ON gwas_snps(chrom,pos);
SQL

gzip -cd data/brain_enhancers.bed.gz | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,""}' > "$WORKDIR/brain.tsv"
gzip -cd data/heart_enhancers.bed.gz | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,""}' > "$WORKDIR/heart.tsv"

sqlite3 "$DB_PATH" <<SQL
BEGIN;
CREATE TEMP TABLE IF NOT EXISTS enh_stage(chrom TEXT,start INT,end INT,score REAL);
DELETE FROM enh_stage;
.mode tabs
.import '$WORKDIR/brain.tsv' enh_stage
INSERT OR IGNORE INTO enhancers(species_id,chrom,start,end,tissue,score,source)
  SELECT '$SPECIES_ID',chrom,start,end,'Brain',score,'$SOURCE_LABEL' FROM enh_stage;
INSERT OR IGNORE INTO enhancer_class(enh_id,class)
  SELECT e.enh_id,'conserved'
  FROM enhancers e LEFT JOIN enhancer_class ec ON ec.enh_id=e.enh_id
  WHERE ec.enh_id IS NULL AND e.species_id='$SPECIES_ID' AND e.tissue='Brain' AND e.source='$SOURCE_LABEL';
COMMIT;

BEGIN;
DELETE FROM enh_stage;
.import '$WORKDIR/heart.tsv' enh_stage
INSERT OR IGNORE INTO enhancers(species_id,chrom,start,end,tissue,score,source)
  SELECT '$SPECIES_ID',chrom,start,end,'Heart',score,'$SOURCE_LABEL' FROM enh_stage;
INSERT OR IGNORE INTO enhancer_class(enh_id,class)
  SELECT e.enh_id,'conserved'
  FROM enhancers e LEFT JOIN enhancer_class ec ON ec.enh_id=e.enh_id
  WHERE ec.enh_id IS NULL AND e.species_id='$SPECIES_ID' AND e.tissue='Heart' AND e.source='$SOURCE_LABEL';
COMMIT;

INSERT OR IGNORE INTO gene_to_enhancer(gene_id,enh_id)
  SELECT g.gene_id,e.enh_id
  FROM genes g JOIN enhancers e ON g.species_id=e.species_id AND g.chrom=e.chrom
  WHERE e.species_id='$SPECIES_ID' AND e.tissue IN ('Brain','Heart')
    AND e.start < g.start + ${LINK_KB}000 AND e.end > g.start - ${LINK_KB}000;

INSERT OR IGNORE INTO snp_to_enhancer(snp_id,enh_id)
  SELECT s.snp_id,e.enh_id
  FROM gwas_snps s JOIN enhancers e ON s.chrom=e.chrom
  WHERE e.species_id='$SPECIES_ID' AND e.tissue IN ('Brain','Heart')
    AND s.pos >= e.start AND s.pos < e.end;
SQL

sqlite3 "$DB_PATH" "SELECT tissue, COUNT(*) AS enhancers FROM enhancers WHERE species_id='$SPECIES_ID' GROUP BY tissue ORDER BY tissue;"
