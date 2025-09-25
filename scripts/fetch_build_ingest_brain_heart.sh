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

need(){ command -v "$1" >/dev/null || { echo "Missing $1"; exit 1; }; }
need curl; need jq; need bedtools; need sqlite3; need gzip; need awk; need sort
command -v bigBedToBed >/dev/null || { echo "Missing bigBedToBed"; exit 1; }

mkdir -p "$BRAIN_H3" "$BRAIN_DN" "$HEART_H3" "$HEART_DN" data "$WORKDIR"

fetch_files(){
  local organ="$1" assay="$2" outdir="$3"
  local url_bed url_bb
  if [[ "$assay" == "H3K27ac" ]]; then
    url_bed="https://www.encodeproject.org/search/?type=File&status=released&assembly=GRCh38&output_type=peaks&file_format=bed&assay_title=Histone%20ChIP-seq&target.label=H3K27ac&biosample_ontology.organ_slims=${organ}&format=json&limit=all"
    url_bb="https://www.encodeproject.org/search/?type=File&status=released&assembly=GRCh38&output_type=peaks&file_format=bigBed&assay_title=Histone%20ChIP-seq&target.label=H3K27ac&biosample_ontology.organ_slims=${organ}&format=json&limit=all"
  else
    url_bed="https://www.encodeproject.org/search/?type=File&status=released&assembly=GRCh38&output_type=peaks&file_format=bed&assay_title=DNase-seq&biosample_ontology.organ_slims=${organ}&format=json&limit=all"
    url_bb="https://www.encodeproject.org/search/?type=File&status=released&assembly=GRCh38&output_type=peaks&file_format=bigBed&assay_title=DNase-seq&biosample_ontology.organ_slims=${organ}&format=json&limit=all"
  fi
  mkdir -p "$outdir"
  curl -sL "$url_bed" | jq -r '.["@graph"][].href' | sed 's#^#https://www.encodeproject.org#' | sed -E 's#/$##' | awk -F/ '{id=$NF; print $0"/@@download/"id".bed.gz"}' | while read -r u; do fn="$outdir/$(basename "$u")"; [[ -f "$fn" ]] || curl -sL -o "$fn" "$u" || true; done
  curl -sL "$url_bb"  | jq -r '.["@graph"][].href' | sed 's#^#https://www.encodeproject.org#' | sed -E 's#/$##' | awk -F/ '{id=$NF; print $0"/@@download/"id".bigBed"}' | while read -r u; do fn="$outdir/$(basename "$u")"; [[ -f "$fn" ]] || curl -sL -o "$fn" "$u" || true; done
  for f in "$outdir"/*.bigBed 2>/dev/null; do [[ -e "$f" ]] || continue; base="${f%.bigBed}"; tmp="${base}.bed"; bigBedToBed "$f" "$tmp" || true; if [[ -s "$tmp" ]]; then gzip -f "$tmp"; fi; done
}

[[ -f "$GENCODE_GTF"   ]] || curl -L -o "$GENCODE_GTF" "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz"
[[ -f "$BLACKLIST_GZ"  ]] || curl -L -o "$BLACKLIST_GZ" "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz"

fetch_files brain H3K27ac "$BRAIN_H3"
fetch_files brain DNase   "$BRAIN_DN"
fetch_files heart H3K27ac "$HEART_H3"
fetch_files heart DNase   "$HEART_DN"

union_dir(){ local indir="$1" out="$2"; if ls "$indir"/* 1>/dev/null 2>&1; then ( gzip -cd "$indir"/*.gz 2>/dev/null || cat "$indir"/*.bed 2>/dev/null || true ) | awk 'BEGIN{OFS="\t"} $1~/^chr[0-9XYM]+$/ {print $1,$2,$3}' | sort -k1,1 -k2,2n | bedtools merge -i - > "$out"; else : > "$out"; fi; }

union_dir "$BRAIN_H3"  "$WORKDIR/brain_h3.bed"
union_dir "$BRAIN_DN"  "$WORKDIR/brain_dnase.bed"
union_dir "$HEART_H3"  "$WORKDIR/heart_h3.bed"
union_dir "$HEART_DN"  "$WORKDIR/heart_dnase.bed"

[[ -s "$WORKDIR/brain_h3.bed"    ]] || cp "$WORKDIR/brain_dnase.bed" "$WORKDIR/brain_h3.bed"
[[ -s "$WORKDIR/brain_dnase.bed" ]] || cp "$WORKDIR/brain_h3.bed"    "$WORKDIR/brain_dnase.bed"
[[ -s "$WORKDIR/heart_h3.bed"    ]] || cp "$WORKDIR/heart_dnase.bed" "$WORKDIR/heart_h3.bed"
[[ -s "$WORKDIR/heart_dnase.bed" ]] || cp "$WORKDIR/heart_h3.bed"    "$WORKDIR/heart_dnase.bed"

prom=$(mktemp)
gzip -cd "$GENCODE_GTF" | awk -v PBP="$PROMOTER_BP" 'BEGIN{OFS="\t"} $0!~/^#/ && $3=="gene" && $1~/^chr[0-9XYM]+$/ {tss0=($7=="+"?$4-1:$5-1); s=tss0-PBP; if(s<0)s=0; e=tss0+PBP; print $1,s,e}' | sort -k1,1 -k2,2n | bedtools merge -i - > "$prom"
bl=$(mktemp)
gzip -cd "$BLACKLIST_GZ" | awk 'BEGIN{OFS="\t"} $1~/^chr/ {print $1,$2,$3}' | sort -k1,1 -k2,2n > "$bl"

bedtools intersect -a "$WORKDIR/brain_h3.bed" -b "$WORKDIR/brain_dnase.bed" -u | bedtools subtract -a - -b "$prom" | bedtools subtract -a - -b "$bl" | sort -k1,1 -k2,2n > data/brain_enhancers.bed
bedtools intersect -a "$WORKDIR/heart_h3.bed" -b "$WORKDIR/heart_dnase.bed" -u | bedtools subtract -a - -b "$prom" | bedtools subtract -a - -b "$bl" | sort -k1,1 -k2,2n > data/heart_enhancers.bed
gzip -f data/brain_enhancers.bed
gzip -f data/heart_enhancers.bed

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

prep_tsv(){ gzip -cd "$1" | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,""}' | sort -k1,1 -k2,2n > "$2"; }
ingest(){
  local tissue="$1" bed="$2" tsv="$3"
  prep_tsv "$bed" "$tsv"
  sqlite3 "$DB_PATH" <<SQL
BEGIN;
CREATE TEMP TABLE IF NOT EXISTS enh_stage(chrom TEXT,start INT,end INT,score REAL);
DELETE FROM enh_stage;
.mode tabs
.import '$tsv' enh_stage
INSERT OR IGNORE INTO enhancers(species_id,chrom,start,end,tissue,score,source)
  SELECT '$SPECIES_ID',chrom,start,end,'$tissue',score,'$SOURCE_LABEL' FROM enh_stage;
INSERT OR IGNORE INTO enhancer_class(enh_id,class)
  SELECT e.enh_id,'conserved'
  FROM enhancers e LEFT JOIN enhancer_class ec ON ec.enh_id=e.enh_id
  WHERE ec.enh_id IS NULL AND e.species_id='$SPECIES_ID' AND e.tissue='$tissue' AND e.source='$SOURCE_LABEL';
COMMIT;
INSERT OR IGNORE INTO gene_to_enhancer(gene_id,enh_id)
  SELECT g.gene_id,e.enh_id
  FROM genes g JOIN enhancers e ON g.species_id=e.species_id AND g.chrom=e.chrom
  WHERE e.species_id='$SPECIES_ID' AND e.tissue='$tissue'
    AND e.start < g.start + ${LINK_KB}000 AND e.end > g.start - ${LINK_KB}000;
INSERT OR IGNORE INTO snp_to_enhancer(snp_id,enh_id)
  SELECT s.snp_id,e.enh_id
  FROM gwas_snps s JOIN enhancers e ON s.chrom=e.chrom
  WHERE e.species_id='$SPECIES_ID' AND e.tissue='$tissue'
    AND s.pos >= e.start AND s.pos < e.end;
SQL
}

ingest "Brain" "data/brain_enhancers.bed.gz" "$WORKDIR/brain.tsv"
ingest "Heart" "data/heart_enhancers.bed.gz" "$WORKDIR/heart.tsv"

sqlite3 "$DB_PATH" "SELECT tissue, COUNT(*) AS enhancers FROM enhancers WHERE species_id='$SPECIES_ID' GROUP BY tissue ORDER BY tissue;"
