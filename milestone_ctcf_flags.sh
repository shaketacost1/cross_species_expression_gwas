#!/usr/bin/env bash
set -euo pipefail

# ==============================
# Milestone 1: CTCF proximity flags (distance + within-window binary)
# ==============================
# What this script does
#   1) Exports per-species, per-tissue enhancer BEDs from your SQLite DB.
#   2) Exports species-matched CTCF BEDs from table `ctcf_sites`.
#   3) Uses bedtools `closest -d` to compute min distance to a CTCF site.
#   4) Writes TSVs with distance and a binary near_ctcf flag (<= WINDOW_BP).
#   5) Imports results back into SQLite as {species}_{tissue}_ctcf_flags.
#   6) Creates/refreshes annotated views {species}_{tissue}_with_ctcf.
#
# Requirements
#   - bedtools, sqlite3
#   - DB schema includes: tables `enhancers`, `species`, `ctcf_sites`
#   - Your species common_names: chicken, dog, pig
#
# Usage
#   PROJ="/Users/shaketaveal/Desktop/Improved Project Alpha" \
#   WINDOW_BP=1000 \
#   SPECIES="chicken,dog,pig" \
#   TISSUES="brain,heart,liver" \
#   ./milestone_ctcf_flags.sh
#
# Notes
#   - WINDOW_BP defaults to 1000 if not set.
#   - Heart method handling includes both H3K27ac:heart* and ATAC:heart*.
#   - If a species lacks CTCF rows, that species is skipped with a warning.
# ==============================

# ---- Config ----
PROJ="${PROJ:-/Users/shaketaveal/Desktop/Improved Project Alpha}"
DB="${DB:-$PROJ/ipa.db}"
TMP="${TMP:-$PROJ/data/tmp}"
OUT="${OUT:-$PROJ/data/processed/ctcf_proximity}"
WINDOW_BP="${WINDOW_BP:-1000}"
SPECIES_LIST="${SPECIES:-chicken,dog,pig}"
TISSUE_LIST="${TISSUES:-brain,heart,liver}"

mkdir -p "$TMP" "$OUT"

echo "[info] PROJ=$PROJ"
echo "[info] DB=$DB"
echo "[info] WINDOW_BP=${WINDOW_BP}bp"
echo "[info] SPECIES=$SPECIES_LIST  TISSUES=$TISSUE_LIST"

require() {
  command -v "$1" >/dev/null 2>&1 || { echo "[fatal] Missing dependency: $1"; exit 1; }
}
require bedtools
require sqlite3

# ---- Helper: build enhancer WHERE clause per tissue ----
# Returns a SQL boolean condition string
tissue_clause() {
  local tissue="$1"
  case "$tissue" in
    brain)
      echo "e.method IN ('H3K27ac:brain','H3K27ac:cortex','H3K27ac:cerebellum')"
      ;;
    heart)
      echo "(e.method LIKE 'H3K27ac:heart%' OR e.method LIKE 'ATAC:heart%')"
      ;;
    liver)
      echo "e.method = 'H3K27ac:liver'"
      ;;
    *)
      echo "1=0"  # unknown tissue
      ;;
  esac
}

# ---- Export species-matched CTCF to BED (BED3) ----
export_ctcf_bed() {
  local species="$1"
  local outbed="$2"
  sqlite3 "$DB
.mode tabs
.headers off
SELECT c.chrom, c.start, c.\"end\"
FROM ctcf_sites c
JOIN species s ON s.species_id = c.species_id
WHERE s.common_name = '$species';
" | awk 'BEGIN{OFS="\t"} NF>=3{print $1,$2,$3}' | sort -k1,1 -k2,2n -k3,3n | uniq > "$outbed" || true

  local n
  n=$(wc -l < "$outbed" 2>/dev/null | awk '{print $1+0}')
  if [[ "$n" -eq 0 ]]; then
    echo "[warn] No CTCF rows for species=$species — skipping CTCF proximity for this species."
    return 1
  else
    echo "[ok] CTCF($species): $n rows → $outbed"
  fi
}

# ---- Export enhancers for a species+tissue to BED6 (with method in col4 for traceability) ----
export_enh_bed() {
  local species="$1"; local tissue="$2"; local outbed="$3"
  local clause
  clause=$(tissue_clause "$tissue")

  sqlite3 "$DB
.mode tabs
.headers off
SELECT e.chrom, e.start, e.\"end\",
       e.method,
       0 AS score,
       '.' AS strand
FROM enhancers e
JOIN species s USING(species_id)
WHERE s.common_name='$species' AND ($clause);
" | awk 'BEGIN{OFS="\t"} NF>=3{print $1,$2,$3,$4,$5,$6}' \
  | sort -k1,1 -k2,2n -k3,3n | uniq > "$outbed"

  local n
  n=$(wc -l < "$outbed" 2>/dev/null | awk '{print $1+0}')
  echo "[ok] Enhancers($species,$tissue): $n rows → $outbed"
}

# ---- Compute closest CTCF + flags ----
closest_and_flags() {
  local enhbed="$1"; local ctcfbed="$2"; local outtsv="$3"; local window="$4"

  # bedtools closest: add min distance in last column (-d)
  # If no CTCF on same chrom, bedtools returns -1 for distance — treat as NA/large
  bedtools closest -d -a "$enhbed" -b "$ctcfbed" 2>/dev/null \
  | awk -v OFS="\t" -v W="$window" '
      BEGIN{NA=999999999}
      {
        chrom=$1; start=$2; end=$3; method=$4;
        d=$NF; if(d<0){d=NA}
        near=(d<=W)?1:0;
        print chrom,start,end,method,d,near
      }' > "$outtsv"

  local n
  n=$(wc -l < "$outtsv" 2>/dev/null | awk '{print $1+0}')
  echo "[ok] closest+flags → $n rows → $outtsv"
}

# ---- Import into SQLite and build annotated views ----
import_sqlite_and_view() {
  local species="$1"; local tissue="$2"; local tsv="$3"; local window="$4"
  local table="${species}_${tissue}_ctcf_flags"
  local view="${species}_${tissue}_with_ctcf"
  local clause
  clause=$(tissue_clause "$tissue")

  sqlite3 "$DB" "
PRAGMA journal_mode=WAL;
CREATE TABLE IF NOT EXISTS ${table}(
  chrom TEXT,
  start INT,
  \"end\" INT,
  method TEXT,
  dist_to_ctcf INT,
  near_ctcf_${window}bp INT
);
DELETE FROM ${table};
"
  sqlite3 "$DB" ".mode tabs" ".import '$tsv' ${table}"
  sqlite3 "$DB" "
CREATE INDEX IF NOT EXISTS idx_${table}_coords ON ${table}(chrom,start,\"end\");

DROP VIEW IF EXISTS ${view};
CREATE VIEW ${view} AS
SELECT e.chrom, e.start, e.\"end\", e.method, e.source,
       f.dist_to_ctcf,
       f.near_ctcf_${window}bp AS near_ctcf
FROM enhancers e
JOIN species s USING(species_id)
LEFT JOIN ${table} f
  ON f.chrom=e.chrom AND f.start=e.start AND f.\"end\"=e.\"end\" AND f.method=e.method
WHERE s.common_name='${species}' AND (${clause});
"
  echo "[sqlite] imported ${table} and refreshed view ${view}"
}

# ---- Summary helper ----
summarize_view() {
  local species="$1"; local tissue="$2"; local window="$3"
  local view="${species}_${tissue}_with_ctcf"
  sqlite3 "$DB" "
SELECT '${species}' AS species,
       '${tissue}'  AS tissue,
       SUM(CASE WHEN near_ctcf=1 THEN 1 ELSE 0 END) AS n_near,
       ROUND(100.0*SUM(CASE WHEN near_ctcf=1 THEN 1 ELSE 0 END)/COUNT(*), 2) AS pct_near,
       ROUND(1.0*AVG(dist_to_ctcf),1) AS mean_dist_bp,
       COUNT(*) AS total
FROM ${view};
"
}

# ==============================
# Main loop
# ==============================
IFS=',' read -r -a SPEC_ARR <<< "$SPECIES_LIST"
IFS=',' read -r -a TISS_ARR <<< "$TISSUE_LIST"

# Export all CTCF once per species
declare -A CTCF_BEDS
for sp in "${SPEC_ARR[@]}"; do
  ctcf_bed="$TMP/${sp}_CTCF.bed"
  if export_ctcf_bed "$sp" "$ctcf_bed"; then
    CTCF_BEDS["$sp"]="$ctcf_bed"
  else
    CTCF_BEDS["$sp"]=""
  fi
done

# Per species+tissue
SUMMARY="$OUT/ctcf_proximity_summary_${WINDOW_BP}bp.tsv"
: > "$SUMMARY"
printf "species\ttissue\tn_near\tpct_near\tmean_dist_bp\ttotal\n" >> "$SUMMARY"

for sp in "${SPEC_ARR[@]}"; do
  ctcf="${CTCF_BEDS[$sp]}"
  if [[ -z "$ctcf" ]]; then
    echo "[warn] Skipping species=$sp (no CTCF exported)."
    continue
  fi
  for tt in "${TISS_ARR[@]}"; do
    enh="$TMP/${sp}_${tt}_enhancers.bed6"
    export_enh_bed "$sp" "$tt" "$enh"

    # If zero enhancers for this combo, skip
    if [[ ! -s "$enh" ]]; then
      echo "[warn] No enhancers for $sp/$tt — skipping."
      continue
    fi

    tsv="$OUT/${sp}_${tt}_ctcf_${WINDOW_BP}bp.tsv"
    closest_and_flags "$enh" "$ctcf" "$tsv" "$WINDOW_BP"
    import_sqlite_and_view "$sp" "$tt" "$tsv" "$WINDOW_BP"

    summarize_view "$sp" "$tt" "$WINDOW_BP" >> "$SUMMARY"
  done
done

echo
echo "---- Summary (${WINDOW_BP}bp) ----"
column -t "$SUMMARY" | sed 's/\t/  /g'
echo "[done] Wrote summary: $SUMMARY"
