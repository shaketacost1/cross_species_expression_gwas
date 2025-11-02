#!/usr/bin/env bash
# dog_pig_heart_atac_to_sqlite.sh
# Merge heart ATAC peaks for dog/pig (if present), load into SQLite, and print QC.
# Usage: bash dog_pig_heart_atac_to_sqlite.sh
# Run from the project root: ~/Desktop/"Improved Project Alpha"

set -euo pipefail

echo "[INFO] Starting dog/pig heart ATAC → SQLite loader"

# --- Helpers ---
have() { command -v "$1" >/dev/null 2>&1; }

# Resolve project root as current dir
ROOT="$(pwd)"
DB="${ROOT}/ipa.db"
BED2SQL="${ROOT}/scripts/etl/bed_to_sqlite.py"

# Sanity checks
if [ ! -f "$DB" ]; then
  echo "[ERROR] Cannot find ipa.db at $DB"
  exit 1
fi
if [ ! -f "$BED2SQL" ]; then
  echo "[ERROR] Cannot find bed_to_sqlite.py at $BED2SQL"
  exit 1
fi
if ! have bedtools; then
  echo "[ERROR] bedtools not found in PATH. Please install bedtools."
  exit 1
fi

# Streaming cat that supports .bed and .bed.gz
cat_bed_stream() {
  local f="$1"
  case "$f" in
    *.gz) gzip -cd "$f" ;;
    *)    cat "$f" ;;
  esac
}

merge_and_load() {
  local species="$1"           # dog or pig
  local raw_dir="$2"           # e.g., data/raw/dog/atac/heart
  local out_dir="data/processed/${species}/enhancers"
  local out_bed="${out_dir}/${species}_heart_ATAC.bed"
  local source="BarkBase_FAANG"
  local method='ATAC:heart (proxy_enhancer)'

  if [ ! -d "$raw_dir" ]; then
    echo "[WARN] No directory: ${raw_dir} — skipping ${species} heart ATAC."
    return 0
  fi

  shopt -s nullglob
  files=( "${raw_dir}"/*.bed "${raw_dir}"/*.bed.gz )
  shopt -u nullglob

  if [ ${#files[@]} -eq 0 ]; then
    echo "[WARN] No peak files found under ${raw_dir} — skipping ${species}."
    return 0
  fi

  echo "[INFO] Merging $(printf '%s ' "${#files[@]}") files for ${species}…"
  mkdir -p "${out_dir}"
  # Merge: take columns 1-3, sort, bedtools merge
  {
    for f in "${files[@]}"; do
      cat_bed_stream "$f" | awk 'BEGIN{OFS="\t"}{print $1,$2,$3}'
    done     | LC_ALL=C sort -k1,1 -k2,2n -k3,3n     | bedtools merge -i -
  } > "${out_bed}"

  echo "[INFO] Created ${out_bed} ($(wc -l < "${out_bed}") rows)"

  echo "[INFO] Loading into SQLite (species=${species}, method='${method}')…"
  python "${BED2SQL}" --db "${DB}" --species "${species}" --type enhancer     --bed "${out_bed}" --source "${source}" --method "${method}"
}

# --- Run for dog and pig ---
merge_and_load "dog" "data/raw/dog/atac/heart"
merge_and_load "pig" "data/raw/pig/atac/heart"

# --- QC ---
echo
echo "[QC] Counts for dog/pig ATAC heart (proxy enhancers):"
sqlite3 "${DB}" <<'SQL'
.headers on
.mode column
SELECT s.common_name AS species, e.method, COUNT(*) AS n
FROM enhancers e JOIN species s USING(species_id)
WHERE s.common_name IN ('dog','pig') AND e.method LIKE 'ATAC:heart%'
GROUP BY 1,2 ORDER BY 1,2;
SQL

echo
echo "[QC] Total heart enhancers by species (H3K27ac heart + ATAC heart proxy):"
sqlite3 "${DB}" <<'SQL'
.headers on
.mode column
SELECT s.common_name AS species, COUNT(*) AS total_heart_enhancers
FROM enhancers e JOIN species s USING(species_id)
WHERE s.common_name IN ('dog','pig')
  AND e.method IN ('H3K27ac:heart','ATAC:heart (proxy_enhancer)')
GROUP BY 1 ORDER BY 1;
SQL

echo "[INFO] Done."
