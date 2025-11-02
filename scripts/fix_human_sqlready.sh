#!/usr/bin/env bash
set -euo pipefail

ROOT="$(pwd)"
ENH_RAW="data/raw/human/enhancers"
CTCF_RAW="data/raw/human/ctcf"
ENH_OUT="data/processed/human/enhancers"
CTCF_OUT="data/processed/human/ctcf"
SQL_OUT="data/processed/sqlready"

mkdir -p "$ENH_OUT" "$CTCF_OUT" "$SQL_OUT"

TISSUES=(brain heart_lv liver)

echo "[1/5] Build per-tissue enhancer BED6+tissue"
for T in "${TISSUES[@]}"; do
  tmp="$ENH_OUT/${T}.tmp"
  : > "$tmp"
  shopt -s nullglob
  for f in "$ENH_RAW/$T"/*.bed "$ENH_RAW/$T"/*.bed.gz; do
    if [[ "$f" == *.gz ]]; then
      gzip -cd "$f"
    else
      cat "$f"
    fi | awk -v OFS='\t' -v tissue="$T" '
      $0 !~ /^#/ {
        chrom=$1; start=$2; end=$3;
        name = (NF>=4 ? $4 : ".");
        score = (NF>=5 ? $5 : 0);
        print chrom, start, end, name, score, tissue;
      }' >> "$tmp"
  done
  if [[ -s "$tmp" ]]; then
    LC_ALL=C sort -k1,1 -k2,2n "$tmp" \
      | bedtools merge -i - -c 5,6 -o max,distinct \
      | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,".",$4,$5}' \
      > "$ENH_OUT/human_${T}_H3K27ac.bed"
  fi
  rm -f "$tmp"
done

echo "[2/5] Build stacked enhancers"
cat "$ENH_OUT"/human_*_H3K27ac.bed 2>/dev/null \
  | LC_ALL=C sort -k1,1 -k2,2n \
  > "$ENH_OUT/human_enhancers.stacked.bed"

echo "[3/5] Rebuild CTCF union"
if [[ -f "$CTCF_OUT/human_ctcf.stacked.bed" ]]; then
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,".",$5,"union"}' \
      "$CTCF_OUT/human_ctcf.stacked.bed" \
    | LC_ALL=C sort -k1,1 -k2,2n \
    | bedtools merge -i - -c 5 -o max \
    | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,".",$4,"union"}' \
    > "$CTCF_OUT/human_ctcf.union.bed"
fi

echo "[4/5] Sanity checks"
for f in "$ENH_OUT"/*.bed "$CTCF_OUT"/*.bed; do
  [[ -e "$f" ]] || continue
  echo "[check] $f"
  awk -F'\t' 'NF!=6{bad++} END{print (bad?"BAD_NCOLS="bad:"OK cols")}' "$f"
  awk -F'\t' '$2<0 || $3<=$2{bad++} END{print (bad?"BAD_COORDS="bad:"OK coords")}' "$f"
done

echo "[5/5] Write SQL-ready TSVs"
for f in "$ENH_OUT"/*.bed "$CTCF_OUT"/*.bed; do
  [[ -e "$f" ]] || continue
  b="$(basename "$f" .bed)"
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6}' "$f" > "$SQL_OUT/${b}.tsv"
done

echo "[done] SQL-ready files in: $SQL_OUT"

