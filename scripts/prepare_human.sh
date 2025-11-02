#!/usr/bin/env bash
set -euo pipefail

# ---- paths & tissues ----
PROJ="$(pwd)"
ENH_RAW="$PROJ/data/raw/human/enhancers"
CTCF_RAW="$PROJ/data/raw/human/ctcf"
ENH_OUT="$PROJ/data/processed/human/enhancers"
CTCF_OUT="$PROJ/data/processed/human/ctcf"
SQL_OUT="$PROJ/data/processed/sqlready"
mkdir -p "$ENH_OUT" "$CTCF_OUT" "$SQL_OUT"

TISSUES=(brain heart_lv liver)

# ---- enhancers: normalize to BED6+tissue; per-tissue merged; stacked ----
: > "$ENH_OUT/human_enhancers.stacked.bed"
for T in "${TISSUES[@]}"; do
  tmp="$ENH_OUT/${T}.tmp"; : > "$tmp"

  # iterate files; don't put "2>/dev/null" in the 'for' list
  for f in "$ENH_RAW/$T"/*.bed*; do
    [[ -e "$f" ]] || continue
    (gzip -cd "$f" 2>/dev/null || cat "$f") \
      | awk -v OFS='\t' -v tissue="$T" '($1!~/^#/){
          c=$1;s=$2;e=$3;n=(NF>=4?$4:".");sc=(NF>=5?$5:0);
          print c,s,e,n,sc,tissue
        }' >> "$tmp"
  done

  out="$ENH_OUT/human_${T}_H3K27ac.bed"
  LC_ALL=C sort -k1,1 -k2,2n "$tmp" \
    | bedtools merge -i - -c 5 -o max > "$out"
  rm -f "$tmp"
  cat "$out" >> "$ENH_OUT/human_enhancers.stacked.bed"
done
LC_ALL=C sort -k1,1 -k2,2n -o "$ENH_OUT/human_enhancers.stacked.bed" "$ENH_OUT/human_enhancers.stacked.bed"

# ---- CTCF: normalize to BED6+tissue; stacked + union ----
tmpctcf="$CTCF_OUT/all.tmp"; : > "$tmpctcf"
for T in "${TISSUES[@]}"; do
  shopt -s nullglob
  for f in "$CTCF_RAW/$T"/*.bed* "$CTCF_RAW/$T"/*.gff*; do
    [[ -e "$f" ]] || continue
    case "$f" in
      *.gff|*.gff3|*.gff.gz|*.gff3.gz)
        (gzip -cd "$f" 2>/dev/null || cat "$f") \
          | awk -v OFS='\t' -v tissue="$T" '($0!~/^#/){print $1,$4-1,$5,".",($6=="."?0:$6),tissue}' >> "$tmpctcf"
        ;;
      *)
        (gzip -cd "$f" 2>/dev/null || cat "$f") \
          | awk -v OFS='\t' -v tissue="$T" '($1!~/^#/){print $1,$2,$3,(NF>=4?$4:"."),(NF>=5?$5:0),tissue}' >> "$tmpctcf"
        ;;
    esac
  done
  shopt -u nullglob
done

LC_ALL=C sort -k1,1 -k2,2n "$tmpctcf" > "$CTCF_OUT/human_ctcf.stacked.bed"
bedtools merge -i "$CTCF_OUT/human_ctcf.stacked.bed" -c 5 -o max > "$CTCF_OUT/human_ctcf.union.bed"
rm -f "$tmpctcf"

# ---- sanity checks ----
for f in "$ENH_OUT"/*.bed "$CTCF_OUT"/*.bed; do
  [[ -e "$f" ]] || continue
  echo "[check] $f"
  awk -F'\t' 'NF!=6{bad++} END{print (bad?"BAD_NCOLS="bad:"OK cols")}' "$f"
  awk -F'\t' '$2<0 || $3<=$2{bad++} END{print (bad?"BAD_COORDS="bad:"OK coords")}' "$f"
done

# ---- SQL-ready TSVs ----
for f in "$ENH_OUT"/*.bed "$CTCF_OUT"/*.bed; do
  [[ -e "$f" ]] || continue
  b="$(basename "$f" .bed)"
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6}' "$f" > "$SQL_OUT/${b}.tsv"
done

echo "[done] SQL-ready files in: $SQL_OUT"
