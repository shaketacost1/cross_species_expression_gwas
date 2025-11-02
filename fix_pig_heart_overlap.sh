#!/usr/bin/env bash
set -euo pipefail
PROJ="/Users/shaketaveal/Desktop/Improved Project Alpha"
OUTDIR="$PROJ/data/processed/hub_conservation"
ENH_DIR="$PROJ/data/processed/pig/enhancers"
TMP="$PROJ/data/tmp"

# 1) choose pig heart enhancer (prefer H3K27ac)
PIG_HEART=""
if [ -s "$ENH_DIR/pig_heart_H3K27ac.bed" ]; then
  PIG_HEART="$ENH_DIR/pig_heart_H3K27ac.bed"
elif [ -s "$ENH_DIR/pig_heart_ATAC.bed" ]; then
  PIG_HEART="$ENH_DIR/pig_heart_ATAC.bed"
else
  echo "ERROR: no pig heart enhancer BED found in $ENH_DIR" >&2
  exit 1
fi
echo "[i] Using pig heart enhancer: $PIG_HEART"

# 2) normalize naming to chr*, sort/unique (safe to run repeatedly)
awk 'BEGIN{OFS="\t"} { $1 = ($1 ~ /^chr/) ? $1 : "chr"$1; print }' "$PIG_HEART" \
  | sort -k1,1 -k2,2n -k3,3n | awk '!seen[$0]++' > "$PIG_HEART.fixed"
mv -f "$PIG_HEART.fixed" "$PIG_HEART"

# 3) confirm sizes
echo -n "[i] pig heart enhancer rows: "
wc -l < "$PIG_HEART" || true
echo -n "[i] heart->pig lifted rows:  "
wc -l < "$TMP/hub_heart_to_pig.bed4" || true

# 4) overlap diagnostics at different thresholds
echo "[i] Intersect diagnostics:"
for frac in 0 0.1 0.5; do
  if [ "$frac" = "0" ]; then
    bedtools intersect -a "$TMP/hub_heart_to_pig.bed4" -b "$PIG_HEART" -wa -u | wc -l | awk -v f="any" '{print "  -f " f ": " $1}'
  else
    bedtools intersect -a "$TMP/hub_heart_to_pig.bed4" -b "$PIG_HEART" -wa -u -f "$frac" | wc -l | awk -v f="$frac" '{print "  -f " f ": " $1}'
  fi
done

# 5) regenerate pig IDs (any overlap; adjust -f if desired)
bedtools intersect -a "$TMP/hub_heart_to_pig.bed4" -b "$PIG_HEART" -wa -u \
  | cut -f4 | awk 'NF' | sort -u > "$OUTDIR/hub_heart_ids_pig.txt"
echo -n "[i] ids_pig count: "
wc -l < "$OUTDIR/hub_heart_ids_pig.txt" || true

# 6) rebuild flags for heart using existing inputs
awk -v FS="\t" -v OFS="\t" \
    -v MOUSE_IDS="$OUTDIR/hub_heart_ids_mouse.txt" \
    -v DOG_IDS="$OUTDIR/hub_heart_ids_dog.txt" \
    -v PIG_IDS="$OUTDIR/hub_heart_ids_pig.txt" \
    -v CHK_IDS="$OUTDIR/hub_heart_ids_chick.txt" '
  function load(file, arr, x){ while((getline x < file)>0){ gsub(/\r/,"",x); if(x!="") arr[x]=1 } close(file) }
  BEGIN{ print "chr","start","end","mouse","dog","pig","chick";
         load(MOUSE_IDS,MM); load(DOG_IDS,DD); load(PIG_IDS,PP); load(CHK_IDS,CC) }
  { id=$1; chr=$2; s=$3; e=$4;
    print chr, s, e, (id in MM?1:0), (id in DD?1:0), (id in PP?1:0), (id in CC?1:0) }
' "$OUTDIR/hub_heart_id2native.tsv" > "$OUTDIR/hub_heart_flags_all.tsv"

echo "[✓] Wrote $OUTDIR/hub_heart_flags_all.tsv"
