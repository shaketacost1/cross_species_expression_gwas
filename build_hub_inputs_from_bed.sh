#!/usr/bin/env bash
set -euo pipefail

PROJ="/Users/shaketaveal/Desktop/Improved Project Alpha"
DB="$PROJ/ipa.db"
TMP="$PROJ/data/tmp";                    mkdir -p "$TMP"
OUTDIR="$PROJ/data/processed/hub_conservation";  mkdir -p "$OUTDIR"


HUMAN_BRAIN="$PROJ/data/processed/human/enhancers/human_brain_H3K27ac.bed"
HUMAN_HEART="$PROJ/data/processed/human/enhancers/human_heart_lv_H3K27ac.bed"
HUMAN_LIVER="$PROJ/data/processed/human/enhancers/human_liver_H3K27ac.bed"


MOUSE_BRAIN="$PROJ/data/processed/mouse/enhancers/mouse_brain_H3K27ac.bed"
MOUSE_HEART="$PROJ/data/processed/mouse/enhancers/mouse_heart_H3K27ac.bed"
MOUSE_LIVER="$PROJ/data/processed/mouse/enhancers/mouse_liver_H3K27ac.bed"

DOG_BRAIN="$PROJ/data/processed/dog/enhancers/dog_brain_H3K27ac.bed"
DOG_HEART="$PROJ/data/processed/dog/enhancers/dog_heart_ATAC.bed"
DOG_LIVER="$PROJ/data/processed/dog/enhancers/dog_liver_H3K27ac.bed"

PIG_BRAIN="$PROJ/data/processed/pig/enhancers/pig_brain_H3K27ac.bed"
PIG_HEART="$PROJ/data/processed/pig/enhancers/pig_heart_ATAC.bed"
PIG_LIVER="$PROJ/data/processed/pig/enhancers/pig_liver_H3K27ac.bed"

CHICK_BRAIN="$PROJ/data/processed/chicken/enhancers/chicken_brain_H3K27ac.bed"
CHICK_HEART="$PROJ/data/processed/chicken/enhancers/chicken_heart_final.chrfix.bed"
CHICK_LIVER="$PROJ/data/processed/chicken/enhancers/chicken_liver_H3K27ac.bed"


REFCHAINS="$PROJ/data/reference/chains"
HG38_TO_MM10="$REFCHAINS/hg38ToMm10.over.chain.gz"
HG38_TO_CANFAM3="$REFCHAINS/hg38ToCanFam3.over.chain.gz"
HG38_TO_SUSSCR11="$REFCHAINS/hg38ToSusScr11.over.chain.gz"
HG38_TO_GALGAL6="$REFCHAINS/hg38ToGalGal6.over.chain.gz"

need() { command -v "$1" >/dev/null || { echo "ERROR: missing \`$1\` in PATH" >&2; exit 1; }; }
need liftOver
need bedtools
need awk
need sort
need cut

mk_id_bed4() { # input BED3 -> output BED4 with id=chr:start-end
  local IN="$1" OUT="$2"
  awk -v OFS='\t' '{print $1,$2,$3,$1":"$2"-"$3}' "$IN" > "$OUT"
}

liftover() { # liftOver BED4 and keep the 4th col (id)
  local IN="$1" CHAIN="$2" OUT="$3" UNM="$4"
  liftOver "$IN" "$CHAIN" "$OUT" "$UNM"
}

intersect_ids() { # intersect lifted BED4 with species enhancer BED, emit IDs (col 4)
  local LIFTED="$1" TARGET="$2" IDS_OUT="$3"
  bedtools intersect -a "$LIFTED" -b "$TARGET" -wa -u \
    | cut -f4 | awk 'NF' | sort -u > "$IDS_OUT"
}

build_one_tissue() {
  local T="$1" HUMAN="$2" MOUSE="$3" DOG="$4" PIG="$5" CHICK="$6"

  echo "[*] Tissue: $T"


  for f in "$HUMAN" "$MOUSE" "$DOG" "$PIG" "$CHICK"; do
    [[ -s "$f" ]] || { echo "  MISSING: $f" >&2; return 1; }
  done


  local NATIVE="$TMP/hub_${T}_human_native.bed4"
  mk_id_bed4 "$HUMAN" "$NATIVE"
  awk -v OFS='\t' '{print $4,$1,$2,$3}' "$NATIVE" > "$OUTDIR/hub_${T}_id2native.tsv"


  local TO_MM="$TMP/hub_${T}_to_mouse.bed4"  UN_MM="$TMP/hub_${T}_to_mouse.unmapped"
  local TO_DG="$TMP/hub_${T}_to_dog.bed4"    UN_DG="$TMP/hub_${T}_to_dog.unmapped"
  local TO_PG="$TMP/hub_${T}_to_pig.bed4"    UN_PG="$TMP/hub_${T}_to_pig.unmapped"
  local TO_CH="$TMP/hub_${T}_to_chick.bed4"  UN_CH="$TMP/hub_${T}_to_chick.unmapped"

  liftover "$NATIVE" "$HG38_TO_MM10"   "$TO_MM" "$UN_MM"
  liftover "$NATIVE" "$HG38_TO_CANFAM3" "$TO_DG" "$UN_DG"
  liftover "$NATIVE" "$HG38_TO_SUSSCR11" "$TO_PG" "$UN_PG"
  liftover "$NATIVE" "$HG38_TO_GALGAL6"  "$TO_CH" "$UN_CH"


  intersect_ids "$TO_MM" "$MOUSE" "$OUTDIR/hub_${T}_ids_mouse.txt"
  intersect_ids "$TO_DG" "$DOG"   "$OUTDIR/hub_${T}_ids_dog.txt"
  intersect_ids "$TO_PG" "$PIG"   "$OUTDIR/hub_${T}_ids_pig.txt"
  intersect_ids "$TO_CH" "$CHICK" "$OUTDIR/hub_${T}_ids_chick.txt"

  awk -v FS="\t" -v OFS="\t" \
      -v MOUSE_IDS="$OUTDIR/hub_${T}_ids_mouse.txt" \
      -v DOG_IDS="$OUTDIR/hub_${T}_ids_dog.txt" \
      -v PIG_IDS="$OUTDIR/hub_${T}_ids_pig.txt" \
      -v CHK_IDS="$OUTDIR/hub_${T}_ids_chick.txt" '
    function load(file, arr, x){ while((getline x < file)>0){ gsub(/\r/,"",x); if(x!="") arr[x]=1 } close(file) }
    BEGIN{ print "chr","start","end","mouse","dog","pig","chick";
           load(MOUSE_IDS,MM); load(DOG_IDS,DD); load(PIG_IDS,PP); load(CHK_IDS,CC) }
    { id=$1; chr=$2; s=$3; e=$4;
      print chr, s, e, (id in MM?1:0), (id in DD?1:0), (id in PP?1:0), (id in CC?1:0) }
  ' "$OUTDIR/hub_${T}_id2native.tsv" > "$OUTDIR/hub_${T}_flags_all.tsv"

  echo "[hub] wrote $OUTDIR/hub_${T}_flags_all.tsv"
}


build_one_tissue heart "$HUMAN_HEART" "$MOUSE_HEART" "$DOG_HEART" "$PIG_HEART" "$CHICK_HEART"
build_one_tissue liver "$HUMAN_LIVER" "$MOUSE_LIVER" "$DOG_LIVER" "$PIG_LIVER" "$CHICK_LIVER"
