#!/usr/bin/env bash
set -euo pipefail
PROJ="/Users/shaketaveal/Desktop/Improved Project Alpha"
OUTDIR="$PROJ/data/processed/hub_conservation"

make_hub_flags() {
  local T="$1"
  awk -v FS="\t" -v OFS="\t" \
      -v MOUSE="$OUTDIR/hub_${T}_ids_mouse.txt" \
      -v DOG="$OUTDIR/hub_${T}_ids_dog.txt" \
      -v PIG="$OUTDIR/hub_${T}_ids_pig.txt" \
      -v CHICK="$OUTDIR/hub_${T}_ids_chick.txt" '
    function load(file, arr, x){ while ((getline x < file) > 0) { gsub(/\r/,"",x); if (x!="") arr[x]=1 } close(file) }
    BEGIN{ print "chr","start","end","mouse","dog","pig","chick";
           load(MOUSE,MM); load(DOG,DD); load(PIG,PP); load(CHICK,CC) }
    { id=$1; chr=$2; s=$3; e=$4;
      print chr, s, e, (id in MM?1:0), (id in DD?1:0), (id in PP?1:0), (id in CC?1:0) }
  ' "$OUTDIR/hub_${T}_id2native.tsv" > "$OUTDIR/hub_${T}_flags_all.tsv"
  echo "[hub] wrote $OUTDIR/hub_${T}_flags_all.tsv"
}

make_hub_flags brain
