#!/usr/bin/env bash
set -euo pipefail

# --- Paths ---
PROJ="/Users/shaketaveal/Desktop/Improved Project Alpha"
DB="$PROJ/ipa.db"
TMP="$PROJ/data/tmp"
OUT="$PROJ/data/processed/hub_conservation"
CHAINS="$PROJ/data/reference/chains"
mkdir -p "$TMP" "$OUT" "$CHAINS"

# --- Inputs: human anchors (hub) ---
H_BRAIN="$PROJ/data/processed/human/enhancers/human_brain_H3K27ac.bed"
H_HEART="$PROJ/data/processed/human/enhancers/human_heart_lv_H3K27ac.bed"
H_LIVER="$PROJ/data/processed/human/enhancers/human_liver_H3K27ac.bed"

# --- Target enhancer sets (same-tissue) ---
# mouse
M_BRAIN="$PROJ/data/processed/mouse/enhancers/mouse_brain_H3K27ac.bed"
M_HEART="$PROJ/data/processed/mouse/enhancers/mouse_heart_H3K27ac.bed"
M_LIVER="$PROJ/data/processed/mouse/enhancers/mouse_liver_H3K27ac.bed"

# dog (BarkBase heart = ATAC, brain/liver = H3K27ac)
D_BRAIN="$PROJ/data/processed/dog/enhancers/dog_brain_H3K27ac.bed"
D_HEART="$PROJ/data/processed/dog/enhancers/dog_heart_ATAC.bed"
D_LIVER="$PROJ/data/processed/dog/enhancers/dog_liver_H3K27ac.bed"

# pig
P_BRAIN="$PROJ/data/processed/pig/enhancers/pig_brain_H3K27ac.bed"
P_HEART="$PROJ/data/processed/pig/enhancers/pig_heart_ATAC.bed"   # use your pig heart set here
P_LIVER="$PROJ/data/processed/pig/enhancers/pig_liver_H3K27ac.bed"

# chicken (heart is RNAproxy; brain/liver are H3K27ac)
C_BRAIN="$PROJ/data/processed/chicken/enhancers/chicken_brain_H3K27ac.bed"
C_HEART="$PROJ/data/processed/chicken/enhancers/chicken_heart_RNAproxy.clean.bed"
C_LIVER="$PROJ/data/processed/chicken/enhancers/chicken_liver_H3K27ac.bed"

# --- Chains: hg38 -> target assemblies ---
# If any are missing: curl -O from UCSC; we’ve used these directions earlier.
HG38_TO_MM10="$CHAINS/hg38ToMm10.over.chain.gz"
HG38_TO_CANFAM3="$CHAINS/hg38ToCanFam3.over.chain.gz"
HG38_TO_SUSSCR11="$CHAINS/hg38ToSusScr11.over.chain.gz"
HG38_TO_GALGAL6="$CHAINS/hg38ToGalGal6.over.chain.gz"

# --- helper: make BED4 with native ID ---
mk4() { in="$1"; out="$2"; awk 'BEGIN{OFS="\t"} NF>=3{print $1,$2,$3,$1":"$2"-"$3}' "$in" > "$out"; }

# --- helper: liftover hub (human) to one target ---
lift_hub() { bed4="$1"; chain="$2"; out="$3"; unm="$4";
  liftOver -bedPlus=4 "$bed4" "$chain" "$out" "$unm"
}

# --- helper: flags by intersecting lifted hub IDs with a target enhancer set ---
mk_flags_one_target() { SPEC="$1"; TISS="$2"; lifted="$3"; target="$4"; LABEL="$5";
  # IDs of lifted hub intervals that overlap target set
  bedtools intersect -u -a "$lifted" -b "$target" | cut -f4 | sort -u > "$OUT/hub_${TISS}_ids_${LABEL}.txt"
}

# --- helper: combine species flags back to human native coords ---
# Arguments: TISSUE  IDMAP  [pairs of label+file ...]
emit_flags() {
  T="$1"; IDMAP="$2"; shift 2
  # Load all ID lists into awk assoc arrays; then print row: chrom start end [per-target flags] [n_targets]
  awk -vFS="\t" -vOFS="\t" "$(
    i=1
    while [ "$#" -gt 1 ]; do
      L="$1"; F="$2"; shift 2
      echo "BEGIN{ while((getline x < \"$F\")>0){ H[L][x]=1 } }"
      labels+=("$L")
    done
    echo '{
      id=$1; chr=$2; s=$3; e=$4; n=0; out=chr FS s FS e;
      '
    # print flags in the same order as provided
    for L in "${labels[@]}"; do
      echo "f = ((H[\"$L\"][id]==1)?1:0); out = out FS f; n += f;"
    done
    echo 'print out FS n
    }'
  )" "$IDMAP"
}

# --- per tissue driver: lift human -> all 4 targets and build a combined flags table ---
run_tissue() {
  T="$1"  # brain|heart|liver
  H_IN="$2"
  # target sets and chains per tissue
  case "$T" in
    brain)
      TG=( mouse "$M_BRAIN" "$HG38_TO_MM10"
           dog   "$D_BRAIN" "$HG38_TO_CANFAM3"
           pig   "$P_BRAIN" "$HG38_TO_SUSSCR11"
           chick "$C_BRAIN" "$HG38_TO_GALGAL6" )
      ;;
    heart)
      TG=( mouse "$M_HEART" "$HG38_TO_MM10"
           dog   "$D_HEART" "$HG38_TO_CANFAM3"
           pig   "$P_HEART" "$HG38_TO_SUSSCR11"
           chick "$C_HEART" "$HG38_TO_GALGAL6" )
      ;;
    liver)
      TG=( mouse "$M_LIVER" "$HG38_TO_MM10"
           dog   "$D_LIVER" "$HG38_TO_CANFAM3"
           pig   "$P_LIVER" "$HG38_TO_SUSSCR11"
           chick "$C_LIVER" "$HG38_TO_GALGAL6" )
      ;;
  esac

  echo ">>> human-hub: $T"

  # 1) BED4 with human-native ID
  H4="$TMP/hub_${T}_human_native.bed4"
  mk4 "$H_IN" "$H4"

  # map id -> native coords
  IDMAP="$OUT/hub_${T}_id2native.tsv"
  awk 'BEGIN{FS=OFS="\t"}{print $4,$1,$2,$3}' "$H4" | sort -u > "$IDMAP"

  # 2) For each target: lift + collect IDs that overlap target
  # Also remember the lifted BED path for potential debugging
  IDS_TO_PASS=()
  i=0
  while [ $i -lt ${#TG[@]} ]; do
    LAB="${TG[$i]}"; TGT="${TG[$((i+1))]}"; CH="${TG[$((i+2))]}"; i=$((i+3))
    LIFT="$TMP/hub_${T}_to_${LAB}.bed4"; UNM="$TMP/hub_${T}_to_${LAB}.unmapped"
    echo "    - liftOver human->$LAB  ($CH)"
    lift_hub "$H4" "$CH" "$LIFT" "$UNM"

    echo "    - intersect lifted->$LAB vs target $LAB ($T)"
    mk_flags_one_target "human" "$T" "$LIFT" "$TGT" "$LAB"

    IDS_TO_PASS+=("$LAB" "$OUT/hub_${T}_ids_${LAB}.txt")
  done

  # 3) Combine per-target ID lists back to human-native coords, writing one row per native enhancer
  OUTFLAGS="$OUT/hub_${T}_flags_all.tsv"
  emit_flags "$T" "$IDMAP" "${IDS_TO_PASS[@]}" > "$OUTFLAGS"

  # columns: chrom  start  end  conserved_in_mouse  conserved_in_dog  conserved_in_pig  conserved_in_chicken  n_targets
  echo "[flags] Wrote $OUTFLAGS"
}

# -------- run all three tissues --------
run_tissue brain "$H_BRAIN"
run_tissue heart "$H_HEART"
run_tissue liver "$H_LIVER"

# -------- Import into SQLite and build views --------
sqlite3 "$DB" "
CREATE TABLE IF NOT EXISTS human_brain_hub_flags(
  chrom TEXT, start INT, \"end\" INT,
  conserved_in_mouse INT, conserved_in_dog INT,
  conserved_in_pig INT, conserved_in_chicken INT,
  n_targets INT
);
CREATE TABLE IF NOT EXISTS human_heart_hub_flags(
  chrom TEXT, start INT, \"end\" INT,
  conserved_in_mouse INT, conserved_in_dog INT,
  conserved_in_pig INT, conserved_in_chicken INT,
  n_targets INT
);
CREATE TABLE IF NOT EXISTS human_liver_hub_flags(
  chrom TEXT, start INT, \"end\" INT,
  conserved_in_mouse INT, conserved_in_dog INT,
  conserved_in_pig INT, conserved_in_chicken INT,
  n_targets INT
);
DELETE FROM human_brain_hub_flags;
DELETE FROM human_heart_hub_flags;
DELETE FROM human_liver_hub_flags;
"
sqlite3 "$DB" ".mode tabs" ".import '$OUT/hub_brain_flags_all.tsv' human_brain_hub_flags"
sqlite3 "$DB" ".mode tabs" ".import '$OUT/hub_heart_flags_all.tsv' human_heart_hub_flags"
sqlite3 "$DB" ".mode tabs" ".import '$OUT/hub_liver_flags_all.tsv' human_liver_hub_flags"

# Annotated views: join back to human enhancers
sqlite3 "$DB" <<'SQL'
DROP VIEW IF EXISTS human_brain_hub_annotated;
CREATE VIEW human_brain_hub_annotated AS
SELECT e.chrom,e.start,e."end",e.method,e.source,
       f.conserved_in_mouse,f.conserved_in_dog,f.conserved_in_pig,f.conserved_in_chicken,f.n_targets
FROM enhancers e JOIN species s USING(species_id)
LEFT JOIN human_brain_hub_flags f
  ON f.chrom=e.chrom AND f.start=e.start AND f."end"=e."end"
WHERE s.common_name='human'
  AND e.method='H3K27ac:brain';

DROP VIEW IF EXISTS human_heart_hub_annotated;
CREATE VIEW human_heart_hub_annotated AS
SELECT e.chrom,e.start,e."end",e.method,e.source,
       f.conserved_in_mouse,f.conserved_in_dog,f.conserved_in_pig,f.conserved_in_chicken,f.n_targets
FROM enhancers e JOIN species s USING(species_id)
LEFT JOIN human_heart_hub_flags f
  ON f.chrom=e.chrom AND f.start=e.start AND f."end"=e."end"
WHERE s.common_name='human'
  AND e.method LIKE 'H3K27ac:heart%';

DROP VIEW IF EXISTS human_liver_hub_annotated;
CREATE VIEW human_liver_hub_annotated AS
SELECT e.chrom,e.start,e."end",e.method,e.source,
       f.conserved_in_mouse,f.conserved_in_dog,f.conserved_in_pig,f.conserved_in_chicken,f.n_targets
FROM enhancers e JOIN species s USING(species_id)
LEFT JOIN human_liver_hub_flags f
  ON f.chrom=e.chrom AND f.start=e.start AND f."end"=e."end"
WHERE s.common_name='human'
  AND e.method='H3K27ac:liver';
SQL

# Quick roll-up summary
sqlite3 "$DB" -cmd ".mode column" -cmd ".headers on" "
WITH U AS (
  SELECT 'brain' AS tissue, n_targets FROM human_brain_hub_flags
  UNION ALL SELECT 'heart', n_targets FROM human_heart_hub_flags
  UNION ALL SELECT 'liver', n_targets FROM human_liver_hub_flags
)
SELECT tissue,
       SUM(CASE WHEN n_targets>=1 THEN 1 ELSE 0 END) AS conserved_any,
       SUM(CASE WHEN n_targets=2 THEN 1 ELSE 0 END) AS conserved_2,
       SUM(CASE WHEN n_targets=3 THEN 1 ELSE 0 END) AS conserved_3,
       SUM(CASE WHEN n_targets=4 THEN 1 ELSE 0 END) AS conserved_all4,
       COUNT(*) AS total,
       ROUND(100.0*SUM(CASE WHEN n_targets>=1 THEN 1 ELSE 0 END)/COUNT(*),2) AS pct_any
FROM U GROUP BY tissue ORDER BY tissue;
"
echo "[done] Human-as-hub conservation complete. Outputs in: $OUT"
