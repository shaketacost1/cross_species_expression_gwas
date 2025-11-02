#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")"   # go to the script's directory

echo "[mouse enhancers]"
find data/raw/mouse -type f \
  \( -iname "*enhancer*.bed" -o -iname "*enhancer*.bed.gz" -o -iname "*H3K27ac*.bed" -o -iname "*H3K27ac*.bed.gz" \) \
  | sort | tee /tmp/mouse_enhancers.list

echo "[mouse CTCF]"
find data/raw/mouse -type f \
  \( -iname "*ctcf*.bed" -o -iname "*ctcf*.bed.gz" -o -iname "*ctcf*.narrowpeak" \) \
  | sort | tee /tmp/mouse_ctcf.list

# Build merged enhancers for each tissue
mkdir -p data/processed/mouse/enhancers data/processed/mouse/ctcf

for tissue in liver heart brain; do
    echo "[processing $tissue]"
    find "data/raw/mouse/enhancers/$tissue" -type f -name "*.bed.gz" -print0 \
      | xargs -0 gzcat \
      | awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' \
      | LC_ALL=C sort -k1,1 -k2,2n \
      | bedtools merge -i - \
      > "data/processed/mouse/enhancers/mouse_${tissue}_H3K27ac.bed"
done

# Load into SQLite
for tissue in brain liver heart; do
    python scripts/etl/bed_to_sqlite.py --db ipa.db --species mouse --type enhancer \
      --bed "data/processed/mouse/enhancers/mouse_${tissue}_H3K27ac.bed" \
      --source ENCODE --method "H3K27ac:${tissue}"
done

# Example for CTCF (brain only here)
gunzip -c data/raw/mouse/ctcf/brain/mouse_brain_CTCF_union.bed.gz \
  > data/processed/mouse/ctcf/mouse_brain_CTCF.bed

python scripts/etl/bed_to_sqlite.py --db ipa.db --species mouse --type ctcf \
  --bed data/processed/mouse/ctcf/mouse_brain_CTCF.bed \
  --source ENCODE --method "CTCF:brain"

