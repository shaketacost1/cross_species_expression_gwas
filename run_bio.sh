#!/bin/bash
set -euo pipefail

# where to store fresh biotype files
OUTDIR="data/biotype_gtf_fresh"
mkdir -p "$OUTDIR"

echo "[1/5] Human"
curl -L -o "$OUTDIR/Homo_sapiens.GRCh38.115.gtf.gz" \
  https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz

echo "[2/5] Mouse"
curl -L -o "$OUTDIR/Mus_musculus.GRCm39.115.gtf.gz" \
  https://ftp.ensembl.org/pub/release-115/gtf/mus_musculus/Mus_musculus.GRCm39.115.gtf.gz

echo "[3/5] Chicken"
curl -L -o "$OUTDIR/Gallus_gallus.GRCg7b.115.gtf.gz" \
  https://ftp.ensembl.org/pub/release-115/gtf/gallus_gallus/Gallus_gallus.GRCg7b.115.gtf.gz

echo "[4/5] Pig"
curl -L -o "$OUTDIR/Sus_scrofa.Sscrofa11.1.115.gtf.gz" \
  https://ftp.ensembl.org/pub/release-115/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.115.gtf.gz

echo "[5/5] Dog"
curl -L -o "$OUTDIR/Canis_lupus_familiaris.ROS_Cfam_1.0.115.gtf.gz" \
  https://ftp.ensembl.org/pub/release-115/gtf/canis_lupus_familiaris/Canis_lupus_familiaris.ROS_Cfam_1.0.115.gtf.gz

echo "[done] Biotype GTFs saved in $OUTDIR"

