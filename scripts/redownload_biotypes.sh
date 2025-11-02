#!/usr/bin/env bash
set -euo pipefail

RELEASE=115
OUTDIR="data/biotype_gtf_fresh"
mkdir -p "$OUTDIR"

# Species -> [ENSEMBL dir, file]
declare -A DIR FILE
DIR[human]="homo_sapiens";               FILE[human]="Homo_sapiens.GRCh38.${RELEASE}.gtf.gz"
DIR[mouse]="mus_musculus";               FILE[mouse]="Mus_musculus.GRCm39.${RELEASE}.gtf.gz"
DIR[chicken]="gallus_gallus";            FILE[chicken]="Gallus_gallus.GRCg7b.${RELEASE}.gtf.gz"
DIR[pig]="sus_scrofa";                   FILE[pig]="Sus_scrofa.Sscrofa11.1.${RELEASE}.gtf.gz"
DIR[dog]="canis_lupus_familiaris";       FILE[dog]="Canis_lupus_familiaris.ROS_Cfam_1.0.${RELEASE}.gtf.gz"

echo "[1/4] Downloading Ensembl GTFs (release ${RELEASE}) into ${OUTDIR}"
for sp in human mouse chicken pig dog; do
  url="https://ftp.ensembl.org/pub/release-${RELEASE}/gtf/${DIR[$sp]}/${FILE[$sp]}"
  dest="${OUTDIR}/${FILE[$sp]}"
  echo "  - ${sp}: $url"
  curl -L --fail --retry 3 -o "$dest" "$url"
done

echo "[2/4] Verifying gzip integrity"
for gz in "${OUTDIR}"/*.gtf.gz; do
  gzip -t "$gz"
done
echo "  OK"

echo "[3/4] Keeping .gtf.gz and also writing uncompressed .gtf for fast parsing"
for gz in "${OUTDIR}"/*.gtf.gz; do
  gtf="${gz%.gz}"
  gunzip -c "$gz" > "$gtf"
done

echo "[4/4] Quick content report (gene rows, missing attrs, top biotypes)"
for gtf in "${OUTDIR}"/*.gtf; do
  echo -e "\n== $(basename "$gtf") =="
  awk -F'\t' '$0!~/^#/ && $3=="gene"{n++} END{print "gene rows:",(n+0)}' "$gtf"

  awk -F'\t' '
    $0!~/^#/ && $3=="gene"{t++; if($9!~/"gene_name"/) mn++; if($9!~/gene_biotype|gene_type/) mb++}
    END{printf "missing gene_name: %d (of %d)\n",(mn+0),(t+0);
        printf "missing biotype  : %d (of %d)\n",(mb+0),(t+0)}' "$gtf"

  echo "top biotypes:"
  awk -F'\t' '
    $0!~/^#/ && $3=="gene"{
      if(match($9,/(gene_biotype|gene_type) "([^"]+)"/,m)) c[m[2]]++; else c["NA"]++
    }
    END{for(b in c) printf "%s\t%d\n",b,c[b]}' "$gtf" | LC_ALL=C sort -k2,2nr | head -10
done

echo -e "\n[done] Fresh GTFs in: ${OUTDIR}"
