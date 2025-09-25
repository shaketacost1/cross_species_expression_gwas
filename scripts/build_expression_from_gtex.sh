#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$ROOT/data"
TMP="$OUT/_tmp_gtex_median.gct.gz"
GCT="$OUT/GTEx_v8_gene_median_tpm.gct.gz"

mkdir -p "$OUT"
rm -f "$TMP" "$GCT"

# Try multiple mirrors; stop at first valid gzip
while read -r url; do
  echo "Trying: $url"
  if curl --fail --location --retry 3 --connect-timeout 20 -o "$TMP" "$url"; then
    if gzip -t "$TMP" 2>/dev/null; then
      mv "$TMP" "$GCT"
      echo "OK: $GCT"
      break
    else
      echo "Not a valid gzip, trying next…"
      rm -f "$TMP"
    fi
  else
    echo "curl failed, trying next…"
  fi
done <<'URLS'
https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz
https://gtexportal.org/static/datasets/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz
https://ftp.broadinstitute.org/pub/gtex/gtx_rnaseq_v8/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz
URLS

if [[ ! -s "$GCT" ]]; then
  echo "Download failed from all mirrors." >&2
  exit 1
fi

# Awk program to convert GCT medians -> symbol,tissue,tpm (Brain/Liver/Heart)
cat > "$OUT/make_expr.awk" <<'AWK'
BEGIN{FS=OFS="\t"}
NR==1{next}                     # GCT version line
NR==2{next}                     # dims line
NR==3{
  for(i=1;i<=NF;i++){
    if($i ~ /^Brain - /){ nb++; brain[nb]=i }
    if($i=="Liver"){ liver=i }
    if($i=="Heart - Atrial Appendage"){ haa=i }
    if($i=="Heart - Left Ventricle"){   hlv=i }
  }
  print "symbol","tissue","tpm"
  next
}
{
  sym=$2                         # Description column = HGNC symbol
  if(sym=="" || sym=="NA") next

  bsum=0; bc=0
  for(k=1;k<=nb;k++){ v=$(brain[k])+0; bsum+=v; bc++ }
  if(bc>0) printf "%s\tBrain\t%.6f\n", sym, (bsum/bc)

  if(liver) printf "%s\tLiver\t%.6f\n", sym, ($(liver)+0)

  hsum=0; hc=0
  if(haa){ hsum+=($(haa)+0); hc++ }
  if(hlv){ hsum+=($(hlv)+0); hc++ }
  if(hc>0) printf "%s\tHeart\t%.6f\n", sym, (hsum/hc)
}
AWK

# Convert to the exact file the app expects
gzip -cd "$GCT" | awk -f "$OUT/make_expr.awk" > "$OUT/expression_tpm.tsv"

# Show quick checks
echo "Lines in data/expression_tpm.tsv:"
wc -l "$OUT/expression_tpm.tsv"
echo "Head:"
head -5 "$OUT/expression_tpm.tsv"
echo "Sample genes (BDNF, PCSK9, TTN):"
grep -E "^(BDNF|PCSK9|TTN)\t" "$OUT/expression_tpm.tsv" | head || true
