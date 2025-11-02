#!/usr/bin/env bash
set -euo pipefail

cd ~/Desktop/"Improved Project Alpha"

OUT="data/processed/chicken/rnaseq/heart"
GTF="data/biotype_gtf_fresh/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.115.gtf"
FA="data/ref/chicken/fasta/Gallus_gallus.GRCg7b.fa"
CHROMS="data/ref/chicken/galGal7b.chrom.sizes"

# tools
command -v bedtools >/dev/null || { echo "ERROR: bedtools not found"; exit 1; }
command -v samtools >/dev/null || { echo "ERROR: samtools not found"; exit 1; }

# chrom sizes
if [[ ! -s "$CHROMS" ]]; then
  [[ -s "${FA}.fai" ]] || samtools faidx "$FA"
  cut -f1,2 "${FA}.fai" > "$CHROMS"
fi

mkdir -p "$OUT/cov" "data/processed/chicken/enhancers"

echo "[1/5] Per-sample RPM coverage…"
for bam in "$OUT"/bam/*.sorted.bam; do
  base=$(basename "$bam" .sorted.bam)
  mapped=$(samtools idxstats "$bam" | awk '{m+=$3} END{print (m>0?m:1)}')
  scale=$(awk -v m="$mapped" 'BEGIN{printf "%.10f", 1000000.0/m}')
  bedtools genomecov -bg -split -scale "$scale" -ibam "$bam" > "$OUT/cov/${base}.rpm.bedgraph"
done

echo "[2/5] Merge bedGraphs with max RPM…"
cat "$OUT"/cov/*.rpm.bedgraph \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - -c 4 -o max > "$OUT/cov/heart_merged_max.rpm.bed"

echo "[3/5] Keep high-coverage tiles…"
THRESH=5
awk -v t="$THRESH" 'BEGIN{OFS="\t"} $4>=t {print $1,$2,$3}' "$OUT/cov/heart_merged_max.rpm.bed" \
  | bedtools sort -g "$CHROMS" \
  | bedtools merge -i - > "$OUT/cov/heart_highcov.bed"

echo "[4/5] Build promoter and exon masks…"
zcat -f "$GTF" \
 | awk 'BEGIN{OFS="\t"} $3=="transcript"{split($9,a,/gene_id "|";/); gid=a[2]; tss=($7=="+")?$4:$5; s=tss-2000; if(s<0)s=0; e=tss+2000; print $1,s,e,gid}' \
 | bedtools sort -g "$CHROMS" | bedtools merge -i - > "$OUT/mask_promoters_2kb.bed"

zcat -f "$GTF" \
 | awk 'BEGIN{OFS="\t"} $3=="exon"{print $1,$4-1,$5}' \
 | bedtools sort -g "$CHROMS" | bedtools merge -i - > "$OUT/mask_exons.bed"

echo "[5/5] Subtract masks → enhancer proxy…"
bedtools subtract -a "$OUT/cov/heart_highcov.bed" -b "$OUT/mask_promoters_2kb.bed" \
| bedtools subtract -a - -b "$OUT/mask_exons.bed" \
| awk 'BEGIN{OFS="\t"} {print $1,$2,$3,".",0,"."}' \
> data/processed/chicken/enhancers/chicken_heart_RNAproxy.bed

echo "Wrote data/processed/chicken/enhancers/chicken_heart_RNAproxy.bed"

