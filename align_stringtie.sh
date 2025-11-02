#!/usr/bin/env bash
set -euo pipefail

cd ~/Desktop/"Improved Project Alpha"

IDX="data/ref/chicken/hisat2_index/galGal7b"
GTF="data/biotype_gtf_fresh/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.115.gtf.gz"
OUT="data/processed/chicken/rnaseq/heart"
RAW="$OUT"   # FASTQs are here

mkdir -p "$OUT/bam" "$OUT/stringtie"

STRAND_HISAT2=""
STRAND_ST=""

while read -r SRR; do
  [[ -z "$SRR" ]] && continue
  echo "=== Aligning $SRR ==="

  # sanity check inputs
  [[ -f "$RAW/${SRR}_1.fastq.gz" && -f "$RAW/${SRR}_2.fastq.gz" ]] || { echo "FASTQs missing for $SRR"; exit 1; }

  hisat2 -p 8 -x "$IDX" $STRAND_HISAT2 \
    -1 "$RAW/${SRR}_1.fastq.gz" \
    -2 "$RAW/${SRR}_2.fastq.gz" |
  samtools sort -@ 8 -o "$OUT/bam/${SRR}.sorted.bam" -

  samtools index "$OUT/bam/${SRR}.sorted.bam"

  stringtie "$OUT/bam/${SRR}.sorted.bam" -p 8 -G "$GTF" -e $STRAND_ST \
    -A "$OUT/stringtie/${SRR}.gene_abund.tsv" \
    -o "$OUT/stringtie/${SRR}.gtf"

  echo "=== Done $SRR ==="
done < "$OUT/run_accessions.txt"
