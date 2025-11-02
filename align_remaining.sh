#!/bin/bash
# align_remaining.sh — align remaining RNA-seq samples for chicken heart

IDX="data/ref/chicken/hisat2_index/galGal7b"
RAW="data/processed/chicken/rnaseq/heart"
OUT="$RAW"
mkdir -p "$OUT/bam"

for SRR in SRR28602663 SRR28602664; do
  echo "=== Aligning $SRR ==="
  hisat2 -p 8 -x "$IDX" \
    -1 "$RAW/${SRR}_1.fastq.gz" \
    -2 "$RAW/${SRR}_2.fastq.gz" \
  | samtools sort -@ 8 -o "$OUT/bam/${SRR}.sorted.bam" -
  samtools index "$OUT/bam/${SRR}.sorted.bam"
  echo "=== Done $SRR ==="
done

