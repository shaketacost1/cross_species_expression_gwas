#!/bin/bash
# stringtie_batch.sh — quantify expression for all chicken heart samples

OUT="data/processed/chicken/rnaseq/heart"
GTF="data/biotype_gtf_fresh/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.115.gtf"
LOG="data/processed/chicken/rnaseq/heart/stringtie_batch.log"

mkdir -p "$OUT/stringtie"
echo "=== StringTie batch run started $(date) ===" > "$LOG"

for BAM in "$OUT"/bam/*.sorted.bam; do
  SRR=$(basename "$BAM" .sorted.bam)
  TSV="$OUT/stringtie/${SRR}.gene_abund.tsv"
  GTF_OUT="$OUT/stringtie/${SRR}.gtf"

  if [[ -s "$TSV" && -s "$GTF_OUT" ]]; then
    echo "[$(date)] Skipping $SRR — already quantified" | tee -a "$LOG"
    continue
  fi

  echo "[$(date)] Quantifying $SRR ..." | tee -a "$LOG"

  stringtie -v "$BAM" -p 8 -G "$GTF" -e \
    -A "$TSV" -o "$GTF_OUT"

  if [[ $? -eq 0 ]]; then
    echo "[$(date)] ✅ Done $SRR" | tee -a "$LOG"
  else
    echo "[$(date)] ❌ Failed $SRR" | tee -a "$LOG"
  fi
done

echo "=== StringTie batch completed $(date) ===" >> "$LOG"

