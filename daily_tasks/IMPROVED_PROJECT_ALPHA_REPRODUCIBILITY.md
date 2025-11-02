
# Improved Project Alpha — End‑to‑End Reproducibility Guide

This single document merges the steps, commands, and data provenance from all provided files into one cohesive, executable-style reference. It is organized by phases (environment, inputs, processing, database loads, liftover/classification, GWAS, and QC/outputs). All paths are relative to the project root unless noted.

---

## 0) Environment & Conventions

**OS / Shell**
- macOS (tested), POSIX shell (bash/zsh). fileciteturn0file7

**Core tools (install via conda/homebrew as needed)**
- `sqlite3`, `bedtools`, `liftOver`, `awk`, `curl`, `bgzip`, `tabix`, `gzip`, UCSC `bigBedToBed`
- Python 3.11+ (some steps used Python 3.12) with `virtualenv`/`conda` as preferred
- Peak calling: `macs3`
- Aligners (for raw ChIP pipelines): `bowtie2`; BAM handling: `samtools`  fileciteturn0file5 fileciteturn0file7

**Python envs (examples)**
```bash
conda create -n macs3 python=3.12 -y
conda activate macs3
pip install macs3
```
fileciteturn0file5

**Working directories**
```
Improved Project Alpha/
  data/
    raw/
    processed/
    reference/
    tmp/
  scripts/
  ipa.db   # SQLite
```
Key subfolders include chain files (`data/reference/chains`), intermediate BEDs (`data/tmp`), processed outputs and cross‑species summaries (`data/processed`, `data/processed/conservation`, `data/processed/reports`). fileciteturn0file7 fileciteturn0file0

---

## 1) Directory Structure (seed)

A practical seed layout reflecting current usage:
```text
Improved Project Alpha/
└── data/
    ├── raw/
    │   ├── dog/
    │   │   ├── enhancers/         # GSE203104 dog H3K27ac nPeak sets
    │   │   └── ctcf/              # dog CTCF peaks (ENA PRJEB2329) and conversions
    │   ├── mouse/
    │   │   └── ctcf/brain/        # bigBed + converted bed.gz; union files
    │   └── chicken/
    │       ├── enhancers/         # brain/liver H3K27ac; heart pending/proxy
    │       └── ctcf/              # native/consensus sets
    ├── processed/
    │   ├── {species}/enhancers/
    │   ├── {species}/ctcf/
    │   └── conservation/
    └── reference/
        ├── chains/
        └── annotations/
```
Populate as data lands (see steps below). fileciteturn0file0 fileciteturn0file6

---

## 2) Data Sources & Provenance

### 2.1 Dog enhancers (H3K27ac)
- Source: GEO **GSE203104** (k27ac narrowPeak sets).  
  Download & extract:
```bash
cd data/raw/dog/enhancers
curl -L -o GSE203104_k27ac.narrowpeak.tar.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE203nnn/GSE203104/suppl/GSE203104_k27ac.narrowpeak.tar.gz"
tar -xvf GSE203104_k27ac.narrowpeak.tar.gz
```
Results: tissue‑specific `*_k27ac_optimal.narrowPeak`. fileciteturn0file0 fileciteturn0file1

### 2.2 Dog CTCF (ChIP‑seq raw → consensus)
- Source: **ENA PRJEB2329** (ERR022286/ERR022304 ChIP; ERR022270/ERR022271 control).  
  Summary of pipeline: FASTQ → Bowtie2 → sorted/indexed BAMs (samtools) → `macs3 callpeak` per replicate → `bedtools` intersection to consensus → merge nearby → `dog_ctcf_consensus_merged.bed`.  
  Outputs under `data/processed/dog/ctcf/`. fileciteturn0file5

### 2.3 Mouse CTCF (Brain)
- Input format: `.bigBed` converted to `.bed.gz` then merged to union.  
```bash
for f in data/raw/mouse/ctcf/brain/*.bigBed; do
  base=$(basename "$f" .bigBed)
  bigBedToBed "$f" stdout | gzip -c > "data/raw/mouse/ctcf/brain/${base}.bed.gz"
done

DIR="data/raw/mouse/ctcf/brain"
gzcat "$DIR"/*.bed.gz | awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' | LC_ALL=C sort -k1,1 -k2,2n | bedtools merge -i - > "$DIR/mouse_brain_CTCF_union.bed"
gzip -f "$DIR/mouse_brain_CTCF_union.bed"
```
Output: `mouse_brain_CTCF_union.bed.gz`. fileciteturn0file0

### 2.4 Chicken enhancers & CTCF
- Enhancers present for brain/liver; heart pending; CTCF sets available. If H3K27ac heart missing, use ATAC or RNA‑derived proxy (see §4.3). fileciteturn0file1 fileciteturn0file6

### 2.5 Pig enhancers & CTCF
- H3K27ac enhancers and CTCF present/ingested; see loader commands below. fileciteturn0file3

### 2.6 Human enhancers & CTCF (reference)
- Human H3K27ac and CTCF are loaded and used for liftover/conservation & GWAS overlaps. fileciteturn0file3 fileciteturn0file4

---

## 3) Database & Annotations

### 3.1 Schema & tables
Database file: `ipa.db`. Core tables include: `enhancers`, `ctcf_sites`, `genes`, `species`, with derived `*_flags`, `*_annotated` views. Example counts by species/biotype were validated. fileciteturn0file4

### 3.2 Load genes/biotypes (per species)
Download GTFs into `reference/annotations/` then load:
```bash
python scripts/etl/gtf_to_genes.py --db ipa.db --species human   --gtf reference/annotations/gencode.vXX.annotation.gtf.gz
python scripts/etl/gtf_to_genes.py --db ipa.db --species mouse   --gtf reference/annotations/gencode.vMXX.annotation.gtf.gz
python scripts/etl/gtf_to_genes.py --db ipa.db --species pig     --gtf reference/annotations/Ensembl_pig.gtf.gz
python scripts/etl/gtf_to_genes.py --db ipa.db --species chicken --gtf reference/annotations/Ensembl_chicken.gtf.gz
python scripts/etl/gtf_to_genes.py --db ipa.db --species dog     --gtf reference/annotations/Ensembl_dog.gtf.gz
```
Verification:
```sql
SELECT s.common_name, COUNT(*) genes, COUNT(DISTINCT biotype) biotypes
FROM genes g JOIN species s USING(species_id) GROUP BY 1;
```
fileciteturn0file1 fileciteturn0file3

---

## 4) Build Enhancer & CTCF Sets (by species)

### 4.1 Dog
**Enhancers (H3K27ac)**: from GEO GSE203104 (see §2.1).  
**CTCF**: replicate peak calling and consensus build (see §2.2).  
Load CTCF consensus to DB:
```bash
python scripts/etl/bed_to_sqlite.py --db ipa.db --species dog --type ctcf   --bed data/processed/dog/ctcf/dog_ctcf_consensus_merged.bed   --source ENA_PRJEB2329 --method MACS3_consensus
```
fileciteturn0file3

### 4.2 Mouse
**CTCF**: load per tissue unions:
```bash
python scripts/etl/bed_to_sqlite.py --db ipa.db --species mouse --type ctcf   --bed data/raw/mouse/ctcf/brain/mouse_brain_CTCF_union.bed.gz --method union
python scripts/etl/bed_to_sqlite.py --db ipa.db --species mouse --type ctcf   --bed data/raw/mouse/ctcf/heart/mouse_heart_CTCF_union.bed.gz --method union
python scripts/etl/bed_to_sqlite.py --db ipa.db --species mouse --type ctcf   --bed data/raw/mouse/ctcf/liver/mouse_liver_CTCF_union.bed.gz --method union
```
**Enhancers (H3K27ac)**: merge liver/heart and load:
```bash
# Liver
find data/raw/mouse/enhancers/liver -type f -name '*.bed.gz' -print0 |   xargs -0 gzcat | awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' |   LC_ALL=C sort -k1,1 -k2,2n | bedtools merge -i - > data/processed/mouse/enhancers/mouse_liver_H3K27ac.bed
python scripts/etl/bed_to_sqlite.py --db ipa.db --species mouse --type enhancer   --bed data/processed/mouse/enhancers/mouse_liver_H3K27ac.bed   --source ENCODE --method 'H3K27ac:liver'

# Heart
find data/raw/mouse/enhancers/heart -type f -name '*.bed.gz' -print0 |   xargs -0 gzcat | awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' |   LC_ALL=C sort -k1,1 -k2,2n | bedtools merge -i - > data/processed/mouse/enhancers/mouse_heart_H3K27ac.bed
python scripts/etl/bed_to_sqlite.py --db ipa.db --species mouse --type enhancer   --bed data/processed/mouse/enhancers/mouse_heart_H3K27ac.bed   --source ENCODE --method 'H3K27ac:heart'
```
fileciteturn0file3

### 4.3 Chicken
If heart H3K27ac missing, construct prioritized heart set using RNA‑proxy and CTCF co‑evidence:
```bash
# RNA‑proxy clean set (example merge step)
bedtools merge -i - -d 150 > chicken_heart_RNAproxy.clean.bed

# CTCF ±1 kb co‑evidence
bedtools window -w 1000   -a chicken_heart_RNAproxy.clean.chr.bed   -b chicken_ctcf_all.chr.bed   > chicken_heart_RNAproxy.ctdf1kb.bed  # (typo fixed to 'ctcf1kb' in downstream)

# Load prioritized view (H3K27ac > RNAproxy+CTCF > RNAproxy)
# See SQL view definition 'chicken_heart_final' in the source document.
```
Materialize final prioritized set and export BED:
```bash
sqlite3 -separator $'\t' ipa.db "SELECT chrom, start, end, 'chicken_heart_final', 0, '.' FROM chicken_heart_final ORDER BY chrom, start" > data/processed/chicken/enhancers/chicken_heart_final.bed
awk 'BEGIN{OFS="\t"}{if($1!~/^chr/)$1="chr"$1;print}' chicken_heart_final.bed | bedtools sort -i - > data/processed/chicken/enhancers/chicken_heart_final.chr.bed
```
fileciteturn0file6 fileciteturn0file1

### 4.4 Pig
Extract/convert, merge, then load (examples):
```bash
# Enhancers (H3K27ac)
python scripts/etl/bed_to_sqlite.py --db ipa.db --species pig --type enhancer   --bed data/raw/pig/enhancers/pig_liver_H3K27ac.bed.gz --source FAANG --method H3K27ac
python scripts/etl/bed_to_sqlite.py --db ipa.db --species pig --type enhancer   --bed data/raw/pig/enhancers/pig_brain_H3K27ac.bed.gz --source FAANG --method H3K27ac
python scripts/etl/bed_to_sqlite.py --db ipa.db --species pig --type enhancer   --bed data/raw/pig/enhancers/pig_heart_H3K27ac.bed.gz --source FAANG --method H3K27ac

# CTCF consensus (example brain)
DIR=data/raw/pig/ctcf/brain
gzcat "$DIR"/*.bed.gz | awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' | sort -k1,1 -k2,2n | bedtools merge -i - > "$DIR/pig_brain_CTCF_union.bed"
gzip -f "$DIR/pig_brain_CTCF_union.bed"
python scripts/etl/bed_to_sqlite.py --db ipa.db --species pig --type ctcf   --bed data/raw/pig/ctcf/pig_brain_CTCF_union.bed.gz --method union
```
fileciteturn0file1

### 4.5 Human
Human enhancer sets (brain/heart/liver) and CTCF are pre‑loaded and serve as anchors for liftover and GWAS integration steps. fileciteturn0file3 fileciteturn0file4

---

## 5) CTCF Proximity (human & mouse example)

Compute nearest CTCF distance per enhancer and summarize multi‑window proximity in SQLite; store row‑level in `ctcf_flags` and roll‑ups in `ctcf_proximity_summary`. Empirically, ≥95–99% enhancers fall within 1kb of CTCF. fileciteturn0file4

Example (conceptual workflow):
```bash
# Export CTCF sites to BED3 (if not already)
# data/tmp/human_CTCF.bed, data/tmp/mouse_CTCF.bed

# For each tissue, compute closest
bedtools closest -d -a human_brain_enhancers.bed -b data/tmp/human_CTCF.bed > data/processed/ctcf_proximity/human_brain_ctcf_closest.tsv

# Load into SQLite and aggregate by windows (100/250/500/1000 bp) → ctcf_proximity_summary
```

---

## 6) Liftover & Cross‑Species Conservation

### 6.1 Chain files
Fetch UCSC chain files into `data/reference/chains/`:
```bash
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/galGal6/liftOver/galGal6ToMm10.over.chain.gz
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/galGal6/liftOver/galGal6ToHg38.over.chain.gz
# Additional: hg38→mm10, hg38→canFam3, hg38→susScr11, etc.
```
fileciteturn0file6 fileciteturn0file4

### 6.2 Run liftover (examples)
```bash
# Human → others
python scripts/etl/liftover_enhancers.py --db ipa.db --src_species human --dst_species mouse   --chain reference/liftover_chains/hg38ToMm10.over.chain.gz
python scripts/etl/liftover_enhancers.py --db ipa.db --src_species human --dst_species pig     --chain reference/liftover_chains/hg38ToSusScr11.over.chain.gz
python scripts/etl/liftover_enhancers.py --db ipa.db --src_species human --dst_species chicken --chain reference/liftover_chains/hg38ToGalGal6.over.chain.gz
python scripts/etl/liftover_enhancers.py --db ipa.db --src_species human --dst_species dog     --chain reference/liftover_chains/hg38ToCanFam3.over.chain.gz
```
Then classify conserved/gained/lost per focal species:
```bash
python scripts/etl/classify_enhancers.py --db ipa.db --focal_species mouse
python scripts/etl/classify_enhancers.py --db ipa.db --focal_species pig
python scripts/etl/classify_enhancers.py --db ipa.db --focal_species chicken
python scripts/etl/classify_enhancers.py --db ipa.db --focal_species dog
```
fileciteturn0file1 fileciteturn0file4 fileciteturn0file7

### 6.3 Chicken heart worked example
- Build RNAproxy → add CTCF co‑evidence → prioritize final set; lift to human/mouse; intersect with matched tissues for conservation; store flags/views; export summaries. See provided SQL templates for lift success rates and reciprocal checks. fileciteturn0file6

---

## 7) GWAS Integration

### 7.1 Export SNPs to BED (hg38)
```bash
DB=ipa.db; BED=data/processed/gwas/bed; mkdir -p "$BED"
sqlite3 -noheader -separator $'\t' "$DB"  'SELECT chrom, pos-1 AS snp_start, pos AS snp_end, rsid, pval, trait_id FROM gwas_snps;'  > "$BED/gwas_snps.hg38.bed"
```
Total SNPs example: `gwas_snps.hg38.bed` with ~1.74M rows. fileciteturn0file2

### 7.2 SNP–enhancer overlaps (per tissue/set)
```bash
# Example (brain, all)
bedtools intersect -a gwas_snps.hg38.bed -b human_brain.bed -wa -wb | awk '{print ...}' > data/processed/gwas/gwas_x_enhancers_brain.tsv
```
Produce counts and trait summaries (`gwas_trait_by_tissue.csv`), per‑trait exports (e.g., CAD heart hubAny), and significance‑filtered lists (p ≤ 5e‑8). fileciteturn0file2

### 7.3 Loader helpers (post‑enhancer load)
```bash
python scripts/etl/load_gwas.py --db ipa.db --trait inflammation --tsv data/raw/human/gwas/inflammation.tsv --species human
python scripts/etl/load_gwas.py --db ipa.db --trait obesity      --tsv data/raw/human/gwas/obesity.tsv      --species human
python scripts/etl/load_gwas.py --db ipa.db --trait alcohol      --tsv data/raw/human/gwas/alcohol.tsv      --species human
```
fileciteturn0file1

---

## 8) Indices & Performance

```bash
sqlite3 ipa.db "PRAGMA journal_mode=WAL; PRAGMA synchronous=NORMAL;"
sqlite3 ipa.db "CREATE INDEX IF NOT EXISTS idx_liftover_src ON liftover_map(src_enhancer_id, dst_species_id, status);"
sqlite3 ipa.db "CREATE INDEX IF NOT EXISTS idx_snp_pos ON gwas_snps(chrom, pos);"
sqlite3 ipa.db "CREATE INDEX IF NOT EXISTS idx_enh_species ON enhancers(species_id, chrom, start, end);"
sqlite3 ipa.db "CREATE INDEX IF NOT EXISTS idx_ctcf_species ON ctcf_sites(species_id, chrom, start, end);"
sqlite3 ipa.db "CREATE INDEX IF NOT EXISTS idx_gene_symbol ON genes(gene_symbol);"
sqlite3 ipa.db "CREATE INDEX IF NOT EXISTS idx_gene_biotype ON genes(species_id, biotype);"
```
fileciteturn0file1

---

## 9) Quality Control, Summaries & Outputs

- **Species/tissue counts present** and confirmed for human, mouse, pig, chicken, dog; see status rollup and next steps. fileciteturn0file3
- **CTCF proximity**: 95–99% within 1 kb; per‑tissue TSVs under `data/processed/ctcf_proximity/`; summaries in `ctcf_proximity_summary`. fileciteturn0file4
- **Conservation summaries**: `data/processed/reports/conservation_summary.csv` with `species,tissue,n_human,pct_human,n_mouse,pct_mouse,total`. Example dog heart ≈ 28% human overlap. fileciteturn0file7
- **GWAS outputs**: overlap tables, trait summaries, and per‑trait gene/SNP exports under `data/processed/reports/`. fileciteturn0file2

---

## 10) Reproducibility Checklist (quick run order)

1. Set up env & tools (§0). fileciteturn0file7
2. Create directory scaffold (§1) and fetch inputs (§2). fileciteturn0file0
3. Load genes/biotypes (§3.2). fileciteturn0file1
4. Build/merge enhancer & CTCF sets per species (§4). fileciteturn0file0 fileciteturn0file5
5. Compute CTCF proximity (§5). fileciteturn0file4
6. Liftover & classify (§6). fileciteturn0file1
7. GWAS integration (§7). fileciteturn0file2
8. Add indices & run QC (§8–9). fileciteturn0file1 fileciteturn0file4

---

## 11) Notes & Gotchas

- Keep assemblies consistent (e.g., galGal6 vs galGal7b validation) and verify chain files via UCSC md5. fileciteturn0file7
- Maintain BED conventions (0‑based start, 1‑based end) across conversions; prefer `bedtools sort` and locale‑stable sorting (`LC_ALL=C`). fileciteturn0file0 fileciteturn0file6
- For mixed modalities (ATAC vs H3K27ac), interpret conservation percentages with caution. fileciteturn0file7
