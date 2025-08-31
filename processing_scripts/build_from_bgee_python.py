import os, tarfile, gzip, csv, sys
from pathlib import Path
import pandas as pd

SPECIES = os.environ.get("BGEE_SPECIES","human,mouse,macaque,pig,chicken").split(",")
SPECIES = [s.strip().lower() for s in SPECIES if s.strip()]
ROOT_IN  = Path("data/raw/bgee")
ROOT_OUT = Path("data/expression")

NEEDED = {
  "geneId":   ["Gene ID","EnsemblGeneId","GeneID","gene_id","ensembl_gene_id","bgeegeneid"],
  "libraryId":["Library ID","RNA-seq library ID","rnaSeqLibraryId","Sample ID","library_id","sample_id"],
  "tpm":      ["TPM","tpm","TPM value","expression TPM"]
}
OPTIONAL = {
  "anatId":   ["Anatomical entity ID","Uberon ID","anat_id","anatomical_entity_id"],
  "anatName": ["Anatomical entity name","Tissue","Organ","anatomical_entity_name","tissue","organ"],
  "stageId":  ["Stage ID","Developmental stage ID","stage_id","dev_stage_id"],
  "stage":    ["Stage name","Stage","Developmental stage","stage_name","dev_stage_name"],
  "sex":      ["Sex","sex"],
  "strain":   ["Strain","strain"],
  "platform": ["Library type","Platform","Library platform","library_platform","platform","library type"]
}

def norm(s:str)->str:
  return ''.join(c.lower() if c.isalnum() else '_' for c in (s or "")).strip('_')

def pick_indices(header, wanted):
  idx = {}; H  = [h.strip() for h in header]; Hn = [norm(h) for h in H]
  for out, options in wanted.items():
    j = None
    for opt in options:
      for i,h in enumerate(H):
        if h.lower() == (opt or "").lower():
          j = i; break
      if j is not None: break
    if j is None:
      On = [norm(opt) for opt in options]
      for i,hh in enumerate(Hn):
        if hh in On:
          j = i; break
    idx[out] = j
  return idx

def ensure_extracted(in_dir: Path):
  if list(in_dir.rglob("*.tsv")) or list(in_dir.rglob("*.tsv.gz")):
    return
  tars = list(in_dir.glob("*.tar")) + list(in_dir.glob("*.tar.gz"))
  if not tars: return
  tarpath = tars[0]
  print(f"   Extracting: {tarpath.name}")
  with tarfile.open(tarpath, "r:*") as tf:
    tf.extractall(in_dir)

def iter_rows_from_file(p: Path):
  opener = gzip.open if p.suffix == ".gz" else open
  with opener(p, 'rt', encoding='utf-8', errors='ignore') as fh:
    r = csv.reader(fh, delimiter='\t')
    try: header = next(r)
    except StopIteration: return [], []
    idx = pick_indices(header, {**NEEDED, **OPTIONAL})
    if None in (idx["geneId"], idx["libraryId"], idx["tpm"]): return [], []
    expr = []; meta_by_lib = {}
    for row in r:
      try:
        gene = row[idx["geneId"]].strip(); lib  = row[idx["libraryId"]].strip(); tpmv = row[idx["tpm"]].strip()
      except Exception:
        continue
      if not gene or not lib: continue
      # fast numeric parse
      try:
        tpm_num = float(tpmv)
      except Exception:
        continue
      expr.append((gene, lib, tpm_num))
      if lib not in meta_by_lib:
        def get(k):
          j = idx.get(k)
          if j is None or j >= len(row): return None
          v = (row[j] or "").strip()
          return v or None
        meta_by_lib[lib] = {
          "libraryId": lib,
          "anatId":   get("anatId"),
          "anatName": get("anatName"),
          "stageId":  get("stageId"),
          "stage":    get("stage"),
          "sex":      get("sex"),
          "strain":   get("strain"),
          "platform": get("platform"),
        }
    return expr, list(meta_by_lib.values())

def process_species(sp: str):
  print(f"== {sp} ==")
  in_dir  = ROOT_IN / sp / "RNA_SEQ"
  out_dir = ROOT_OUT / sp
  out_dir.mkdir(parents=True, exist_ok=True)

  ensure_extracted(in_dir)

  files = [p for p in in_dir.rglob("*") if (p.suffix == ".tsv" or p.suffixes[-2:] == [".tsv",".gz"])]
  if not files:
    print(f"   No TSV files found under {in_dir}")
    return 1

  expr_chunks = []; meta_chunks = []
  for p in files:
    e, m = iter_rows_from_file(p)
    if e: expr_chunks.extend(e)
    if m: meta_chunks.extend(m)

  if not expr_chunks:
    print(f"   Parsed 0 expression rows for {sp}.")
    return 2

  expr = pd.DataFrame(expr_chunks, columns=["geneId","libraryId","tpm"])
  agg  = expr.groupby(["geneId","libraryId"], as_index=False)["tpm"].mean()
  mat  = agg.pivot(index="geneId", columns="libraryId", values="tpm").fillna(0).sort_index()

  cols = ["libraryId","anatId","anatName","stageId","stage","sex","strain","platform"]
  meta = (pd.DataFrame(meta_chunks).drop_duplicates().sort_values("libraryId")
          if meta_chunks else pd.DataFrame(columns=cols))

  mat.to_csv(out_dir / "tpm_matrix.csv")
  meta.to_csv(out_dir / "sample_metadata.csv", index=False)
  print(f"   Wrote: {out_dir/'tpm_matrix.csv'}")
  print(f"   Wrote: {out_dir/'sample_metadata.csv'}")
  return 0

def main():
  codes = [process_species(sp) for sp in SPECIES]
  sys.exit(0 if all(c==0 for c in codes) else 1)

if __name__ == "__main__":
  main()
