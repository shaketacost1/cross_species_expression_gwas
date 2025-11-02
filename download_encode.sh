#!/usr/bin/env bash
set -euo pipefail

enc() {
  local query="$1"   # experiment-level filters
  local outdir="$2"
  mkdir -p "$outdir"
  echo "Query(experiments): $query"
  echo "Out: $outdir"

  # Pull experiments (human only) with embedded file objects
  exp_url="https://www.encodeproject.org/search/?type=Experiment&status=released&limit=all&frame=embedded&replicates.library.biosample.donor.organism.scientific_name=Homo%20sapiens&${query}"
  exp_json="$(curl -s -H 'Accept: application/json' "$exp_url")"

  exp_count="$(printf '%s' "$exp_json" | jq -r '.["@graph"] | length')"
  echo "Experiments found: $exp_count"
  if [ "$exp_count" -eq 0 ]; then
    echo "No experiments for: $query"
    return 0
  fi

  # Collect downloadable PEAK files in BED/bigBed for GRCh38
  hrefs="$(printf '%s' "$exp_json" | jq -r '
      .["@graph"][]?.files[]?
      | select(.status=="released")
      | select(.href and (.no_file_available|not))
      | select(.assembly=="GRCh38")
      | select((.output_category//"") == "peaks" or (.output_type//"" | test("peaks"; "i")))
      | select((.file_format//"") == "bed" or (.file_format//"") == "bigBed")
      | .href
    ' | sort -u)"

  echo "Files matched: $(printf '%s\n' "$hrefs" | sed '/^$/d' | wc -l | tr -d ' ')"

  if [ -z "$hrefs" ]; then
    echo "No matching files after filtering (try loosening filters)."
    return 0
  fi

  printf '%s\n' "$hrefs" | while read -r href; do
    [ -z "${href:-}" ] && continue
    fname="$(basename "$href")"
    url="https://www.encodeproject.org${href}?download=true"
    echo "→ Downloading $fname"
    curl -L -s -o "${outdir}/${fname}" "$url"
  done

  echo "Downloaded to $outdir"
}

# ---- H3K27ac (enhancers proxy) into separate subfolders
enc "assay_title=Histone%20ChIP-seq&target.label=H3K27ac&biosample_ontology.term_name=liver"                          data/raw/human/enhancers/liver
enc "assay_title=Histone%20ChIP-seq&target.label=H3K27ac&biosample_ontology.term_name=heart%20left%20ventricle"       data/raw/human/enhancers/heart_lv
enc "assay_title=Histone%20ChIP-seq&target.label=H3K27ac&biosample_ontology.organ_slims=brain"                         data/raw/human/enhancers/brain

# ---- CTCF (TF ChIP-seq) into separate subfolders
enc "assay_title=TF%20ChIP-seq&target.label=CTCF&biosample_ontology.term_name=liver"                                   data/raw/human/ctcf/liver
enc "assay_title=TF%20ChIP-seq&target.label=CTCF&biosample_ontology.term_name=heart%20left%20ventricle"                data/raw/human/ctcf/heart_lv
enc "assay_title=TF%20ChIP-seq&target.label=CTCF&biosample_ontology.organ_slims=brain"                                 data/raw/human/ctcf/brain

# =========================
# Mouse (Mus musculus, mm10)
# =========================
enc_mouse() {
  local query="$1"   # experiment-level filters
  local outdir="$2"
  mkdir -p "$outdir"
  echo "Query(experiments, mouse): $query"
  echo "Out: $outdir"

  # Experiments (mouse) with embedded files
  exp_url="https://www.encodeproject.org/search/?type=Experiment&status=released&limit=all&frame=embedded&replicates.library.biosample.donor.organism.scientific_name=Mus%20musculus&${query}"
  exp_json="$(curl -s -H 'Accept: application/json' "$exp_url")"

  exp_count="$(printf '%s' "$exp_json" | jq -r '.["@graph"] | length')"
  echo "Experiments found: $exp_count"
  [ "$exp_count" -eq 0 ] && { echo "No experiments for: $query"; return 0; }

  # Collect downloadable peak files (mm10) in BED/bigBed
  hrefs="$(printf '%s' "$exp_json" | jq -r '
      .["@graph"][]?.files[]?
      | select(.status=="released")
      | select(.href and (.no_file_available|not))
      | select(.assembly=="mm10")
      | select((.output_category//"") == "peaks" or (.output_type//"" | test("peaks"; "i")))
      | select((.file_format//"") == "bed" or (.file_format//"") == "bigBed")
      | .href
    ' | sort -u)"

  echo "Files matched: $(printf '%s\n' "$hrefs" | sed '/^$/d' | wc -l | tr -d ' ')"
  [ -z "$hrefs" ] && { echo "No matching files after filtering."; return 0; }

  printf '%s\n' "$hrefs" | while read -r href; do
    [ -z "${href:-}" ] && continue
    fname="$(basename "$href")"
    url="https://www.encodeproject.org${href}?download=true"
    echo "→ Downloading $fname"
    curl -L -s -o "${outdir}/${fname}" "$url"
  done

  echo "Downloaded to $outdir"
}

# --- Mouse H3K27ac (enhancers proxy) ---
enc_mouse "assay_title=Histone%20ChIP-seq&target.label=H3K27ac&biosample_ontology.organ_slims=liver" data/raw/mouse/enhancers/liver
enc_mouse "assay_title=Histone%20ChIP-seq&target.label=H3K27ac&biosample_ontology.organ_slims=heart" data/raw/mouse/enhancers/heart
enc_mouse "assay_title=Histone%20ChIP-seq&target.label=H3K27ac&biosample_ontology.organ_slims=brain" data/raw/mouse/enhancers/brain

# --- Mouse CTCF ---
enc_mouse "assay_title=TF%20ChIP-seq&target.label=CTCF&biosample_ontology.organ_slims=liver" data/raw/mouse/ctcf/liver
enc_mouse "assay_title=TF%20ChIP-seq&target.label=CTCF&biosample_ontology.organ_slims=heart" data/raw/mouse/ctcf/heart
enc_mouse "assay_title=TF%20ChIP-seq&target.label=CTCF&biosample_ontology.organ_slims=brain" data/raw/mouse/ctcf/brain
