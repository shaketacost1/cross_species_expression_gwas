#!/usr/bin/env bash
set -euo pipefail

# 1) Root
ROOT="$HOME/Desktop/Project Alpha"
cd "$ROOT"

# 2) Create folders
mkdir -p bin \
         data/chains \
         data/peaks/{human,mouse,pig} \
         data/liftover/{mouse_to_hg38,pig_to_hg38} \
         logs tmp scripts

# 3) Download chain files if missing
cd "$ROOT/data/chains"
[[ -f mm39ToHg38.over.chain.gz ]] || curl -fsSLO https://hgdownload.soe.ucsc.edu/goldenPath/mm39/liftOver/mm39ToHg38.over.chain.gz
[[ -f susScr11ToHg38.over.chain.gz ]] || curl -fsSLO https://hgdownload.soe.ucsc.edu/goldenPath/susScr11/liftOver/susScr11ToHg38.over.chain.gz

# 4) Move existing files into place (no overwrite)
cd "$ROOT"
mv -n human.qc.bed                 data/peaks/human/ 2>/dev/null || true
mv -n human.qc.sorted.bed          data/peaks/human/ 2>/dev/null || true
mv -n human.sorted.bed             data/peaks/human/ 2>/dev/null || true
mv -n human.hg38.bed               data/peaks/human/ 2>/dev/null || true

mv -n mouse.qc.bed                 data/peaks/mouse/ 2>/dev/null || true
mv -n mouse.qc.sorted.bed          data/peaks/mouse/ 2>/dev/null || true
mv -n mouse.sorted.bed             data/peaks/mouse/ 2>/dev/null || true
mv -n mouse.liftover.hg38.bed      data/liftover/mouse_to_hg38/ 2>/dev/null || true
mv -n mouse.unmapped.bed           data/liftover/mouse_to_hg38/ 2>/dev/null || true

mv -n pig.qc.bed                   data/peaks/pig/ 2>/dev/null || true
mv -n pig.qc.sorted.bed            data/peaks/pig/ 2>/dev/null || true
mv -n pig.sorted.bed               data/peaks/pig/ 2>/dev/null || true
mv -n pig.liftover.hg38.bed        data/liftover/pig_to_hg38/ 2>/dev/null || true
mv -n pig.unmapped.bed             data/liftover/pig_to_hg38/ 2>/dev/null || true

# 5) Project-local shims for tools
cd "$ROOT/bin"
if command -v /opt/homebrew/bin/liftOver >/dev/null 2>&1; then
  ln -sf /opt/homebrew/bin/liftOver liftOver
elif command -v liftOver >/dev/null 2>&1; then
  ln -sf "$(command -v liftOver)" liftOver
else
  echo "WARNING: liftOver not found. Try: brew install ucsc-liftover" >&2
fi

if command -v /opt/homebrew/bin/bedtools >/dev/null 2>&1; then
  ln -sf /opt/homebrew/bin/bedtools bedtools
elif command -v bedtools >/dev/null 2>&1; then
  ln -sf "$(command -v bedtools)" bedtools
else
  echo "WARNING: bedtools not found. Try: brew install bedtools" >&2
fi

# 6) Add project bin to PATH (once)
cd "$ROOT"
if ! grep -q '### Project Alpha PATH' "$HOME/.zshrc" 2>/dev/null; then
  cat >> "$HOME/.zshrc" <<'ZRC'

### Project Alpha PATH
if [ -d "$HOME/Desktop/Project Alpha/bin" ]; then
  export PATH="$HOME/Desktop/Project Alpha/bin:$PATH"
fi
ZRC
  echo "Added Project Alpha bin to PATH in ~/.zshrc"
fi

# 7) Makefile with common tasks
cat > "$ROOT/Makefile" <<'MK'
SHELL := /bin/bash
ROOT  := $(HOME)/Desktop/Project\ Alpha
BIN   := $(ROOT)/bin
DATA  := $(ROOT)/data
PEAKS := $(DATA)/peaks
CHAINS:= $(DATA)/chains
LIFT  := $(DATA)/liftover
SORT := LC_ALL=C sort -k1,1 -k2,2n

.PHONY: all qc liftover app db vacuum tree

all: tree qc liftover

tree:
	@echo "Project tree:" && find . -maxdepth 3 -type d | sed 's|^\./||'

qc:
	$(SORT) $(PEAKS)/human/../../human.qc.bed | uniq > $(PEAKS)/human/human.qc.sorted.bed
	$(SORT) $(PEAKS)/mouse/../../mouse.qc.bed | uniq > $(PEAKS)/mouse/mouse.qc.sorted.bed
	$(SORT) $(PEAKS)/pig/../../pig.qc.bed     | uniq > $(PEAKS)/pig/pig.qc.sorted.bed
	cp -f $(PEAKS)/human/human.qc.sorted.bed $(PEAKS)/human/human.hg38.bed

liftover:
	$(BIN)/liftOver $(PEAKS)/mouse/mouse.qc.sorted.bed $(CHAINS)/mm39ToHg38.over.chain.gz $(LIFT)/mouse_to_hg38/mouse.liftover.hg38.bed $(LIFT)/mouse_to_hg38/mouse.unmapped.bed
	$(BIN)/liftOver $(PEAKS)/pig/pig.qc.sorted.bed     $(CHAINS)/susScr11ToHg38.over.chain.gz $(LIFT)/pig_to_hg38/pig.liftover.hg38.bed   $(LIFT)/pig_to_hg38/pig.unmapped.bed

db:
	@Rscript -e 'DBI::dbConnect(RSQLite::SQLite(), "data/regland.sqlite") |> DBI::dbDisconnect()'

vacuum:
	@Rscript -e 'con<-DBI::dbConnect(RSQLite::SQLite(),"data/regland.sqlite"); DBI::dbExecute(con,"VACUUM"); DBI::dbExecute(con,"ANALYZE"); DBI::dbDisconnect(con)'
MK

echo "✅ Tidy complete."
