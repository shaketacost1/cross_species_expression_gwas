# Chain files for liftover
1. Download:
   - mm10→hg38: http://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToHg38.over.chain.gz
   - hg38→mm10: http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToMm10.over.chain.gz
2. Uncompress (and strip CR if needed):

## 2. Uncompress (strip CR if needed)
gunzip -c resources/chains/mm10ToHg38.over.chain.gz | tr -d '\r' > resources/chains/mm10ToHg38.over.chain
gunzip -c resources/chains/hg38ToMm10.over.chain.gz | tr -d '\r' > resources/chains/hg38ToMm10.over.chain

## 3. Quick checks
# first line should begin with 'chain'
head -1 resources/chains/mm10ToHg38.over.chain
# count of chain records (should be tens of thousands)
grep -c '^chain ' resources/chains/*over.chain

## 4. Set environment variables (recommended via .Renviron)
# In project .Renviron (or ~/.Renviron):
MM10_TO_HG38_CHAIN=resources/chains/mm10ToHg38.over.chain
HG38_TO_MM10_CHAIN=resources/chains/hg38ToMm10.over.chain

# Verify R sees them:
# Rscript -e 'cat(Sys.getenv("MM10_TO_HG38_CHAIN"), "\n", Sys.getenv("HG38_TO_MM10_CHAIN"), "\n")'

## 5. Keep large chain files out of git
# Add this to .gitignore (already set in this repo):
resources/chains/
