# Sprint 1 — Human enhancer conservation vs mouse (RBH)

- **Total human enhancers:** 4,067,523
- **Conserved (RBH):** 34,101
- **Human-specific vs mouse (RBH):** 4,033,422

![Counts](figs/human_vs_mouse_enhancer_conservation_RBH_counts.png)

## How to view conserved enhancers (RBH)
**UCSC Genome Browser (hg38):**
1. Go to genome.ucsc.edu → *My Data* → *Custom Tracks*.
2. Choose **GRCh38/hg38**.
3. Upload `docs/share/human_enhancers_conserved_with_mouse_RBH.bed.gz`.
4. Track name: *Conserved enhancers (RBH)*, Type: BED.
   
**IGV (hg38):**
1. *Genomes* → *Load Genome from Server* → **Human hg38**.
2. *File* → *Load from File…* and select `docs/share/human_enhancers_conserved_with_mouse_RBH.bed.gz`.

