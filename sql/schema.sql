PRAGMA foreign_keys = ON;

-- --- core dimension tables ---
CREATE TABLE species (
  species_id INTEGER PRIMARY KEY,
  common_name TEXT UNIQUE NOT NULL,
  scientific_name TEXT NOT NULL,
  assembly TEXT NOT NULL
);

CREATE TABLE genes (
  gene_id INTEGER PRIMARY KEY,
  species_id INTEGER NOT NULL REFERENCES species(species_id),
  gene_symbol TEXT,
  ensembl_gene_id TEXT,
  chrom TEXT,
  start INT,
  end INT,
  strand TEXT CHECK(strand IN ('+','-')),
  biotype TEXT,              -- e.g., protein_coding, lincRNA
  UNIQUE(species_id, ensembl_gene_id)
);

-- --- regulatory features ---
CREATE TABLE enhancers (
  enhancer_id INTEGER PRIMARY KEY,
  species_id INTEGER NOT NULL REFERENCES species(species_id),
  chrom TEXT, start INT, end INT,
  source TEXT,              -- e.g., EnsemblRegulatory, VISTA, ATAC-H3K27ac
  score REAL,               -- optional peak score
  method TEXT,              -- peak caller or annotation source
  accession TEXT            -- optional external ID
);

CREATE INDEX idx_enhancers_region ON enhancers(species_id, chrom, start, end);

CREATE TABLE ctcf_sites (
  ctcf_id INTEGER PRIMARY KEY,
  species_id INTEGER NOT NULL REFERENCES species(species_id),
  chrom TEXT, start INT, end INT,
  score REAL,
  method TEXT
);
CREATE INDEX idx_ctcf_region ON ctcf_sites(species_id, chrom, start, end);

-- --- enhancer ⇄ gene proximity/links ---
CREATE TABLE enhancer_gene_links (
  egl_id INTEGER PRIMARY KEY,
  enhancer_id INTEGER NOT NULL REFERENCES enhancers(enhancer_id) ON DELETE CASCADE,
  gene_id INTEGER NOT NULL REFERENCES genes(gene_id) ON DELETE CASCADE,
  link_type TEXT,           -- nearest_tss, ABC, correlation, etc.
  distance_bp INT
);
CREATE INDEX idx_egl_gene ON enhancer_gene_links(gene_id);

-- --- liftover & cross-species mapping ---
CREATE TABLE liftover_map (
  map_id INTEGER PRIMARY KEY,
  src_enhancer_id INTEGER NOT NULL REFERENCES enhancers(enhancer_id) ON DELETE CASCADE,
  src_species_id INTEGER NOT NULL REFERENCES species(species_id),
  dst_species_id INTEGER NOT NULL REFERENCES species(species_id),
  dst_chrom TEXT, dst_start INT, dst_end INT,
  status TEXT CHECK(status IN ('mapped','unmapped','ambiguous')),
  chain TEXT                 -- which chain file used
);

-- --- conservation labels after liftover/classification ---
CREATE TABLE enhancer_classification (
  class_id INTEGER PRIMARY KEY,
  enhancer_id INTEGER NOT NULL REFERENCES enhancers(enhancer_id) ON DELETE CASCADE,
  label TEXT NOT NULL CHECK(label IN ('conserved','gained','lost','unknown')),
  basis TEXT,               -- rule or evidence (e.g., 3/5 species present)
  notes TEXT
);

-- --- expression (minimal, flexible) ---
CREATE TABLE expression (
  expr_id INTEGER PRIMARY KEY,
  gene_id INTEGER NOT NULL REFERENCES genes(gene_id) ON DELETE CASCADE,
  species_id INTEGER NOT NULL REFERENCES species(species_id),
  tissue TEXT,
  condition TEXT,
  tpm REAL
);
CREATE INDEX idx_expr_gene ON expression(gene_id);

-- --- GWAS (coarse, per your request) ---
CREATE TABLE gwas_traits (
  trait_id INTEGER PRIMARY KEY,
  trait_group TEXT NOT NULL,         -- inflammation, obesity, alcohol
  description TEXT
);

CREATE TABLE gwas_snps (
  snp_id INTEGER PRIMARY KEY,
  rsid TEXT UNIQUE,
  chrom TEXT, pos INT,
  trait_id INTEGER NOT NULL REFERENCES gwas_traits(trait_id),
  ancestry TEXT,                     -- optional
  reference TEXT                     -- optional study/source
);

-- Overlap table to attach SNPs to enhancers (post-coordinate harmonization)
CREATE TABLE snp_enhancer_overlap (
  overlap_id INTEGER PRIMARY KEY,
  snp_id INTEGER NOT NULL REFERENCES gwas_snps(snp_id) ON DELETE CASCADE,
  species_id INTEGER NOT NULL REFERENCES species(species_id),
  enhancer_id INTEGER REFERENCES enhancers(enhancer_id) ON DELETE CASCADE
);

-- --- handy views ---
-- Enhancers near CTCF (<= 5kb by default). Adjust as needed.
CREATE VIEW v_enhancers_near_ctcf AS
SELECT e.enhancer_id, e.species_id, e.chrom, e.start, e.end,
       MIN(ABS(((e.start+e.end)/2) - ((c.start+c.end)/2))) AS mid_dist_to_ctcf_bp
FROM enhancers e
JOIN ctcf_sites c
  ON e.species_id = c.species_id
 AND e.chrom = c.chrom
 AND NOT (e.end < c.start - 5000 OR e.start > c.end + 5000)
GROUP BY e.enhancer_id;

-- Trait enrichment counts by label (for quick summaries)
CREATE VIEW v_trait_by_label AS
SELECT t.trait_group, ec.label, COUNT(*) AS n_snps_in_enhancers
FROM snp_enhancer_overlap seo
JOIN enhancers e ON e.enhancer_id = seo.enhancer_id
LEFT JOIN enhancer_classification ec ON ec.enhancer_id = e.enhancer_id
JOIN gwas_snps s ON s.snp_id = seo.snp_id
JOIN gwas_traits t ON t.trait_id = s.trait_id
GROUP BY t.trait_group, ec.label;

