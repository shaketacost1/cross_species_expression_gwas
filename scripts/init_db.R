
library(DBI); library(RSQLite)
db <- dbConnect(SQLite(), file.path("data", "regland.sqlite"))

dbExecute(db, "CREATE TABLE IF NOT EXISTS species(
  species_id TEXT PRIMARY KEY, name TEXT, genome_build TEXT
)")
dbExecute(db, "CREATE TABLE IF NOT EXISTS genes(
  gene_id INTEGER PRIMARY KEY, symbol TEXT, species_id TEXT,
  chrom TEXT, start INTEGER, end INTEGER
)")
dbExecute(db, "CREATE TABLE IF NOT EXISTS enhancers(
  enh_id INTEGER PRIMARY KEY, species_id TEXT, chrom TEXT, start INTEGER, end INTEGER, score REAL, source TEXT
)")
dbExecute(db, "CREATE TABLE IF NOT EXISTS enhancer_class(
  enh_id INTEGER, class TEXT, support_n INTEGER, notes TEXT
)")
dbExecute(db, "CREATE TABLE IF NOT EXISTS gene_to_enhancer(
  gene_id INTEGER, enh_id INTEGER, method TEXT, distance_bp INTEGER
)")
dbExecute(db, "CREATE TABLE IF NOT EXISTS gwas_snps(
  snp_id INTEGER PRIMARY KEY, chrom TEXT, pos INTEGER, rsid TEXT, trait TEXT, pval REAL, source TEXT
)")
dbExecute(db, "CREATE TABLE IF NOT EXISTS snp_to_enhancer(
  snp_id INTEGER, enh_id INTEGER, overlap_bp INTEGER
)")

# seed minimal demo
dbExecute(db, "INSERT OR IGNORE INTO species(species_id,name,genome_build) VALUES
  ('human_hg38', 'Human', 'hg38'), ('mouse_mm39', 'Mouse', 'mm39')")
dbExecute(db, "INSERT INTO genes(symbol,species_id,chrom,start,end) VALUES
  ('BDNF', 'human_hg38', 'chr11', 27600000, 27750000)")
dbExecute(db, "INSERT INTO enhancers(species_id,chrom,start,end,score,source) VALUES
  ('human_hg38', 'chr11', 27680000, 27680500, 25.1, 'demo'),
  ('human_hg38', 'chr11', 27720000, 27720600, 19.4, 'demo')")
dbExecute(db, "INSERT INTO enhancer_class(enh_id,class,support_n,notes) VALUES
  (1, 'conserved', 2, 'demo'), (2, 'human_gained', 1, 'demo')")
dbExecute(db, "INSERT INTO gene_to_enhancer(gene_id,enh_id,method,distance_bp) VALUES
  (1,1,'nearest', 0), (1,2,'nearest', 0)")
dbExecute(db, "INSERT INTO gwas_snps(chrom,pos,rsid,trait,pval,source) VALUES
  ('chr11', 27680250, 'rs999', 'Addiction', 1e-7, 'demo')")
dbExecute(db, "INSERT INTO snp_to_enhancer(snp_id,enh_id,overlap_bp) VALUES
  (1,1,50)")

dbDisconnect(db)
message("SQLite demo created at data/regland.sqlite")

