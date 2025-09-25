Shiny (keep your current UI)
Needed in repo (separate from Django):
shiny/
app.R
R/ helper scripts (if any)
renv.lock (freeze R package versions)
README.md with run instructions
Runtime data source (point to):
A downloadable regland.sqlite (from GitHub Release), not committed
NOT needed: large .bed files/chains in Git
Middle Layer (API / “middleware”)
This is basically your Django REST layer.
Needed in repo:
backend/
manage.py
backend/settings.py (use env vars), urls.py
apps/regland/
models.py (schema for enhancers, enhancer_class, gwas_snps, snp_to_enhancer)
serializers.py
views.py / viewsets.py
urls.py
migrations/
management/commands/ (e.g., import_release_artifacts.py, seed_small_demo.py)
requirements.txt / pyproject.toml
README.md (setup, env vars, import steps)
API contract:
openapi.yaml (or DRF Spectacular/DRF-YASG to generate it)
Tiny test/fixture data:
fixtures/small_demo.json (a few rows for each table)
NOT needed: large .bed / .chain.gz or full SQLite—ship those in Releases and import via a command.
Typical endpoints (example):
GET /api/enhancers?species_id=mouse_mm39&chrom=chr1&start=..&end=..
GET /api/labels/summary (counts per class)
GET /api/gwas/overlap?enh_id=...
POST /api/admin/import-release (admin-only; triggers management command)
Backend Data Pipeline (repro & heavy data)
Keep pipeline code, not the big outputs.
Needed in repo:
pipeline/
Makefile (your QC + liftOver rules)
scripts/ (bash/R/Python ETL)
README.md (how to reproduce; expected inputs/outputs)
.gitignore rules for big files:
data/**/*.bed, work/**/*.bed, *.bed, *.chain.gz, *.tgz
Release artifacts (uploaded to GitHub Releases, not committed):
artifacts/human_hg38_and_sorted.tgz
artifacts/mouse_to_hg38.tgz
artifacts/pig_to_hg38.tgz
Optionally: regland.sqlite built snapshot
Nice-to-have:
scripts/build_sqlite_from_artifacts.py (or R script) to reconstruct regland.sqlite from the .tgz bundles.
Minimal file checklist (copy/paste)
Repo root
.gitignore (ignores big artifacts)
README.md
openapi.yaml (or docs/api.md)
For Django (API)
backend/ (Django project)
backend/apps/regland/ with models.py, serializers.py, views.py, migrations/
backend/requirements.txt
backend/management/commands/import_release_artifacts.py
fixtures/small_demo.json (optional)
For Web frontend (if you build it)
frontend/ with package.json, src/, .env.example
For Shiny (if you keep it)
shiny/app.R, renv.lock, README.md
For pipeline/repro
pipeline/Makefile
scripts/ (QC, liftOver, clustering, DB build)
artifacts/ (local staging; upload contents to Releases)
TL;DR
Frontend: web app (React) or Shiny, not both in the same deploy; keep Shiny separate if you still use it.
Middle layer (API): Django app + models/serializers/views + tiny fixtures + import command; no big data in git.
Backend/pipeline: Makefile + scripts in repo; publish heavy .bed and .chain.gz as GitHub Releases and load them when needed.
