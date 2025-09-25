

---

## Project Structure

This project is split into three layers: **frontend**, **middleware (API)**, and **backend/pipeline**. Large genomic files (`.bed`, `.chain.gz`, etc.) are **not stored in Git**. Instead, they are published as GitHub Releases and imported as needed.

### 1. Frontend

 **Shiny app**.


* `shiny/`

  * `app.R`
  * `R/` helper scripts
  * `renv.lock` (freeze R dependencies)
  * `README.md` (instructions to run the app)
* Runtime points to a SQLite database (`regland.sqlite`) published in Releases

---

### 2. Middleware (API Layer)

Implemented with **Django REST Framework**.

* `backend/`

  * `manage.py`
  * `backend/settings.py`, `urls.py`
  * `apps/regland/`

    * `models.py` (schema: `enhancers`, `enhancer_class`, `gwas_snps`, `snp_to_enhancer`)
    * `serializers.py`
    * `views.py` / `viewsets.py`
    * `urls.py`
    * `migrations/`
    * `management/commands/` (e.g. `import_release_artifacts.py`)
  * `requirements.txt` or `pyproject.toml`
* `fixtures/small_demo.json` – tiny test dataset (few rows per table)
* `openapi.yaml` – API spec (or generate via DRF Spectacular/YASG)

**Example Endpoints:**

* `GET /api/enhancers?species_id=mouse_mm39&chrom=chr1&start=..&end=..`
* `GET /api/labels/summary`
* `GET /api/gwas/overlap?enh_id=...`

---

### 3. Backend / Data Pipeline

Contains reproducible ETL code and orchestration (but not large data).

* `pipeline/`

  * `Makefile` – QC + LiftOver + clustering rules
  * `scripts/` – Bash/R/Python ETL scripts
  * `README.md` – instructions to rerun pipeline
* `.gitignore` excludes large artifacts (`*.bed`, `*.chain.gz`, `work/**`)
* `artifacts/` – local staging of compressed bundles; upload to **GitHub Releases**

**Release Artifacts Example:**

* `human_hg38_and_sorted.tgz`
* `mouse_to_hg38.tgz`
* `pig_to_hg38.tgz`
* `regland.sqlite` (optional prebuilt database)

---

### TL;DR

* **Frontend** → Web client (React) *or* Shiny app (not both in same deploy)
* **Middleware** → Django REST API + models + serializers + views + small demo fixtures
* **Backend** → Pipeline code (Makefile + scripts), with **large data shipped via Releases**

---
project-root/
├── frontend/                # Web client (React) OR Shiny app
│   ├── src/                 # Components, pages, API client
│   ├── public/
│   ├── package.json
│   ├── .env.example
│   └── app.R                # (if using Shiny instead of React)
│
├── backend/                 # Django REST API
│   ├── manage.py
│   ├── backend/
│   │   ├── settings.py
│   │   └── urls.py
│   └── apps/
│       └── regland/
│           ├── models.py
│           ├── serializers.py
│           ├── views.py
│           ├── urls.py
│           └── migrations/
│
├── pipeline/                # Data pipeline (reproducible ETL)
│   ├── Makefile
│   ├── scripts/
│   │   ├── qc.sh
│   │   ├── liftover.sh
│   │   └── cluster.sh
│   └── README.md
│
├── artifacts/               # Bundled outputs (upload to Releases)
│   ├── human_hg38_and_sorted.tgz
│   ├── mouse_to_hg38.tgz
│   └── pig_to_hg38.tgz
│
├── fixtures/                # Tiny JSON/SQL demo datasets
│   └── small_demo.json
│
├── docs/                    # Documentation & API contract
│   └── openapi.yaml
│
├── .gitignore
├── requirements.txt         # Django / Python deps
├── renv.lock                # R deps (if Shiny used)
└── README.md

