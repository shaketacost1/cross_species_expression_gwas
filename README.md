# Project Alpha - Python Shiny

This folder contains a Python port of the Project Alpha Shiny app and Playwright tests used for CI.

Quick local run (recommended inside a virtualenv):

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python -m playwright install chromium
# start the app
python python_app.py --port 8003
```

Run smoke Playwright test:

```bash
REGLAND_TEST_URL=http://127.0.0.1:8003 pytest -q tests/test_smoke_playwright.py::test_smoke_all_charts -q -s
```

CI: GitHub Actions workflow `.github/workflows/ci.yml` starts the app on port 8000 and runs the smoke test.
