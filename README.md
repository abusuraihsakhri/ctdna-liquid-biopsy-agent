# CTDNA Liquid Biopsy Agent

A research prototype for synthetic ctDNA audit examples, CHIP filtering heuristics, serial variant allele frequency (VAF) analysis, and cross-platform concordance exploration.

**Research and education only—not for clinical use.** Thresholds, CHIP probabilities, molecular response categories, assumed detection limits, relapse risk scores, and time-to-threshold extrapolations are illustrative and **not clinically validated**. They are not RECIST response criteria. Do not use these outputs to diagnose MRD, predict relapse, or direct treatment.

## Browser demonstrator

The standalone [browser interface](web/index.html) evaluates three example rules inside the current browser tab. It is static HTML/JavaScript, with no Python runtime, backend requests, analytics, persistent input storage, or signed audit records. **Use synthetic identifiers only.**

GitHub Pages deployment automation is included. The repository must be configured for GitHub Pages **GitHub Actions** publishing before a live URL can be verified.

## Features and architecture

- **agents/** — Pydantic/FastAPI example supervisor, three heuristic threshold/descriptor workers, an in-memory HMAC-SHA256 audit chain, and a deterministic mock chat provider.
- **ctdna_liquid_biopsy_agent/** — CHIP/WBC matching, serial single-variant VAF trend analysis, exploratory molecular response labels, multi-platform concordance, and a separate CLI/API.
- **ctdna_sentinel.py** — legacy standalone CLI.
- **web/index.html** — independent, browser-only implementation of the three example enterprise workers.
- **enrichment.py** and **simulator.py** — experimental and benchmark workflows.

These Python APIs have different schemas. The browser implements only illustrative enterprise rules: primary metric >25, secondary metric >12, a critical flag, and selected discordance descriptors. These are not assay-specific clinical limits.

## Install and test

Python 3.10 or later is required.

```bash
git clone https://github.com/abusuraihsakhri/ctdna-liquid-biopsy-agent.git
cd ctdna-liquid-biopsy-agent
python -m pip install -e ".[server,test]"
python -m pytest -q
```

Core dependencies include Pydantic; the **server** extra supplies FastAPI/Uvicorn and the **test** extra supplies Pytest/HTTPX. Only the deterministic **mock** model provider is implemented; external model providers are not connected.

## Command-line workflows

```bash
# Generic synthetic-task supervisor
ctdna-liquid-biopsy-engine audit --task-id TASK-001 --target SYNTH-01 --primary 28.5 --secondary 14.2 --critical --status DISCORDANT

# Clinical-package example
ctdna-clinical audit --case-id SYNTH-CASE --primary 26 --secondary 12 --status DISCORDANT

# Batch with the provided clinical sample dataset
ctdna-clinical batch -i sample.csv -o results.csv

# Start the enterprise API locally
ctdna-liquid-biopsy-engine serve --host 127.0.0.1 --port 8000
```

The enterprise batch CLI expects task_id, target_identifier, primary_metric, secondary_metric, is_critical_flag, and status_descriptor. The clinical batch CLI accepts case_id, patient_synthetic_id, metric_primary, metric_secondary, is_stat, and status_flag. Boolean CSV values such as True/False, 1/0, and yes/no are parsed explicitly; unrecognized values are rejected. CSV exports neutralize potentially executable spreadsheet formulas.

## REST endpoints

The enterprise service exposes:

| Method | Route | Purpose |
| --- | --- | --- |
| GET | /health | Health check |
| GET | /metrics | In-process counters |
| POST | /api/audit | Evaluate a synthetic task |
| POST | /api/chat | Deterministic mock reply |
| GET | /api/audit/logs | Inspect in-memory audit chain |

The clinical package has a separate create_app() in its server.py with /health, /api/audit, and /api/chat. Its request schema uses case_id, primary_metric, secondary_metric, status_flag, and optional synthetic metadata.

## Docker Compose

```bash
export AUDIT_SECRET_KEY="$(python -c 'import secrets; print(secrets.token_hex(32))')"
docker compose up --build
```

Compose binds the service to **127.0.0.1:8000** on the host and requires a persistent signing key. Audit records are in memory and are lost on process restart. No database, production-grade access control, or persistent clinical audit storage is implemented.

## Privacy and limitations

Do not submit protected health information (PHI) to Python APIs, CLI batch files, or the browser UI. The identifier-pattern guard is **not** comprehensive de-identification, HIPAA compliance certification, or a guarantee against PHI disclosure. The enterprise audit log endpoint has no authentication; do not expose the server to the public Internet. The browser console only uses in-tab memory for submitted example values and does not send them to a remote service.

CHIP filtering, concordance scores, platform profiles, and molecular relapse heuristics require analytical validation, measured error characteristics, calibrated uncertainty and assay-specific cutoffs. The static demonstrator does not run Python or Pyodide; complete calculations require a local Python environment or a managed backend.

## Browser support and license

The static demo uses standard HTML/CSS/JavaScript for recent Chrome, Firefox, Safari and Edge and has no external CDN or JavaScript dependencies. [MIT License](LICENSE).
