# CTDNA Liquid Biopsy Agent

> **Domain:** Clinical Decision Support & Biomedical Computing
> **Reference Guidelines & Standards:** CAP / CLSI / ISO / AMP / ASCO Standards

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

**CTDNA Liquid Biopsy Agent** is an advanced analytical and computational platform implementing Plasma ctDNA Variant Allele Kinetics & Molecular Relapse Tracker. It provides:

- **Multi-agent clinical decision support** with specialized sub-agents for VAF kinetics tracking, MRD status classification, and clonal evolution detection
- **CHIP filtering** to distinguish true somatic tumor-derived ctDNA from clonal hematopoiesis background noise
- **MRD tracking** with exponential decay modeling, molecular response classification, and relapse prediction
- **Multi-platform concordance analysis** across Guardant360, FoundationOne Liquid CDx, and Signatera
- **Zero-PHI outbound protection** with AST and regex inspection blocking sensitive identifiers
- **Tamper-evident HMAC-SHA256 audit trail** with cryptographic chain verification

---

## 🏗️ Architecture

The project contains two parallel implementations:

### Enterprise Agent System (`agents/`)
- FastAPI REST API with OpenAPI 3.1 specification
- Pydantic v2 data models with validation
- HMAC-SHA256 tamper-evident audit trail
- PHI outbound guard with pattern detection
- Prometheus metrics exporter
- WebSocket telemetry streamer
- Active learning Bayesian calibration engine

### Clinical Domain Package (`ctdna_liquid_biopsy_agent/`)
- CHIP (Clonal Hematopoiesis of Indeterminate Potential) filtering module
- MRD (Molecular Residual Disease) tracking with VAF trend analysis
- Multi-platform concordance analyzer
- Clinical domain engine with guideline rules
- Specialized sub-agents for clinical audit

### Legacy Module (`ctdna_sentinel.py`)
- Standalone CLI with coordinator pattern
- Domain knowledge registry
- Basic sub-agent evaluation

---

## 🚀 Installation

```bash
# Clone the repository
git clone https://github.com/abusuraihsakhri/ctdna-liquid-biopsy-agent.git
cd ctdna-liquid-biopsy-agent

# Install dependencies
pip install fastapi uvicorn pydantic pytest

# Or install as a package
pip install -e .
```

---

## 💻 CLI Quickstart & Usage

### 1. Run Single Task Evaluation
```bash
python cli.py audit --task-id TASK-001 --target KEY-01 --primary 28.5 --secondary 14.2 --critical --status DISCORDANT
```

### 2. Query Supervisory Chat
```bash
python cli.py chat "What is the system status?"
```

### 3. Batch Process CSV Records
```bash
python cli.py batch -i sample.csv -o results.csv
```

### 4. Verify Audit Trail Integrity
```bash
python cli.py verify-audit
```

### 5. Launch FastAPI REST Server
```bash
python cli.py serve --host 127.0.0.1 --port 8000
```

### Parameter Reference
| Argument | Description | Default |
|:---------|:------------|:--------|
| `--task-id` | Unique task/case identifier | TASK-2026-001 |
| `--target` | Target identifier (de-identified) | KEY-TARGET-01 |
| `--primary` | Primary measurement value | 28.5 |
| `--secondary` | Secondary kinetic/confidence score | 14.2 |
| `--critical` | Trigger emergency escalation | False |
| `--status` | Status/phenotype descriptor | DISCORDANT |

### Input Data Schema

| Field | Description | Requirement |
|:------|:------------|:------------|
| `case_id` | Unique case identifier | Required |
| `patient_synthetic_id` | De-identified patient token | Required |
| `metric_primary` | Primary VAF/kinetic measurement | Required |
| `metric_secondary` | Secondary confidence score | Required |
| `is_stat` | STAT emergency escalation flag | Required |
| `status_flag` | Phenotype/discordance descriptor | Required |

---

## 🌐 REST API Endpoints

When running the FastAPI server:

| Method | Endpoint | Description |
|:-------|:---------|:------------|
| GET | `/health` | Service health check |
| GET | `/metrics` | Operational metrics |
| POST | `/api/audit` | Submit task for evaluation |
| POST | `/api/chat` | Query supervisory assistant |
| GET | `/api/audit/logs` | Retrieve audit trail with integrity verification |

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active AST and regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition with signature verification.
* **Air-Gapped LLM Reasoning Adapter:** Agnostic integration for local Ollama instances (`llama3`, `mistral`), Claude 3.5 Sonnet, GPT-4o, and deterministic test mocks.
* **Active Learning Bayesian Calibration:** Dynamic tracker updating worker reliability weights and monitoring Brier calibration drift.
* **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`).

### Security Configuration

Set the audit secret key via environment variable:
```bash
export AUDIT_SECRET_KEY="your-secure-random-key"
```

If not set, a random key is generated at runtime (recommended for production to set explicitly).

---

## 🧪 Testing & Verification

Run the full automated test suite:
```bash
pytest -v
```

Run with coverage:
```bash
pytest -v --cov=agents --cov=ctdna_liquid_biopsy_agent
```

Execute high-throughput batch simulation benchmarks:
```bash
python simulator.py 1000
```

---

## 🐳 Container Deployment

### Docker
```bash
docker build -t ctdna-liquid-biopsy-agent .
docker run -p 8000:8000 -e AUDIT_SECRET_KEY="your-secure-key" ctdna-liquid-biopsy-agent
```

### Docker Compose
```bash
export AUDIT_SECRET_KEY="your-secure-key"
docker-compose up -d
```

---

## 📁 Project Structure

```
ctdna-liquid-biopsy-agent/
├── agents/                          # Enterprise agent system
│   ├── api.py                       # FastAPI REST server
│   ├── base.py                      # Security, PHI guard, audit trail
│   ├── models.py                    # Pydantic v2 data models
│   ├── supervisor.py                # Master orchestrator
│   ├── workers.py                   # Specialized evaluation workers
│   ├── llm_factory.py               # LLM provider factory
│   ├── metrics.py                   # Prometheus metrics
│   ├── learning.py                  # Bayesian calibration engine
│   └── streamer.py                  # WebSocket telemetry
├── ctdna_liquid_biopsy_agent/       # Clinical domain package
│   ├── agents.py                    # Clinical sub-agents
│   ├── chip_filter.py               # CHIP variant filtering
│   ├── mrd_tracker.py               # MRD tracking with VAF trends
│   ├── concordance.py               # Multi-platform concordance
│   ├── engine.py                    # Clinical domain rules
│   ├── models.py                    # Clinical data models
│   ├── server.py                    # Clinical FastAPI server
│   └── cli.py                       # Clinical CLI
├── tests/                           # Test suite
│   ├── test_clinical_modules.py     # Clinical module tests
│   ├── test_ctdna_liquid_biopsy_agent.py  # Enterprise agent tests
│   └── test_enrichment.py           # Enrichment engine tests
├── cli.py                           # Main CLI entry point
├── ctdna_sentinel.py                # Legacy standalone module
├── enrichment.py                    # Enrichment feature engines
├── simulator.py                     # High-throughput simulator
├── web/index.html                   # Operations console UI
├── Dockerfile                       # Container definition
├── docker-compose.yml               # Multi-service orchestration
├── pyproject.toml                   # Project configuration
└── README.md                        # This file
```

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
