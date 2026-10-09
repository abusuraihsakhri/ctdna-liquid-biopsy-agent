"""Regression coverage for import, CSV, security, and ctDNA edge cases."""
import csv
import math

import pytest
from pydantic import ValidationError

from agents.base import AuditTrail
from agents.llm_factory import LLMFactory
from agents.models import SystemTaskPayload
from cli import main as enterprise_main
from ctdna_liquid_biopsy_agent.cli import main as clinical_main
from ctdna_liquid_biopsy_agent.mrd_tracker import MRDTracker, SerialVAFMeasurement
from ctdna_liquid_biopsy_agent.chip_filter import CHIPFilter, PlasmaVariant, WBCSequenceData
from input_validation import parse_csv_bool, safe_csv_cell


@pytest.mark.parametrize("source,expected", [
    ("True", True), ("true", True), ("1", True), ("yes", True),
    ("False", False), ("false", False), ("0", False), ("no", False),
    ("", False), (None, False),
])
def test_csv_boolean(source, expected):
    assert parse_csv_bool(source) is expected


def test_invalid_csv_boolean_rejected():
    with pytest.raises(ValueError):
        parse_csv_bool("sometimes")


@pytest.mark.parametrize("value", ["=1+2", "+SUM(A1:A2)", "@cmd", " -1+2", "\t=1+2"])
def test_csv_formula_escaped(value):
    assert safe_csv_cell(value).startswith("'")


def test_both_batch_clis_parse_false_and_escape_formulas(tmp_path):
    enterprise_in = tmp_path / "enterprise.csv"
    enterprise_in.write_text("task_id,target_identifier,primary_metric,secondary_metric,is_critical_flag,status_descriptor\n"
                             "T1,=1+2,10,1,False,NOMINAL\n")
    enterprise_out = tmp_path / "enterprise-out.csv"
    assert enterprise_main(["batch", "-i", str(enterprise_in), "-o", str(enterprise_out)]) == 0
    with enterprise_out.open(newline="") as source:
        result = next(csv.DictReader(source))
    assert result["overall_urgency"] == "ROUTINE"
    assert result["target_identifier"] == "'=1+2"

    clinical_in = tmp_path / "clinical.csv"
    clinical_in.write_text("case_id,patient_synthetic_id,metric_primary,metric_secondary,is_stat,status_flag\n"
                           "T2,=2+2,10,1,False,NORMAL\n")
    clinical_out = tmp_path / "clinical-out.csv"
    assert clinical_main(["batch", "-i", str(clinical_in), "-o", str(clinical_out)]) == 0
    with clinical_out.open(newline="") as source:
        result = next(csv.DictReader(source))
    assert result["stat_critical_alerts"] == "0"
    assert result["patient_synthetic_id"] == "'=2+2"


def test_nonfinite_metrics_rejected():
    for bad in (math.inf, -math.inf, math.nan):
        with pytest.raises(ValidationError):
            SystemTaskPayload(task_id="TEST", target_identifier="KEY", primary_metric=bad)


def test_mock_does_not_claim_clinical_verification():
    response = LLMFactory.create("mock").invoke("How does it work?")
    assert "no clinical verification" in response
    with pytest.raises(ValueError):
        LLMFactory.create("openai")


def test_audit_snapshot_is_not_mutable():
    audit = AuditTrail("test-key")
    audit.log("tester", "test", "action", {"ok": True})
    snapshot = audit.get_trail()
    snapshot[0]["current_hash"] = "tampered"
    assert audit.verify_integrity()


def _point(day, vaf, sample, variant="EGFR-L858R"):
    return SerialVAFMeasurement(sample_id=sample, date="2026-01-01", time_from_treatment_start_days=day,
                                vaf_percent=vaf, variant_id=variant, depth=1000, alt_reads=10)


def test_mrd_rebound_prioritized_over_baseline_decline():
    report = MRDTracker().track_mrd("CASE", "SYNTH", [
        _point(0, 10, "A"), _point(30, 1, "B"), _point(60, 4, "C")
    ])
    assert report.molecular_response == "PMD"


def test_mrd_undetectable_to_rising_is_not_declining():
    report = MRDTracker().track_mrd("CASE", "SYNTH", [
        _point(0, 0, "A"), _point(30, 0.5, "B")
    ])
    assert report.vaf_trend == "RISING"
    assert report.molecular_response == "PMD"


def test_mrd_rejects_mixed_variants_and_nonfinite_vaf():
    with pytest.raises(ValueError):
        MRDTracker().track_mrd("C", "S", [
            _point(0, 1, "A"), _point(30, 2, "B", "TP53")
        ])
    with pytest.raises(ValueError):
        MRDTracker().track_mrd("C", "S", [_point(0, float("nan"), "A")])


def test_wbc_positive_discordant_variant_is_not_included():
    plasma = [PlasmaVariant("V1", "DNMT3A", 1.0, 1000, 10)]
    wbc = [WBCSequenceData("W1", "V1", "DNMT3A", 9.0, 1000, 90)]
    result = CHIPFilter().filter_variants("C", "S", plasma, wbc)
    assert result.uncertain_variants == 1
    assert result.tumor_variants_retained == 0


def test_enterprise_api_health_and_audit():
    from fastapi.testclient import TestClient
    from agents.api import app

    client = TestClient(app)
    assert client.get("/health").status_code == 200
    audit = client.post("/api/audit", json={
        "task_id": "SYNTH-API-01", "target_identifier": "SYNTH-SPEC-01",
        "primary_metric": 10, "secondary_metric": 2,
        "is_critical_flag": False, "status_descriptor": "NOMINAL"
    })
    assert audit.status_code == 200
    assert audit.json()["overall_urgency"] == "ROUTINE"
    assert client.get("/api/audit/logs").json()["verified"] is True


def test_clinical_api_health_and_audit():
    from fastapi.testclient import TestClient
    from ctdna_liquid_biopsy_agent.server import create_app

    app = create_app()
    assert app is not None
    client = TestClient(app)
    assert client.get("/health").status_code == 200
    response = client.post("/api/audit", json={
        "case_id": "SYNTH-01", "patient_synthetic_id": "SYNTH-PT",
        "primary_metric": 10, "secondary_metric": 2,
        "is_stat": False, "status_flag": "NORMAL"
    })
    assert response.status_code == 200
    assert response.json()["overall_status"] == "CONCORDANT_NORMAL"


def test_concordance_all_negative_is_agreement():
    from ctdna_liquid_biopsy_agent.concordance import ConcordanceAnalyzer, PlatformCtDNAResult
    records = [
        PlatformCtDNAResult("GUARDANT360", "V1", "EGFR", 0.0, 1000, False, 0.1),
        PlatformCtDNAResult("SIGNATERA", "V1", "EGFR", 0.0, 1000, False, 0.01),
    ]
    result = ConcordanceAnalyzer().analyze_concordance("C", "S", records)
    assert result.overall_concordance_score == 100.0
    assert result.variant_concordance[0].consensus_call == "NOT_DETECTED"


def test_concordance_does_not_mark_unassayed_platform_as_negative():
    from ctdna_liquid_biopsy_agent.concordance import ConcordanceAnalyzer, PlatformCtDNAResult
    records = [
        PlatformCtDNAResult("GUARDANT360", "V1", "EGFR", 2.0, 1000, True, 0.1),
        PlatformCtDNAResult("SIGNATERA", "V2", "KRAS", 1.0, 1000, True, 0.01),
    ]
    result = ConcordanceAnalyzer().analyze_concordance("C", "S", records)
    for variant in result.variant_concordance:
        assert variant.platforms_missed == []


def test_concordance_rejects_duplicate_platform_variant():
    from ctdna_liquid_biopsy_agent.concordance import ConcordanceAnalyzer, PlatformCtDNAResult
    v = PlatformCtDNAResult("GUARDANT360", "V1", "EGFR", 2.0, 1000, True, 0.1)
    with pytest.raises(ValueError):
        ConcordanceAnalyzer().analyze_concordance("C", "S", [v, v])
