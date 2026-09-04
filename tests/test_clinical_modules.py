"""
Automated Pytest for ctdna-liquid-biopsy-agent Clinical Domain Modules.
Tests CHIP Filter, MRD Tracker, Concordance Analyzer, and Clinical Engine.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import pytest
from ctdna_liquid_biopsy_agent.chip_filter import (
    CHIPFilter, PlasmaVariant, WBCSequenceData, CHIPFilterResult,
)
from ctdna_liquid_biopsy_agent.mrd_tracker import (
    MRDTracker, SerialVAFMeasurement, MRDTrackingReport,
)
from ctdna_liquid_biopsy_agent.concordance import (
    ConcordanceAnalyzer, PlatformCtDNAResult, ConcordanceAnalysis,
)
from ctdna_liquid_biopsy_agent.engine import ClinicalDomainEngine
from ctdna_liquid_biopsy_agent.models import ClinicalCasePayload, UrgencyLevel
from ctdna_liquid_biopsy_agent.agents import (
    VAFKineticsTrackerAgent, MRDStatusClassifierAgent,
    ClonalEvolutionDetectorAgent, LiquidBiopsyCoordinator,
)


# =============================================================================
# CHIP FILTER TESTS
# =============================================================================

class TestCHIPFilter:
    def setup_method(self):
        self.filter = CHIPFilter()

    def test_tumor_variant_not_in_wbc(self):
        """Variant absent from WBC and not a CHIP gene should be INCLUDED."""
        plasma = [PlasmaVariant(variant_id="V1", gene="BRCA1", vaf_percent=5.0, depth=1000, alt_reads=50)]
        result = self.filter.filter_variants("CASE-001", "SYNTH-01", plasma, [], patient_age=55)
        assert result.tumor_variants_retained == 1
        assert result.chip_variants_excluded == 0
        assert result.filtered_variants[0].filtering_action == "INCLUDE"

    def test_chip_variant_in_wbc_similar_vaf(self):
        """Variant in WBC with similar VAF should be EXCLUDED as CHIP."""
        plasma = [PlasmaVariant(variant_id="V1", gene="DNMT3A", vaf_percent=8.0, depth=1000, alt_reads=80)]
        wbc = [WBCSequenceData(sample_id="W1", variant_id="V1", gene="DNMT3A", vaf_percent=7.5, depth=500, alt_reads=38)]
        result = self.filter.filter_variants("CASE-001", "SYNTH-01", plasma, wbc, patient_age=70)
        assert result.chip_variants_excluded == 1
        assert result.filtered_variants[0].is_chip is True

    def test_tumor_variant_higher_vaf_than_wbc(self):
        """Variant with plasma VAF >> WBC VAF should be INCLUDED as tumor."""
        plasma = [PlasmaVariant(variant_id="V1", gene="EGFR", vaf_percent=15.0, depth=1000, alt_reads=150)]
        wbc = [WBCSequenceData(sample_id="W1", variant_id="V1", gene="EGFR", vaf_percent=2.0, depth=500, alt_reads=10)]
        result = self.filter.filter_variants("CASE-001", "SYNTH-01", plasma, wbc, patient_age=60)
        assert result.tumor_variants_retained == 1
        assert result.filtered_variants[0].filtering_action == "INCLUDE"

    def test_chip_gene_absent_from_wbc_uncertain(self):
        """Known CHIP gene absent from WBC within typical range should be UNCERTAIN."""
        plasma = [PlasmaVariant(variant_id="V1", gene="TET2", vaf_percent=5.0, depth=1000, alt_reads=50)]
        result = self.filter.filter_variants("CASE-001", "SYNTH-01", plasma, [], patient_age=75)
        assert result.uncertain_variants == 1
        assert result.filtered_variants[0].filtering_action == "UNCERTAIN"

    def test_empty_variant_list(self):
        """Empty plasma variant list should return empty result."""
        result = self.filter.filter_variants("CASE-001", "SYNTH-01", [], [], patient_age=50)
        assert result.total_variants == 0
        assert result.tumor_variants_retained == 0

    def test_chip_gene_info_lookup(self):
        """Should return CHIP gene info for known genes."""
        info = self.filter.get_chip_gene_info("DNMT3A")
        assert info is not None
        assert info["prevalence_rank"] == 1

    def test_chip_gene_info_unknown(self):
        """Should return None for unknown genes."""
        info = self.filter.get_chip_gene_info("UNKNOWN_GENE")
        assert info is None


# =============================================================================
# MRD TRACKER TESTS
# =============================================================================

class TestMRDTracker:
    def setup_method(self):
        self.tracker = MRDTracker()

    def test_declining_vaf_trend(self):
        """Declining VAF measurements should produce DECLINING trend."""
        measurements = [
            SerialVAFMeasurement(sample_id="S1", date="2026-01-01", time_from_treatment_start_days=0, vaf_percent=10.0, variant_id="V1", depth=1000, alt_reads=100),
            SerialVAFMeasurement(sample_id="S2", date="2026-02-01", time_from_treatment_start_days=30, vaf_percent=5.0, variant_id="V1", depth=1000, alt_reads=50),
            SerialVAFMeasurement(sample_id="S3", date="2026-03-01", time_from_treatment_start_days=60, vaf_percent=1.0, variant_id="V1", depth=1000, alt_reads=10),
        ]
        report = self.tracker.track_mrd("CASE-001", "SYNTH-01", measurements)
        assert report.vaf_trend == "DECLINING"
        assert report.relapse_risk_score < 50

    def test_rising_vaf_trend(self):
        """Rising VAF measurements should produce RISING trend."""
        measurements = [
            SerialVAFMeasurement(sample_id="S1", date="2026-01-01", time_from_treatment_start_days=0, vaf_percent=0.5, variant_id="V1", depth=1000, alt_reads=5),
            SerialVAFMeasurement(sample_id="S2", date="2026-02-01", time_from_treatment_start_days=30, vaf_percent=2.0, variant_id="V1", depth=1000, alt_reads=20),
            SerialVAFMeasurement(sample_id="S3", date="2026-03-01", time_from_treatment_start_days=60, vaf_percent=5.0, variant_id="V1", depth=1000, alt_reads=50),
        ]
        report = self.tracker.track_mrd("CASE-001", "SYNTH-01", measurements)
        assert report.vaf_trend == "RISING"
        assert report.relapse_risk_score > 50

    def test_cmr_achieved(self):
        """Undetectable VAF should indicate CMR."""
        measurements = [
            SerialVAFMeasurement(sample_id="S1", date="2026-01-01", time_from_treatment_start_days=0, vaf_percent=5.0, variant_id="V1", depth=1000, alt_reads=50),
            SerialVAFMeasurement(sample_id="S2", date="2026-02-01", time_from_treatment_start_days=30, vaf_percent=0.005, variant_id="V1", depth=1000, alt_reads=0),
        ]
        report = self.tracker.track_mrd("CASE-001", "SYNTH-01", measurements)
        assert report.cmr_achieved is True
        assert report.molecular_response == "CMR"

    def test_empty_measurements(self):
        """Empty measurements should return empty report."""
        report = self.tracker.track_mrd("CASE-001", "SYNTH-01", [])
        assert report.total_variants == 0 if hasattr(report, 'total_variants') else len(report.serial_measurements) == 0
        assert report.vaf_trend == "STABLE"

    def test_waterfall_data_computed(self):
        """Waterfall plot data should be computed for all measurements."""
        measurements = [
            SerialVAFMeasurement(sample_id="S1", date="2026-01-01", time_from_treatment_start_days=0, vaf_percent=10.0, variant_id="V1", depth=1000, alt_reads=100),
            SerialVAFMeasurement(sample_id="S2", date="2026-02-01", time_from_treatment_start_days=30, vaf_percent=5.0, variant_id="V1", depth=1000, alt_reads=50),
        ]
        report = self.tracker.track_mrd("CASE-001", "SYNTH-01", measurements)
        assert len(report.waterfall_plot_data) == 2
        assert report.waterfall_plot_data[0]["percent_change_from_baseline"] == 0.0

    def test_response_definition_lookup(self):
        """Should return response definition for known categories."""
        definition = self.tracker.get_response_definition("CMR")
        assert definition is not None
        assert definition.category == "CMR"


# =============================================================================
# CONCORDANCE ANALYZER TESTS
# =============================================================================

class TestConcordanceAnalyzer:
    def setup_method(self):
        self.analyzer = ConcordanceAnalyzer()

    def test_concordant_variant(self):
        """Variant detected on multiple platforms should be concordant."""
        results = [
            PlatformCtDNAResult(platform="GUARDANT360", variant_id="V1", gene="EGFR", vaf_percent=2.0, depth=1000, detected=True, sensitivity=0.1),
            PlatformCtDNAResult(platform="SIGNATERA", variant_id="V1", gene="EGFR", vaf_percent=1.8, depth=2000, detected=True, sensitivity=0.01),
        ]
        analysis = self.analyzer.analyze_concordance("CASE-001", "SYNTH-01", results)
        assert len(analysis.consensus_variants) == 1
        assert analysis.variant_concordance[0].concordant is True

    def test_discordant_variant(self):
        """Variant detected on only one platform should be uncertain."""
        results = [
            PlatformCtDNAResult(platform="GUARDANT360", variant_id="V1", gene="EGFR", vaf_percent=2.0, depth=1000, detected=True, sensitivity=0.1),
            PlatformCtDNAResult(platform="SIGNATERA", variant_id="V1", gene="EGFR", vaf_percent=0.0, depth=2000, detected=False, sensitivity=0.01),
        ]
        analysis = self.analyzer.analyze_concordance("CASE-001", "SYNTH-01", results)
        assert analysis.variant_concordance[0].consensus_call == "UNCERTAIN"

    def test_empty_results(self):
        """Empty results should return empty analysis."""
        analysis = self.analyzer.analyze_concordance("CASE-001", "SYNTH-01", [])
        assert analysis.overall_concordance_score == 0.0
        assert len(analysis.variant_concordance) == 0

    def test_platform_profile_lookup(self):
        """Should return platform profile for known platforms."""
        profile = self.analyzer.get_platform_profile("GUARDANT360")
        assert profile is not None
        assert profile["panel_genes"] == 73

    def test_platform_recommendation(self):
        """Should recommend a platform for monitoring."""
        results = [
            PlatformCtDNAResult(platform="GUARDANT360", variant_id="V1", gene="EGFR", vaf_percent=2.0, depth=1000, detected=True, sensitivity=0.1),
            PlatformCtDNAResult(platform="SIGNATERA", variant_id="V1", gene="EGFR", vaf_percent=1.8, depth=2000, detected=True, sensitivity=0.01),
        ]
        analysis = self.analyzer.analyze_concordance("CASE-001", "SYNTH-01", results)
        assert analysis.recommended_platform_for_monitoring != "UNKNOWN"


# =============================================================================
# CLINICAL ENGINE TESTS
# =============================================================================

class TestClinicalDomainEngine:
    def test_primary_index_exceeded(self):
        """Value above baseline should return alert dict."""
        result = ClinicalDomainEngine.evaluate_primary_index(25.0)
        assert result is not None
        assert "exceeds" in result["finding"].lower()

    def test_primary_index_normal(self):
        """Value below baseline should return None."""
        result = ClinicalDomainEngine.evaluate_primary_index(15.0)
        assert result is None

    def test_secondary_kinetics_stat(self):
        """STAT flag should trigger escalation."""
        result = ClinicalDomainEngine.evaluate_secondary_kinetics(5.0, is_stat=True)
        assert result is not None
        assert "STAT" in result["title"]

    def test_secondary_kinetics_above_limit(self):
        """Value above limit should trigger escalation."""
        result = ClinicalDomainEngine.evaluate_secondary_kinetics(15.0, is_stat=False)
        assert result is not None

    def test_secondary_kinetics_normal(self):
        """Value below limit without STAT should return None."""
        result = ClinicalDomainEngine.evaluate_secondary_kinetics(5.0, is_stat=False)
        assert result is None

    def test_biomarker_discordance(self):
        """Discordant status should return alert dict."""
        result = ClinicalDomainEngine.evaluate_biomarker_concordance("DISCORDANT", {})
        assert result is not None
        assert "divergence" in result["finding"].lower() or "discordance" in result["finding"].lower()

    def test_biomarker_normal(self):
        """Normal status should return None."""
        result = ClinicalDomainEngine.evaluate_biomarker_concordance("NORMAL", {})
        assert result is None


# =============================================================================
# CLINICAL AGENT TESTS
# =============================================================================

class TestClinicalAgents:
    def test_vaf_kinetics_tracker_alert(self):
        """Primary metric above threshold should produce alert."""
        agent = VAFKineticsTrackerAgent()
        case = ClinicalCasePayload(
            case_id="CASE-001", patient_synthetic_id="SYNTH-01",
            primary_metric=25.0, secondary_metric=5.0,
            status_flag="NORMAL", is_stat=False,
        )
        alerts = agent.audit(case)
        assert len(alerts) == 1
        assert alerts[0].urgency == UrgencyLevel.WARNING

    def test_vaf_kinetics_tracker_normal(self):
        """Primary metric below threshold should produce no alerts."""
        agent = VAFKineticsTrackerAgent()
        case = ClinicalCasePayload(
            case_id="CASE-001", patient_synthetic_id="SYNTH-01",
            primary_metric=15.0, secondary_metric=5.0,
            status_flag="NORMAL", is_stat=False,
        )
        alerts = agent.audit(case)
        assert len(alerts) == 0

    def test_mrd_classifier_stat(self):
        """STAT flag should produce critical alert."""
        agent = MRDStatusClassifierAgent()
        case = ClinicalCasePayload(
            case_id="CASE-001", patient_synthetic_id="SYNTH-01",
            primary_metric=10.0, secondary_metric=5.0,
            status_flag="NORMAL", is_stat=True,
        )
        alerts = agent.audit(case)
        assert len(alerts) == 1
        assert alerts[0].urgency == UrgencyLevel.STAT_CRITICAL

    def test_clonal_evolution_discordant(self):
        """Discordant status should produce advisory alert."""
        agent = ClonalEvolutionDetectorAgent()
        case = ClinicalCasePayload(
            case_id="CASE-001", patient_synthetic_id="SYNTH-01",
            primary_metric=10.0, secondary_metric=5.0,
            status_flag="DISCORDANT", is_stat=False,
        )
        alerts = agent.audit(case)
        assert len(alerts) == 1
        assert alerts[0].urgency == UrgencyLevel.ADVISORY

    def test_coordinator_process_case(self):
        """Coordinator should process case and return dossier."""
        coord = LiquidBiopsyCoordinator()
        case = ClinicalCasePayload(
            case_id="CASE-001", patient_synthetic_id="SYNTH-01",
            primary_metric=10.0, secondary_metric=5.0,
            status_flag="NORMAL", is_stat=False,
        )
        dossier = coord.process_case(case)
        assert dossier["case_id"] == "CASE-001"
        assert dossier["overall_status"] == "CONCORDANT_NORMAL"

    def test_coordinator_critical_case(self):
        """Coordinator should detect critical status."""
        coord = LiquidBiopsyCoordinator()
        case = ClinicalCasePayload(
            case_id="CASE-001", patient_synthetic_id="SYNTH-01",
            primary_metric=10.0, secondary_metric=15.0,
            status_flag="NORMAL", is_stat=True,
        )
        dossier = coord.process_case(case)
        assert dossier["overall_status"] == "CRITICAL_ACTION_REQUIRED"
        assert dossier["stat_critical_alerts"] > 0

    def test_coordinator_chat_status(self):
        """Chat query about status should return case count."""
        coord = LiquidBiopsyCoordinator()
        response = coord.query_supervisory_chat("What is the status?")
        assert "tracking" in response.lower() or "cases" in response.lower()

    def test_coordinator_chat_guidelines(self):
        """Chat query about guidelines should return standards info."""
        coord = LiquidBiopsyCoordinator()
        response = coord.query_supervisory_chat("What are the guidelines?")
        assert "AMP" in response or "CAP" in response or "standard" in response.lower()
