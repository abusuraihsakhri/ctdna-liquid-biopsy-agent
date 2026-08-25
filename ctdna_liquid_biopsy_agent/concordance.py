"""
Multi-Platform ctDNA Concordance Analysis Module.
Cross-references ctDNA results from different platforms (Guardant360, FoundationOne Liquid CDx,
Signatera) with confidence-weighted consensus variant calls.
Domain: Precision Oncology — Liquid Biopsy
"""
import datetime
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any


# Platform-specific characteristics
PLATFORM_PROFILES = {
    "GUARDANT360": {
        "panel_genes": 73,
        "sensitivity_vaf": 0.1,  # ~0.1% VAF LOD
        "type": "plasma_only",
        "reporting_threshold": 0.1,
    },
    "FOUNDATIONONE_LIQ_CDX": {
        "panel_genes": 324,
        "sensitivity_vaf": 0.2,  # ~0.2% VAF LOD
        "type": "plasma_only",
        "reporting_threshold": 0.25,
    },
    "SIGNATERA": {
        "panel_genes": 16,  # personalized
        "sensitivity_vaf": 0.01,  # ~0.01% VAF LOD (tumor-informed)
        "type": "tumor_informed",
        "reporting_threshold": 0.01,
    },
    "CUSTOM": {
        "panel_genes": 0,
        "sensitivity_vaf": 0.5,
        "type": "variable",
        "reporting_threshold": 0.5,
    },
}


@dataclass
class PlatformCtDNAResult:
    platform: str  # "GUARDANT360", "FOUNDATIONONE_LIQ_CDX", "SIGNATERA", "CUSTOM"
    variant_id: str
    gene: str
    vaf_percent: float
    depth: int
    detected: bool
    sensitivity: float  # platform-specific detection sensitivity


@dataclass
class VariantConcordance:
    variant_id: str
    gene: str
    platforms_detected: List[str]
    platforms_missed: List[str]
    concordant: bool
    confidence_weighted_vaf: float
    consensus_call: str  # "DETECTED", "NOT_DETECTED", "UNCERTAIN"


@dataclass
class PlatformDiscrepancy:
    variant_id: str
    platform_a: str
    platform_b: str
    vaf_a: float
    vaf_b: float
    possible_reason: str  # "SENSITIVITY_DIFFERENCE", "SAMPLING_VARIATION", "TECHNICAL_ARTIFACT"


@dataclass
class ConcordanceAnalysis:
    case_id: str
    patient_synthetic_id: str
    platforms_analyzed: List[str]
    variant_concordance: List[VariantConcordance]
    overall_concordance_score: float  # 0-100
    consensus_variants: List[str]
    platform_discrepancies: List[PlatformDiscrepancy]
    recommended_platform_for_monitoring: str
    timestamp: str = field(default_factory=lambda: datetime.datetime.now(datetime.timezone.utc).isoformat())


class ConcordanceAnalyzer:
    """
    Analyzes concordance of ctDNA results across multiple testing platforms.
    Produces confidence-weighted consensus variant calls.
    """

    def __init__(self):
        self._platform_profiles = PLATFORM_PROFILES

    def analyze_concordance(
        self,
        case_id: str,
        patient_synthetic_id: str,
        platform_results: List[PlatformCtDNAResult],
    ) -> ConcordanceAnalysis:
        """
        Analyze concordance across multiple platform results.

        Args:
            case_id: Unique case identifier
            patient_synthetic_id: De-identified patient ID
            platform_results: Results from different ctDNA platforms

        Returns:
            ConcordanceAnalysis with consensus calls and discrepancy analysis
        """
        if not platform_results:
            return self._empty_analysis(case_id, patient_synthetic_id)

        # Get unique platforms
        platforms = list(set(r.platform for r in platform_results))

        # Group results by variant_id
        variant_groups: Dict[str, List[PlatformCtDNAResult]] = {}
        for r in platform_results:
            variant_groups.setdefault(r.variant_id, []).append(r)

        # Analyze each variant
        concordances: List[VariantConcordance] = []
        discrepancies: List[PlatformDiscrepancy] = []

        for variant_id, results in variant_groups.items():
            vc, disc = self._analyze_variant(variant_id, results, platforms)
            concordances.append(vc)
            discrepancies.extend(disc)

        # Compute overall concordance score
        overall_score = self._compute_overall_concordance(concordances)

        # Identify consensus variants
        consensus = [vc.variant_id for vc in concordances if vc.consensus_call == "DETECTED"]

        # Recommend platform for monitoring
        recommended = self._recommend_platform(platform_results, concordances)

        return ConcordanceAnalysis(
            case_id=case_id,
            patient_synthetic_id=patient_synthetic_id,
            platforms_analyzed=platforms,
            variant_concordance=concordances,
            overall_concordance_score=overall_score,
            consensus_variants=consensus,
            platform_discrepancies=discrepancies,
            recommended_platform_for_monitoring=recommended,
        )

    def _analyze_variant(
        self,
        variant_id: str,
        results: List[PlatformCtDNAResult],
        all_platforms: List[str],
    ) -> tuple:
        """Analyze concordance for a single variant across platforms."""
        detected_platforms = [r.platform for r in results if r.detected]
        detected_set = set(detected_platforms)
        missed_platforms = [p for p in all_platforms if p not in detected_set]

        n_detected = len(detected_platforms)
        n_total = len(all_platforms)

        # Confidence-weighted VAF
        if n_detected > 0:
            total_weight = 0.0
            weighted_vaf = 0.0
            for r in results:
                if r.detected:
                    weight = r.sensitivity  # higher sensitivity = higher weight
                    weighted_vaf += r.vaf_percent * weight
                    total_weight += weight
            cw_vaf = weighted_vaf / total_weight if total_weight > 0 else 0.0
        else:
            cw_vaf = 0.0

        # Concordance determination
        concordant = n_detected >= 2 or (n_detected == n_total and n_total >= 2)

        # Consensus call
        if n_detected >= 2:
            consensus = "DETECTED"
        elif n_detected == 1 and n_total > 1:
            consensus = "UNCERTAIN"
        elif n_detected == 0:
            consensus = "NOT_DETECTED"
        else:
            consensus = "DETECTED"

        # Identify discrepancies
        disc = self._find_discrepancies(variant_id, results)

        return (
            VariantConcordance(
                variant_id=variant_id,
                gene=results[0].gene if results else "",
                platforms_detected=detected_platforms,
                platforms_missed=missed_platforms,
                concordant=concordant,
                confidence_weighted_vaf=round(cw_vaf, 4),
                consensus_call=consensus,
            ),
            disc,
        )

    def _find_discrepancies(
        self, variant_id: str, results: List[PlatformCtDNAResult]
    ) -> List[PlatformDiscrepancy]:
        """Find discrepancies between platform results for a variant."""
        discrepancies = []
        detected = [r for r in results if r.detected]
        not_detected = [r for r in results if not r.detected]

        for d in detected:
            for nd in not_detected:
                # Determine likely reason
                d_profile = self._platform_profiles.get(d.platform, {})
                nd_profile = self._platform_profiles.get(nd.platform, {})
                d_sens = d_profile.get("sensitivity_vaf", 0.5)
                nd_sens = nd_profile.get("sensitivity_vaf", 0.5)

                if d.vaf_percent < nd_sens:
                    reason = "SENSITIVITY_DIFFERENCE"
                elif abs(d_sens - nd_sens) > 0.1:
                    reason = "SENSITIVITY_DIFFERENCE"
                else:
                    reason = "SAMPLING_VARIATION"

                discrepancies.append(PlatformDiscrepancy(
                    variant_id=variant_id,
                    platform_a=d.platform,
                    platform_b=nd.platform,
                    vaf_a=d.vaf_percent,
                    vaf_b=0.0,
                    possible_reason=reason,
                ))

        return discrepancies

    def _compute_overall_concordance(self, concordances: List[VariantConcordance]) -> float:
        """Compute overall concordance score (0-100)."""
        if not concordances:
            return 0.0

        concordant_count = sum(1 for vc in concordances if vc.concordant)
        return round((concordant_count / len(concordances)) * 100, 1)

    def _recommend_platform(
        self,
        results: List[PlatformCtDNAResult],
        concordances: List[VariantConcordance],
    ) -> str:
        """Recommend the best platform for ongoing monitoring."""
        # Score platforms by detection rate and sensitivity
        platform_scores: Dict[str, float] = {}
        platform_counts: Dict[str, Dict[str, int]] = {}

        for r in results:
            if r.platform not in platform_counts:
                platform_counts[r.platform] = {"detected": 0, "total": 0}
            platform_counts[r.platform]["total"] += 1
            if r.detected:
                platform_counts[r.platform]["detected"] += 1

        for platform, counts in platform_counts.items():
            detection_rate = counts["detected"] / counts["total"] if counts["total"] > 0 else 0
            profile = self._platform_profiles.get(platform, {})
            sensitivity = 1.0 / profile.get("sensitivity_vaf", 0.5)  # lower LOD = better
            platform_scores[platform] = detection_rate * 0.6 + (sensitivity / 100) * 0.4

        if not platform_scores:
            return "UNKNOWN"

        return max(platform_scores, key=platform_scores.get)

    def _empty_analysis(self, case_id: str, patient_synthetic_id: str) -> ConcordanceAnalysis:
        """Return empty analysis when no results provided."""
        return ConcordanceAnalysis(
            case_id=case_id,
            patient_synthetic_id=patient_synthetic_id,
            platforms_analyzed=[],
            variant_concordance=[],
            overall_concordance_score=0.0,
            consensus_variants=[],
            platform_discrepancies=[],
            recommended_platform_for_monitoring="UNKNOWN",
        )

    def get_platform_profile(self, platform: str) -> Optional[Dict[str, Any]]:
        """Look up platform-specific characteristics."""
        return self._platform_profiles.get(platform)
