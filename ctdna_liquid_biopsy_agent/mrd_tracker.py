"""
MRD (Molecular Residual Disease) Tracking Module with VAF Trend Analysis.
Implements exponential decay modeling, molecular response classification,
and relapse prediction for ctDNA longitudinal monitoring.
Domain: Precision Oncology — Liquid Biopsy
"""
import math
import datetime
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any


@dataclass
class SerialVAFMeasurement:
    sample_id: str
    date: str
    time_from_treatment_start_days: int
    vaf_percent: float
    variant_id: str
    depth: int
    alt_reads: int


@dataclass
class MolecularResponseClassification:
    category: str  # "CMR", "PMR", "SMD", "PMD"
    definition: str
    vaf_threshold: str
    imaging_correlation: str


@dataclass
class MRDTrackingReport:
    case_id: str
    patient_synthetic_id: str
    serial_measurements: List[SerialVAFMeasurement]
    vaf_trend: str  # "DECLINING", "STABLE", "RISING"
    vaf_half_life_days: Optional[float]
    molecular_response: str  # "CMR", "PMR", "SMD", "PMD"
    cmr_achieved: bool
    cmr_date: Optional[str]
    relapse_risk_score: float  # 0-100
    predicted_time_to_relapse_days: Optional[int]
    waterfall_plot_data: List[Dict[str, Any]]
    timestamp: str = field(default_factory=lambda: datetime.datetime.now(datetime.timezone.utc).isoformat())


# Molecular response thresholds
LOD_THRESHOLD = 0.01  # Limit of detection: 0.01% VAF
CMR_THRESHOLD = LOD_THRESHOLD
PMR_DECLINE_THRESHOLD = 0.50  # 50% decline from baseline
PMD_INCREASE_THRESHOLD = 0.50  # 50% increase from nadir

# Molecular response definitions
RESPONSE_DEFINITIONS = {
    "CMR": MolecularResponseClassification(
        category="CMR",
        definition="Complete Molecular Response — ctDNA undetectable below limit of detection",
        vaf_threshold=f"< {LOD_THRESHOLD}%",
        imaging_correlation="Expected CR or PR on imaging",
    ),
    "PMR": MolecularResponseClassification(
        category="PMR",
        definition="Partial Molecular Response — ctDNA declining ≥50% from baseline",
        vaf_threshold="≥50% decline from baseline",
        imaging_correlation="Likely PR on imaging",
    ),
    "SMD": MolecularResponseClassification(
        category="SMD",
        definition="Stable Molecular Disease — ctDNA <50% decline, no significant increase",
        vaf_threshold="<50% decline, <50% increase",
        imaging_correlation="Stable disease on imaging",
    ),
    "PMD": MolecularResponseClassification(
        category="PMD",
        definition="Progressive Molecular Disease — ctDNA increase ≥50% or new variant detected",
        vaf_threshold="≥50% increase from nadir",
        imaging_correlation="Likely PD on imaging",
    ),
}


class MRDTracker:
    """
    Tracks ctDNA VAF trends over time using exponential decay modeling.
    Classifies molecular response and predicts relapse risk.
    """

    def __init__(self):
        self._response_defs = RESPONSE_DEFINITIONS

    def track_mrd(
        self,
        case_id: str,
        patient_synthetic_id: str,
        measurements: List[SerialVAFMeasurement],
    ) -> MRDTrackingReport:
        """
        Analyze serial VAF measurements to produce an MRD tracking report.

        Args:
            case_id: Unique case identifier
            patient_synthetic_id: De-identified patient ID
            measurements: List of serial VAF measurements sorted by date

        Returns:
            MRDTrackingReport with trend analysis, molecular response, and relapse risk
        """
        if not measurements:
            return self._empty_report(case_id, patient_synthetic_id)

        sorted_measurements = sorted(measurements, key=lambda m: m.time_from_treatment_start_days)

        # Compute waterfall plot data (percent change from baseline)
        baseline_vaf = sorted_measurements[0].vaf_percent
        waterfall_data = self._compute_waterfall(sorted_measurements, baseline_vaf)

        # Fit exponential decay model
        half_life = self._fit_exponential_decay(sorted_measurements)

        # Determine VAF trend
        vaf_trend = self._classify_trend(sorted_measurements)

        # Classify molecular response
        molecular_response = self._classify_molecular_response(sorted_measurements, baseline_vaf)

        # Check for CMR
        cmr_achieved, cmr_date = self._check_cmr(sorted_measurements)

        # Compute relapse risk score
        relapse_risk = self._compute_relapse_risk(
            sorted_measurements, vaf_trend, molecular_response, cmr_achieved
        )

        # Predict time to relapse if rising
        predicted_relapse_days = self._predict_relapse_time(
            sorted_measurements, vaf_trend, half_life
        )

        return MRDTrackingReport(
            case_id=case_id,
            patient_synthetic_id=patient_synthetic_id,
            serial_measurements=sorted_measurements,
            vaf_trend=vaf_trend,
            vaf_half_life_days=half_life,
            molecular_response=molecular_response,
            cmr_achieved=cmr_achieved,
            cmr_date=cmr_date,
            relapse_risk_score=relapse_risk,
            predicted_time_to_relapse_days=predicted_relapse_days,
            waterfall_plot_data=waterfall_data,
        )

    def _fit_exponential_decay(self, measurements: List[SerialVAFMeasurement]) -> Optional[float]:
        """
        Fit VAF(t) = VAF_0 * exp(-lambda * t) using least-squares on log-transformed data.
        Returns half-life in days, or None if insufficient data.
        """
        # Filter to detectable VAFs only
        detectable = [m for m in measurements if m.vaf_percent > LOD_THRESHOLD]
        if len(detectable) < 2:
            return None

        # Log-linear regression: ln(VAF) = ln(VAF_0) - lambda * t
        n = len(detectable)
        sum_t = sum(m.time_from_treatment_start_days for m in detectable)
        sum_log_v = sum(math.log(m.vaf_percent) for m in detectable)
        sum_t2 = sum(m.time_from_treatment_start_days ** 2 for m in detectable)
        sum_tv = sum(
            m.time_from_treatment_start_days * math.log(m.vaf_percent)
            for m in detectable
        )

        denom = n * sum_t2 - sum_t ** 2
        if denom == 0:
            return None

        # lambda = -(n * sum_tv - sum_t * sum_log_v) / denom
        lam = -(n * sum_tv - sum_t * sum_log_v) / denom

        if lam <= 0:
            return None  # Not decaying

        half_life = math.log(2) / lam
        return round(half_life, 1)

    def _classify_trend(self, measurements: List[SerialVAFMeasurement]) -> str:
        """Classify overall VAF trend as DECLINING, STABLE, or RISING."""
        if len(measurements) < 2:
            return "STABLE"

        detectable = [m for m in measurements if m.vaf_percent > LOD_THRESHOLD]
        if len(detectable) < 2:
            return "DECLINING"

        first_half = detectable[: len(detectable) // 2 + 1]
        second_half = detectable[len(detectable) // 2 :]

        avg_first = sum(m.vaf_percent for m in first_half) / len(first_half)
        avg_second = sum(m.vaf_percent for m in second_half) / len(second_half)

        change_pct = (avg_second - avg_first) / avg_first if avg_first > 0 else 0

        if change_pct < -0.20:
            return "DECLINING"
        elif change_pct > 0.20:
            return "RISING"
        return "STABLE"

    def _classify_molecular_response(
        self, measurements: List[SerialVAFMeasurement], baseline_vaf: float
    ) -> str:
        """Classify molecular response per RECIST 1.1 molecular criteria."""
        if not measurements:
            return "SMD"

        latest = measurements[-1]
        latest_vaf = latest.vaf_percent

        # Check for CMR: undetectable
        if latest_vaf < CMR_THRESHOLD:
            return "CMR"

        # Check for PMR: ≥50% decline from baseline
        if baseline_vaf > 0:
            decline = (baseline_vaf - latest_vaf) / baseline_vaf
            if decline >= PMR_DECLINE_THRESHOLD:
                return "PMR"

        # Check for PMD: ≥50% increase from nadir
        nadir_vaf = min(m.vaf_percent for m in measurements)
        if nadir_vaf > 0:
            increase = (latest_vaf - nadir_vaf) / nadir_vaf
            if increase >= PMD_INCREASE_THRESHOLD:
                return "PMD"

        return "SMD"

    def _check_cmr(self, measurements: List[SerialVAFMeasurement]) -> tuple:
        """Check if CMR was achieved and when."""
        for m in measurements:
            if m.vaf_percent < CMR_THRESHOLD:
                return True, m.date
        return False, None

    def _compute_relapse_risk(
        self,
        measurements: List[SerialVAFMeasurement],
        trend: str,
        response: str,
        cmr_achieved: bool,
    ) -> float:
        """
        Compute relapse risk score (0-100).
        Rising VAF after CMR → high risk; sustained CMR → low risk.
        """
        score = 50.0  # baseline

        if trend == "RISING":
            score += 30
        elif trend == "DECLINING":
            score -= 20

        if response == "CMR":
            score -= 25
        elif response == "PMD":
            score += 25
        elif response == "PMR":
            score -= 10

        # Rising after CMR is especially concerning
        if cmr_achieved and trend == "RISING":
            score += 20

        # Clamp to 0-100
        return max(0.0, min(100.0, round(score, 1)))

    def _predict_relapse_time(
        self,
        measurements: List[SerialVAFMeasurement],
        trend: str,
        half_life: Optional[float],
    ) -> Optional[int]:
        """Predict days to clinical relapse based on VAF trajectory."""
        if trend != "RISING" or not measurements:
            return None

        detectable = [m for m in measurements if m.vaf_percent > LOD_THRESHOLD]
        if len(detectable) < 2:
            return None

        # Simple linear extrapolation to a clinical threshold (e.g., 1.0% VAF)
        clinical_threshold = 1.0
        latest = detectable[-1]
        if latest.vaf_percent >= clinical_threshold:
            return 0  # Already at threshold

        prev = detectable[-2]
        days_between = latest.time_from_treatment_start_days - prev.time_from_treatment_start_days
        if days_between <= 0:
            return None

        rate_per_day = (latest.vaf_percent - prev.vaf_percent) / days_between
        if rate_per_day <= 0:
            return None

        days_to_threshold = (clinical_threshold - latest.vaf_percent) / rate_per_day
        return max(0, int(round(days_to_threshold)))

    def _compute_waterfall(
        self, measurements: List[SerialVAFMeasurement], baseline_vaf: float
    ) -> List[Dict[str, Any]]:
        """Compute waterfall plot data (percent change from baseline)."""
        data = []
        for m in measurements:
            if baseline_vaf > 0:
                pct_change = ((m.vaf_percent - baseline_vaf) / baseline_vaf) * 100
            else:
                pct_change = 0.0
            data.append({
                "sample_id": m.sample_id,
                "date": m.date,
                "day": m.time_from_treatment_start_days,
                "vaf_percent": m.vaf_percent,
                "percent_change_from_baseline": round(pct_change, 2),
            })
        return data

    def _empty_report(self, case_id: str, patient_synthetic_id: str) -> MRDTrackingReport:
        """Return an empty report when no measurements are available."""
        return MRDTrackingReport(
            case_id=case_id,
            patient_synthetic_id=patient_synthetic_id,
            serial_measurements=[],
            vaf_trend="STABLE",
            vaf_half_life_days=None,
            molecular_response="SMD",
            cmr_achieved=False,
            cmr_date=None,
            relapse_risk_score=50.0,
            predicted_time_to_relapse_days=None,
            waterfall_plot_data=[],
        )

    def get_response_definition(self, category: str) -> Optional[MolecularResponseClassification]:
        """Look up the definition for a molecular response category."""
        return self._response_defs.get(category)
