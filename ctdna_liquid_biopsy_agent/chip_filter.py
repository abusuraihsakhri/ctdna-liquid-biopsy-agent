"""
CHIP (Clonal Hematopoiesis of Indeterminate Potential) Filtering Module.
Distinguishes true somatic tumor-derived ctDNA from CHIP-derived background noise
using matched WBC sequencing data.
Domain: Precision Oncology — Liquid Biopsy
"""
import datetime
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any


# Known CHIP-associated genes with typical characteristics
CHIP_GENES = {
    "DNMT3A": {"prevalence_rank": 1, "typical_vaf_range": (2.0, 20.0), "progression_risk": "LOW"},
    "TET2": {"prevalence_rank": 2, "typical_vaf_range": (2.0, 15.0), "progression_risk": "LOW"},
    "ASXL1": {"prevalence_rank": 3, "typical_vaf_range": (2.0, 10.0), "progression_risk": "MODERATE"},
    "TP53": {"prevalence_rank": 4, "typical_vaf_range": (1.0, 10.0), "progression_risk": "HIGH"},
    "JAK2": {"prevalence_rank": 5, "typical_vaf_range": (1.0, 5.0), "progression_risk": "MODERATE"},
    "SF3B1": {"prevalence_rank": 6, "typical_vaf_range": (2.0, 12.0), "progression_risk": "LOW"},
    "SRSF2": {"prevalence_rank": 7, "typical_vaf_range": (2.0, 10.0), "progression_risk": "MODERATE"},
    "U2AF1": {"prevalence_rank": 8, "typical_vaf_range": (2.0, 8.0), "progression_risk": "MODERATE"},
    "PPM1D": {"prevalence_rank": 9, "typical_vaf_range": (2.0, 8.0), "progression_risk": "LOW"},
    "IDH1": {"prevalence_rank": 10, "typical_vaf_range": (2.0, 10.0), "progression_risk": "MODERATE"},
    "IDH2": {"prevalence_rank": 11, "typical_vaf_range": (2.0, 10.0), "progression_risk": "MODERATE"},
    "CBL": {"prevalence_rank": 12, "typical_vaf_range": (2.0, 8.0), "progression_risk": "LOW"},
}

# Age-dependent CHIP prevalence estimates
CHIP_PREVALENCE_BY_AGE = {
    (0, 40): 0.01,    # <1%
    (40, 50): 0.02,   # ~2%
    (50, 60): 0.05,   # ~5%
    (60, 70): 0.08,   # ~8%
    (70, 80): 0.12,   # ~12%
    (80, 120): 0.20,  # ~20%
}


@dataclass
class WBCSequenceData:
    sample_id: str
    variant_id: str
    gene: str
    vaf_percent: float
    depth: int
    alt_reads: int


@dataclass
class PlasmaVariant:
    variant_id: str
    gene: str
    vaf_percent: float
    depth: int
    alt_reads: int
    mutation_type: str = "SNV"  # SNV, INDEL, CNV


@dataclass
class CHIPFilterResult:
    variant_id: str
    gene: str
    vaf_percent: float
    is_chip: bool
    chip_gene: Optional[str]
    chip_probability: float  # 0-1
    filtering_action: str  # "EXCLUDE", "INCLUDE", "UNCERTAIN"
    rationale: str


@dataclass
class FilteredVariantList:
    case_id: str
    patient_synthetic_id: str
    total_variants: int
    chip_variants_excluded: int
    tumor_variants_retained: int
    uncertain_variants: int
    filtered_variants: List[CHIPFilterResult]
    chip_prevalence: str  # "HIGH", "MODERATE", "LOW", "NONE"
    timestamp: str = field(default_factory=lambda: datetime.datetime.now(datetime.timezone.utc).isoformat())


class CHIPFilter:
    """
    Filters CHIP-derived variants from ctDNA variant lists using matched
    WBC sequencing data and known CHIP gene characteristics.
    """

    def __init__(self):
        self._chip_genes = CHIP_GENES
        self._prevalence_by_age = CHIP_PREVALENCE_BY_AGE

    def filter_variants(
        self,
        case_id: str,
        patient_synthetic_id: str,
        plasma_variants: List[PlasmaVariant],
        wbc_variants: List[WBCSequenceData],
        patient_age: Optional[int] = None,
    ) -> FilteredVariantList:
        """
        Filter plasma variants to exclude CHIP-derived background.

        Args:
            case_id: Unique case identifier
            patient_synthetic_id: De-identified patient ID
            plasma_variants: Variants detected in plasma
            wbc_variants: Variants detected in matched WBC sample
            patient_age: Patient age for CHIP prevalence estimation

        Returns:
            FilteredVariantList with CHIP classification for each variant
        """
        # Build WBC lookup by variant_id and gene
        wbc_by_id = {v.variant_id: v for v in wbc_variants}
        wbc_by_gene: Dict[str, List[WBCSequenceData]] = {}
        for v in wbc_variants:
            wbc_by_gene.setdefault(v.gene, []).append(v)

        # Estimate expected CHIP prevalence
        chip_prevalence = self._estimate_chip_prevalence(patient_age)

        results: List[CHIPFilterResult] = []
        for pv in plasma_variants:
            result = self._classify_variant(pv, wbc_by_id, wbc_by_gene, chip_prevalence)
            results.append(result)

        chip_excluded = sum(1 for r in results if r.filtering_action == "EXCLUDE")
        tumor_retained = sum(1 for r in results if r.filtering_action == "INCLUDE")
        uncertain = sum(1 for r in results if r.filtering_action == "UNCERTAIN")

        prevalence_label = self._prevalence_label(chip_prevalence)

        return FilteredVariantList(
            case_id=case_id,
            patient_synthetic_id=patient_synthetic_id,
            total_variants=len(results),
            chip_variants_excluded=chip_excluded,
            tumor_variants_retained=tumor_retained,
            uncertain_variants=uncertain,
            filtered_variants=results,
            chip_prevalence=prevalence_label,
        )

    def _classify_variant(
        self,
        plasma: PlasmaVariant,
        wbc_by_id: Dict[str, WBCSequenceData],
        wbc_by_gene: Dict[str, List[WBCSequenceData]],
        chip_prevalence: float,
    ) -> CHIPFilterResult:
        """Classify a single plasma variant as CHIP or tumor-derived."""
        gene = plasma.gene
        is_chip_gene = gene in self._chip_genes

        # Check if variant is present in WBC by exact ID match
        wbc_match = wbc_by_id.get(plasma.variant_id)

        if wbc_match is not None:
            # Variant found in WBC — compare VAFs
            wbc_vaf = wbc_match.vaf_percent
            plasma_vaf = plasma.vaf_percent

            # If WBC VAF is similar to plasma VAF (within 50% relative), likely CHIP
            if wbc_vaf > 0:
                ratio = plasma_vaf / wbc_vaf
                if 0.5 <= ratio <= 2.0:
                    # Similar VAFs → CHIP
                    chip_prob = 0.9 if is_chip_gene else 0.7
                    return CHIPFilterResult(
                        variant_id=plasma.variant_id,
                        gene=gene,
                        vaf_percent=plasma.vaf_percent,
                        is_chip=True,
                        chip_gene=gene if is_chip_gene else None,
                        chip_probability=chip_prob,
                        filtering_action="EXCLUDE",
                        rationale=f"Variant present in WBC at similar VAF ({wbc_vaf:.2f}%) — consistent with CHIP origin",
                    )
                elif ratio > 2.0:
                    # Plasma VAF much higher → likely tumor with some WBC contamination
                    chip_prob = 0.2
                    return CHIPFilterResult(
                        variant_id=plasma.variant_id,
                        gene=gene,
                        vaf_percent=plasma.vaf_percent,
                        is_chip=False,
                        chip_gene=None,
                        chip_probability=chip_prob,
                        filtering_action="INCLUDE",
                        rationale=f"Plasma VAF ({plasma_vaf:.2f}%) significantly higher than WBC ({wbc_vaf:.2f}%) — likely tumor-derived",
                    )

        # Variant NOT in WBC — likely tumor
        if wbc_match is None:
            # Check if it's a known CHIP gene (could be low-level CHIP below WBC detection)
            if is_chip_gene:
                chip_info = self._chip_genes[gene]
                typical_low, typical_high = chip_info["typical_vaf_range"]
                if typical_low <= plasma.vaf_percent <= typical_high:
                    chip_prob = 0.3 + chip_prevalence
                    action = "UNCERTAIN"
                    rationale = (
                        f"Known CHIP gene ({gene}) but absent from WBC — could be tumor or low-level CHIP. "
                        f"Recommend WBC re-sequencing at higher depth"
                    )
                else:
                    chip_prob = 0.1
                    action = "INCLUDE"
                    rationale = f"VAF ({plasma.vaf_percent:.2f}%) outside typical CHIP range for {gene} — likely tumor"
            else:
                chip_prob = 0.05
                action = "INCLUDE"
                rationale = "Variant absent from WBC and not a known CHIP gene — likely tumor-derived"

            return CHIPFilterResult(
                variant_id=plasma.variant_id,
                gene=gene,
                vaf_percent=plasma.vaf_percent,
                is_chip=chip_prob > 0.5,
                chip_gene=gene if chip_prob > 0.5 and is_chip_gene else None,
                chip_probability=min(1.0, chip_prob),
                filtering_action=action,
                rationale=rationale,
            )

        # Fallback
        return CHIPFilterResult(
            variant_id=plasma.variant_id,
            gene=gene,
            vaf_percent=plasma.vaf_percent,
            is_chip=False,
            chip_gene=None,
            chip_probability=0.1,
            filtering_action="INCLUDE",
            rationale="Unable to classify — defaulting to include",
        )

    def _estimate_chip_prevalence(self, age: Optional[int]) -> float:
        """Estimate age-dependent CHIP prevalence."""
        if age is None:
            return 0.05  # default ~5%
        for (low, high), prevalence in self._prevalence_by_age.items():
            if low <= age < high:
                return prevalence
        return 0.20  # very elderly

    def _prevalence_label(self, prevalence: float) -> str:
        """Convert prevalence fraction to label."""
        if prevalence >= 0.10:
            return "HIGH"
        elif prevalence >= 0.05:
            return "MODERATE"
        elif prevalence >= 0.01:
            return "LOW"
        return "NONE"

    def get_chip_gene_info(self, gene: str) -> Optional[Dict[str, Any]]:
        """Look up CHIP characteristics for a gene."""
        return self._chip_genes.get(gene)
