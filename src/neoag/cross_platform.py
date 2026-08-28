from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Mapping

from .utils import read_tsv, to_float


DNA_SOURCES = {"SNV", "INDEL"}

STATUS_POLICY = {
    "CROSS_PLATFORM_PASS_CONCORDANT": ("HIGH", "common_multiplier", 1.0, "", "no"),
    "ALT_PRESENT_BELOW_PASS_OR_CALLER_DIFFERENCE": ("MODERATE", "low_level_support_multiplier", 1.0, "", "yes"),
    "COVERED_NO_ALT_SAMPLE_OR_ASSAY_DIFFERENCE": ("MODERATE", "sample_specific_multiplier", 0.85, "C_CAUTION", "yes"),
    "OTHER_COVERED_BUT_LIMITED_POWER_AT_SOURCE_VAF": ("LOW_UNASSESSED", "power_limited_multiplier", 1.0, "", "yes"),
    "OTHER_LOW_OR_NO_COVERAGE": ("LOW_UNASSESSED", "coverage_unassessed_multiplier", 1.0, "", "yes"),
    "WEAK_OR_ABSENT_ALT_REVIEW": ("LOW", "weak_other_support_multiplier", 0.9, "C_CAUTION", "yes"),
    "SOURCE_INDEL_NOT_REPRODUCED_REASSEMBLY_REQUIRED": ("LOW", "source_unreproduced_multiplier", 0.5, "C_CAUTION", "yes"),
    "SOURCE_PASS_NOT_REPRODUCED_BY_PILEUP": ("LOW", "source_unreproduced_multiplier", 0.5, "C_CAUTION", "yes"),
    "SOURCE_WEAK_EXACT_PILEUP_SUPPORT": ("LOW", "source_weak_multiplier", 0.75, "C_CAUTION", "yes"),
    "COMPLEX_REQUIRES_HAPLOTYPE_REVIEW": ("LOW", "source_unreproduced_multiplier", 0.5, "C_CAUTION", "yes"),
    "NORMAL_SUPPORT_REVIEW": ("HIGH_CAUTION", "normal_support_multiplier", 0.25, "D", "yes"),
    "PHASED_COMPONENTS_MIXED_REVIEW": ("LOW", "source_weak_multiplier", 0.75, "C_CAUTION", "yes"),
}

EVIDENCE_FIELDS = (
    "comparison_status", "cross_platform_status", "source_vcf_tumor_ad",
    "source_vcf_tumor_af", "wes_tumor_depth", "wes_tumor_alt_count",
    "wes_tumor_alt_vaf", "wgs_tumor_depth", "wgs_tumor_alt_count",
    "wgs_tumor_alt_vaf", "normal_depth", "normal_alt_count", "normal_alt_vaf",
    "other_zero_alt_probability_at_source_pileup_vaf",
)


def canonical_variant_key(chrom: Any, pos: Any, ref: Any, alt: Any) -> str:
    chrom_text = str(chrom or "")
    if chrom_text.startswith("chr"):
        chrom_text = chrom_text[3:]
    ref_text = str(ref or "").upper()
    alt_text = str(alt or "").upper()
    try:
        position = int(float(str(pos)))
    except (TypeError, ValueError):
        return ""
    while len(ref_text) > 1 and len(alt_text) > 1 and ref_text[-1] == alt_text[-1]:
        ref_text, alt_text = ref_text[:-1], alt_text[:-1]
    while len(ref_text) > 1 and len(alt_text) > 1 and ref_text[0] == alt_text[0]:
        ref_text, alt_text, position = ref_text[1:], alt_text[1:], position + 1
    if not chrom_text or not ref_text or not alt_text:
        return ""
    return f"{chrom_text}:{position}:{ref_text}:{alt_text}"


def load_cross_platform_evidence(path: str | Path | None) -> dict[str, dict[str, str]]:
    if not path or not Path(path).is_file():
        return {}
    return {
        str(row.get("variant_key") or ""): row
        for row in read_tsv(path)
        if str(row.get("variant_key") or "")
    }


def annotate_event_cross_platform(
    event: dict[str, Any],
    evidence: Mapping[str, Mapping[str, Any]],
    profile: Mapping[str, Any],
) -> dict[str, Any]:
    source = str(event.get("mutation_source") or event.get("event_type") or "").upper()
    if source not in DNA_SOURCES:
        event.setdefault("cross_platform_status", "NOT_APPLICABLE")
        event.setdefault("cross_platform_confidence", "NOT_APPLICABLE")
        event.setdefault("cross_platform_multiplier", "1.0000")
        event.setdefault("cross_platform_review_required", "no")
        return event

    key = canonical_variant_key(event.get("chrom"), event.get("pos"), event.get("ref"), event.get("alt"))
    event["cross_platform_variant_key"] = key
    row = evidence.get(key)
    component_rows = []
    if not row and event.get("component_event_ids"):
        for match in re.finditer(r"chr([^:;|]+):(\d+)([A-Za-z]+)>([A-Za-z]+)", str(event.get("component_event_ids"))):
            component_key = canonical_variant_key(match.group(1), match.group(2), match.group(3), match.group(4))
            if evidence.get(component_key):
                component_rows.append(evidence[component_key])
        if component_rows and len(component_rows) == len(str(event.get("component_event_ids")).split(";")):
            statuses = {str(item.get("cross_platform_status") or "") for item in component_rows}
            if statuses == {"CROSS_PLATFORM_PASS_CONCORDANT"}:
                row = {
                    "comparison_status": "COMMON",
                    "cross_platform_status": "CROSS_PLATFORM_PASS_CONCORDANT",
                }
                event["cross_platform_variant_key"] = ";".join(
                    str(item.get("variant_key") or "") for item in component_rows
                )
            else:
                row = {
                    "comparison_status": "MIXED_COMPONENTS",
                    "cross_platform_status": "PHASED_COMPONENTS_MIXED_REVIEW",
                }
    if not row:
        event["cross_platform_status"] = "UNASSESSED_NOT_IN_COMPARISON"
        event["cross_platform_confidence"] = "UNASSESSED"
        event["cross_platform_multiplier"] = "1.0000"
        event["cross_platform_review_required"] = "yes"
        return event

    for field in EVIDENCE_FIELDS:
        if row.get(field) not in {None, ""}:
            event[field] = row.get(field)
    status = str(row.get("cross_platform_status") or "UNASSESSED")
    confidence, config_key, default_multiplier, default_cap, review = STATUS_POLICY.get(
        status, ("UNASSESSED", "unassessed_multiplier", 1.0, "", "yes")
    )
    config = profile.get("cross_platform", {})
    event["cross_platform_confidence"] = confidence
    event["cross_platform_multiplier"] = f"{to_float(config.get(config_key), default_multiplier):.4f}"
    event["cross_platform_priority_cap"] = str(config.get(f"{config_key}_priority_cap", default_cap) or "")
    event["cross_platform_review_required"] = review
    return event


def propagate_cross_platform_to_peptide(peptide: dict[str, Any], event: Mapping[str, Any]) -> dict[str, Any]:
    for field in (
        "cross_platform_variant_key", "comparison_status", "cross_platform_status",
        "cross_platform_confidence", "cross_platform_multiplier",
        "cross_platform_priority_cap", "cross_platform_review_required",
        "source_vcf_tumor_ad", "source_vcf_tumor_af", "wes_tumor_depth",
        "wes_tumor_alt_count", "wes_tumor_alt_vaf", "wgs_tumor_depth",
        "wgs_tumor_alt_count", "wgs_tumor_alt_vaf", "normal_depth",
        "normal_alt_count", "normal_alt_vaf",
        "other_zero_alt_probability_at_source_pileup_vaf",
    ):
        peptide[field] = event.get(field, peptide.get(field, ""))
    return peptide


def cross_platform_recommendation(status: str) -> str:
    return {
        "CROSS_PLATFORM_PASS_CONCORDANT": "WES/WGS concordant DNA event",
        "ALT_PRESENT_BELOW_PASS_OR_CALLER_DIFFERENCE": "cross-platform low-level ALT support; review caller filters",
        "COVERED_NO_ALT_SAMPLE_OR_ASSAY_DIFFERENCE": "sample/time-point-specific DNA evidence; validate the intended specimen",
        "OTHER_COVERED_BUT_LIMITED_POWER_AT_SOURCE_VAF": "cross-platform assessment is power-limited at the observed VAF",
        "SOURCE_INDEL_NOT_REPRODUCED_REASSEMBLY_REQUIRED": "source InDel requires local reassembly or IGV haplotype review",
        "SOURCE_PASS_NOT_REPRODUCED_BY_PILEUP": "source PASS call was not reproduced by targeted pileup",
        "SOURCE_WEAK_EXACT_PILEUP_SUPPORT": "source call has weak exact pileup support",
        "NORMAL_SUPPORT_REVIEW": "ALT support in matched normal; exclude until artifact/germline review",
    }.get(str(status or ""), "")
