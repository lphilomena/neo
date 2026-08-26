from __future__ import annotations

"""Event-specific candidate source-chain confidence (C1-C4).

This module answers a narrower question than the peptide-HLA recommendation grade:
can the reported candidate peptide be traced back to a biologically plausible and
technically supported SNV, InDel, Fusion or Splice event?

C1-C4 are *source-chain confidence tiers* and must not be confused with R1-R4:

* C1: complete source chain with independent/cross-modal orthogonal support.
* C2: complete source chain with strong computational evidence but no orthogonal confirmation.
* C3: plausible source chain with one or more applicable unresolved/low-power requirements.
* C4: refuted or invalid source chain (hard failure).

The classifier is deterministic, reason-code driven and explicitly separates
NOT_APPLICABLE, UNASSESSED, INDETERMINATE_LOW_POWER, NEGATIVE and CONFLICT.
"""

from dataclasses import dataclass, asdict
import json
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence

from .evidence_states import (
    derive_event_authenticity,
    derive_mutant_specificity,
    derive_rna_support,
)
from .utils import read_tsv, write_tsv


SOURCE_CHAIN_RULE_VERSION = "source-chain-v1.0"

SUPPORTED = "SUPPORTED"
NOT_APPLICABLE = "NOT_APPLICABLE"
UNASSESSED = "UNASSESSED"
INDETERMINATE_LOW_POWER = "INDETERMINATE_LOW_POWER"
NEGATIVE = "NEGATIVE"
CONFLICT = "CONFLICT"

VALID_REQUIREMENT_STATUSES = {
    SUPPORTED,
    NOT_APPLICABLE,
    UNASSESSED,
    INDETERMINATE_LOW_POWER,
    NEGATIVE,
    CONFLICT,
}

SOURCE_CHAIN_TIERS = {
    "C1": {
        "grade": 3,
        "label": "ORTHOGONALLY_CONFIRMED_COMPLETE_SOURCE_CHAIN",
        "meaning": "Complete event-to-transcript/ORF-to-peptide chain with independent or cross-modal confirmation.",
    },
    "C2": {
        "grade": 2,
        "label": "STRONG_COMPUTATIONAL_COMPLETE_SOURCE_CHAIN",
        "meaning": "Complete source chain supported by primary computational/read evidence; orthogonal confirmation not yet available.",
    },
    "C3": {
        "grade": 1,
        "label": "PLAUSIBLE_INCOMPLETE_SOURCE_CHAIN",
        "meaning": "Plausible source chain with one or more applicable unresolved, low-power or negative-but-not-refuting requirements.",
    },
    "C4": {
        "grade": 0,
        "label": "REFUTED_OR_INVALID_SOURCE_CHAIN",
        "meaning": "Source event, transcript/ORF, peptide reconstruction or traceability is refuted or invalid.",
    },
}


@dataclass(frozen=True)
class RequirementAssessment:
    name: str
    label: str
    status: str
    reason_code: str
    reason: str
    source_fields: tuple[str, ...] = ()
    core: bool = True
    fatal_if_negative: bool = False
    fatal_if_conflict: bool = True

    def as_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["source_fields"] = list(self.source_fields)
        return payload


@dataclass(frozen=True)
class SourceChainResult:
    track: str
    tier: str
    label: str
    grade: int
    rule_version: str
    orthogonal_status: str
    orthogonal_sources: tuple[str, ...]
    requirements: tuple[RequirementAssessment, ...]
    reason_codes: tuple[str, ...]
    hard_failure: bool
    hard_failure_codes: tuple[str, ...]
    missing_requirements: tuple[str, ...]
    low_power_requirements: tuple[str, ...]
    negative_requirements: tuple[str, ...]
    conflict_requirements: tuple[str, ...]
    supported_requirements: tuple[str, ...]
    not_applicable_requirements: tuple[str, ...]

    def as_row(self) -> dict[str, str]:
        counts = {
            status: sum(req.status == status for req in self.requirements)
            for status in VALID_REQUIREMENT_STATUSES
        }
        requirement_map = {
            req.name: {
                "label": req.label,
                "status": req.status,
                "reason_code": req.reason_code,
                "reason": req.reason,
                "core": req.core,
                "fatal_if_negative": req.fatal_if_negative,
                "fatal_if_conflict": req.fatal_if_conflict,
                "source_fields": list(req.source_fields),
            }
            for req in self.requirements
        }
        row = {
            "source_chain_track": self.track,
            "source_chain_confidence_tier": self.tier,
            "source_chain_confidence_label": self.label,
            "source_chain_confidence_grade": str(self.grade),
            "source_chain_rule_version": self.rule_version,
            "source_chain_orthogonal_status": self.orthogonal_status,
            "source_chain_orthogonal_sources": ",".join(self.orthogonal_sources),
            "source_chain_hard_failure": "yes" if self.hard_failure else "no",
            "source_chain_hard_failure_codes": ",".join(self.hard_failure_codes),
            "source_chain_reason_codes": ",".join(self.reason_codes),
            "source_chain_confidence_reason_codes": ",".join(self.reason_codes),
            "source_chain_supported_requirements": ",".join(self.supported_requirements),
            "source_chain_missing_requirements": ",".join(self.missing_requirements),
            "source_chain_low_power_requirements": ",".join(self.low_power_requirements),
            "source_chain_negative_requirements": ",".join(self.negative_requirements),
            "source_chain_conflict_requirements": ",".join(self.conflict_requirements),
            "source_chain_not_applicable_requirements": ",".join(self.not_applicable_requirements),
            "source_chain_requirement_count": str(len(self.requirements)),
            "source_chain_supported_count": str(counts[SUPPORTED]),
            "source_chain_unassessed_count": str(counts[UNASSESSED]),
            "source_chain_low_power_count": str(counts[INDETERMINATE_LOW_POWER]),
            "source_chain_negative_count": str(counts[NEGATIVE]),
            "source_chain_conflict_count": str(counts[CONFLICT]),
            "source_chain_not_applicable_count": str(counts[NOT_APPLICABLE]),
            "source_chain_requirement_statuses": ";".join(
                f"{req.name}:{req.status}" for req in self.requirements
            ),
            "source_chain_requirement_details": json.dumps(
                requirement_map, sort_keys=True, ensure_ascii=True, separators=(",", ":")
            ),
        }
        requirement_aliases = {
            "event_authenticity_status": ("event_read_qc",),
            "orthogonal_confirmation_status": ("orthogonal_confirmation",),
            "transcript_orf_status": (
                "transcript_codon_protein", "orf_reconstruction",
                "fusion_transcript_orf", "splice_transcript_orf",
            ),
            "novel_sequence_status": ("novel_sequence", "novel_tail_nmd"),
            "phasing_status": ("read_backed_phasing",),
            "normal_background_status": ("normal_background",),
        }
        by_name = {req.name: req for req in self.requirements}
        for field, names in requirement_aliases.items():
            req = next((by_name[name] for name in names if name in by_name), None)
            row[field] = req.status if req else NOT_APPLICABLE
        for req in self.requirements:
            prefix = f"source_chain_{req.name}"
            row[f"{prefix}_applicability"] = (
                NOT_APPLICABLE if req.status == NOT_APPLICABLE else "APPLICABLE"
            )
            row[f"{prefix}_status"] = req.status
            row[f"{prefix}_value"] = req.reason
            row[f"{prefix}_reason_code"] = req.reason_code
            row[f"{prefix}_source"] = ",".join(req.source_fields)
            row[f"{prefix}_conflict"] = "yes" if req.status == CONFLICT else "no"
        return row


# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------

_EMPTY = {"", "NA", "N/A", "NONE", "NULL", ".", "NAN", "UNASSESSED", "UNRESOLVED", "NOT_AVAILABLE"}
_TRUE = {"TRUE", "YES", "Y", "1", "PASS", "SUPPORTED", "CONFIRMED", "VALID", "COMPLETE", "PRESENT"}
_FALSE = {"FALSE", "NO", "N", "0", "FAIL", "FAILED", "INVALID", "ABSENT", "REJECT", "REFUTED"}
_CONFLICT_TOKENS = ("CONFLICT", "DISCORDANT", "INCONSISTENT", "MISMATCH")
_LOW_POWER_TOKENS = ("LOW_COVERAGE", "LOW COVERAGE", "LOW_POWER", "LOW POWER", "INSUFFICIENT_COVERAGE", "INSUFFICIENT READ", "LOW_CONFIDENCE")
_NEGATIVE_TOKENS = ("NOT_DETECTED", "NO_SUPPORT", "NEGATIVE", "ABSENT", "UNSUPPORTED")
_INVALID_TOKENS = ("INVALID", "ARTIFACT", "WRONG", "IMPOSSIBLE", "REFUTED", "FAILED", "REJECT")
_PASS_TOKENS = ("PASS", "SUPPORTED", "CONFIRMED", "VALID", "COMPLETE", "CONCORDANT", "HIGH_CONFIDENCE", "HIGH CONFIDENCE", "IN_FRAME", "IN-FRAME")


def _section(rules: Mapping[str, Any], name: str) -> Mapping[str, Any]:
    value = rules.get(name, {}) if isinstance(rules, Mapping) else {}
    return value if isinstance(value, Mapping) else {}


def _raw(row: Mapping[str, Any], field: str) -> str:
    value = row.get(field, "")
    return "" if value is None else str(value).strip()


def _present(row: Mapping[str, Any], *fields: str) -> tuple[str, ...]:
    return tuple(field for field in fields if _raw(row, field).upper() not in _EMPTY)


def _text(row: Mapping[str, Any], *fields: str) -> str:
    return " ".join(_raw(row, field).upper() for field in fields if _raw(row, field).upper() not in _EMPTY)


def _number(row: Mapping[str, Any], *fields: str) -> float | None:
    for field in fields:
        value = _raw(row, field)
        if value.upper() in _EMPTY:
            continue
        try:
            return float(value)
        except ValueError:
            continue
    return None


def _truthy(row: Mapping[str, Any], *fields: str) -> bool | None:
    for field in fields:
        value = _raw(row, field).upper()
        if value in _EMPTY:
            continue
        if value in _TRUE or any(token in value for token in ("SUPPORTED", "CONFIRMED", "PASS")):
            return True
        if value in _FALSE or any(token in value for token in ("NOT_SUPPORTED", "NOT SUPPORTED", "FAILED", "INVALID")):
            return False
    return None


def _assessment(
    name: str,
    label: str,
    status: str,
    reason_code: str,
    reason: str,
    fields: Sequence[str] = (),
    *,
    core: bool = True,
    fatal_if_negative: bool = False,
    fatal_if_conflict: bool = True,
) -> RequirementAssessment:
    if status not in VALID_REQUIREMENT_STATUSES:
        raise ValueError(f"Unsupported source-chain requirement status: {status}")
    return RequirementAssessment(
        name=name,
        label=label,
        status=status,
        reason_code=reason_code,
        reason=reason,
        source_fields=tuple(fields),
        core=core,
        fatal_if_negative=fatal_if_negative,
        fatal_if_conflict=fatal_if_conflict,
    )


def _status_from_text(
    text: str,
    *,
    supported_code: str,
    negative_code: str,
    unassessed_code: str,
    conflict_code: str,
    low_power_code: str,
    invalid_is_negative: bool = True,
) -> tuple[str, str]:
    upper = text.upper().strip()
    if not upper:
        return UNASSESSED, unassessed_code
    if any(token in upper for token in _CONFLICT_TOKENS):
        return CONFLICT, conflict_code
    if any(token in upper for token in _LOW_POWER_TOKENS):
        return INDETERMINATE_LOW_POWER, low_power_code
    if any(token in upper for token in _INVALID_TOKENS):
        return (NEGATIVE if invalid_is_negative else CONFLICT), negative_code
    if any(token in upper for token in _NEGATIVE_TOKENS):
        return NEGATIVE, negative_code
    if any(token in upper for token in _PASS_TOKENS):
        return SUPPORTED, supported_code
    return UNASSESSED, unassessed_code


def source_chain_track(row: Mapping[str, Any]) -> str:
    """Return the upstream event source class, not the downstream peptide consequence."""

    explicit = _text(row, "source_chain_track", "candidate_source_track")
    for track in ("SNV", "INDEL", "FUSION", "SPLICE"):
        if explicit == track:
            return track

    source = _text(row, "mutation_source", "source_event_type", "variant_type")
    event = _text(row, "event_type", "consequence", "peptide_consequence")
    if "FUSION" in source or "FUSION" in event:
        return "FUSION"
    # A DNA SNV/InDel that causes a splice consequence remains a DNA-source candidate.
    if any(token in source for token in ("INDEL", "INSERTION", "DELETION", "DUPLICATION", "FRAMESHIFT")):
        return "INDEL"
    if any(token in source for token in ("SNV", "SNP", "MISSENSE", "SUBSTITUTION")):
        return "SNV"
    if "SPLICE" in source or any(token in source for token in ("JUNCTION", "EXON_SKIP", "INTRON_RETENTION")):
        return "SPLICE"
    if any(token in event for token in ("FRAMESHIFT", "INFRAME", "IN_FRAME", "INSERTION", "DELETION", "DUPLICATION", "INDEL")):
        return "INDEL"
    if any(token in event for token in ("SNV", "MISSENSE", "STOP_GAINED", "NONSENSE")):
        return "SNV"
    if "SPLICE" in event or "JUNCTION" in event:
        return "SPLICE"
    return "OTHER"


def _rule_threshold(rules: Mapping[str, Any], key: str, default: float) -> float:
    return float(_section(rules, "source_chain").get(key, default))


# ---------------------------------------------------------------------------
# Shared requirement evaluators
# ---------------------------------------------------------------------------


def _event_qc_requirement(row: Mapping[str, Any], rules: Mapping[str, Any], *, track: str) -> RequirementAssessment:
    state = derive_event_authenticity(row, rules)
    text = _text(row, "event_authenticity_state", "filter_status", "cross_platform_status", "comparison_status")
    fields = _present(row, "event_authenticity_state", "filter_status", "cross_platform_status", "comparison_status", "tumor_depth", "tumor_alt_count", "tumor_vaf")
    if state.get("hard_fail") or state["state"] == "EVENT_ARTIFACT_RISK":
        return _assessment("event_read_qc", "Event/read-level QC", NEGATIVE, str(state.get("hard_code") or "SC_EVENT_ARTIFACT"), str(state["reason"]), fields, fatal_if_negative=True)
    if state.get("conflict") or state["state"] == "EVENT_CONFLICT":
        return _assessment("event_read_qc", "Event/read-level QC", CONFLICT, "SC_EVENT_EVIDENCE_CONFLICT", str(state["reason"]), fields, fatal_if_conflict=True)
    if state["state"] in {"EVENT_CONFIRMED", "EVENT_STRONG"}:
        return _assessment("event_read_qc", "Event/read-level QC", SUPPORTED, f"SC_{track}_EVENT_QC_SUPPORTED", str(state["reason"]), fields, fatal_if_negative=True)
    if state["state"] in {"EVENT_PARTIAL", "EVENT_SAMPLE_SPECIFIC"}:
        return _assessment("event_read_qc", "Event/read-level QC", INDETERMINATE_LOW_POWER, f"SC_{track}_EVENT_QC_PARTIAL", str(state["reason"]), fields, fatal_if_negative=True)

    # Raw DNA fallback for SNV/InDel when no pre-derived event state is available.
    if track in {"SNV", "INDEL"}:
        depth = _number(row, "tumor_depth", "wes_tumor_depth", "wgs_tumor_depth")
        alt = _number(row, "tumor_alt_count", "wes_tumor_alt_count", "wgs_tumor_alt_count")
        filter_text = _text(row, "filter_status", "vcf_filter")
        min_depth = _rule_threshold(rules, "dna_min_depth", 20)
        min_alt = _rule_threshold(rules, "dna_min_alt_reads", 3)
        if depth is not None and depth < min_depth:
            return _assessment("event_read_qc", "Event/read-level QC", INDETERMINATE_LOW_POWER, "SC_DNA_LOW_DEPTH", f"tumor depth={depth:g} below {min_depth:g}", fields, fatal_if_negative=True)
        if alt is not None and alt < min_alt:
            return _assessment("event_read_qc", "Event/read-level QC", INDETERMINATE_LOW_POWER, "SC_DNA_LOW_ALT_READS", f"tumor ALT reads={alt:g} below {min_alt:g}", fields, fatal_if_negative=True)
        if depth is not None and alt is not None and depth >= min_depth and alt >= min_alt and (not filter_text or "PASS" in filter_text):
            return _assessment("event_read_qc", "Event/read-level QC", SUPPORTED, "SC_DNA_READ_QC_SUPPORTED", f"depth={depth:g}; ALT reads={alt:g}; filter={filter_text or 'not supplied'}", fields, fatal_if_negative=True)
    return _assessment("event_read_qc", "Event/read-level QC", UNASSESSED, f"SC_{track}_EVENT_QC_UNASSESSED", "No sufficient event/read-level QC evidence", fields, fatal_if_negative=True)


def _somatic_status_requirement(row: Mapping[str, Any], rules: Mapping[str, Any]) -> RequirementAssessment:
    del rules
    fields = _present(row, "matched_normal_status", "normal_support_status", "normal_depth", "normal_alt_count", "normal_alt_vaf", "germline_route_status")
    route = _text(row, "germline_route_status", "candidate_origin_route")
    if "GERMLINE" in route and any(token in route for token in ("REVIEWED", "SEPARATE", "APPROVED")):
        return _assessment("somatic_or_germline_route", "Somatic status or explicit germline route", SUPPORTED, "SC_GERMLINE_ROUTE_EXPLICIT", route, fields)
    status = _text(row, "matched_normal_status", "normal_support_status")
    normal_depth = _number(row, "normal_depth")
    normal_alt = _number(row, "normal_alt_count")
    normal_vaf = _number(row, "normal_alt_vaf")
    if "ALT_SUPPORTED" in status or (normal_alt is not None and normal_alt >= 2 and (normal_vaf or 0) >= 0.01):
        return _assessment("somatic_or_germline_route", "Somatic status or explicit germline route", NEGATIVE, "SC_MATCHED_NORMAL_ALT_SUPPORTED", f"matched-normal support: status={status}; depth={normal_depth}; ALT={normal_alt}; VAF={normal_vaf}", fields, fatal_if_negative=True)
    if any(token in status for token in _CONFLICT_TOKENS):
        return _assessment("somatic_or_germline_route", "Somatic status or explicit germline route", CONFLICT, "SC_NORMAL_STATUS_CONFLICT", status, fields, fatal_if_conflict=True)
    if normal_depth is not None:
        min_depth = 10
        if normal_depth < min_depth:
            return _assessment("somatic_or_germline_route", "Somatic status or explicit germline route", INDETERMINATE_LOW_POWER, "SC_NORMAL_LOW_COVERAGE", f"normal depth={normal_depth:g} below {min_depth}", fields)
        if (normal_alt or 0) == 0:
            return _assessment("somatic_or_germline_route", "Somatic status or explicit germline route", SUPPORTED, "SC_MATCHED_NORMAL_NEGATIVE", f"normal depth={normal_depth:g}; ALT={normal_alt or 0:g}", fields)
    if any(token in status for token in ("SOMATIC", "NORMAL_NEGATIVE", "TUMOR_ONLY_PASS", "PASS")):
        return _assessment("somatic_or_germline_route", "Somatic status or explicit germline route", SUPPORTED, "SC_SOMATIC_STATUS_SUPPORTED", status, fields)
    return _assessment("somatic_or_germline_route", "Somatic status or explicit germline route", UNASSESSED, "SC_SOMATIC_STATUS_UNASSESSED", "Matched-normal or germline-route evidence unavailable", fields)


def _rna_or_direct_requirement(row: Mapping[str, Any], rules: Mapping[str, Any], *, track: str) -> RequirementAssessment:
    fields = _present(
        row,
        "rna_support_status", "rna_alt_reads", "rna_vaf", "rna_depth", "rna_junction_reads", "junction_reads",
        "protein_evidence_status", "proteomics_status", "ligandome_evidence_status", "peptide_mass_spec_status",
    )
    protein = _text(row, "protein_evidence_status", "proteomics_status", "ligandome_evidence_status", "peptide_mass_spec_status")
    if any(token in protein for token in _PASS_TOKENS):
        return _assessment("rna_or_direct_evidence", "RNA mutant/junction or direct protein/peptide evidence", SUPPORTED, "SC_DIRECT_PROTEIN_OR_PEPTIDE_SUPPORTED", protein, fields)
    if any(token in protein for token in _CONFLICT_TOKENS):
        return _assessment("rna_or_direct_evidence", "RNA mutant/junction or direct protein/peptide evidence", CONFLICT, "SC_DIRECT_EVIDENCE_CONFLICT", protein, fields)
    state = derive_rna_support(row, rules)
    if state["state"] == "RNA_CONFIRMED":
        return _assessment("rna_or_direct_evidence", "RNA mutant/junction or direct protein/peptide evidence", SUPPORTED, f"SC_{track}_RNA_SUPPORTED", str(state["reason"]), fields)
    if state["state"] == "RNA_LOW_SUPPORT":
        return _assessment("rna_or_direct_evidence", "RNA mutant/junction or direct protein/peptide evidence", INDETERMINATE_LOW_POWER, f"SC_{track}_RNA_LOW_SUPPORT", str(state["reason"]), fields)
    if state["state"] == "GENE_EXPRESSION_ONLY":
        return _assessment("rna_or_direct_evidence", "RNA mutant/junction or direct protein/peptide evidence", UNASSESSED, "SC_GENE_EXPRESSION_ONLY", str(state["reason"]), fields)
    if state["state"] == "RNA_NEGATIVE":
        depth = _number(row, "rna_depth", "junction_depth", "rna_total_reads")
        min_depth = _rule_threshold(rules, "rna_min_evaluable_depth", 10)
        if depth is not None and depth < min_depth:
            return _assessment("rna_or_direct_evidence", "RNA mutant/junction or direct protein/peptide evidence", INDETERMINATE_LOW_POWER, "SC_RNA_NEGATIVE_LOW_POWER", f"RNA not detected with depth={depth:g} below {min_depth:g}", fields)
        return _assessment("rna_or_direct_evidence", "RNA mutant/junction or direct protein/peptide evidence", NEGATIVE, "SC_RNA_NOT_DETECTED", str(state["reason"]), fields, fatal_if_negative=False)
    return _assessment("rna_or_direct_evidence", "RNA mutant/junction or direct protein/peptide evidence", UNASSESSED, "SC_RNA_OR_DIRECT_EVIDENCE_UNASSESSED", str(state["reason"]), fields)


def _phasing_requirement(row: Mapping[str, Any]) -> RequirementAssessment:
    fields = _present(row, "phasing_required", "proximal_variant_count", "phase_group_id", "component_event_ids", "haplotype_status", "phase_confidence", "phase_support_reads")
    required = _truthy(row, "phasing_required")
    proximal = _number(row, "proximal_variant_count")
    components = [value for value in _raw(row, "component_event_ids").replace(";", ",").split(",") if value.strip()]
    if required is None:
        required = bool((proximal or 0) > 0 or len(components) > 1 or _raw(row, "phase_group_id"))
    if not required:
        return _assessment("read_backed_phasing", "Read-backed phasing for proximal variants", NOT_APPLICABLE, "SC_PHASING_NOT_APPLICABLE", "No proximal variant requiring phasing", fields, core=False, fatal_if_negative=True)
    text = _text(row, "haplotype_status", "phase_confidence", "phasing_status")
    if any(token in text for token in ("IMPOSSIBLE", "DISPROVED", "TRANS", "INCOMPATIBLE")):
        return _assessment("read_backed_phasing", "Read-backed phasing for proximal variants", NEGATIVE, "SC_PHASING_DISPROVES_PEPTIDE", text, fields, fatal_if_negative=True)
    if any(token in text for token in _CONFLICT_TOKENS):
        return _assessment("read_backed_phasing", "Read-backed phasing for proximal variants", CONFLICT, "SC_PHASING_CONFLICT", text, fields, fatal_if_conflict=True)
    if any(token in text for token in ("PHASED", "CIS", "HIGH_CONFIDENCE", "SUPPORTED")):
        return _assessment("read_backed_phasing", "Read-backed phasing for proximal variants", SUPPORTED, "SC_READ_BACKED_PHASING_SUPPORTED", text, fields, fatal_if_negative=True)
    if any(token in text for token in _LOW_POWER_TOKENS):
        return _assessment("read_backed_phasing", "Read-backed phasing for proximal variants", INDETERMINATE_LOW_POWER, "SC_PHASING_LOW_POWER", text, fields, fatal_if_negative=True)
    return _assessment("read_backed_phasing", "Read-backed phasing for proximal variants", UNASSESSED, "SC_PHASING_REQUIRED_UNASSESSED", text or "Proximal variants require read-backed phasing", fields, fatal_if_negative=True)


def _peptide_hla_traceability_requirement(row: Mapping[str, Any], *, require_orf: bool = False) -> RequirementAssessment:
    fields = _present(row, "peptide_id", "peptide", "mutant_peptide", "hla_allele", "event_id", "transcript_id", "transcript_hypothesis_id", "orf_id", "origin_peptide_id", "source_record_id", "generation_status")
    peptide = _raw(row, "peptide") or _raw(row, "mutant_peptide")
    hla = _raw(row, "hla_allele") or _raw(row, "hla") or _raw(row, "restricting_hla")
    event_id = _raw(row, "event_id")
    orf_link = _raw(row, "orf_id") or _raw(row, "transcript_hypothesis_id") or _raw(row, "transcript_id") or _raw(row, "origin_peptide_id")
    status = _text(row, "peptide_hla_traceability_status", "provenance_status", "generation_status", "junction_match_status")
    if any(token in status for token in _INVALID_TOKENS):
        return _assessment("peptide_hla_traceability", "Peptide-HLA traceability to event/transcript/ORF", NEGATIVE, "SC_PEPTIDE_HLA_TRACEABILITY_INVALID", status, fields, fatal_if_negative=True)
    if any(token in status for token in _CONFLICT_TOKENS):
        return _assessment("peptide_hla_traceability", "Peptide-HLA traceability to event/transcript/ORF", CONFLICT, "SC_PEPTIDE_HLA_TRACEABILITY_CONFLICT", status, fields, fatal_if_conflict=True)
    if peptide and hla and event_id and (orf_link or not require_orf):
        return _assessment("peptide_hla_traceability", "Peptide-HLA traceability to event/transcript/ORF", SUPPORTED, "SC_PEPTIDE_HLA_TRACEABLE", f"event={event_id}; peptide={peptide}; HLA={hla}; ORF/transcript={orf_link or 'not required'}", fields, fatal_if_negative=True)
    missing = [name for name, value in (("peptide", peptide), ("HLA", hla), ("event_id", event_id), ("ORF/transcript", orf_link if require_orf else "not_required")) if not value]
    return _assessment("peptide_hla_traceability", "Peptide-HLA traceability to event/transcript/ORF", UNASSESSED, "SC_PEPTIDE_HLA_TRACEABILITY_INCOMPLETE", "missing=" + ",".join(missing), fields, fatal_if_negative=True)


def _novel_sequence_requirement(row: Mapping[str, Any], *, track: str) -> RequirementAssessment:
    fields = _present(row, "mutant_specificity_status", "mutant_specificity_gate_status", "contains_novel_aa", "crosses_junction", "mutation_positions_in_peptide", "novel_tail_status", "peptide", "wildtype_peptide")
    specificity = derive_mutant_specificity(row, {})
    if specificity.get("hard_fail") or specificity["state"] == "NON_MUTANT_SEQUENCE":
        return _assessment("novel_sequence", "Candidate peptide contains the event-derived novel sequence", NEGATIVE, "SC_PEPTIDE_NOT_EVENT_DERIVED", str(specificity["reason"]), fields, fatal_if_negative=True)
    if specificity["state"] == "WT_BETTER" and track in {"SNV", "INDEL"}:
        # WT-better is a recommendation-level problem; the peptide can still be correctly reconstructed.
        mutation_positions = _raw(row, "mutation_positions_in_peptide")
        if mutation_positions:
            return _assessment("novel_sequence", "Candidate peptide contains the event-derived novel sequence", SUPPORTED, "SC_PEPTIDE_CONTAINS_MUTATION_WT_BETTER", f"mutation positions={mutation_positions}; WT presentation may be better", fields, fatal_if_negative=True)
    explicit_novel = _truthy(row, "contains_novel_aa", "crosses_junction")
    mutation_positions = _raw(row, "mutation_positions_in_peptide")
    if explicit_novel is True and track == "FUSION":
        boundary_position = _raw(row, "junction_position_in_peptide_1based") or _raw(row, "junction_offset_in_peptide")
        left_peptide = _raw(row, "fusion_left_peptide")
        right_peptide = _raw(row, "fusion_right_peptide")
        if not (boundary_position and left_peptide and right_peptide):
            return _assessment("novel_sequence", "Candidate peptide contains the event-derived novel sequence", NEGATIVE, "SC_FUSION_BOUNDARY_MAPPING_INCOMPLETE", "Fusion peptide was labelled junction-spanning but residue-level left/right boundary mapping is incomplete", fields, fatal_if_negative=True)
    if explicit_novel is True or mutation_positions:
        return _assessment("novel_sequence", "Candidate peptide contains the event-derived novel sequence", SUPPORTED, "SC_NOVEL_SEQUENCE_PRESENT", f"contains_novel={explicit_novel}; mutation_positions={mutation_positions}", fields, fatal_if_negative=True)
    if explicit_novel is False and track in {"FUSION", "SPLICE"}:
        return _assessment("novel_sequence", "Candidate peptide contains the event-derived novel sequence", NEGATIVE, "SC_JUNCTION_PEPTIDE_DOES_NOT_CROSS_BOUNDARY", "Peptide is explicitly mapped to only one side of the junction", fields, fatal_if_negative=True)
    if specificity["state"] in {"MT_SPECIFIC", "MARGINAL_MT_ADVANTAGE", "MT_WT_SIMILAR", "WT_BETTER"}:
        return _assessment("novel_sequence", "Candidate peptide contains the event-derived novel sequence", SUPPORTED, "SC_MUTANT_PEPTIDE_RECONSTRUCTED", str(specificity["reason"]), fields, fatal_if_negative=True)
    return _assessment("novel_sequence", "Candidate peptide contains the event-derived novel sequence", UNASSESSED, "SC_NOVEL_SEQUENCE_UNASSESSED", "No explicit mutation/junction/novel-tail mapping", fields, fatal_if_negative=True)


def _orthogonal_confirmation(row: Mapping[str, Any], *, track: str, rna_requirement: RequirementAssessment) -> tuple[str, tuple[str, ...], str]:
    sources: list[str] = []
    status_text = _text(
        row,
        "orthogonal_confirmation_status", "targeted_validation_status", "rt_pcr_status", "sanger_status",
        "long_read_status", "independent_library_status", "dna_sv_confirmation_status", "protein_evidence_status",
        "proteomics_status", "ligandome_evidence_status", "cross_platform_status",
        "cross_site_status", "cross_site_exact_support", "sample_identity_status",
    )
    if any(token in status_text for token in ("REFUTED", "FAILED_WITH_ADEQUATE_POWER", "NEGATIVE_WITH_ADEQUATE")):
        return NEGATIVE, tuple(sources), "SC_ORTHOGONAL_REFUTED"

    # Explicit non-confirmation / not-performed states must be handled before
    # positive token matching.  Otherwise strings such as ``NOT_CONFIRMED``
    # would contain the substring ``CONFIRMED`` and be misclassified as
    # orthogonal support.
    explicit_not_performed_tokens = (
        "NOT_CONFIRMED",
        "NOT CONFIRMED",
        "NOT_PERFORMED",
        "NOT PERFORMED",
        "NOT_TESTED",
        "NOT TESTED",
        "NO_ORTHOGONAL_CONFIRMATION",
        "NO ORTHOGONAL CONFIRMATION",
    )
    explicit_not_performed = any(token in status_text for token in explicit_not_performed_tokens)

    fields_to_labels = {
        "rt_pcr_status": "RT-PCR",
        "sanger_status": "Sanger",
        "long_read_status": "long-read",
        "independent_library_status": "independent-library",
        "dna_sv_confirmation_status": "DNA-SV",
        "protein_evidence_status": "proteomics",
        "proteomics_status": "proteomics",
        "ligandome_evidence_status": "ligandome",
        "targeted_validation_status": "targeted-validation",
        "orthogonal_confirmation_status": "orthogonal-assay",
    }
    for field, label in fields_to_labels.items():
        field_text = _raw(row, field).upper()
        if not field_text:
            continue
        if any(token in field_text for token in explicit_not_performed_tokens):
            continue
        if any(token in field_text for token in _PASS_TOKENS):
            sources.append(label)

    cross = _text(row, "cross_platform_status", "comparison_status")
    cross_site = _text(row, "cross_site_status")
    if "EXACT_SHARED" in cross_site and _truthy(row, "cross_site_exact_support") is True:
        sources.append("independent-tumor-site-RNA")
    if track in {"SNV", "INDEL"}:
        if "CROSS_PLATFORM_PASS_CONCORDANT" in cross:
            sources.append("WES/WGS-cross-platform")
        # DNA event plus RNA allele/protein evidence is cross-modal confirmation.
        if rna_requirement.status == SUPPORTED:
            sources.append("DNA/RNA-or-protein-cross-modal")
    # Multiple callers on the same RNA BAM are computational consistency, not orthogonal.

    if sources:
        return SUPPORTED, tuple(dict.fromkeys(sources)), "SC_ORTHOGONAL_SUPPORTED"
    if explicit_not_performed:
        return UNASSESSED, tuple(), "SC_ORTHOGONAL_NOT_PERFORMED"
    if any(token in status_text for token in _CONFLICT_TOKENS):
        return CONFLICT, tuple(), "SC_ORTHOGONAL_CONFLICT"
    if any(token in status_text for token in _LOW_POWER_TOKENS):
        return INDETERMINATE_LOW_POWER, tuple(), "SC_ORTHOGONAL_LOW_POWER"
    return UNASSESSED, tuple(), "SC_ORTHOGONAL_NOT_PERFORMED"


# ---------------------------------------------------------------------------
# Track-specific evaluators
# ---------------------------------------------------------------------------


def _snv_requirements(row: Mapping[str, Any], rules: Mapping[str, Any]) -> list[RequirementAssessment]:
    event_qc = _event_qc_requirement(row, rules, track="SNV")
    somatic = _somatic_status_requirement(row, rules)
    rna = _rna_or_direct_requirement(row, rules, track="SNV")
    transcript_fields = _present(row, "transcript_id", "codon_change", "protein_change", "combined_protein_change", "event_name", "consequence")
    transcript = _raw(row, "transcript_id")
    protein = _raw(row, "protein_change") or _raw(row, "combined_protein_change") or _raw(row, "event_name")
    consequence = _text(row, "consequence", "peptide_consequence")
    transcript_req = _assessment(
        "transcript_codon_protein",
        "Affected transcript, codon and amino-acid change are resolved",
        SUPPORTED if transcript and (protein or "MISSENSE" in consequence or "STOP" in consequence) else UNASSESSED,
        "SC_SNV_TRANSCRIPT_CODON_RESOLVED" if transcript and (protein or consequence) else "SC_SNV_TRANSCRIPT_CODON_UNASSESSED",
        f"transcript={transcript or 'NA'}; protein/codon={protein or consequence or 'NA'}",
        transcript_fields,
        fatal_if_negative=True,
    )
    phasing = _phasing_requirement(row)
    novel = _novel_sequence_requirement(row, track="SNV")
    trace = _peptide_hla_traceability_requirement(row, require_orf=False)
    sequence_qc = _snv_sequence_qc_requirement(row, rules)
    return [event_qc, sequence_qc, somatic, rna, transcript_req, phasing, novel, trace]


def _snv_sequence_qc_requirement(
    row: Mapping[str, Any], rules: Mapping[str, Any]
) -> RequirementAssessment:
    fields = _present(
        row,
        "base_quality", "mean_alt_base_quality", "mapping_quality", "mean_alt_mapping_quality",
        "strand_bias_status", "ffpe_artifact_status", "low_complexity_status",
        "paralogous_region_status", "mapping_artifact_status",
    )
    text = _text(
        row, "strand_bias_status", "ffpe_artifact_status", "low_complexity_status",
        "paralogous_region_status", "mapping_artifact_status",
    )
    if any(token in text for token in (
        "FFPE_ARTIFACT", "STRAND_BIAS_FAIL", "LOW_COMPLEXITY_ARTIFACT",
        "PARALOG_ARTIFACT", "MAPPING_ARTIFACT", "PSEUDOGENE_ARTIFACT", "INVALID",
    )):
        return _assessment(
            "snv_sequence_qc", "Base/MAP quality, strand bias, FFPE and mappability QC",
            NEGATIVE, "SC_SNV_SEQUENCE_QC_ARTIFACT", text, fields, fatal_if_negative=True,
        )
    if any(token in text for token in _CONFLICT_TOKENS):
        return _assessment(
            "snv_sequence_qc", "Base/MAP quality, strand bias, FFPE and mappability QC",
            CONFLICT, "SC_SNV_SEQUENCE_QC_CONFLICT", text, fields, fatal_if_conflict=True,
        )
    baseq = _number(row, "mean_alt_base_quality", "base_quality")
    mapq = _number(row, "mean_alt_mapping_quality", "mapping_quality")
    min_baseq = _rule_threshold(rules, "snv_min_base_quality", 20)
    min_mapq = _rule_threshold(rules, "snv_min_mapping_quality", 20)
    if (
        any(token in text for token in _LOW_POWER_TOKENS)
        or (baseq is not None and baseq < min_baseq)
        or (mapq is not None and mapq < min_mapq)
    ):
        return _assessment(
            "snv_sequence_qc", "Base/MAP quality, strand bias, FFPE and mappability QC",
            INDETERMINATE_LOW_POWER, "SC_SNV_SEQUENCE_QC_LOW_POWER",
            f"baseQ={baseq}; MAPQ={mapq}; {text}", fields, fatal_if_negative=True,
        )
    explicit_pass = any(token in text for token in _PASS_TOKENS)
    if (baseq is not None and baseq >= min_baseq and mapq is not None and mapq >= min_mapq) or explicit_pass:
        return _assessment(
            "snv_sequence_qc", "Base/MAP quality, strand bias, FFPE and mappability QC",
            SUPPORTED, "SC_SNV_SEQUENCE_QC_PASS",
            f"baseQ={baseq}; MAPQ={mapq}; {text or 'quality thresholds passed'}", fields,
            fatal_if_negative=True,
        )
    return _assessment(
        "snv_sequence_qc", "Base/MAP quality, strand bias, FFPE and mappability QC",
        UNASSESSED, "SC_SNV_SEQUENCE_QC_UNASSESSED",
        "Base/MAP quality, FFPE and mappability evidence unavailable", fields,
        fatal_if_negative=True,
    )


def _indel_normalization_requirement(row: Mapping[str, Any]) -> RequirementAssessment:
    fields = _present(
        row, "variant_normalization_status", "left_normalized", "minimal_representation_status",
        "normalized_variant_key", "variant_key", "biological_event_id_before_normalization",
        "biological_event_id_after_normalization", "normalization_equivalence_status",
    )
    text = _text(row, "variant_normalization_status", "minimal_representation_status", "left_normalized")
    before = _raw(row, "biological_event_id_before_normalization")
    after = _raw(row, "biological_event_id_after_normalization")
    equivalence = _text(row, "normalization_equivalence_status")
    if (before and after and before != after) or any(token in equivalence for token in ("NOT_EQUIVALENT", "IDENTITY_CHANGED", "DIFFERENT_EVENT")):
        return _assessment("normalized_representation", "Left-normalized minimal InDel representation", NEGATIVE, "SC_INDEL_NORMALIZATION_CHANGED_EVENT_IDENTITY", f"before={before}; after={after}; {equivalence}", fields, fatal_if_negative=True)
    if any(token in text for token in ("WRONG_EVENT", "IDENTITY_CHANGED", "INVALID")):
        return _assessment("normalized_representation", "Left-normalized minimal InDel representation", NEGATIVE, "SC_INDEL_REPRESENTATION_INVALID", text, fields, fatal_if_negative=True)
    if any(token in text for token in _CONFLICT_TOKENS):
        return _assessment("normalized_representation", "Left-normalized minimal InDel representation", CONFLICT, "SC_INDEL_REPRESENTATION_CONFLICT", text, fields, fatal_if_conflict=True)
    if (before and after and before == after) or "EQUIVALENT" in equivalence or _raw(row, "normalized_variant_key") or any(token in text for token in ("NORMALIZED", "MINIMAL", "PASS", "TRUE")):
        return _assessment("normalized_representation", "Left-normalized minimal InDel representation", SUPPORTED, "SC_INDEL_REPRESENTATION_NORMALIZED", text or _raw(row, "normalized_variant_key"), fields, fatal_if_negative=True)
    return _assessment("normalized_representation", "Left-normalized minimal InDel representation", UNASSESSED, "SC_INDEL_REPRESENTATION_UNASSESSED", "No normalization audit field", fields, fatal_if_negative=True)


def _indel_local_context_requirement(row: Mapping[str, Any]) -> RequirementAssessment:
    fields = _present(row, "local_realign_status", "repeat_context_status", "homopolymer_status", "microhomology_status", "local_assembly_status")
    text = _text(row, "local_realign_status", "repeat_context_status", "homopolymer_status", "microhomology_status", "local_assembly_status")
    if any(token in text for token in ("ARTIFACT", "FAILED_REALIGN", "NOT_SUPPORTED_AFTER_REALIGN", "INVALID")):
        return _assessment("local_context_qc", "Local realignment/repeat/homopolymer/microhomology QC", NEGATIVE, "SC_INDEL_LOCAL_REALIGNMENT_REFUTES", text, fields, fatal_if_negative=True)
    if any(token in text for token in _CONFLICT_TOKENS):
        return _assessment("local_context_qc", "Local realignment/repeat/homopolymer/microhomology QC", CONFLICT, "SC_INDEL_LOCAL_CONTEXT_CONFLICT", text, fields, fatal_if_conflict=True)
    if any(token in text for token in _LOW_POWER_TOKENS) or any(token in text for token in ("REPEAT_RISK", "HOMOPOLYMER_RISK", "MICROHOMOLOGY_RISK")):
        return _assessment("local_context_qc", "Local realignment/repeat/homopolymer/microhomology QC", INDETERMINATE_LOW_POWER, "SC_INDEL_LOCAL_CONTEXT_CAUTION", text, fields, fatal_if_negative=True)
    if any(token in text for token in _PASS_TOKENS):
        return _assessment("local_context_qc", "Local realignment/repeat/homopolymer/microhomology QC", SUPPORTED, "SC_INDEL_LOCAL_CONTEXT_QC_PASS", text, fields, fatal_if_negative=True)
    return _assessment("local_context_qc", "Local realignment/repeat/homopolymer/microhomology QC", UNASSESSED, "SC_INDEL_LOCAL_CONTEXT_UNASSESSED", "Local realignment/repeat-context audit unavailable", fields, fatal_if_negative=True)


def _indel_orf_requirement(row: Mapping[str, Any]) -> RequirementAssessment:
    fields = _present(row, "frame_status", "rna_frame_status", "reading_frame", "orf_status", "orf_id", "novel_tail_status", "stop_position", "nmd_risk_status", "protein_change", "combined_protein_change")
    text = _text(row, "frame_status", "rna_frame_status", "reading_frame", "orf_status", "novel_tail_status")
    if any(token in text for token in ("FRAME_ERROR", "WRONG_FRAME", "INVALID_ORF", "ORF_INVALID", "TRANSLATION_ERROR")):
        return _assessment("orf_reconstruction", "In-frame/frameshift ORF reconstruction", NEGATIVE, "SC_INDEL_ORF_INVALID", text, fields, fatal_if_negative=True)
    consequence = _text(row, "event_type", "consequence", "peptide_consequence", "protein_change", "combined_protein_change")
    if any(token in text for token in ("VALID", "CONFIRMED", "IN_FRAME", "FRAME_OK")) or _raw(row, "orf_id") or any(token in consequence for token in ("FRAMESHIFT", "INFRAME", "IN_FRAME")):
        return _assessment("orf_reconstruction", "In-frame/frameshift ORF reconstruction", SUPPORTED, "SC_INDEL_ORF_RECONSTRUCTED", text or consequence, fields, fatal_if_negative=True)
    return _assessment("orf_reconstruction", "In-frame/frameshift ORF reconstruction", UNASSESSED, "SC_INDEL_ORF_UNASSESSED", "Reading frame/ORF reconstruction unavailable", fields, fatal_if_negative=True)


def _novel_tail_nmd_requirement(row: Mapping[str, Any]) -> RequirementAssessment:
    fields = _present(row, "novel_tail_status", "contains_novel_aa", "stop_position", "premature_stop_position", "nmd_risk_status", "nmd_prediction")
    text = _text(row, "novel_tail_status", "nmd_risk_status", "nmd_prediction")
    contains = _truthy(row, "contains_novel_aa")
    if any(token in text for token in ("INVALID", "WRONG_TAIL", "NOT_IN_ORF")):
        return _assessment("novel_tail_nmd", "Novel tail, stop position and NMD annotation", NEGATIVE, "SC_INDEL_NOVEL_TAIL_INVALID", text, fields, fatal_if_negative=True)
    if contains is True and (_raw(row, "nmd_risk_status") or _raw(row, "nmd_prediction") or _raw(row, "stop_position") or _raw(row, "premature_stop_position")):
        return _assessment("novel_tail_nmd", "Novel tail, stop position and NMD annotation", SUPPORTED, "SC_INDEL_NOVEL_TAIL_NMD_ANNOTATED", text or "novel tail and stop/NMD annotated", fields)
    if contains is True:
        return _assessment("novel_tail_nmd", "Novel tail, stop position and NMD annotation", INDETERMINATE_LOW_POWER, "SC_INDEL_NOVEL_TAIL_NMD_INCOMPLETE", "Novel tail present; stop/NMD annotation incomplete", fields)
    # In-frame InDels do not require a novel tail/NMD analysis.
    consequence = _text(row, "event_type", "consequence", "peptide_consequence")
    if "FRAMESHIFT" not in consequence:
        return _assessment("novel_tail_nmd", "Novel tail, stop position and NMD annotation", NOT_APPLICABLE, "SC_INDEL_NMD_NOT_APPLICABLE", "Non-frameshift InDel", fields, core=False)
    return _assessment("novel_tail_nmd", "Novel tail, stop position and NMD annotation", UNASSESSED, "SC_INDEL_NOVEL_TAIL_NMD_UNASSESSED", "Frameshift novel-tail/NMD evidence unavailable", fields)


def _indel_requirements(row: Mapping[str, Any], rules: Mapping[str, Any]) -> list[RequirementAssessment]:
    return [
        _indel_normalization_requirement(row),
        _indel_local_context_requirement(row),
        _event_qc_requirement(row, rules, track="INDEL"),
        _somatic_status_requirement(row, rules),
        _rna_or_direct_requirement(row, rules, track="INDEL"),
        _indel_orf_requirement(row),
        _novel_tail_nmd_requirement(row),
        _phasing_requirement(row),
        _novel_sequence_requirement(row, track="INDEL"),
        _peptide_hla_traceability_requirement(row, require_orf=True),
    ]


def _fusion_breakpoint_requirement(row: Mapping[str, Any]) -> RequirementAssessment:
    fields = _present(row, "breakpoint1", "breakpoint2", "gene5", "gene3", "fusion_gene", "fusion_name", "strand1", "strand2", "orientation_status", "breakpoint_status")
    text = _text(row, "breakpoint_status", "orientation_status", "fusion_direction_status")
    gene_pair = _raw(row, "fusion_gene") or _raw(row, "fusion_name") or (f"{_raw(row, 'gene5')}::{_raw(row, 'gene3')}" if _raw(row, "gene5") and _raw(row, "gene3") else "")
    breakpoints = bool(_raw(row, "breakpoint1") and _raw(row, "breakpoint2")) or bool(_raw(row, "junction_chrom") and _raw(row, "junction_start") and _raw(row, "junction_end"))
    if any(token in text for token in ("WRONG", "INVALID", "REVERSED", "IMPOSSIBLE")):
        return _assessment("breakpoint_definition", "Fusion partners, exact breakpoints, orientation and strand", NEGATIVE, "SC_FUSION_BREAKPOINT_ORIENTATION_INVALID", text, fields, fatal_if_negative=True)
    if any(token in text for token in _CONFLICT_TOKENS):
        return _assessment("breakpoint_definition", "Fusion partners, exact breakpoints, orientation and strand", CONFLICT, "SC_FUSION_BREAKPOINT_CONFLICT", text, fields, fatal_if_conflict=True)
    if gene_pair and breakpoints:
        return _assessment("breakpoint_definition", "Fusion partners, exact breakpoints, orientation and strand", SUPPORTED, "SC_FUSION_BREAKPOINT_DEFINED", f"fusion={gene_pair}; breakpoints defined", fields, fatal_if_negative=True)
    if gene_pair:
        return _assessment("breakpoint_definition", "Fusion partners, exact breakpoints, orientation and strand", UNASSESSED, "SC_FUSION_GENE_PAIR_ONLY", f"gene pair={gene_pair}; exact breakpoint unavailable", fields, fatal_if_negative=True)
    return _assessment("breakpoint_definition", "Fusion partners, exact breakpoints, orientation and strand", UNASSESSED, "SC_FUSION_BREAKPOINT_UNASSESSED", "Fusion pair/breakpoint unavailable", fields, fatal_if_negative=True)


def _fusion_split_requirement(row: Mapping[str, Any], rules: Mapping[str, Any]) -> RequirementAssessment:
    fields = _present(
        row, "split_reads", "rna_junction_reads", "junction_reads", "unique_split_reads",
        "junction_support_status", "anchor_size", "duplicate_removed_status",
        "junction_sequence_uniqueness_status", "independent_start_sites", "unique_molecules",
    )
    text = _text(
        row, "junction_support_status", "breakpoint_support_status",
        "duplicate_removed_status", "junction_sequence_uniqueness_status",
    )
    if any(token in text for token in ("ARTIFACT", "NOT_SUPPORTED", "INVALID", "NON_UNIQUE", "DUPLICATE_ONLY")):
        return _assessment("split_read_support", "Unique split-read support for the fusion breakpoint", NEGATIVE, "SC_FUSION_SPLIT_READ_REFUTED", text, fields, fatal_if_negative=True)
    reads = _number(row, "unique_split_reads", "split_reads", "rna_junction_reads", "junction_reads")
    minimum = _rule_threshold(rules, "junction_min_reads", _rule_threshold(rules, "fusion_min_split_reads", 3))
    strong = _rule_threshold(rules, "junction_strong_reads", _rule_threshold(rules, "fusion_strong_split_reads", 10))
    if reads is None:
        return _assessment("split_read_support", "Unique split-read support for the fusion breakpoint", UNASSESSED, "SC_FUSION_SPLIT_READ_UNASSESSED", "Split-read count unavailable", fields, fatal_if_negative=True)
    if reads >= minimum and not any(token in text for token in _LOW_POWER_TOKENS):
        return _assessment("split_read_support", "Unique split-read support for the fusion breakpoint", SUPPORTED, "SC_FUSION_SPLIT_READ_SUPPORTED", f"split/junction reads={reads:g}; strong threshold={strong:g}", fields, fatal_if_negative=True)
    if reads > 0:
        return _assessment("split_read_support", "Unique split-read support for the fusion breakpoint", INDETERMINATE_LOW_POWER, "SC_FUSION_SPLIT_READ_LOW", f"split/junction reads={reads:g} below {minimum:g}", fields, fatal_if_negative=True)
    return _assessment("split_read_support", "Unique split-read support for the fusion breakpoint", NEGATIVE, "SC_FUSION_SPLIT_READ_ZERO", "No split/junction read support", fields, fatal_if_negative=True)


def _fusion_supplementary_structure_requirement(row: Mapping[str, Any], rules: Mapping[str, Any]) -> RequirementAssessment:
    fields = _present(row, "rna_spanning_reads", "spanning_pairs", "tools_detected", "caller_support_frac", "independent_start_sites", "dna_sv_confirmation_status")
    spanning = _number(row, "rna_spanning_reads", "spanning_pairs")
    callers = _number(row, "tools_detected", "caller_count")
    independent = _number(row, "independent_start_sites", "unique_molecules")
    split = _number(row, "unique_split_reads", "split_reads", "rna_junction_reads", "junction_reads")
    strong = _rule_threshold(rules, "fusion_strong_split_reads", 10)
    if (spanning or 0) > 0 or (callers or 0) >= 2 or (independent or 0) >= 2 or (split or 0) >= strong or any(token in _text(row, "dna_sv_confirmation_status") for token in _PASS_TOKENS):
        return _assessment("supplementary_structure_support", "Spanning pairs/independent reads or supplementary structural evidence", SUPPORTED, "SC_FUSION_SUPPLEMENTARY_STRUCTURE_SUPPORTED", f"spanning={spanning}; callers={callers}; independent={independent}; split={split}", fields)
    if spanning is not None or callers is not None or independent is not None or split is not None:
        return _assessment("supplementary_structure_support", "Spanning pairs/independent reads or supplementary structural evidence", INDETERMINATE_LOW_POWER, "SC_FUSION_SUPPLEMENTARY_STRUCTURE_LOW", f"spanning={spanning}; callers={callers}; independent={independent}; split={split}", fields)
    return _assessment("supplementary_structure_support", "Spanning pairs/independent reads or supplementary structural evidence", UNASSESSED, "SC_FUSION_SUPPLEMENTARY_STRUCTURE_UNASSESSED", "No spanning/caller/structural support metrics", fields)


def _normal_background_requirement(row: Mapping[str, Any], *, track: str, mandatory: bool = True) -> RequirementAssessment:
    fields = _present(
        row, "normal_junction_assessment_status", "normal_junction_status", "readthrough_status",
        "mapping_artifact_status", "homology_artifact_status", "normal_junction_seen",
        "normal_junction_depth", "normal_junction_coverage", "normal_gene_coverage",
        "adjacent_gene_status", "liver_expression_context",
    )
    text = _text(row, "normal_junction_assessment_status", "normal_junction_status", "readthrough_status", "mapping_artifact_status", "homology_artifact_status")
    explicit_artifact = any(token in text for token in ("READTHROUGH_CONFIRMED", "MAPPING_ARTIFACT", "HOMOLOGY_ARTIFACT"))
    normal_detected = any(token in text for token in ("SUPPORTED_IN_NORMAL", "NORMAL_JUNCTION_DETECTED", "EXACT_MATCH")) and "NOT_DETECTED" not in text
    if explicit_artifact or normal_detected:
        return _assessment("normal_background", "Normal-junction/read-through/mapping-artifact exclusion", NEGATIVE, f"SC_{track}_NORMAL_BACKGROUND_REFUTES", text, fields, fatal_if_negative=True)
    if any(token in text for token in _CONFLICT_TOKENS):
        return _assessment("normal_background", "Normal-junction/read-through/mapping-artifact exclusion", CONFLICT, f"SC_{track}_NORMAL_BACKGROUND_CONFLICT", text, fields, fatal_if_conflict=True)
    if any(token in text for token in _LOW_POWER_TOKENS) or "NOT_DETECTED_LOW_COVERAGE" in text:
        return _assessment("normal_background", "Normal-junction/read-through/mapping-artifact exclusion", INDETERMINATE_LOW_POWER, f"SC_{track}_NORMAL_BACKGROUND_LOW_POWER", text, fields, core=mandatory, fatal_if_negative=True)
    adequate = any(token in text for token in ("NOT_DETECTED_ADEQUATE_COVERAGE", "ADEQUATE_COVERAGE_NEGATIVE", "EXCLUDED_WITH_ADEQUATE_COVERAGE"))
    coverage = _number(row, "normal_junction_depth", "normal_junction_coverage", "normal_gene_coverage")
    if adequate or (coverage is not None and coverage >= 10 and any(token in text for token in ("NOT_DETECTED", "ABSENT", "NO_MATCH", "EXCLUDED"))):
        return _assessment("normal_background", "Normal-junction/read-through/mapping-artifact exclusion", SUPPORTED, f"SC_{track}_NORMAL_BACKGROUND_CLEAR", text, fields, fatal_if_negative=True)
    if any(token in text for token in ("NOT_DETECTED", "ABSENT", "NO_MATCH", "EXCLUDED", "PASS")):
        return _assessment("normal_background", "Normal-junction/read-through/mapping-artifact exclusion", INDETERMINATE_LOW_POWER, f"SC_{track}_NORMAL_BACKGROUND_COVERAGE_UNASSESSED", f"{text}; normal coverage={coverage}", fields, core=mandatory, fatal_if_negative=True)
    return _assessment("normal_background", "Normal-junction/read-through/mapping-artifact exclusion", UNASSESSED, f"SC_{track}_NORMAL_BACKGROUND_UNASSESSED", "Normal-junction/read-through/background filter unavailable", fields, core=mandatory, fatal_if_negative=True)


def _fusion_frame_orf_requirement(row: Mapping[str, Any]) -> RequirementAssessment:
    fields = _present(row, "frame_status", "rna_frame_status", "bp1_frame", "bp2_frame", "gene5", "gene3", "exon_boundary", "fusion_protein_sequence", "orf_id", "transcript_hypothesis_id", "orientation_status")
    text = _text(row, "frame_status", "rna_frame_status", "orientation_status", "orf_status")
    if any(token in text for token in ("WRONG", "INVALID", "IMPOSSIBLE", "FRAME_ERROR", "ORF_INVALID")):
        return _assessment("fusion_transcript_orf", "5'/3' orientation, exon connection, reading frame and fusion ORF", NEGATIVE, "SC_FUSION_ORF_INVALID", text, fields, fatal_if_negative=True)
    if any(token in text for token in _CONFLICT_TOKENS):
        return _assessment("fusion_transcript_orf", "5'/3' orientation, exon connection, reading frame and fusion ORF", CONFLICT, "SC_FUSION_ORF_CONFLICT", text, fields, fatal_if_conflict=True)
    if any(token in text for token in ("IN_FRAME", "IN-FRAME", "INFRAME", "FRAME_CONFIRMED", "VALID")) and (_raw(row, "fusion_protein_sequence") or _raw(row, "orf_id") or _raw(row, "transcript_hypothesis_id")):
        return _assessment("fusion_transcript_orf", "5'/3' orientation, exon connection, reading frame and fusion ORF", SUPPORTED, "SC_FUSION_ORF_RECONSTRUCTED", text, fields, fatal_if_negative=True)
    if any(token in text for token in ("IN_FRAME", "IN-FRAME", "INFRAME", "FRAME_CONFIRMED")):
        return _assessment("fusion_transcript_orf", "5'/3' orientation, exon connection, reading frame and fusion ORF", INDETERMINATE_LOW_POWER, "SC_FUSION_FRAME_ONLY_ORF_INCOMPLETE", text, fields, fatal_if_negative=True)
    return _assessment("fusion_transcript_orf", "5'/3' orientation, exon connection, reading frame and fusion ORF", UNASSESSED, "SC_FUSION_ORF_UNASSESSED", "Fusion orientation/frame/ORF incomplete", fields, fatal_if_negative=True)


def _fusion_requirements(row: Mapping[str, Any], rules: Mapping[str, Any]) -> list[RequirementAssessment]:
    return [
        _fusion_breakpoint_requirement(row),
        _fusion_split_requirement(row, rules),
        _fusion_supplementary_structure_requirement(row, rules),
        _normal_background_requirement(row, track="FUSION", mandatory=True),
        _fusion_frame_orf_requirement(row),
        _novel_sequence_requirement(row, track="FUSION"),
        _peptide_hla_traceability_requirement(row, require_orf=True),
    ]


def _splice_junction_definition_requirement(row: Mapping[str, Any]) -> RequirementAssessment:
    fields = _present(row, "genome_build", "junction_chrom", "junction_start", "junction_end", "junction_strand", "canonical_junction_id", "junction_resolution_status", "junction_coordinate_system")
    text = _text(row, "junction_resolution_status", "junction_coordinate_system", "coordinate_validation_status")
    if any(token in text for token in ("INVALID", "COORDINATE_ERROR", "BUILD_MISMATCH", "FAILED")):
        return _assessment("junction_definition", "Exact splice junction with build/chromosome/intron coordinates/strand", NEGATIVE, "SC_SPLICE_JUNCTION_COORDINATE_INVALID", text, fields, fatal_if_negative=True)
    exact = bool(
        (_raw(row, "canonical_junction_id") or (_raw(row, "junction_chrom") and _raw(row, "junction_start") and _raw(row, "junction_end")))
        and (_raw(row, "junction_strand") or "STRAND_RESOLVED" in text)
    )
    if exact:
        return _assessment("junction_definition", "Exact splice junction with build/chromosome/intron coordinates/strand", SUPPORTED, "SC_SPLICE_JUNCTION_DEFINED", _raw(row, "canonical_junction_id") or text, fields, fatal_if_negative=True)
    if any(token in text for token in _CONFLICT_TOKENS):
        return _assessment("junction_definition", "Exact splice junction with build/chromosome/intron coordinates/strand", CONFLICT, "SC_SPLICE_JUNCTION_CONFLICT", text, fields, fatal_if_conflict=True)
    return _assessment("junction_definition", "Exact splice junction with build/chromosome/intron coordinates/strand", UNASSESSED, "SC_SPLICE_JUNCTION_UNRESOLVED", "Exact junction coordinates/strand unavailable", fields, fatal_if_negative=True)


def _splice_unique_read_requirement(row: Mapping[str, Any], rules: Mapping[str, Any]) -> RequirementAssessment:
    fields = _present(row, "unique_split_reads", "rna_junction_reads", "junction_reads", "junction_support_status", "junction_support_conflict")
    conflict = _text(row, "junction_support_conflict", "junction_support_status")
    if any(token in conflict for token in _CONFLICT_TOKENS):
        return _assessment("unique_split_read_support", "Exact unique split-read support", CONFLICT, "SC_SPLICE_JUNCTION_READ_CONFLICT", conflict, fields, fatal_if_conflict=True)
    reads = _number(row, "unique_split_reads", "rna_junction_reads", "junction_reads")
    minimum = _rule_threshold(rules, "junction_min_reads", 3)
    if reads is None:
        return _assessment("unique_split_read_support", "Exact unique split-read support", UNASSESSED, "SC_SPLICE_JUNCTION_READ_UNASSESSED", "Exact junction read count unavailable", fields, fatal_if_negative=True)
    if reads >= minimum:
        return _assessment("unique_split_read_support", "Exact unique split-read support", SUPPORTED, "SC_SPLICE_JUNCTION_READ_SUPPORTED", f"junction reads={reads:g}", fields, fatal_if_negative=True)
    if reads > 0:
        return _assessment("unique_split_read_support", "Exact unique split-read support", INDETERMINATE_LOW_POWER, "SC_SPLICE_JUNCTION_READ_LOW", f"junction reads={reads:g} below {minimum:g}", fields, fatal_if_negative=True)
    return _assessment("unique_split_read_support", "Exact unique split-read support", NEGATIVE, "SC_SPLICE_JUNCTION_READ_ZERO", "junction reads=0", fields, fatal_if_negative=True)


def _splice_mapping_qc_requirement(row: Mapping[str, Any]) -> RequirementAssessment:
    fields = _present(row, "anchor_size", "overhang", "mapping_quality", "junction_mapq", "repeat_context_status", "pseudogene_risk_status", "junction_support_status")
    text = _text(row, "repeat_context_status", "pseudogene_risk_status", "junction_support_status", "mapping_artifact_status")
    if any(token in text for token in ("MAPPING_ARTIFACT", "PSEUDOGENE_ARTIFACT", "INVALID", "WRONG_JUNCTION")):
        return _assessment("junction_mapping_qc", "Anchor/overhang/MAPQ/repeat/pseudogene QC", NEGATIVE, "SC_SPLICE_MAPPING_ARTIFACT", text, fields, fatal_if_negative=True)
    if any(token in text for token in _CONFLICT_TOKENS):
        return _assessment("junction_mapping_qc", "Anchor/overhang/MAPQ/repeat/pseudogene QC", CONFLICT, "SC_SPLICE_MAPPING_CONFLICT", text, fields, fatal_if_conflict=True)
    anchor = _number(row, "anchor_size", "overhang")
    mapq = _number(row, "mapping_quality", "junction_mapq")
    if any(token in text for token in _LOW_POWER_TOKENS) or (anchor is not None and anchor < 8) or (mapq is not None and mapq < 20):
        return _assessment("junction_mapping_qc", "Anchor/overhang/MAPQ/repeat/pseudogene QC", INDETERMINATE_LOW_POWER, "SC_SPLICE_MAPPING_LOW_POWER", f"anchor={anchor}; MAPQ={mapq}; {text}", fields, fatal_if_negative=True)
    if any(token in text for token in _PASS_TOKENS) or (anchor is not None and anchor >= 8 and (mapq is None or mapq >= 20)):
        return _assessment("junction_mapping_qc", "Anchor/overhang/MAPQ/repeat/pseudogene QC", SUPPORTED, "SC_SPLICE_MAPPING_QC_PASS", f"anchor={anchor}; MAPQ={mapq}; {text}", fields, fatal_if_negative=True)
    return _assessment("junction_mapping_qc", "Anchor/overhang/MAPQ/repeat/pseudogene QC", UNASSESSED, "SC_SPLICE_MAPPING_QC_UNASSESSED", "Anchor/MAPQ/repeat-risk audit unavailable", fields, fatal_if_negative=True)


def _splice_event_structure_requirement(row: Mapping[str, Any]) -> RequirementAssessment:
    fields = _present(row, "splice_event_type", "event_type", "event_name", "reference_path_status", "alternative_path_status", "junction_ids")
    event_type = _text(row, "splice_event_type", "event_type", "event_name")
    known = ("SE", "A3SS", "A5SS", "RI", "MXE", "CRYPTIC", "EXON_SKIP", "INTRON_RETENTION", "ALT_DONOR", "ALT_ACCEPTOR", "SPLICE")
    if any(token in event_type for token in ("INVALID", "UNRESOLVABLE", "WRONG_EVENT")):
        return _assessment("splice_event_structure", "SE/A3SS/A5SS/RI/cryptic-exon event structure", NEGATIVE, "SC_SPLICE_EVENT_STRUCTURE_INVALID", event_type, fields, fatal_if_negative=True)
    if any(token in event_type for token in known):
        return _assessment("splice_event_structure", "SE/A3SS/A5SS/RI/cryptic-exon event structure", SUPPORTED, "SC_SPLICE_EVENT_STRUCTURE_RESOLVED", event_type, fields, fatal_if_negative=True)
    return _assessment("splice_event_structure", "SE/A3SS/A5SS/RI/cryptic-exon event structure", UNASSESSED, "SC_SPLICE_EVENT_STRUCTURE_UNASSESSED", "Splice-event structure unavailable", fields, fatal_if_negative=True)


def _splice_paths_requirement(row: Mapping[str, Any]) -> RequirementAssessment:
    fields = _present(row, "reference_path_status", "alternative_path_status", "junction_match_status", "junction_match_method", "reference_junction_id", "alternative_junction_id")
    text = _text(row, "reference_path_status", "alternative_path_status", "junction_match_status")
    if any(token in text for token in ("WRONG_JUNCTION_FILL", "INVALID", "MISMATCH")):
        return _assessment("reference_alternative_paths", "Reference and alternative splice paths are linked", NEGATIVE, "SC_SPLICE_PATH_LINK_INVALID", text, fields, fatal_if_negative=True)
    if any(token in text for token in _CONFLICT_TOKENS):
        return _assessment("reference_alternative_paths", "Reference and alternative splice paths are linked", CONFLICT, "SC_SPLICE_PATH_CONFLICT", text, fields, fatal_if_conflict=True)
    if any(token in text for token in ("MATCHED", "RESOLVED", "ANNOTATED", "SUPPORTED", "PASS")) or (_raw(row, "reference_junction_id") and _raw(row, "alternative_junction_id")):
        return _assessment("reference_alternative_paths", "Reference and alternative splice paths are linked", SUPPORTED, "SC_SPLICE_REFERENCE_ALT_PATH_RESOLVED", text, fields, fatal_if_negative=True)
    return _assessment("reference_alternative_paths", "Reference and alternative splice paths are linked", UNASSESSED, "SC_SPLICE_REFERENCE_ALT_PATH_UNASSESSED", "Reference/alternative path relation unavailable", fields, fatal_if_negative=True)


def _splice_orf_requirement(row: Mapping[str, Any]) -> RequirementAssessment:
    fields = _present(row, "transcript_hypothesis_id", "orf_id", "orf_evidence_grade", "frame_status", "rna_frame_status", "translation_status", "protein_sequence")
    text = _text(row, "orf_evidence_grade", "frame_status", "rna_frame_status", "translation_status")
    if any(token in text for token in ("INVALID", "NO_ORF", "UNTRANSLATABLE", "WRONG_FRAME", "IMPOSSIBLE")):
        return _assessment("splice_transcript_orf", "Splice transcript hypothesis and translatable ORF", NEGATIVE, "SC_SPLICE_ORF_INVALID", text, fields, fatal_if_negative=True)
    if any(token in text for token in _CONFLICT_TOKENS):
        return _assessment("splice_transcript_orf", "Splice transcript hypothesis and translatable ORF", CONFLICT, "SC_SPLICE_ORF_CONFLICT", text, fields, fatal_if_conflict=True)
    if _raw(row, "transcript_hypothesis_id") and _raw(row, "orf_id") and (any(token in text for token in _PASS_TOKENS) or not text):
        return _assessment("splice_transcript_orf", "Splice transcript hypothesis and translatable ORF", SUPPORTED, "SC_SPLICE_ORF_RECONSTRUCTED", text or "transcript hypothesis + ORF ID", fields, fatal_if_negative=True)
    if any(token in text for token in _PASS_TOKENS):
        return _assessment("splice_transcript_orf", "Splice transcript hypothesis and translatable ORF", SUPPORTED, "SC_SPLICE_ORF_RECONSTRUCTED", text, fields, fatal_if_negative=True)
    return _assessment("splice_transcript_orf", "Splice transcript hypothesis and translatable ORF", UNASSESSED, "SC_SPLICE_ORF_UNASSESSED", "Complete splice transcript/ORF unavailable", fields, fatal_if_negative=True)


def _splice_usage_requirement(row: Mapping[str, Any]) -> RequirementAssessment:
    fields = _present(
        row, "psi", "delta_psi", "junction_usage_status", "tumor_junction_frequency",
        "normal_junction_frequency", "known_normal_isoform_status", "normal_isoform_status",
    )
    text = _text(
        row, "junction_usage_status", "known_normal_isoform_status", "normal_isoform_status",
    )
    if any(token in text for token in ("KNOWN_NORMAL_ISOFORM_CONFIRMED", "BROAD_NORMAL_ISOFORM", "NORMAL_ISOFORM_MATCH")):
        return _assessment(
            "splice_usage", "PSI/junction usage and known-normal-isoform context",
            NEGATIVE, "SC_SPLICE_KNOWN_NORMAL_ISOFORM", text, fields, fatal_if_negative=True,
        )
    if any(token in text for token in _CONFLICT_TOKENS):
        return _assessment(
            "splice_usage", "PSI/junction usage and known-normal-isoform context",
            CONFLICT, "SC_SPLICE_USAGE_CONFLICT", text, fields, fatal_if_conflict=True,
        )
    psi = _number(row, "psi")
    delta_psi = _number(row, "delta_psi")
    tumor_freq = _number(row, "tumor_junction_frequency")
    normal_freq = _number(row, "normal_junction_frequency")
    if any(token in text for token in _LOW_POWER_TOKENS):
        return _assessment(
            "splice_usage", "PSI/junction usage and known-normal-isoform context",
            INDETERMINATE_LOW_POWER, "SC_SPLICE_USAGE_LOW_POWER",
            f"PSI={psi}; delta_PSI={delta_psi}; tumor_frequency={tumor_freq}; normal_frequency={normal_freq}; {text}",
            fields, fatal_if_negative=True,
        )
    if psi is not None or delta_psi is not None or tumor_freq is not None or any(token in text for token in _PASS_TOKENS):
        return _assessment(
            "splice_usage", "PSI/junction usage and known-normal-isoform context",
            SUPPORTED, "SC_SPLICE_USAGE_ASSESSED",
            f"PSI={psi}; delta_PSI={delta_psi}; tumor_frequency={tumor_freq}; normal_frequency={normal_freq}; {text}",
            fields, fatal_if_negative=True,
        )
    return _assessment(
        "splice_usage", "PSI/junction usage and known-normal-isoform context",
        UNASSESSED, "SC_SPLICE_USAGE_UNASSESSED",
        "PSI/junction usage and known-normal-isoform evidence unavailable", fields,
        fatal_if_negative=True,
    )


def _splice_translation_safety_requirement(row: Mapping[str, Any]) -> RequirementAssessment:
    fields = _present(
        row, "translation_direction_status", "protein_translation_direction",
        "nmd_risk_status", "nmd_risk", "nmd_reason", "frame_status", "rna_frame_status",
    )
    text = _text(
        row, "translation_direction_status", "protein_translation_direction",
        "nmd_risk_status", "nmd_risk", "frame_status", "rna_frame_status",
    )
    if any(token in text for token in ("WRONG_DIRECTION", "ANTISENSE_TRANSLATION", "UNTRANSLATABLE", "INVALID")):
        return _assessment(
            "splice_translation_safety", "Translation direction, frame and NMD risk",
            NEGATIVE, "SC_SPLICE_TRANSLATION_INVALID", text, fields, fatal_if_negative=True,
        )
    if any(token in text for token in _CONFLICT_TOKENS):
        return _assessment(
            "splice_translation_safety", "Translation direction, frame and NMD risk",
            CONFLICT, "SC_SPLICE_TRANSLATION_CONFLICT", text, fields, fatal_if_conflict=True,
        )
    if any(token in text for token in ("HIGH_NMD_RISK", "NMD_LIKELY", "PREMATURE_STOP")):
        return _assessment(
            "splice_translation_safety", "Translation direction, frame and NMD risk",
            INDETERMINATE_LOW_POWER, "SC_SPLICE_NMD_CAUTION", text, fields,
            fatal_if_negative=True,
        )
    direction_ok = any(token in text for token in ("SENSE", "CORRECT_DIRECTION", "TRANSLATION_VALID"))
    frame_ok = any(token in text for token in ("IN_FRAME", "IN-FRAME", "INFRAME", "FRAME_CONFIRMED"))
    nmd_assessed = bool(_raw(row, "nmd_risk_status") or _raw(row, "nmd_risk"))
    if direction_ok and frame_ok and nmd_assessed:
        return _assessment(
            "splice_translation_safety", "Translation direction, frame and NMD risk",
            SUPPORTED, "SC_SPLICE_TRANSLATION_NMD_ASSESSED", text, fields,
            fatal_if_negative=True,
        )
    return _assessment(
        "splice_translation_safety", "Translation direction, frame and NMD risk",
        UNASSESSED, "SC_SPLICE_TRANSLATION_NMD_UNASSESSED",
        "Translation direction/frame/NMD assessment incomplete", fields,
        fatal_if_negative=True,
    )


def _splice_requirements(row: Mapping[str, Any], rules: Mapping[str, Any]) -> list[RequirementAssessment]:
    return [
        _splice_junction_definition_requirement(row),
        _splice_unique_read_requirement(row, rules),
        _splice_mapping_qc_requirement(row),
        _splice_event_structure_requirement(row),
        _splice_paths_requirement(row),
        _splice_usage_requirement(row),
        _splice_orf_requirement(row),
        _splice_translation_safety_requirement(row),
        _novel_sequence_requirement(row, track="SPLICE"),
        _peptide_hla_traceability_requirement(row, require_orf=True),
        _normal_background_requirement(row, track="SPLICE", mandatory=True),
    ]


_REQUIREMENT_BUILDERS: dict[str, Callable[[Mapping[str, Any], Mapping[str, Any]], list[RequirementAssessment]]] = {
    "SNV": _snv_requirements,
    "INDEL": _indel_requirements,
    "FUSION": _fusion_requirements,
    "SPLICE": _splice_requirements,
}


# ---------------------------------------------------------------------------
# Public classification API
# ---------------------------------------------------------------------------


def derive_source_chain_confidence(
    row: Mapping[str, Any],
    rules: Mapping[str, Any] | None = None,
) -> SourceChainResult:
    effective_rules = rules or {}
    track = source_chain_track(row)
    if track not in _REQUIREMENT_BUILDERS:
        req = _assessment(
            "source_chain_track",
            "Source-chain event class",
            NOT_APPLICABLE,
            "SC_TRACK_NOT_APPLICABLE",
            f"Unsupported source-chain track={track}",
            (),
            core=False,
        )
        return SourceChainResult(
            track=track,
            tier="C3",
            label=SOURCE_CHAIN_TIERS["C3"]["label"],
            grade=SOURCE_CHAIN_TIERS["C3"]["grade"],
            rule_version=str(_section(effective_rules, "source_chain").get("rule_version", SOURCE_CHAIN_RULE_VERSION)),
            orthogonal_status=NOT_APPLICABLE,
            orthogonal_sources=(),
            requirements=(req,),
            reason_codes=(req.reason_code,),
            hard_failure=False,
            hard_failure_codes=(),
            missing_requirements=(),
            low_power_requirements=(),
            negative_requirements=(),
            conflict_requirements=(),
            supported_requirements=(),
            not_applicable_requirements=(req.name,),
        )

    requirements = _REQUIREMENT_BUILDERS[track](row, effective_rules)
    rna_requirement = next((req for req in requirements if req.name == "rna_or_direct_evidence"), _assessment("rna_or_direct_evidence", "RNA/direct evidence", NOT_APPLICABLE, "SC_RNA_NOT_APPLICABLE", "Not required for this track", (), core=False))
    orthogonal_status, orthogonal_sources, orthogonal_code = _orthogonal_confirmation(
        row, track=track, rna_requirement=rna_requirement
    )
    orthogonal_req = _assessment(
        "orthogonal_confirmation",
        "Independent/cross-modal orthogonal confirmation",
        orthogonal_status,
        orthogonal_code,
        ",".join(orthogonal_sources) or "No qualifying independent/cross-modal confirmation",
        _present(row, "orthogonal_confirmation_status", "targeted_validation_status", "rt_pcr_status", "sanger_status", "long_read_status", "independent_library_status", "dna_sv_confirmation_status", "protein_evidence_status", "cross_platform_status"),
        core=False,
        fatal_if_negative=True,
    )
    requirements.append(orthogonal_req)

    # Explicit global invalidity and sample/coordinate failures are always C4.
    global_text = _text(
        row,
        "sample_identity_status", "sample_mix_status", "coordinate_validation_status",
        "reference_alt_orientation_status", "genome_build_status", "source_chain_override_status",
    )
    global_hard_codes: list[str] = []
    if any(token in global_text for token in ("SAMPLE_MIX", "SAMPLE_MISMATCH", "WRONG_SAMPLE")):
        global_hard_codes.append("SC_SAMPLE_IDENTITY_INVALID")
    if any(token in global_text for token in ("COORDINATE_ERROR", "BUILD_MISMATCH", "REF_ALT_REVERSED", "WRONG_REF_ALT")):
        global_hard_codes.append("SC_COORDINATE_OR_ALLELE_INVALID")
    if "FORCE_C4" in global_text:
        global_hard_codes.append("SC_FORCE_C4_OVERRIDE")

    hard_codes = list(global_hard_codes)
    for req in requirements:
        if req.status == NEGATIVE and req.fatal_if_negative:
            hard_codes.append(req.reason_code)
        if req.status == CONFLICT and req.fatal_if_conflict:
            hard_codes.append(req.reason_code)
    hard_codes = list(dict.fromkeys(hard_codes))

    core_requirements = [req for req in requirements if req.core and req.status != NOT_APPLICABLE]
    all_core_supported = bool(core_requirements) and all(req.status == SUPPORTED for req in core_requirements)
    if hard_codes:
        tier = "C4"
    elif all_core_supported and orthogonal_status == SUPPORTED:
        tier = "C1"
    elif all_core_supported:
        tier = "C2"
    else:
        tier = "C3"

    supported = tuple(req.name for req in requirements if req.status == SUPPORTED)
    missing = tuple(req.name for req in requirements if req.status == UNASSESSED)
    low_power = tuple(req.name for req in requirements if req.status == INDETERMINATE_LOW_POWER)
    negative = tuple(req.name for req in requirements if req.status == NEGATIVE)
    conflicts = tuple(req.name for req in requirements if req.status == CONFLICT)
    not_applicable = tuple(req.name for req in requirements if req.status == NOT_APPLICABLE)
    reason_codes = tuple(dict.fromkeys(req.reason_code for req in requirements if req.status != SUPPORTED or req.name == "orthogonal_confirmation"))

    return SourceChainResult(
        track=track,
        tier=tier,
        label=SOURCE_CHAIN_TIERS[tier]["label"],
        grade=SOURCE_CHAIN_TIERS[tier]["grade"],
        rule_version=str(_section(effective_rules, "source_chain").get("rule_version", SOURCE_CHAIN_RULE_VERSION)),
        orthogonal_status=orthogonal_status,
        orthogonal_sources=orthogonal_sources,
        requirements=tuple(requirements),
        reason_codes=reason_codes,
        hard_failure=bool(hard_codes),
        hard_failure_codes=tuple(hard_codes),
        missing_requirements=missing,
        low_power_requirements=low_power,
        negative_requirements=negative,
        conflict_requirements=conflicts,
        supported_requirements=supported,
        not_applicable_requirements=not_applicable,
    )


def source_chain_consensus_fields(
    row: Mapping[str, Any], rules: Mapping[str, Any] | None = None
) -> dict[str, str]:
    return derive_source_chain_confidence(row, rules).as_row()


def build_source_chain_table(
    input_tsv: str | Path,
    output_tsv: str | Path,
    requirements_tsv: str | Path | None = None,
    rules: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    rows = read_tsv(input_tsv)
    output: list[dict[str, str]] = []
    long_rows: list[dict[str, str]] = []
    tier_counts: dict[str, int] = {tier: 0 for tier in SOURCE_CHAIN_TIERS}
    track_counts: dict[str, int] = {}
    for row in rows:
        result = derive_source_chain_confidence(row, rules)
        merged = dict(row)
        merged.update(result.as_row())
        output.append(merged)
        tier_counts[result.tier] = tier_counts.get(result.tier, 0) + 1
        track_counts[result.track] = track_counts.get(result.track, 0) + 1
        for req in result.requirements:
            long_rows.append({
                "sample_id": _raw(row, "sample_id"),
                "event_id": _raw(row, "event_id"),
                "peptide_id": _raw(row, "peptide_id"),
                "gene": _raw(row, "gene"),
                "source_chain_track": result.track,
                "source_chain_confidence_tier": result.tier,
                "requirement_name": req.name,
                "requirement_label": req.label,
                "requirement_applicability": NOT_APPLICABLE if req.status == NOT_APPLICABLE else "APPLICABLE",
                "requirement_status": req.status,
                "requirement_value": req.reason,
                "requirement_core": "yes" if req.core else "no",
                "fatal_if_negative": "yes" if req.fatal_if_negative else "no",
                "fatal_if_conflict": "yes" if req.fatal_if_conflict else "no",
                "reason_code": req.reason_code,
                "reason": req.reason,
                "source_fields": ",".join(req.source_fields),
                "requirement_conflict": "yes" if req.status == CONFLICT else "no",
                "rule_version": result.rule_version,
            })
    fields = list(rows[0]) if rows else []
    source_fields = list(output[0]) if output else []
    fields.extend(field for field in source_fields if field not in fields)
    write_tsv(output_tsv, output, fields)
    if requirements_tsv:
        write_tsv(requirements_tsv, long_rows)
    return {
        "rows": len(output),
        "requirements": len(long_rows),
        "output": str(output_tsv),
        "requirements_output": str(requirements_tsv or ""),
        "tier_counts": tier_counts,
        "track_counts": track_counts,
        "rule_version": str(_section(rules or {}, "source_chain").get("rule_version", SOURCE_CHAIN_RULE_VERSION)),
    }


__all__ = [
    "SOURCE_CHAIN_RULE_VERSION",
    "SOURCE_CHAIN_TIERS",
    "SUPPORTED",
    "NOT_APPLICABLE",
    "UNASSESSED",
    "INDETERMINATE_LOW_POWER",
    "NEGATIVE",
    "CONFLICT",
    "RequirementAssessment",
    "SourceChainResult",
    "source_chain_track",
    "derive_source_chain_confidence",
    "source_chain_consensus_fields",
    "build_source_chain_table",
]
