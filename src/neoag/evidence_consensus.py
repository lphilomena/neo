from __future__ import annotations

from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import tomllib
from typing import Any, Mapping

from .evidence_states import DERIVERS, derive_all_states, event_track
from .pareto import nondominated_fronts
from .utils import read_tsv, write_tsv


DEFAULT_RULES: dict[str, Any] = {
    "metadata": {
        "name": "sarcoma_evidence_consensus_v1",
        "version": "1.0",
        "status": "PROVISIONAL_RESEARCH_ONLY",
        "mode": "parallel",
        "replace_weighted_ranking": False,
        "unassessed_is_negative": False,
    },
    "presentation": {
        "netmhcpan_el_strong": 0.5,
        "netmhcpan_el_supported": 2.0,
        "mhcflurry_presentation_strong": 0.70,
        "mhcflurry_presentation_supported": 0.50,
        "netmhcstabpan_supported": 1.4,
        "prime_rank_strong": 2.0,
        "prime_rank_supported": 5.0,
        "bigmhc_im_strong": 0.70,
        "bigmhc_im_supported": 0.50,
        "deepimmuno_supported": 0.50,
    },
    "rna": {
        "snv_min_alt_reads": 3,
        "snv_min_vaf": 0.02,
        "junction_min_reads": 3,
        "junction_strong_reads": 10,
    },
    "completeness": {
        "allow_r1_with_hla_loh_unassessed": False,
        "allow_r1_with_safety_partial": False,
        "allow_r1_with_ccf_low_confidence": False,
    },
    "manual_review": {"genes": ["KRAS", "TP53"], "events": ["EWSR1::WT1"]},
    "output": {"peptide_id_tie_break": True, "event_deduplicate": True},
}

BASE_DIMENSIONS = [
    "event_authenticity_grade",
    "rna_support_grade",
    "presentation_consensus_grade",
    "mutant_specificity_grade",
    "safety_grade",
    "hla_appm_grade",
    "clonality_grade",
    "evidence_completeness_grade",
]
PARETO_DERIVED_FIELDS = [
    "novel_tail_evidence_grade",
    "junction_authenticity_grade",
    "junction_reads_grade",
    "frame_evidence_grade",
    "normal_junction_safety_grade",
]
PARETO_DIMENSIONS_BY_TRACK = {
    "MISSENSE": [
        "event_authenticity_grade", "rna_support_grade", "presentation_consensus_grade",
        "mutant_specificity_grade", "safety_grade", "hla_appm_grade", "clonality_grade",
        "evidence_completeness_grade",
    ],
    "FRAMESHIFT": [
        "event_authenticity_grade", "rna_support_grade", "novel_tail_evidence_grade",
        "presentation_consensus_grade", "safety_grade", "hla_appm_grade", "clonality_grade",
    ],
    "FUSION": [
        "junction_authenticity_grade", "junction_reads_grade", "frame_evidence_grade",
        "presentation_consensus_grade", "normal_junction_safety_grade", "hla_appm_grade",
        "evidence_completeness_grade",
    ],
    "SPLICE": [
        "junction_authenticity_grade", "junction_reads_grade", "frame_evidence_grade",
        "presentation_consensus_grade", "normal_junction_safety_grade", "hla_appm_grade",
        "evidence_completeness_grade",
    ],
    "DNA_SV": [
        "event_authenticity_grade", "rna_support_grade", "novel_tail_evidence_grade",
        "presentation_consensus_grade", "safety_grade", "hla_appm_grade", "clonality_grade",
        "evidence_completeness_grade",
    ],
    "MANUAL_REVIEW": BASE_DIMENSIONS,
    "OTHER": BASE_DIMENSIONS,
}
DIMENSIONS = BASE_DIMENSIONS + PARETO_DERIVED_FIELDS
STATE_NAMES = (
    "event_authenticity", "rna_support", "presentation_consensus",
    "mutant_specificity", "clonality", "hla_appm", "safety",
    "evidence_completeness",
)
STATE_OUTPUT_FIELDS = tuple(f"{name}_state" for name in STATE_NAMES)
GRADE_ORDER = {"R1": 1, "R2": 2, "R3": 3, "R4": 4}
CAP_TO_GRADE = {
    "A": "R1", "NONE": "R1", "B": "R2", "B_CAUTION": "R2",
    "C": "R3", "C_CAUTION": "R3", "D": "R4",
}
CAP_FIELDS = (
    "priority_cap", "safety_priority_cap", "cross_platform_priority_cap",
    "mutant_specificity_priority_cap",
)
CONSENSUS_FIELDS = (
    "legacy_weighted_rank", "biological_event_track", "evidence_track", "pareto_dimensions",
    "hard_failure", "hard_failure_reasons",
    "legacy_priority_cap", "consensus_priority_cap", "evidence_grade_cap", "evidence_grade_cap_reasons",
    "manual_review_required", "manual_review_reason",
    *STATE_OUTPUT_FIELDS, *DIMENSIONS, "safety_completeness_grade",
    "ccf_confidence_state", "ccf_confidence_grade",
    "netmhcpan_tiebreak_rank", "mhcflurry_tiebreak_score",
    "evidence_grade_uncapped", "evidence_grade", "pareto_front", "track_rank",
    "evidence_rank", "evidence_rank_key", "evidence_consensus_score", "evidence_completeness_score",
    "evidence_assessed_layers", "evidence_missing_layers", "evidence_conflict_layers",
    "evidence_layer_states", "consensus_action", "recommended_next_steps", "consensus_trace",
)
CONFLICT_FIELDS = (
    "peptide_id", "event_id", "gene", "evidence_track", "layer", "state",
    "field", "selected_source", "selected_value", "other_source", "other_value",
    "precedence_version", "conflict_type", "reason", "recommended_action",
)


def load_consensus_rules(path: str | Path | None = None) -> dict[str, Any]:
    if path is None:
        return json.loads(json.dumps(DEFAULT_RULES))
    with Path(path).open("rb") as handle:
        loaded = tomllib.load(handle)
    merged = json.loads(json.dumps(DEFAULT_RULES))
    for section, values in loaded.items():
        if isinstance(values, Mapping) and isinstance(merged.get(section), Mapping):
            merged[section].update(values)
        else:
            merged[section] = values
    return merged


def _number(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def _strictest_cap(row: Mapping[str, Any]) -> str:
    caps = [str(row.get(field, "")).strip().upper() for field in CAP_FIELDS]
    caps = [cap for cap in caps if cap in CAP_TO_GRADE]
    return max(caps, key=lambda cap: GRADE_ORDER[CAP_TO_GRADE[cap]], default="")


def _manual_review(row: Mapping[str, Any], rules: Mapping[str, Any]) -> tuple[bool, str]:
    cfg = rules.get("manual_review", {})
    gene = str(row.get("gene", "")).strip().upper()
    event_text = " ".join(str(row.get(field, "")) for field in ("gene", "event_id", "fusion_name")).upper()
    genes = {str(value).upper() for value in cfg.get("genes", [])}
    events = {str(value).upper() for value in cfg.get("events", [])}
    reasons = []
    if gene in genes:
        reasons.append(f"manual-review gene={gene}")
    for event in events:
        if event and event in event_text:
            reasons.append(f"manual-review event={event}")
    return bool(reasons), ";".join(reasons)


def _first_number(row: Mapping[str, Any], *fields: str) -> float | None:
    for field in fields:
        value = row.get(field)
        if value in (None, "", "NA", "N/A", "."):
            continue
        try:
            return float(value)
        except (TypeError, ValueError):
            continue
    return None


def _pareto_derived_grades(
    source: Mapping[str, Any],
    states: Mapping[str, Mapping[str, Any]],
    biological_track: str,
    rules: Mapping[str, Any],
) -> dict[str, int]:
    text = " ".join(str(source.get(field, "")).upper() for field in (
        "event_type", "peptide_consequence", "consequence", "frame_status", "rna_frame_status",
        "contains_novel_aa", "novel_tail_status",
    ))
    novel_tail = 3 if _is_clear_novel_sequence(source) or any(token in text for token in (
        "FRAMESHIFT", "FRAME_SHIFT", "NOVEL_TAIL", "NOVEL AA", "CONTAINS_NOVEL_AA TRUE",
    )) else 0

    junction_reads = _first_number(source, "rna_junction_reads", "junction_reads", "split_reads")
    rna_cfg = rules.get("rna", {})
    minimum = float(rna_cfg.get("junction_min_reads", 3))
    strong = float(rna_cfg.get("junction_strong_reads", 10))
    if junction_reads is None:
        junction_grade = 0
    elif junction_reads >= strong:
        junction_grade = 3
    elif junction_reads >= minimum:
        junction_grade = 2
    elif junction_reads > 0:
        junction_grade = 1
    else:
        junction_grade = 0

    frame_text = " ".join(str(source.get(field, "")).upper() for field in (
        "frame_status", "rna_frame_status", "in_frame", "reading_frame",
    ))
    if any(token in frame_text for token in ("IN_FRAME", "IN-FRAME", "INFRAME", "FRAME_CONFIRMED")):
        frame_grade = 3
    elif any(token in frame_text for token in ("OUT_OF_FRAME", "OUT-OF-FRAME", "FRAMESHIFT")):
        frame_grade = 2 if novel_tail else 1
    elif frame_text.strip():
        frame_grade = 1
    else:
        frame_grade = 0

    normal_text = " ".join(str(source.get(field, "")).upper() for field in (
        "normal_junction_status", "normal_junction_assessment_status", "normal_junction_exact_match",
    ))
    if any(token in normal_text for token in ("DETECTED", "EXACT_MATCH", "SUPPORTED_IN_NORMAL")) and "NOT_DETECTED" not in normal_text:
        normal_junction_grade = 0
    elif any(token in normal_text for token in ("NOT_DETECTED", "ABSENT", "NO_MATCH")):
        normal_junction_grade = 3
    elif normal_text.strip() and not any(token in normal_text for token in ("UNASSESSED", "MISSING", "NOT_AVAILABLE")):
        normal_junction_grade = 2
    else:
        normal_junction_grade = 1

    junction_authenticity = int(states["event_authenticity"]["grade"])
    if biological_track not in {"FUSION", "SPLICE"}:
        junction_authenticity = 0
        junction_grade = 0
        frame_grade = 0
        normal_junction_grade = 0
    return {
        "novel_tail_evidence_grade": novel_tail,
        "junction_authenticity_grade": junction_authenticity,
        "junction_reads_grade": junction_grade,
        "frame_evidence_grade": frame_grade,
        "normal_junction_safety_grade": normal_junction_grade,
    }


def _safety_completeness_grade(
    source: Mapping[str, Any], safety_state: Mapping[str, Any],
) -> int:
    text = " ".join(str(source.get(field, "")).strip().upper() for field in (
        "safety_evidence_completeness", "event_safety_evidence_completeness",
        "reference_proteome_status", "normal_ligandome_status", "normal_junction_assessment_status",
    ))
    if any(token in text for token in ("UNASSESSED", "MISSING", "NOT_AVAILABLE", "INCOMPLETE")):
        return 0
    if any(token in text for token in ("COMPLETE", "FULL", "HIGH")):
        return 3
    if any(token in text for token in ("PARTIAL", "MEDIUM")):
        return 2
    if "LOW" in text:
        return 1
    state = str(safety_state.get("state", ""))
    return 3 if state == "SAFETY_PASS" else 1 if safety_state.get("assessed") else 0


def _ccf_confidence(
    source: Mapping[str, Any], clonality_state: Mapping[str, Any],
) -> tuple[str, int]:
    text = " ".join(str(source.get(field, "")).strip().upper() for field in (
        "ccf_confidence", "ccf_resolution", "ccf_status", "clonality_status",
    ))
    if any(token in text for token in ("UNASSESSED", "UNRESOLVED", "MISSING", "NOT_AVAILABLE")):
        return "CCF_UNASSESSED", 0
    if "LOW" in text:
        return "CCF_LOW_CONFIDENCE", 1
    if any(token in text for token in ("MEDIUM", "MODERATE", "PARTIAL")):
        return "CCF_MEDIUM_CONFIDENCE", 2
    if any(token in text for token in ("HIGH", "CONFIDENT", "CLONAL_LIKE")):
        return "CCF_HIGH_CONFIDENCE", 3
    if not clonality_state.get("assessed"):
        return "CCF_UNASSESSED", 0
    return "CCF_CONFIDENCE_UNSPECIFIED", 1


def _format_tiebreak_number(value: float | None) -> str:
    return "NA" if value is None else f"{value:.8g}"


def _is_clear_novel_sequence(source: Mapping[str, Any]) -> bool:
    novel = str(source.get("contains_novel_aa", "")).strip().lower() in {"true", "yes", "1"}
    junction = str(source.get("crosses_junction", "")).strip().lower() in {"true", "yes", "1"}
    consequence = " ".join(str(source.get(field, "")).upper() for field in ("peptide_consequence", "event_type"))
    return novel or junction or any(token in consequence for token in ("NOVEL_TAIL", "NOVEL JUNCTION", "FRAMESHIFT_NOVEL"))


def _uncapped_grade(
    source: Mapping[str, Any],
    states: Mapping[str, Mapping[str, Any]],
    track: str,
) -> str:
    event = states["event_authenticity"]["state"]
    rna = states["rna_support"]["state"]
    presentation = states["presentation_consensus"]["state"]
    specificity = states["mutant_specificity"]["state"]
    hla = states["hla_appm"]["state"]
    safety = states["safety"]["state"]
    clonality = states["clonality"]
    completeness = states["evidence_completeness"]["grade"]
    cross = " ".join(str(source.get(field, "")).upper() for field in ("cross_platform_status", "comparison_status"))

    if (
        any(state["hard_fail"] for state in states.values())
        or event == "EVENT_ARTIFACT_RISK"
        or presentation == "PRESENTATION_WEAK"
        or specificity in {"WT_BETTER", "NON_MUTANT_SEQUENCE"}
        or hla in {"RESTRICTING_HLA_LOST", "MAJOR_APPM_DEFECT"}
        or safety in {"SAFETY_HIGH_RISK", "SAFETY_REJECT"}
        or (event == "EVENT_CONFLICT" and "SOURCE_PASS_NOT_REPRODUCED_BY_PILEUP" not in cross)
    ):
        return "R4"

    mutant_specific = specificity == "MT_SPECIFIC" or _is_clear_novel_sequence(source)
    ccf_supported_or_not_required = clonality["grade"] >= 2 or track in {"FUSION", "SPLICE"}
    safety_complete = safety == "SAFETY_PASS" and not str(source.get("safety_missing_layers", "")).strip()
    r1 = (
        event in {"EVENT_STRONG", "EVENT_CONFIRMED"}
        and rna == "RNA_CONFIRMED"
        and presentation == "PRESENTATION_CONSISTENT_STRONG"
        and mutant_specific
        and hla == "HLA_APPM_RETAINED"
        and safety_complete
        and ccf_supported_or_not_required
        and completeness >= 3
    )
    if r1:
        return "R1"

    core_established = (
        event in {"EVENT_STRONG", "EVENT_CONFIRMED"}
        and rna == "RNA_CONFIRMED"
        and presentation in {"PRESENTATION_CONSISTENT_STRONG", "PRESENTATION_MODERATE"}
        and specificity not in {"WT_BETTER", "NON_MUTANT_SEQUENCE", "UNASSESSED"}
        and hla not in {"RESTRICTING_HLA_LOST", "MAJOR_APPM_DEFECT", "HLA_LOH_UNASSESSED"}
        and safety not in {"SAFETY_HIGH_RISK", "SAFETY_REJECT"}
    )
    cautions = 0
    cautions += clonality["grade"] < 2 and track not in {"FUSION", "SPLICE"}
    cautions += hla == "HLA_APPM_CAUTION" or "LOW" in str(source.get("appm_evidence_completeness", "")).upper()
    cautions += safety in {"SAFETY_PARTIAL", "SAFETY_REVIEW"}
    cautions += specificity == "MARGINAL_MT_ADVANTAGE"
    cautions += completeness < 3
    if core_established and cautions <= 1:
        return "R2"
    return "R3"


def _next_steps(
    grade: str,
    source: Mapping[str, Any],
    states: Mapping[str, Mapping[str, Any]],
) -> tuple[str, str]:
    if grade == "R1":
        return "FIRST_BATCH_EXPERIMENTAL_PRIORITY", "confirm assay design; retain matched WT control; proceed to first-batch experimental validation"
    if grade == "R2":
        return "ADVANCE_WITH_CAUTION", "resolve the single caution where feasible; retain matched WT control; then proceed to validation"
    if grade == "R4":
        return "DO_NOT_ADVANCE", "retain only for audit or mechanism-driven manual review; do not advance automatically"
    steps: list[str] = []
    if states["rna_support"]["state"] != "RNA_CONFIRMED":
        steps.extend(["targeted RNA", "IGV", "RT-PCR/Sanger"])
    if event_track(source) == "FUSION" and not states["clonality"]["assessed"]:
        steps.extend(["RT-PCR/Sanger", "second fusion caller", "orthogonal breakpoint review"])
    if states["presentation_consensus"]["state"] in {"PRESENTATION_SINGLE_TOOL", "PRESENTATION_DISCORDANT", "PRESENTATION_UNASSESSED"}:
        steps.append("second presentation tool/group")
    normal_junction = " ".join(str(source.get(field, "")).upper() for field in (
        "normal_junction_assessment_status", "event_normal_junction_assessment_status",
    ))
    if not normal_junction or "UNASSESSED" in normal_junction:
        steps.append("normal tissue junction check")
    phasing = " ".join(str(source.get(field, "")).upper() for field in ("haplotype_status", "phase_confidence"))
    if any(token in phasing for token in ("REQUIRED", "UNPHASED", "LOW_CONFIDENCE")):
        steps.extend(["IGV phasing review", "read-backed phasing"])
    if states["mutant_specificity"]["state"] != "MT_SPECIFIC" and not _is_clear_novel_sequence(source):
        steps.append("complete MT/WT paired prediction")
    if "LOW" in str(source.get("appm_evidence_completeness", "")).upper() or states["hla_appm"]["state"] == "HLA_LOH_UNASSESSED":
        steps.append("complete APPM/HLA LOH evidence")
    if "LOW" in " ".join(str(source.get(field, "")).upper() for field in ("ccf_confidence", "ccf_resolution")):
        steps.append("review purity/CNV/CCF evidence")
    if states["safety"]["state"] == "SAFETY_PARTIAL":
        steps.append("complete normal proteome/ligandome/junction safety references")
    if states["event_authenticity"]["state"] not in {"EVENT_STRONG", "EVENT_CONFIRMED"}:
        steps.append("second caller or targeted event validation")
    unique = list(dict.fromkeys(steps))
    return "EVIDENCE_COMPLETION_FIRST", "; ".join(unique or ["complete missing evidence before experimental progression"])


def _constrained_grade(uncapped: str, cap: str, hard_failure: bool) -> str:
    minimum = 4 if hard_failure else GRADE_ORDER.get(CAP_TO_GRADE.get(cap, "R1"), 1)
    return f"R{max(GRADE_ORDER[uncapped], minimum)}"


def _derived_grade_caps(
    source: Mapping[str, Any],
    states: Mapping[str, Mapping[str, Any]],
    track: str,
) -> list[tuple[str, str]]:
    caps: list[tuple[str, str]] = []
    mutation_source = str(source.get("mutation_source", "")).upper()
    ccf_resolution = str(source.get("ccf_resolution", "")).upper()
    if track == "FUSION" and ("RNA_ONLY" in mutation_source or "RNA_ONLY" in ccf_resolution or not states["clonality"]["assessed"]):
        caps.append(("R3", "CAP_RNA_ONLY_FUSION"))
    normal_junction = " ".join(str(source.get(field, "")).upper() for field in (
        "normal_junction_assessment_status", "event_normal_junction_assessment_status",
    ))
    if not normal_junction or any(token in normal_junction for token in ("UNASSESSED", "MISSING", "NOT_AVAILABLE")):
        if track in {"FUSION", "SPLICE"}:
            caps.append(("R3", "CAP_NORMAL_JUNCTION_UNASSESSED"))
    if states["safety"]["state"] == "SAFETY_PARTIAL":
        caps.append(("R3", "CAP_SAFETY_PARTIAL"))
    phasing = " ".join(str(source.get(field, "")).upper() for field in ("haplotype_status", "phase_confidence"))
    if any(token in phasing for token in ("REQUIRED", "UNPHASED", "LOW_CONFIDENCE")):
        caps.append(("R3", "CAP_PHASING_REQUIRED"))
    capture = " ".join(str(source.get(field, "")).upper() for field in ("evidence_scope", "capture_limited", "wes_confidence_tier"))
    if track == "DNA_SV" and ("CAPTURE" in capture or "TRUE" in capture or "WES" in capture):
        caps.append(("R3", "CAP_WES_CAPTURE_LIMITED_SV"))
    specificity = states["mutant_specificity"]["state"]
    if specificity == "MT_WT_SIMILAR":
        caps.append(("R3", "CAP_MT_WT_SIMILAR"))
    elif specificity == "MARGINAL_MT_ADVANTAGE":
        caps.append(("R2", "CAP_MARGINAL_MT_ADVANTAGE"))
    elif specificity == "WT_BETTER":
        caps.append(("R4", "CAP_WT_BETTER"))
    elif specificity == "UNASSESSED" and not _is_clear_novel_sequence(source):
        caps.append(("R2", "CAP_MUTANT_SPECIFICITY_UNASSESSED"))
    hla_state = states["hla_appm"]["state"]
    if hla_state == "HLA_LOH_UNASSESSED":
        caps.append(("R2", "CAP_HLA_LOH_UNASSESSED"))
    elif hla_state == "MAJOR_APPM_DEFECT":
        caps.append(("R4", "CAP_MAJOR_APPM_DEFECT"))
    appm_completeness = str(source.get("appm_evidence_completeness", "")).upper()
    if "LOW" in appm_completeness:
        caps.append(("R2", "CAP_APPM_EVIDENCE_LOW"))
    ccf_confidence = " ".join(str(source.get(field, "")).upper() for field in ("ccf_confidence", "ccf_resolution"))
    if "LOW" in ccf_confidence:
        caps.append(("R2", "CAP_CCF_LOW_CONFIDENCE"))
    presentation = states["presentation_consensus"]["state"]
    if presentation == "PRESENTATION_SINGLE_TOOL":
        caps.append(("R3", "CAP_SINGLE_PRESENTATION_TOOL"))
    elif presentation == "PRESENTATION_DISCORDANT":
        caps.append(("R3", "CAP_PRESENTATION_DISCORDANT"))
    cross = " ".join(str(source.get(field, "")).upper() for field in ("cross_platform_status", "comparison_status"))
    if "SOURCE_PASS_NOT_REPRODUCED_BY_PILEUP" in cross:
        caps.append(("R3", "CAP_SOURCE_UNREPRODUCED"))
    if states["event_authenticity"]["state"] == "EVENT_SAMPLE_SPECIFIC":
        caps.append(("R3", "CAP_EVENT_SAMPLE_SPECIFIC"))
    if states["rna_support"]["state"] in {"RNA_NEGATIVE", "GENE_EXPRESSION_ONLY"}:
        caps.append(("R3", f"CAP_{states['rna_support']['state']}"))
    return caps


def _normalized_row(
    source: Mapping[str, Any],
    rules: Mapping[str, Any],
    legacy_rank: int,
) -> tuple[dict[str, str], list[dict[str, str]]]:
    states = derive_all_states(source, rules)
    biological_track = event_track(source)
    review, review_reason = _manual_review(source, rules)
    pareto_track = "MANUAL_REVIEW" if review else biological_track
    row = {key: str(value or "") for key, value in source.items()}
    row["legacy_weighted_rank"] = str(legacy_rank)
    row["biological_event_track"] = biological_track
    row["evidence_track"] = pareto_track
    row["pareto_dimensions"] = ",".join(PARETO_DIMENSIONS_BY_TRACK.get(pareto_track, BASE_DIMENSIONS))
    assessed = [name for name, state in states.items() if state["assessed"]]
    missing = [name for name, state in states.items() if not state["assessed"]]
    conflicts = [name for name, state in states.items() if state["conflict"]]
    source_conflict_fields = [value for value in str(source.get("evidence_conflict_fields", "")).split(",") if value]
    if source_conflict_fields:
        conflicts.append("source_precedence")
    hard = [
        f"{state.get('hard_code') or 'HARD_UNSPECIFIED'}:{name}:{state['reason']}"
        for name, state in states.items() if state["hard_fail"]
    ]
    for name, state in states.items():
        row[f"{name}_state"] = str(state["state"])
        row[f"{name}_grade"] = str(state["grade"])
    row.update({key: str(value) for key, value in _pareto_derived_grades(
        source, states, biological_track, rules,
    ).items()})
    netmhcpan_rank = _first_number(
        source, "netmhcpan_mt_rank_el", "netmhcpan_el_rank", "el_rank", "binding_rank",
    )
    mhcflurry_score = _first_number(source, "mhcflurry_presentation_score")
    row["safety_completeness_grade"] = str(_safety_completeness_grade(source, states["safety"]))
    ccf_confidence_state, ccf_confidence_grade = _ccf_confidence(source, states["clonality"])
    row["ccf_confidence_state"] = ccf_confidence_state
    row["ccf_confidence_grade"] = str(ccf_confidence_grade)
    row["netmhcpan_tiebreak_rank"] = _format_tiebreak_number(netmhcpan_rank)
    row["mhcflurry_tiebreak_score"] = _format_tiebreak_number(mhcflurry_score)
    cap = _strictest_cap(source)
    uncapped = _uncapped_grade(source, states, biological_track)
    derived_caps = _derived_grade_caps(source, states, biological_track)
    if cap:
        derived_caps.append((CAP_TO_GRADE[cap], f"CAP_EXISTING_PRIORITY_{cap}"))
    grade_cap = max((value for value, _ in derived_caps), key=lambda value: GRADE_ORDER[value], default="R1")
    grade = _constrained_grade(uncapped, next((name for name, value in CAP_TO_GRADE.items() if value == grade_cap), ""), bool(hard))
    action, next_steps = _next_steps(grade, source, states)
    if source_conflict_fields:
        review = True
        source_reason = "source conflicts=" + ",".join(source_conflict_fields)
        review_reason = ";".join(value for value in (review_reason, source_reason) if value)
    row.update({
        "hard_failure": "yes" if hard else "no",
        "hard_failure_reasons": ";".join(hard),
        "legacy_priority_cap": cap,
        "consensus_priority_cap": grade_cap,
        "evidence_grade_cap": grade_cap,
        "evidence_grade_cap_reasons": ",".join(reason for _, reason in derived_caps),
        "manual_review_required": "yes" if review else "no",
        "manual_review_reason": review_reason,
        "evidence_grade_uncapped": uncapped,
        "evidence_grade": grade,
        "evidence_consensus_score": f"{sum(state['grade'] for state in states.values()) / (3 * len(states)):.4f}",
        "evidence_completeness_score": f"{len(assessed) / len(states):.4f}",
        "evidence_assessed_layers": ",".join(assessed),
        "evidence_missing_layers": ",".join(missing),
        "evidence_conflict_layers": ",".join(conflicts),
        "evidence_layer_states": ";".join(f"{name}:{state['state']}" for name, state in states.items()),
        "consensus_action": action,
        "recommended_next_steps": next_steps,
        "pareto_front": "",
        "track_rank": "",
        "evidence_rank": "",
        "evidence_rank_key": "",
    })
    trace = []
    if hard:
        trace.append("hard_failure")
    if cap:
        trace.append(f"priority_cap={cap}")
    if derived_caps:
        trace.append("grade_caps=" + ",".join(reason for _, reason in derived_caps))
    if missing:
        trace.append("missing=" + ",".join(missing))
    if conflicts:
        trace.append("conflict=" + ",".join(conflicts))
    if review:
        trace.append("manual_review")
    if grade != uncapped:
        trace.append(f"grade_capped={uncapped}->{grade}")
    row["consensus_trace"] = ";".join(trace or ["evidence_grade_only"])

    state_row = {
        "peptide_id": str(source.get("peptide_id", "")),
        "event_id": str(source.get("event_id", "")),
        "sample_id": str(source.get("sample_id", "")),
        "gene": str(source.get("gene", "")),
        "event_type": str(source.get("event_type", "")),
        "biological_event_track": biological_track,
        "evidence_track": row["evidence_track"],
        "pareto_dimensions": row["pareto_dimensions"],
        "legacy_weighted_rank": str(legacy_rank),
        "legacy_efficacy_score": str(source.get("efficacy_score", "")),
        "legacy_final_priority": str(source.get("final_priority", "")),
    }
    for name, state in states.items():
        state_row[f"{name}_state"] = str(state["state"])
        state_row[f"{name}_grade"] = str(state["grade"])
        state_row[f"{name}_assessed"] = "yes" if state["assessed"] else "no"
        state_row[f"{name}_reason"] = str(state["reason"])
    state_row.update({field: row[field] for field in (
        "hard_failure", "hard_failure_reasons", "legacy_priority_cap", "consensus_priority_cap",
        "evidence_grade_cap", "evidence_grade_cap_reasons",
        "manual_review_required", "manual_review_reason", "evidence_grade_uncapped",
        "evidence_grade", "evidence_consensus_score", "evidence_completeness_score",
        "evidence_assessed_layers", "evidence_missing_layers", "evidence_conflict_layers",
        "consensus_action", "recommended_next_steps", "consensus_trace",
    )})
    conflict_rows = [{
        "peptide_id": str(source.get("peptide_id", "")),
        "event_id": str(source.get("event_id", "")),
        "gene": str(source.get("gene", "")),
        "evidence_track": row["evidence_track"],
        "layer": name,
        "state": str(states[name]["state"]),
        "field": "",
        "selected_source": "",
        "selected_value": "",
        "other_source": "",
        "other_value": "",
        "precedence_version": str(source.get("evidence_source_precedence_version", "")),
        "conflict_type": "DERIVED_STATE_CONFLICT",
        "reason": str(states[name]["reason"]),
        "recommended_action": "manual evidence reconciliation before candidate progression",
    } for name in conflicts if name != "source_precedence"]
    try:
        source_details = json.loads(str(source.get("evidence_conflict_details", "[]")) or "[]")
    except json.JSONDecodeError:
        source_details = []
    for detail in source_details if isinstance(source_details, list) else []:
        conflict_rows.append({
            "peptide_id": str(source.get("peptide_id", "")),
            "event_id": str(source.get("event_id", "")),
            "gene": str(source.get("gene", "")),
            "evidence_track": row["evidence_track"],
            "layer": "source_precedence",
            "state": "CONFLICT",
            "field": str(detail.get("field", "")),
            "selected_source": str(detail.get("selected_source", "")),
            "selected_value": str(detail.get("selected_value", "")),
            "other_source": str(detail.get("other_source", "")),
            "other_value": str(detail.get("other_value", "")),
            "precedence_version": str(detail.get("precedence_version", source.get("evidence_source_precedence_version", ""))),
            "conflict_type": str(detail.get("conflict_type", "NONEMPTY_SOURCE_DISAGREEMENT")),
            "reason": f"authoritative source selected for field {detail.get('field', '')}",
            "recommended_action": "review source disagreement; retain precedence-selected value unless evidence provenance is wrong",
        })
    return row, [state_row, *conflict_rows]


def score_evidence_consensus(row: Mapping[str, Any], rules: Mapping[str, Any] | None = None) -> dict[str, str]:
    """Compatibility helper returning normalized consensus fields for one row."""

    normalized, _ = _normalized_row(row, rules or DEFAULT_RULES, 1)
    return {field: normalized[field] for field in CONSENSUS_FIELDS if field in normalized}


def _assign_pareto(rows: list[dict[str, str]]) -> None:
    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for index, row in enumerate(rows):
        row["_pareto_id"] = str(index)
        groups[(row["evidence_track"], row["evidence_grade"])].append(row)
    for (track, _grade), group in groups.items():
        dimensions = PARETO_DIMENSIONS_BY_TRACK.get(track, BASE_DIMENSIONS)
        fronts = nondominated_fronts(group, dimensions)
        for row in group:
            row["pareto_front"] = str(fronts[row["_pareto_id"]])
    for row in rows:
        row.pop("_pareto_id", None)


def _build_evidence_rank_key(row: Mapping[str, Any]) -> str:
    netmhcpan = _first_number(row, "netmhcpan_tiebreak_rank")
    mhcflurry = _first_number(row, "mhcflurry_tiebreak_score")
    return "|".join((
        str(row.get("evidence_grade", "R4")),
        str(row.get("evidence_track", "OTHER")),
        f"F{int(_number(row.get('pareto_front')) or 999999)}",
        f"{row.get('safety_state', 'UNASSESSED')}_COMPLETENESS_{int(_number(row.get('safety_completeness_grade')))}",
        str(row.get("event_authenticity_state", "EVENT_UNASSESSED")),
        str(row.get("rna_support_state", "RNA_UNASSESSED")),
        str(row.get("presentation_consensus_state", "PRESENTATION_UNASSESSED")),
        str(row.get("mutant_specificity_state", "UNASSESSED")),
        str(row.get("hla_appm_state", "HLA_LOH_UNASSESSED")),
        f"{row.get('clonality_state', 'UNASSESSED')}_{row.get('ccf_confidence_state', 'CCF_UNASSESSED')}",
        str(row.get("evidence_completeness_state", "LOW")),
        f"NETMHCPAN_EL={_format_tiebreak_number(netmhcpan)}",
        f"MHCFLURRY={_format_tiebreak_number(mhcflurry)}",
        str(row.get("peptide_id", "")),
    ))


def _rank_key(row: Mapping[str, Any], rules: Mapping[str, Any]) -> tuple[Any, ...]:
    peptide_tie_break = bool(rules.get("output", {}).get("peptide_id_tie_break", True))
    netmhcpan = _first_number(row, "netmhcpan_tiebreak_rank")
    mhcflurry = _first_number(row, "mhcflurry_tiebreak_score")
    return (
        GRADE_ORDER.get(str(row.get("evidence_grade", "R4")), 4),
        str(row.get("evidence_track", "")),
        int(_number(row.get("pareto_front")) or 999999),
        -int(_number(row.get("safety_completeness_grade"))),
        -int(_number(row.get("event_authenticity_grade"))),
        -int(_number(row.get("rna_support_grade"))),
        -int(_number(row.get("presentation_consensus_grade"))),
        -int(_number(row.get("mutant_specificity_grade"))),
        -int(_number(row.get("hla_appm_grade"))),
        -int(_number(row.get("ccf_confidence_grade"))),
        -int(_number(row.get("evidence_completeness_grade"))),
        netmhcpan if netmhcpan is not None else float("inf"),
        -mhcflurry if mhcflurry is not None else float("inf"),
        str(row.get("peptide_id", "")) if peptide_tie_break else "",
    )


def _row_text(row: Mapping[str, Any], *fields: str) -> str:
    for field in fields:
        value = str(row.get(field, "")).strip()
        if value and value.upper() not in {"NA", "N/A", "NONE", "."}:
            return value
    return ""


def _representative_fields(row: Mapping[str, Any], index: int) -> dict[str, str]:
    prefix = f"representative_{index}_"
    return {
        f"{prefix}event_id": _row_text(row, "event_id"),
        f"{prefix}peptide_id": _row_text(row, "peptide_id"),
        f"{prefix}peptide": _row_text(row, "peptide", "mutant_peptide", "mt_peptide"),
        f"{prefix}hla_allele": _row_text(row, "hla_allele", "hla", "allele", "restricting_hla"),
        f"{prefix}evidence_rank": _row_text(row, "evidence_rank"),
        f"{prefix}evidence_grade": _row_text(row, "evidence_grade"),
        f"{prefix}pareto_front": _row_text(row, "pareto_front"),
        f"{prefix}redundancy_group": _row_text(row, "redundancy_group", "peptide_redundancy_group", "overlap_group"),
        f"{prefix}evidence_rank_key": _row_text(row, "evidence_rank_key"),
    }


def _event_output(peptides: list[dict[str, str]], deduplicate: bool) -> list[dict[str, str]]:
    if not deduplicate:
        return []
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in peptides:
        event_id = _row_text(row, "event_id") or f"NO_EVENT:{row.get('peptide_id', '')}"
        phase_group = _row_text(row, "phase_group_id", "haplotype_group_id")
        group_id = f"PHASE:{phase_group}" if phase_group else f"EVENT:{event_id}"
        groups[group_id].append(row)
    events = []
    for group_id, group in groups.items():
        ordered = sorted(group, key=lambda item: int(item["evidence_rank"]))
        best_by_event_hla: list[dict[str, str]] = []
        seen_event_hla: set[tuple[str, str]] = set()
        for row in ordered:
            event_id = _row_text(row, "event_id") or f"NO_EVENT:{row.get('peptide_id', '')}"
            hla = _row_text(row, "hla_allele", "hla", "allele", "restricting_hla")
            key = (event_id, hla)
            if key in seen_event_hla:
                continue
            seen_event_hla.add(key)
            best_by_event_hla.append(row)

        representatives: list[dict[str, str]] = []
        seen_redundancy: set[str] = set()
        for row in best_by_event_hla:
            redundancy = _row_text(row, "redundancy_group", "peptide_redundancy_group", "overlap_group")
            if redundancy and redundancy in seen_redundancy:
                continue
            if redundancy:
                seen_redundancy.add(redundancy)
            representatives.append(row)
            if len(representatives) == 2:
                break
        if not representatives:
            representatives = ordered[:1]

        best = representatives[0]
        member_event_ids = sorted({_row_text(row, "event_id") for row in group if _row_text(row, "event_id")})
        phase_groups = sorted({_row_text(row, "phase_group_id", "haplotype_group_id") for row in group if _row_text(row, "phase_group_id", "haplotype_group_id")})
        event_row = {
            "event_group_id": group_id,
            "event_id": _row_text(best, "event_id"),
            "member_event_ids": ",".join(member_event_ids),
            "member_event_count": str(len(member_event_ids)),
            "phase_group_id": ",".join(phase_groups),
            "sample_id": str(best.get("sample_id", "")),
            "gene": str(best.get("gene", "")),
            "event_type": str(best.get("event_type", "")),
            "biological_event_track": str(best.get("biological_event_track", "")),
            "evidence_track": best["evidence_track"],
            "best_peptide_id": str(best.get("peptide_id", "")),
            "best_peptide": _row_text(best, "peptide", "mutant_peptide", "mt_peptide"),
            "best_hla_allele": _row_text(best, "hla_allele", "hla", "allele", "restricting_hla"),
            "best_peptide_evidence_rank": best["evidence_rank"],
            "best_evidence_grade": best["evidence_grade"],
            "best_pareto_front": best["pareto_front"],
            "best_evidence_rank_key": best["evidence_rank_key"],
            "peptide_count": str(len(group)),
            "event_hla_candidate_count": str(len(best_by_event_hla)),
            "representative_count": str(len(representatives)),
            "representative_selection_rule": "best per event_id+HLA; deduplicate redundancy_group; maximum 2 per event/phase group",
            "hard_failure_peptide_count": str(sum(row["hard_failure"] == "yes" for row in group)),
            "manual_review_required": "yes" if any(row["manual_review_required"] == "yes" for row in group) else "no",
            "consensus_action": best["consensus_action"],
            "recommended_next_steps": best["recommended_next_steps"],
            "event_consensus_trace": best["consensus_trace"],
        }
        for index, representative in enumerate(representatives, 1):
            event_row.update(_representative_fields(representative, index))
        if len(representatives) < 2:
            event_row.update(_representative_fields({}, 2))
        events.append(event_row)
    events.sort(key=lambda row: int(row["best_peptide_evidence_rank"]))
    for rank, row in enumerate(events, 1):
        row["event_evidence_rank"] = str(rank)
    return events


def _comparison_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    output = []
    for row in rows:
        legacy = int(row["legacy_weighted_rank"])
        consensus = int(row["evidence_rank"])
        shift = legacy - consensus
        direction = "improved" if shift > 0 else "decreased" if shift < 0 else "unchanged"
        output.append({
            "peptide_id": str(row.get("peptide_id", "")),
            "event_id": str(row.get("event_id", "")),
            "gene": str(row.get("gene", "")),
            "evidence_track": row["evidence_track"],
            "legacy_weighted_rank": str(legacy),
            "legacy_efficacy_score": str(row.get("efficacy_score", "")),
            "legacy_final_priority": str(row.get("final_priority", "")),
            "evidence_rank": str(consensus),
            "rank_shift_weighted_minus_consensus": str(shift),
            "rank_shift_direction": direction,
            "evidence_grade": row["evidence_grade"],
            "pareto_front": row["pareto_front"],
            "evidence_rank_key": row["evidence_rank_key"],
            "consensus_action": row["consensus_action"],
            "recommended_next_steps": row["recommended_next_steps"],
            "difference_reason": f"{row['consensus_trace']};rank_{direction};grade={row['evidence_grade']};pareto_front={row['pareto_front']}",
        })
    return output


def _write_provenance(path: str | Path, rules: Mapping[str, Any], result: Mapping[str, Any]) -> None:
    target = Path(path)
    payload: dict[str, Any] = {}
    if target.is_file():
        try:
            payload = json.loads(target.read_text())
        except (OSError, json.JSONDecodeError):
            payload = {}
    rules_json = json.dumps(rules, ensure_ascii=True, sort_keys=True, separators=(",", ":"))
    payload["evidence_consensus"] = {
        "algorithm": "discrete_state_grade_track_pareto_v2",
        "status": rules.get("metadata", {}).get("status", "PROVISIONAL_RESEARCH_ONLY"),
        "rules_name": rules.get("metadata", {}).get("name", ""),
        "rules_version": rules.get("metadata", {}).get("version", ""),
        "rules_sha256": hashlib.sha256(rules_json.encode()).hexdigest(),
        "replace_weighted_ranking": False,
        "outputs": {key: value for key, value in result.items() if key.startswith("output_")},
    }
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n")


def rank_evidence_consensus(
    comprehensive_tsv: str | Path,
    output_peptides_tsv: str | Path,
    output_events_tsv: str | Path,
    output_states_tsv: str | Path,
    output_conflicts_tsv: str | Path,
    rules: Mapping[str, Any],
    provenance_json: str | Path | None = None,
) -> dict[str, Any]:
    """Run the independent provisional evidence-consensus ranking."""

    source_rows = read_tsv(comprehensive_tsv)
    original_fields = list(source_rows[0]) if source_rows else []
    rows: list[dict[str, str]] = []
    states: list[dict[str, str]] = []
    conflicts: list[dict[str, str]] = []
    for legacy_rank, source in enumerate(source_rows, 1):
        normalized, records = _normalized_row(source, rules, legacy_rank)
        rows.append(normalized)
        states.append(records[0])
        conflicts.extend(records[1:])
    _assign_pareto(rows)
    for row in rows:
        row["evidence_rank_key"] = _build_evidence_rank_key(row)
    rows.sort(key=lambda row: _rank_key(row, rules))
    tracks: Counter[str] = Counter()
    for rank, row in enumerate(rows, 1):
        row["evidence_rank"] = str(rank)
        tracks[row["evidence_track"]] += 1
        row["track_rank"] = str(tracks[row["evidence_track"]])

    peptide_fields = original_fields + [field for field in CONSENSUS_FIELDS if field not in original_fields]
    write_tsv(output_peptides_tsv, rows, peptide_fields)
    write_tsv(output_states_tsv, states)
    write_tsv(output_conflicts_tsv, conflicts, CONFLICT_FIELDS)
    event_rows = _event_output(rows, bool(rules.get("output", {}).get("event_deduplicate", True)))
    write_tsv(output_events_tsv, event_rows)
    comparison_path = Path(output_peptides_tsv).with_name("weighted_vs_consensus_comparison.tsv")
    write_tsv(comparison_path, _comparison_rows(rows))
    result = {
        "rows": len(rows),
        "events": len(event_rows),
        "conflicts": len(conflicts),
        "output_peptides": str(output_peptides_tsv),
        "output_events": str(output_events_tsv),
        "output_states": str(output_states_tsv),
        "output_conflicts": str(output_conflicts_tsv),
        "output_comparison": str(comparison_path),
        "grade_counts": dict(sorted(Counter(row["evidence_grade"] for row in rows).items())),
        "track_counts": dict(sorted(Counter(row["evidence_track"] for row in rows).items())),
        "legacy_ranking_modified": False,
    }
    if provenance_json:
        _write_provenance(provenance_json, rules, result)
    return result


def build_evidence_consensus(
    comprehensive_tsv: str | Path,
    outdir: str | Path,
    rules: Mapping[str, Any] | None = None,
    *,
    peptide_output: str | Path | None = None,
    provenance_json: str | Path | None = None,
) -> dict[str, Any]:
    output_dir = Path(outdir)
    output_dir.mkdir(parents=True, exist_ok=True)
    result = rank_evidence_consensus(
        comprehensive_tsv,
        peptide_output or output_dir / "ranked_peptides.evidence_consensus.tsv",
        output_dir / "ranked_events.evidence_consensus.tsv",
        output_dir / "evidence_states.tsv",
        output_dir / "evidence_conflicts.tsv",
        rules or load_consensus_rules(),
        provenance_json,
    )
    return {
        "rows": result["rows"],
        "ranked_peptides": result["output_peptides"],
        "ranked_events": result["output_events"],
        "evidence_states": result["output_states"],
        "evidence_conflicts": result["output_conflicts"],
        "comparison": result["output_comparison"],
        "grade_counts": result["grade_counts"],
        "track_counts": result["track_counts"],
        "legacy_ranking_modified": False,
    }


def rank_by_evidence_consensus(
    input_tsv: str | Path,
    output_tsv: str | Path,
    rules: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Compatibility API that writes the complete bundle beside output_tsv."""

    result = build_evidence_consensus(input_tsv, Path(output_tsv).parent, rules, peptide_output=output_tsv)
    return {
        "rows": result["rows"],
        "output": result["ranked_peptides"],
        **{key: value for key, value in result.items() if key not in {"rows", "ranked_peptides"}},
    }
