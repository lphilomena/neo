"""Three-layer neoantigen data model.

Layer 1 — mutation_source: DNA-level event origin (SNV / InDel / SV / CNV / TAA).
Layer 2 — peptide_consequence: transcript/protein neoantigen class.
Layer 3 — immunology evidence: unified scoring domain (binding, presentation, …).
"""

from __future__ import annotations

from typing import Any, Mapping

from .utils import clamp, norm_tpm, to_float, truthy

# Layer 1: upstream mutation sources
MUTATION_SOURCES = frozenset({"SNV", "InDel", "SV", "CNV", "TAA", "Other"})

# Layer 2: peptide / neoantigen consequence classes
PEPTIDE_CONSEQUENCES = frozenset({
    "missense",
    "frameshift",
    "fusion",
    "splice_junction",
    "exon_deletion_junction",
    "insertion",
    "other",
    "noncoding",
})

# Layer 3: immunology evaluation dimensions (Project B unified scoring)
IMMUNOLOGY_DIMENSIONS = (
    "hla_binding",
    "hla_presentation",
    "expression",
    "rna_junction_support",
    "clonality",
    "tumor_specificity",
    "normal_tissue_safety",
    "apm_integrity",
    "immunogenicity",
)

L3_FIELD_PREFIX = "l3_"


def _norm(s: Any) -> str:
    return str(s or "").strip()


def infer_mutation_source(
    *,
    event_type: str = "",
    tool: str = "",
    consequence: str = "",
) -> str:
    """Infer Layer-1 mutation_source from legacy event_type / tool / consequence."""
    et = _norm(event_type)
    et_u = et.upper()
    cons = _norm(consequence).lower()
    tool_l = _norm(tool).lower()

    if et.startswith("SV_") or et_u == "SV":
        return "SV"
    if et_u in {"SNV", "SNP"}:
        return "SNV"
    if et_u in {"INDEL", "INS", "DEL", "INSERTION", "DELETION"}:
        return "InDel"
    if et_u == "TAA":
        return "TAA"
    if et_u == "CNV" or "copy_number" in cons:
        return "CNV"
    if "frameshift" in cons or "inframe_indel" in cons:
        return "InDel"
    if et_u == "FUSION" or tool_l in {"pvacfuse", "easyfuse", "svphase1"} or "svphase1" in tool_l:
        return "SV"
    if tool_l == "pvacseq":
        if "indel" in cons or et_u in {"INS", "DEL"}:
            return "InDel"
        return "SNV"
    return "Other"


def infer_peptide_consequence(
    *,
    event_type: str = "",
    consequence: str = "",
    tool: str = "",
) -> str:
    """Infer Layer-2 peptide_consequence from event metadata."""
    et = _norm(event_type)
    cons = _norm(consequence).lower()
    tool_l = _norm(tool).lower()

    sv_map = {
        "SV_Fusion": "fusion",
        "SV_Frameshift": "frameshift",
        "SV_Junction": "splice_junction",
        "SV_Insertion": "insertion",
        "SV_Noncoding": "noncoding",
    }
    if et in sv_map:
        return sv_map[et]

    if et == "Fusion" or tool_l == "pvacfuse" or "fusion" in cons:
        return "fusion"
    if et == "Splice" or tool_l == "pvacsplice" or "splice" in cons:
        return "splice_junction"
    if et == "TAA":
        return "other"
    if "missense" in cons:
        return "missense"
    if "frameshift" in cons:
        return "frameshift"
    if "inframe_deletion" in cons or "disruptive_inframe" in cons:
        return "exon_deletion_junction"
    if "insertion" in cons or "duplication" in cons:
        return "insertion"
    if et in {"SNV", "InDel", "INDEL"}:
        if "frameshift" in cons:
            return "frameshift"
        if et != "SNV" or "missense" in cons:
            return "missense" if et == "SNV" else "frameshift"
        return "missense"
    if cons:
        return "other"
    return "other"


def enrich_event_layers(event: Mapping[str, Any]) -> dict[str, str]:
    """Fill Layer-1/2 fields on an event row (non-destructive for explicit values)."""
    out = dict(event)
    if not _norm(out.get("mutation_source")):
        out["mutation_source"] = infer_mutation_source(
            event_type=_norm(out.get("event_type")),
            tool=_norm(out.get("source")),
            consequence=_norm(out.get("consequence")),
        )
    if not _norm(out.get("peptide_consequence")):
        out["peptide_consequence"] = infer_peptide_consequence(
            event_type=_norm(out.get("event_type")),
            consequence=_norm(out.get("consequence")),
            tool=_norm(out.get("source")),
        )
    return out


def enrich_peptide_layers(
    peptide: Mapping[str, Any],
    event: Mapping[str, Any] | None = None,
) -> dict[str, str]:
    """Fill Layer-1/2 fields on a peptide row, optionally from parent event."""
    out = dict(peptide)
    ev = event or {}
    if not _norm(out.get("mutation_source")):
        out["mutation_source"] = _norm(ev.get("mutation_source")) or infer_mutation_source(
            event_type=_norm(out.get("event_type") or ev.get("event_type")),
            tool=_norm(out.get("source_tool") or ev.get("source")),
            consequence=_norm(ev.get("consequence")),
        )
    if not _norm(out.get("peptide_consequence")):
        out["peptide_consequence"] = _norm(ev.get("peptide_consequence")) or infer_peptide_consequence(
            event_type=_norm(out.get("event_type") or ev.get("event_type")),
            consequence=_norm(ev.get("consequence")),
            tool=_norm(out.get("source_tool") or ev.get("source")),
        )
    if not _norm(out.get("rna_junction_reads")) and ev:
        out["rna_junction_reads"] = _norm(ev.get("rna_junction_reads"))
    return out


def l3_dimension_field(name: str) -> str:
    return f"{L3_FIELD_PREFIX}{name}_score"


def default_l3_weights(profile: Mapping[str, Any]) -> dict[str, float]:
    """Layer-3 weights; falls back to legacy score_weights mapping."""
    explicit = profile.get("l3_weights")
    if explicit:
        return {k: float(v) for k, v in explicit.items()}
    vw = profile.get("score_weights", {})
    ew = profile.get("event_weights", {})
    return {
        "event_confidence": float(ew.get("event_confidence", 0.10)),
        "expression": float(ew.get("event_expression", 0.08)),
        "clonality": float(ew.get("clonality", 0.07)),
        "tumor_specificity": float(ew.get("tumor_specificity", 0.05)),
        "hla_binding": float(vw.get("binding_evidence", 0.20)),
        "hla_presentation": float(vw.get("presentation_evidence", 0.25)),
        "rna_junction_support": 0.05,
        "normal_tissue_safety": 0.05,
        "apm_integrity": 0.05,
        "immunogenicity": float(vw.get("immunogenicity", 0.15)),
    }


def safety_dimension_score(status: str) -> float:
    s = _norm(status).upper()
    if s == "FAIL":
        return 0.0
    if s == "CAUTION":
        return 0.5
    return 1.0


def rna_junction_dimension_score(reads: Any) -> float:
    n = to_float(reads, 0.0)
    if n <= 0:
        return 0.0
    return clamp(n / 10.0)


def compute_l3_dimension_scores(
    peptide: Mapping[str, Any],
    event: Mapping[str, Any],
    presentation: Mapping[str, Any],
    profile: Mapping[str, Any],
    *,
    appm: float = 1.0,
    ccf: float = 1.0,
    immuno: float = 0.0,
    high_expression_tpm: float = 20.0,
) -> dict[str, str]:
    """Compute explicit Layer-3 dimension scores (0–1) for unified immunology scoring."""
    safety = peptide.get("safety_status") or event.get("safety_status", "PASS")
    scores = {
        "event_confidence": clamp(to_float(event.get("event_confidence"), 0.0)),
        "expression": norm_tpm(event.get("event_expression"), high_expression_tpm),
        "clonality": clamp(to_float(event.get("clonality"), 0.0)),
        "tumor_specificity": clamp(to_float(event.get("tumor_specificity"), 0.0)),
        "hla_binding": clamp(to_float(presentation.get("binding_evidence_score"), 0.0)),
        "hla_presentation": clamp(to_float(presentation.get("presentation_evidence_score"), 0.0)),
        "rna_junction_support": rna_junction_dimension_score(
            peptide.get("rna_junction_reads") or event.get("rna_junction_reads")
        ),
        "normal_tissue_safety": safety_dimension_score(str(safety)),
        "apm_integrity": clamp(appm),
        "immunogenicity": clamp(immuno),
    }
    weights = default_l3_weights(profile)
    total_w = sum(weights.get(k, 0.0) for k in scores)
    composite = 0.0
    if total_w > 0:
        composite = sum(scores[k] * weights.get(k, 0.0) for k in scores) / total_w
    out: dict[str, str] = {l3_dimension_field(k): f"{v:.4f}" for k, v in scores.items()}
    out["immunology_composite_score"] = f"{clamp(composite):.4f}"
    return out
