"""TESLA-style presentation gates before immunogenicity-weighted ranking."""

from __future__ import annotations

from typing import Any, Mapping

from .utils import to_float


def evaluate_presentation_gate(
    peptide: Mapping[str, Any],
    event: Mapping[str, Any],
    presentation: Mapping[str, Any],
    profile: Mapping[str, Any],
) -> dict[str, str]:
    """Return gate status fields applied as a multiplier in peptide scoring."""
    cfg = profile.get("gates", {})
    if not cfg.get("enabled", True):
        return {
            "presentation_gate_status": "PASS",
            "presentation_gate_reason": "gates_disabled",
            "presentation_gate_multiplier": "1.0000",
        }

    reasons: list[str] = []
    passed = True

    grades = cfg.get("require_presentation_grades") or []
    if grades:
        grade = presentation.get("presentation_evidence_grade", "")
        if grade not in grades:
            passed = False
            reasons.append(f"grade={grade or 'missing'}")

    max_el = float(cfg.get("max_el_rank", 2.0))
    if max_el > 0:
        el_raw = presentation.get("netmhcpan_el_rank") or peptide.get("netmhcpan_el_rank", "")
        if str(el_raw).strip() not in {"", "99", "99.0"}:
            el = to_float(el_raw, 99.0)
            if el > max_el:
                passed = False
                reasons.append(f"el_rank={el:.2f}")

    min_stab = float(cfg.get("min_stabpan_score", 0.0))
    if min_stab > 0:
        stab_raw = presentation.get("netmhcstabpan_score") or peptide.get("netmhcstabpan_score", "")
        if str(stab_raw).strip():
            stab = to_float(stab_raw, 0.0)
            if stab < min_stab:
                passed = False
                reasons.append(f"stabpan={stab:.2f}")

    min_tpm = float(cfg.get("min_event_expression_tpm", 0.0))
    if min_tpm > 0:
        tpm = to_float(event.get("event_expression"), 0.0)
        if tpm < min_tpm:
            passed = False
            reasons.append(f"tpm={tpm:.2f}")

    min_vaf = float(cfg.get("min_tumor_vaf", 0.0))
    if min_vaf > 0:
        vaf = to_float(event.get("tumor_vaf"), 0.0)
        if vaf < min_vaf:
            passed = False
            reasons.append(f"vaf={vaf:.4f}")

    consequence = str(event.get("peptide_consequence") or "").lower()
    mutation_source = str(event.get("mutation_source") or "").upper()
    is_junction = consequence in {"fusion", "splice_junction"} and mutation_source != "INDEL"
    if is_junction:
        min_junction = float(cfg.get("min_rna_junction_reads", 0.0))
        if min_junction > 0:
            reads = to_float(event.get("rna_junction_reads"), 0.0)
            if reads < min_junction:
                passed = False
                reasons.append(f"rna_junction_reads={reads:.0f}")
    else:
        min_alt = float(cfg.get("min_rna_alt_reads", 0.0))
        if min_alt > 0:
            alt_reads = to_float(event.get("rna_alt_reads"), 0.0)
            if alt_reads < min_alt:
                passed = False
                reasons.append(f"rna_alt_reads={alt_reads:.0f}")
        min_rna_vaf = float(cfg.get("min_rna_vaf", 0.0))
        if min_rna_vaf > 0:
            rna_vaf = to_float(event.get("rna_vaf"), 0.0)
            if rna_vaf < min_rna_vaf:
                passed = False
                reasons.append(f"rna_vaf={rna_vaf:.4f}")
        min_allele_expression = float(cfg.get("min_allele_expression", 0.0))
        if min_allele_expression > 0:
            gene_tpm = to_float(
                event.get("gene_expression_tpm") or event.get("event_expression"), 0.0
            )
            rna_vaf = to_float(event.get("rna_vaf"), 0.0)
            allele_expression = gene_tpm * rna_vaf
            if allele_expression < min_allele_expression:
                passed = False
                reasons.append(f"allele_expression={allele_expression:.4f}")

    if passed:
        return {
            "presentation_gate_status": "PASS",
            "presentation_gate_reason": "ok",
            "presentation_gate_multiplier": "1.0000",
        }

    mode = str(cfg.get("failure_mode", "multiply")).lower()
    fail_mult = float(cfg.get("failure_multiplier", 0.25))
    if mode == "hard_zero":
        fail_mult = 0.0

    return {
        "presentation_gate_status": "FAIL",
        "presentation_gate_reason": ";".join(reasons) if reasons else "gate_failed",
        "presentation_gate_multiplier": f"{fail_mult:.4f}",
    }
