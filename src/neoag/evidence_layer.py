"""Materialize standard intermediate evidence TSVs for multi-entry Project B scoring."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from .config import load_profile
from .model_layers import enrich_event_layers, enrich_peptide_layers
from .safety import apply_event_safety, apply_peptide_safety, load_normal_expression, load_normal_hla_ligands
from .schemas import (
    EXPRESSION_EVIDENCE_FIELDS,
    RNA_JUNCTION_EVIDENCE_FIELDS,
    SAFETY_EVIDENCE_FIELDS,
)
from .evidence_provenance import ProvenanceRecord, provenance_derived, attach_provenance
from .utils import first, read_tsv, safe_id, to_float, write_tsv


def _expression_source(expression_path: str | Path | None) -> str:
    if expression_path and Path(expression_path).is_file():
        return str(expression_path)
    return "raw_events.event_expression"


def build_expression_evidence(
    raw_events: str | Path,
    out_path: str | Path,
    *,
    expression_path: str | Path | None = None,
    sample_id: str = "",
    provenance: ProvenanceRecord | None = None,
) -> list[dict[str, str]]:
    """Event-level expression evidence from raw events (+ optional TPM table join)."""
    tpm_map: dict[str, float] = {}
    if expression_path and Path(expression_path).is_file():
        from .sv.evidence import load_expression

        tpm_map = load_expression(expression_path)

    rows: list[dict[str, str]] = []
    for ev in read_tsv(raw_events):
        ev = enrich_event_layers(ev)
        gene = str(ev.get("gene") or "")
        tpm = to_float(ev.get("event_expression"), 0.0)
        if gene in tpm_map:
            tpm = max(tpm, tpm_map[gene])
        rows.append({
            "event_id": ev.get("event_id", ""),
            "sample_id": ev.get("sample_id", sample_id),
            "gene": gene,
            "event_expression": f"{tpm:.4f}",
            "expression_tpm": f"{tpm:.4f}",
            "expression_source": _expression_source(expression_path),
            "mutation_source": ev.get("mutation_source", ""),
            "peptide_consequence": ev.get("peptide_consequence", ""),
        })
    prov = provenance or provenance_derived("expression_evidence", out_path, upstream=_expression_source(expression_path))
    write_tsv(out_path, attach_provenance(rows, prov), EXPRESSION_EVIDENCE_FIELDS)
    return rows


def _rna_support_status(alt_reads: str, depth: str, vaf: str, junction_reads: int) -> str:
    alt = to_float(alt_reads, 0.0)
    dp = to_float(depth, 0.0)
    vf = to_float(vaf, 0.0)
    if alt > 0 and (dp > 0 or vf > 0):
        return "RNA_ALT_SUPPORTED"
    if junction_reads > 0:
        return "RNA_JUNCTION_SUPPORTED"
    if dp > 0:
        return "RNA_ALT_NOT_DETECTED"
    return "UNASSESSED"


def _targeted_validation(reads: int, source: str, mutation_source: str, consequence: str) -> dict[str, str]:
    method = "junction"
    text = f"{mutation_source} {consequence}".lower()
    if "fusion" in text:
        method = "fusion_targeted_rna"
    elif "splice" in text or "junction" in text:
        method = "splice_junction_targeted_rna"
    if not source:
        return {
            "targeted_validation_status": "UNASSESSED",
            "targeted_validation_source": "",
            "targeted_validation_method": method,
        }
    return {
        "targeted_validation_status": "SUPPORTED" if reads > 0 else "NO_TARGETED_SUPPORT",
        "targeted_validation_source": source,
        "targeted_validation_method": method,
    }


def build_rna_junction_evidence(
    raw_events: str | Path,
    raw_peptides: str | Path,
    out_path: str | Path,
    *,
    junction_path: str | Path | None = None,
    fusion_evidence_path: str | Path | None = None,
    rna_vaf_path: str | Path | None = None,
    sample_id: str = "",
    provenance: ProvenanceRecord | None = None,
) -> list[dict[str, str]]:
    """RNA allele and junction support at event and peptide level."""
    from .adapters.rna_vaf import choose_rna_vaf_support, load_rna_vaf_support

    extra: dict[str, int] = {}
    if junction_path and Path(junction_path).is_file():
        from .sv.evidence import load_junction_reads

        extra = load_junction_reads(junction_path)

    fusion_reads: dict[str, int] = {}
    if fusion_evidence_path and Path(fusion_evidence_path).is_file():
        for row in read_tsv(fusion_evidence_path):
            if row.get("filter_status") != "pass":
                continue
            eid = row.get("event_id", "")
            reads = int(to_float(row.get("rna_junction_reads"), 0.0))
            if eid:
                fusion_reads[eid] = max(fusion_reads.get(eid, 0), reads)

    rna_vaf = load_rna_vaf_support(rna_vaf_path)

    def build_row(source_row: Mapping[str, Any], *, peptide_id: str, row_source: str) -> dict[str, str]:
        enriched = enrich_peptide_layers(source_row) if peptide_id else enrich_event_layers(source_row)
        eid = enriched.get("event_id", "")
        gene = str(enriched.get("gene") or "")
        reads = int(to_float(enriched.get("rna_junction_reads"), 0.0))
        for key in (eid, gene, gene.replace("::", "_")):
            if key and key in extra:
                reads = max(reads, extra[key])
        if eid in fusion_reads:
            reads = max(reads, fusion_reads[eid])
        source = ""
        if eid in fusion_reads:
            source = str(fusion_evidence_path)
        elif any(k in extra for k in (eid, gene, gene.replace("::", "_")) if k):
            source = str(junction_path) if junction_path else "rna_junction"
        elif reads > 0:
            source = row_source
        rna = choose_rna_vaf_support(enriched, rna_vaf, row_source)
        rna_fields = rna.as_fields()
        mutation_source = enriched.get("mutation_source", "")
        consequence = enriched.get("peptide_consequence", "")
        return {
            "evidence_id": safe_id(f"RNAJ_{peptide_id or eid}"),
            "event_id": eid,
            "peptide_id": peptide_id,
            "sample_id": enriched.get("sample_id", sample_id),
            "gene": gene,
            "gene_pair": gene if "::" in gene else "",
            "junction_reads": str(reads),
            "junction_source": source or row_source,
            "mutation_source": mutation_source,
            "peptide_consequence": consequence,
            **rna_fields,
            "rna_support_status": _rna_support_status(
                rna_fields.get("rna_alt_reads", ""),
                rna_fields.get("rna_depth", ""),
                rna_fields.get("rna_vaf", ""),
                reads,
            ),
            **_targeted_validation(reads, source, mutation_source, consequence),
        }

    rows: list[dict[str, str]] = []
    seen: set[str] = set()

    for ev in read_tsv(raw_events):
        row = build_row(ev, peptide_id="", row_source="raw_events")
        if row["evidence_id"] not in seen:
            seen.add(row["evidence_id"])
            rows.append(row)

    for pep in read_tsv(raw_peptides):
        pid = pep.get("peptide_id", "")
        row = build_row(pep, peptide_id=pid, row_source="raw_peptides")
        if row["evidence_id"] in seen:
            continue
        seen.add(row["evidence_id"])
        rows.append(row)

    prov = provenance or provenance_derived(
        "rna_junction_evidence",
        out_path,
        upstream=str(rna_vaf_path or junction_path or fusion_evidence_path or "raw_events"),
    )
    write_tsv(out_path, attach_provenance(rows, prov), RNA_JUNCTION_EVIDENCE_FIELDS)
    return rows

def build_safety_evidence(
    raw_events: str | Path,
    raw_peptides: str | Path,
    out_path: str | Path,
    profile: Mapping[str, Any],
    *,
    normal_expression: str | Path | None = None,
    normal_hla_ligands: str | Path | None = None,
    provenance: ProvenanceRecord | None = None,
) -> list[dict[str, str]]:
    """Pre-score safety evidence for events and peptides."""
    norm_expr = load_normal_expression(normal_expression)
    norm_lig = load_normal_hla_ligands(normal_hla_ligands)
    event_map: dict[str, dict[str, str]] = {}
    rows: list[dict[str, str]] = []

    for ev in read_tsv(raw_events):
        ev = enrich_event_layers(dict(ev))
        ev = apply_event_safety(ev, profile, norm_expr)
        event_map[ev["event_id"]] = ev
        rows.append({
            "evidence_id": safe_id(f"SAFE_EVT_{ev['event_id']}"),
            "level": "event",
            "event_id": ev.get("event_id", ""),
            "peptide_id": "",
            "sample_id": ev.get("sample_id", ""),
            "gene": ev.get("gene", ""),
            "peptide": "",
            "safety_status": ev.get("safety_status", ""),
            "safety_reason": ev.get("safety_reason", ""),
            "normal_tissue_max_tpm": ev.get("normal_tissue_max_tpm", ""),
            "normal_hspc_tpm": ev.get("normal_hspc_tpm", ""),
            "critical_tissue_hit": ev.get("critical_tissue_hit", ""),
            "normal_hla_ligand_overlap": "",
        })

    for pep in read_tsv(raw_peptides):
        pep = enrich_peptide_layers(dict(pep))
        ev = event_map.get(pep.get("event_id", ""), {})
        if not ev:
            continue
        pep = apply_peptide_safety(pep, ev, profile, norm_lig)
        rows.append({
            "evidence_id": safe_id(f"SAFE_PEP_{pep.get('peptide_id', '')}"),
            "level": "peptide",
            "event_id": pep.get("event_id", ""),
            "peptide_id": pep.get("peptide_id", ""),
            "sample_id": pep.get("sample_id", ""),
            "gene": pep.get("gene", ""),
            "peptide": pep.get("peptide", ""),
            "safety_status": pep.get("safety_status", ""),
            "safety_reason": pep.get("safety_reason", ""),
            "normal_tissue_max_tpm": ev.get("normal_tissue_max_tpm", ""),
            "normal_hspc_tpm": ev.get("normal_hspc_tpm", ""),
            "critical_tissue_hit": ev.get("critical_tissue_hit", ""),
            "normal_hla_ligand_overlap": pep.get("normal_hla_ligand_overlap", ""),
        })

    prov = provenance or provenance_derived(
        "safety_evidence",
        out_path,
        upstream=f"normal_expression:{normal_expression};normal_hla_ligands:{normal_hla_ligands}",
    )
    write_tsv(out_path, attach_provenance(rows, prov), SAFETY_EVIDENCE_FIELDS)
    return rows


def build_standard_evidence_layer(
    outdir: str | Path,
    profile: Mapping[str, Any] | str,
    *,
    raw_events: str | Path | None = None,
    raw_peptides: str | Path | None = None,
    expression: str | Path | None = None,
    rna_junction: str | Path | None = None,
    fusion_evidence: str | Path | None = None,
    rna_vaf: str | Path | None = None,
    normal_expression: str | Path | None = None,
    normal_hla_ligands: str | Path | None = None,
    sample_id: str = "",
) -> dict[str, str]:
    """Write expression / RNA junction / safety evidence under the standard layout."""
    outdir = Path(outdir)
    if isinstance(profile, str):
        profile = load_profile(profile)

    parsed = outdir / "parsed"
    safety_dir = outdir / "safety"
    parsed.mkdir(parents=True, exist_ok=True)
    safety_dir.mkdir(parents=True, exist_ok=True)

    events_path = Path(raw_events) if raw_events else parsed / "raw_events.tsv"
    peptides_path = Path(raw_peptides) if raw_peptides else parsed / "raw_peptides.tsv"
    fusion_evidence_path = Path(fusion_evidence) if fusion_evidence else parsed / "fusion_evidence.tsv"
    if not fusion_evidence_path.is_file():
        fusion_evidence_path = None

    expr_out = parsed / "expression_evidence.tsv"
    rna_out = parsed / "rna_junction_evidence.tsv"
    safe_out = safety_dir / "safety_evidence.tsv"

    build_expression_evidence(events_path, expr_out, expression_path=expression, sample_id=sample_id)
    build_rna_junction_evidence(
        events_path,
        peptides_path,
        rna_out,
        junction_path=rna_junction,
        fusion_evidence_path=fusion_evidence_path,
        rna_vaf_path=rna_vaf,
        sample_id=sample_id,
    )
    build_safety_evidence(
        events_path,
        peptides_path,
        safe_out,
        profile,
        normal_expression=normal_expression,
        normal_hla_ligands=normal_hla_ligands,
    )

    return {
        "expression_evidence": str(expr_out),
        "rna_junction_evidence": str(rna_out),
        "safety_evidence": str(safe_out),
        **(
            {"fusion_evidence": str(fusion_evidence_path)}
            if fusion_evidence_path
            else {}
        ),
    }
