"""Run v0.3 presentation + scoring on SV Phase 1 raw_events/raw_peptides."""

from __future__ import annotations

import datetime
import shutil
from pathlib import Path

from ..adapters.mhcflurry import parse_mhcflurry, write_mhcflurry_evidence
from ..adapters.netmhcpan import parse_netmhcpan, write_netmhcpan_evidence
from ..appm_lite import build_appm_lite
from ..ccf_v2 import build_ccf_2
from ..config import load_profile
from ..immunogenicity_composite import apply_immunogenicity_evidence, run_immunogenicity_predictors
from ..presentation import build_presentation_evidence
from ..reports_dual import load_report_bundle, make_dual_reports
from ..schemas import PRESENTATION_FIELDS
from ..scoring import score
from ..tools.registry import RunContext
from ..tools.runner import run_mhcflurry, run_netmhcpan
from ..utils import copy_if_different, read_tsv, write_json, write_tsv
from ..validation import make_validation_plan
from ..evidence_provenance import ProvenanceRegistry
from ..evidence_layer import build_standard_evidence_layer
from ..peptide_safety_gate import build_peptide_safety_gate
from ..immune_escape import build_immune_escape_evidence


def _resolve_optional(path: str | Path | None) -> Path | None:
    if not path:
        return None
    p = Path(path)
    return p if p.is_file() else None


def run_sv_score(
    *,
    outdir: str | Path,
    profile_name_or_path: str,
    sample_id: str,
    raw_events: str | Path,
    raw_peptides: str | Path,
    netmhcpan: str | Path | None = None,
    mhcflurry: str | Path | None = None,
    vep_appm: str | Path | None = None,
    expression: str | Path | None = None,
    hla_loh: str | Path | None = None,
    purity: str | Path | None = None,
    cnv: str | Path | None = None,
    normal_expression: str | Path | None = None,
    normal_hla_ligands: str | Path | None = None,
    binding_stub: bool = False,
    immunogenicity_stub: bool = False,
    run_binding: bool = True,
    tool_executables: dict[str, str] | None = None,
    reference_proteome: str | Path | None = None,
    normal_junctions: str | Path | None = None,
) -> dict[str, str]:
    """NetMHCpan/MHCflurry → presentation → APPM/CCF → score → report."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    profile = load_profile(profile_name_or_path)

    parsed = outdir / "parsed"
    parsed.mkdir(exist_ok=True)
    pres = outdir / "presentation"
    pres.mkdir(exist_ok=True)
    appm_dir = outdir / "appm"
    clon = outdir / "clonality"
    scoring = outdir / "scoring"
    safety_dir = outdir / "safety"
    immune_dir = outdir / "immune_escape"
    reports = outdir / "reports"
    for d in (appm_dir, clon, scoring, safety_dir, immune_dir, reports):
        d.mkdir(exist_ok=True)

    events_path = parsed / "raw_events.tsv"
    peptides_path = parsed / "raw_peptides.tsv"
    copy_if_different(raw_events, events_path)
    copy_if_different(raw_peptides, peptides_path)

    tools_dir = outdir / "tools"
    tools_dir.mkdir(exist_ok=True)
    net_raw = _resolve_optional(netmhcpan)
    mhc_raw = _resolve_optional(mhcflurry)
    provenance_registry = ProvenanceRegistry()
    ctx = RunContext(
        sample_id=sample_id,
        outdir=outdir,
        stub=binding_stub,
        raw_peptides=peptides_path,
        executables=tool_executables or {},
    )
    if run_binding:
        if not net_raw:
            net_raw = tools_dir / "netmhcpan.xls"
            run_netmhcpan(ctx, net_raw)
        net_prov = ctx.tool_provenance.get("netmhcpan")
        if net_prov:
            provenance_registry.set(net_prov)
        elif net_raw:
            provenance_registry.register_passthrough("netmhcpan", net_raw)
        if not mhc_raw:
            mhc_raw = tools_dir / "mhcflurry.csv"
            run_mhcflurry(ctx, mhc_raw)
        if mhc_raw and not provenance_registry.has("mhcflurry"):
            provenance_registry.register_real("mhcflurry", mhc_raw) if not binding_stub else provenance_registry.register_stub("mhcflurry")

    net_path = pres / "netmhcpan_evidence.tsv"
    mhc_path = pres / "mhcflurry_evidence.tsv"
    write_netmhcpan_evidence(
        net_path,
        parse_netmhcpan(net_raw, sample_id),
        ctx.tool_provenance.get("netmhcpan") if net_raw else None,
    )
    write_mhcflurry_evidence(mhc_path, parse_mhcflurry(mhc_raw, sample_id))

    pres_path = pres / "presentation_evidence.tsv"
    build_presentation_evidence(
        peptides_path, net_path, mhc_path, profile, pres_path, provenance_registry=provenance_registry,
    )
    pres_rows = read_tsv(pres_path)
    immuno_ctx = RunContext(
        sample_id=sample_id,
        outdir=outdir,
        stub=immunogenicity_stub,
        raw_peptides=peptides_path,
        executables=tool_executables or {},
    )
    immuno_paths = run_immunogenicity_predictors(
        peptides_path, outdir, profile, immuno_ctx, provenance_registry=provenance_registry,
    )
    apply_immunogenicity_evidence(pres_rows, immuno_paths, profile)
    write_tsv(pres_path, pres_rows, PRESENTATION_FIELDS)

    _, appm_summary = build_appm_lite(
        sample_id,
        _resolve_optional(vep_appm),
        _resolve_optional(expression),
        _resolve_optional(hla_loh),
        profile,
        appm_dir,
        cnv_tsv=_resolve_optional(cnv),
        raw_peptides=peptides_path,
        tumor_purity_tsv=_resolve_optional(purity),
    )
    ccf_path = clon / "ccf_2.tsv"
    ccf_lite_alias = clon / "ccf_lite.tsv"
    build_ccf_2(
        events_path,
        _resolve_optional(purity),
        _resolve_optional(cnv),
        profile,
        ccf_path,
    )
    shutil.copy2(ccf_path, ccf_lite_alias)

    evidence_paths = build_standard_evidence_layer(
        outdir,
        profile,
        raw_events=events_path,
        raw_peptides=peptides_path,
        expression=_resolve_optional(expression),
        normal_expression=_resolve_optional(normal_expression),
        normal_hla_ligands=_resolve_optional(normal_hla_ligands),
        sample_id=sample_id,
    )
    peptide_safety_path = safety_dir / "peptide_safety.tsv"
    event_safety_path = safety_dir / "event_safety.tsv"
    build_peptide_safety_gate(raw_events=events_path, raw_peptides=peptides_path, out_peptide_safety=peptide_safety_path, out_event_safety=event_safety_path, profile=profile, normal_expression=_resolve_optional(normal_expression), normal_hla_ligands=_resolve_optional(normal_hla_ligands), reference_proteome=_resolve_optional(reference_proteome), normal_junctions=_resolve_optional(normal_junctions))
    immune_paths = build_immune_escape_evidence(sample_id=sample_id, raw_peptides=peptides_path, outdir=immune_dir, vep_tsv=_resolve_optional(vep_appm), expression_tsv=_resolve_optional(expression), cnv_tsv=_resolve_optional(cnv), hla_loh_tsv=_resolve_optional(hla_loh), profile=profile, appm_gene_status=appm_dir / "appm_gene_status.tsv", appm_pathway_status=appm_dir / "appm_pathway_status.tsv", ccf_tsv=ccf_path)

    ranked_events = scoring / "ranked_events.tsv"
    ranked_peptides = scoring / "ranked_peptides.tsv"
    evs, peps = score(
        events_path,
        peptides_path,
        pres_path,
        appm_dir / "appm_summary.tsv",
        ccf_path,
        normal_expression,
        normal_hla_ligands,
        profile,
        ranked_events,
        ranked_peptides,
        peptide_safety_tsv=peptide_safety_path,
        peptide_escape_flags_tsv=immune_paths["peptide_escape_flags"],
        appm_peptide_modifiers_tsv=appm_dir / "appm_peptide_modifiers.tsv",
    )
    from ..validation import validation_plan_fieldnames
    val_rows = make_validation_plan(peps, outdir=outdir)
    val_path = scoring / "validation_plan.tsv"
    write_tsv(val_path, val_rows, validation_plan_fieldnames())
    prov_payload = {
        "created_at": datetime.datetime.now(datetime.UTC).isoformat().replace("+00:00", "Z"),
        "sample_id": sample_id,
        "profile": profile["_profile_name"],
        "raw_events": str(events_path),
        "raw_peptides": str(peptides_path),
        "binding_stub": binding_stub,
        "tools": provenance_registry.to_json(),
        "warning": "SV Phase 1 scoring prototype; validate SVs and binding experimentally.",
    }
    report_bundle = load_report_bundle(
        profile=profile,
        events=evs,
        peptides=peps,
        appm_summary=appm_summary,
        validation_rows=val_rows,
        outdir=outdir,
        provenance=prov_payload,
        sample_id=sample_id,
    )
    report_paths = make_dual_reports(reports, report_bundle)
    report_path = report_paths["evidence_report"]
    prov = outdir / "provenance.sv_score.json"
    write_json(prov, prov_payload)
    return {
        "raw_events": str(events_path),
        "raw_peptides": str(peptides_path),
        "presentation_evidence": str(pres_path),
        "expression_evidence": evidence_paths["expression_evidence"],
        "rna_junction_evidence": evidence_paths["rna_junction_evidence"],
        "safety_evidence": evidence_paths["safety_evidence"],
        "peptide_safety": str(peptide_safety_path),
        "event_safety": str(event_safety_path),
        "immune_escape_summary": immune_paths["immune_escape_summary"],
        "peptide_escape_flags": immune_paths["peptide_escape_flags"],
        "appm_summary": str(appm_dir / "appm_summary.tsv"),
        "ccf_2": str(ccf_path),
        "ccf_lite": str(ccf_lite_alias),
        "ranked_events": str(ranked_events),
        "ranked_peptides": str(ranked_peptides),
        "validation_plan": str(val_path),
        "evidence_report": str(report_path),
        "evidence_report_patient": report_paths["evidence_report_patient"],
        "evidence_report_technical": report_paths["evidence_report_technical"],
        "provenance": str(prov),
    }
