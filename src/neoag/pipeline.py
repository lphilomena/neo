from __future__ import annotations
from pathlib import Path
import datetime
from .config import load_profile
from .adapters.pvactools_parser import parse_pvactools_outputs
from .adapters.netmhcpan import parse_netmhcpan, write_netmhcpan_evidence
from .adapters.mhcflurry import parse_mhcflurry, write_mhcflurry_evidence
from .adapters.netmhcstabpan import parse_netmhcstabpan, write_netmhcstabpan_evidence
from .immunogenicity_composite import apply_immunogenicity_evidence, run_immunogenicity_predictors
from .presentation import build_presentation_evidence
from .appm_lite import build_appm_lite
from .ccf_v2 import build_ccf_2
from .scoring import score
from .validation import make_validation_plan
from .reports_dual import load_report_bundle, make_dual_reports
from .evidence_provenance import ProvenanceRegistry
from .evidence_layer import build_standard_evidence_layer
from .peptide_safety_gate import build_peptide_safety_gate
from .immune_escape import build_immune_escape_evidence
from .schemas import PRESENTATION_FIELDS
from .comprehensive_evidence import build_comprehensive_peptide_evidence
from .evidence_consensus import build_evidence_consensus, load_consensus_rules
from .tools.registry import RunContext
from .utils import copy_if_different, read_tsv, write_tsv, write_json
import shutil


def run(
    outdir,
    profile_name_or_path,
    sample_id,
    pvac_paths=None,
    netmhcpan=None,
    mhcflurry=None,
    netmhcstabpan=None,
    vep_appm=None,
    expression=None,
    hla_loh=None,
    purity=None,
    cnv=None,
    normal_expression=None,
    normal_hla_ligands=None,
    immunogenicity_stub=False,
    tool_executables=None,
    *,
    transcript_expression=None,
    rna_vaf=None,
    raw_events=None,
    raw_peptides=None,
    rna_junction=None,
    entry_mode=None,
    reference_proteome=None,
    normal_junctions=None,
    cancer_gene_list=None,
):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    profile = load_profile(profile_name_or_path)
    parsed = outdir / "parsed"
    parsed.mkdir(exist_ok=True)
    pres = outdir / "presentation"
    pres.mkdir(exist_ok=True)
    appm = outdir / "appm"
    appm.mkdir(exist_ok=True)
    clon = outdir / "clonality"
    clon.mkdir(exist_ok=True)
    scoring = outdir / "scoring"
    scoring.mkdir(exist_ok=True)
    safety_dir = outdir / "safety"
    immune_dir = outdir / "immune_escape"
    safety_dir.mkdir(exist_ok=True)
    immune_dir.mkdir(exist_ok=True)
    reports = outdir / "reports"
    reports.mkdir(exist_ok=True)

    raw_events_path = parsed / "raw_events.tsv"
    raw_peptides_path = parsed / "raw_peptides.tsv"

    if raw_events and raw_peptides:
        copy_if_different(raw_events, raw_events_path)
        copy_if_different(raw_peptides, raw_peptides_path)
    elif pvac_paths:
        parse_pvactools_outputs(pvac_paths, sample_id, profile["_profile_name"], raw_events_path, raw_peptides_path)
    else:
        raise ValueError("run requires pvac_paths or pre-built raw_events + raw_peptides")

    provenance_registry = ProvenanceRegistry()
    if vep_appm:
        provenance_registry.register_passthrough("vep", vep_appm)
    if hla_loh:
        hla_loh_tool = "hla_loh"
        try:
            rows = read_tsv(hla_loh)
            if rows:
                hla_loh_tool = first(rows[0], ["evidence_tool"], hla_loh_tool) or hla_loh_tool
        except Exception:
            pass
        path_l = str(hla_loh).lower()
        if hla_loh_tool == "hla_loh":
            if "spechla" in path_l:
                hla_loh_tool = "spechla"
            elif "lohhla" in path_l:
                hla_loh_tool = "lohhla"
        provenance_registry.register_converted(hla_loh_tool, hla_loh)
    if purity:
        provenance_registry.register_converted("facets", purity)
    elif cnv:
        provenance_registry.register_converted("facets", cnv)
    if cancer_gene_list:
        provenance_registry.register_passthrough("cancer_gene_list", cancer_gene_list)

    net_path = None
    if netmhcpan:
        net_path = pres / "netmhcpan_evidence.tsv"
        net_rec = provenance_registry.register_passthrough("netmhcpan", netmhcpan)
        write_netmhcpan_evidence(net_path, parse_netmhcpan(netmhcpan, sample_id), net_rec)
    mhc_path = None
    if mhcflurry:
        mhc_path = pres / "mhcflurry_evidence.tsv"
        mhc_rec = provenance_registry.register_passthrough("mhcflurry", mhcflurry)
        write_mhcflurry_evidence(mhc_path, parse_mhcflurry(mhcflurry, sample_id), mhc_rec)
    stab_path = None
    if netmhcstabpan:
        stab_path = pres / "netmhcstabpan_evidence.tsv"
        stab_rec = provenance_registry.register_passthrough("netmhcstabpan", netmhcstabpan)
        write_netmhcstabpan_evidence(stab_path, parse_netmhcstabpan(netmhcstabpan, sample_id), stab_rec)

    pres_path = pres / "presentation_evidence.tsv"
    build_presentation_evidence(
        raw_peptides_path, net_path, mhc_path, profile, pres_path, stab_path,
        provenance_registry=provenance_registry,
    )
    pres_rows = read_tsv(pres_path)
    immuno_ctx = RunContext(
        sample_id=sample_id,
        outdir=outdir,
        stub=immunogenicity_stub,
        raw_peptides=raw_peptides_path,
        executables=tool_executables or {},
    )
    immuno_paths = run_immunogenicity_predictors(
        raw_peptides_path, outdir, profile, immuno_ctx, provenance_registry=provenance_registry,
    )
    apply_immunogenicity_evidence(pres_rows, immuno_paths, profile)
    write_tsv(pres_path, pres_rows, PRESENTATION_FIELDS)

    appm_rows, appm_summary = build_appm_lite(sample_id, vep_appm, expression, hla_loh, profile, appm, cnv_tsv=cnv, raw_peptides=raw_peptides_path)
    ccf_path = clon / "ccf_2.tsv"
    ccf_lite_alias = clon / "ccf_lite.tsv"
    build_ccf_2(raw_events_path, purity, cnv, profile, ccf_path)
    shutil.copy2(ccf_path, ccf_lite_alias)

    evidence_paths = build_standard_evidence_layer(
        outdir,
        profile,
        raw_events=raw_events_path,
        raw_peptides=raw_peptides_path,
        expression=expression,
        transcript_expression=transcript_expression,
        rna_junction=rna_junction,
        rna_vaf=rna_vaf,
        normal_expression=normal_expression,
        normal_hla_ligands=normal_hla_ligands,
        sample_id=sample_id,
    )
    peptide_safety_path = safety_dir / "peptide_safety.tsv"
    event_safety_path = safety_dir / "event_safety.tsv"
    build_peptide_safety_gate(raw_events=raw_events_path, raw_peptides=raw_peptides_path, out_peptide_safety=peptide_safety_path, out_event_safety=event_safety_path, profile=profile, normal_expression=normal_expression, normal_hla_ligands=normal_hla_ligands, reference_proteome=reference_proteome, normal_junctions=normal_junctions)
    immune_paths = build_immune_escape_evidence(
        sample_id=sample_id,
        raw_peptides=raw_peptides_path,
        outdir=immune_dir,
        vep_tsv=vep_appm,
        expression_tsv=expression,
        cnv_tsv=cnv,
        hla_loh_tsv=hla_loh,
        profile=profile,
        appm_gene_status=appm / "appm_gene_status.tsv",
        appm_pathway_status=appm / "appm_pathway_status.tsv",
        ccf_tsv=ccf_path,
    )

    ranked_events = scoring / "ranked_events.tsv"
    ranked_peptides = scoring / "ranked_peptides.tsv"
    evs, peps = score(
        raw_events_path,
        raw_peptides_path,
        pres_path,
        appm / "appm_summary.tsv",
        ccf_path,
        normal_expression,
        normal_hla_ligands,
        profile,
        ranked_events,
        ranked_peptides,
        peptide_safety_tsv=peptide_safety_path,
        event_safety_tsv=event_safety_path,
        peptide_escape_flags_tsv=immune_paths["peptide_escape_flags"],
        appm_peptide_modifiers_tsv=appm / "appm_peptide_modifiers.tsv",
        cancer_gene_list_tsv=cancer_gene_list,
    )
    from .validation import validation_plan_fieldnames
    val_rows = make_validation_plan(peps, outdir=outdir)
    val_path = scoring / "validation_plan.tsv"
    write_tsv(val_path, val_rows, validation_plan_fieldnames())
    comprehensive_path = scoring / "comprehensive_peptide_evidence.tsv"
    comprehensive_summary = build_comprehensive_peptide_evidence(
        output_tsv=comprehensive_path,
        ranked_peptides=ranked_peptides,
        raw_peptides=raw_peptides_path,
        raw_events=raw_events_path,
        presentation_evidence=pres_path,
        appm_peptide_modifiers=appm / "appm_peptide_modifiers.tsv",
        ccf_2=ccf_path,
        expression_evidence=evidence_paths["expression_evidence"],
        rna_junction_evidence=evidence_paths["rna_junction_evidence"],
        peptide_safety=peptide_safety_path,
        event_safety=event_safety_path,
        peptide_escape_flags=immune_paths["peptide_escape_flags"],
        validation_plan=val_path,
        conflicts_tsv=scoring / "evidence_source_conflicts.tsv",
    )
    consensus_rules_path = Path(__file__).resolve().parents[2] / "configs/ranking/sarcoma_evidence_consensus_v1.toml"
    consensus_rules = load_consensus_rules(consensus_rules_path if consensus_rules_path.is_file() else None)
    evidence_consensus_summary = build_evidence_consensus(
        comprehensive_path,
        scoring,
        consensus_rules,
        weighted_baseline_tsv=ranked_peptides,
    )
    evidence_consensus_path = Path(evidence_consensus_summary["ranked_peptides"])
    prov_payload = {
        "created_at": datetime.datetime.now(datetime.UTC).isoformat().replace("+00:00", "Z"),
        "sample_id": sample_id,
        "profile": profile["_profile_name"],
        "entry_mode": entry_mode,
        "tools": provenance_registry.to_json(),
        "warning": "Computational prototype only.",
        "parallel_rankings": {
            "legacy_weighted": str(ranked_peptides),
            "evidence_consensus": str(evidence_consensus_path),
            "evidence_states": evidence_consensus_summary["evidence_states"],
            "evidence_conflicts": evidence_consensus_summary["evidence_conflicts"],
            "evidence_source_conflicts": comprehensive_summary["conflicts_tsv"],
            "event_consensus": evidence_consensus_summary["ranked_events"],
            "comparison": evidence_consensus_summary["comparison"],
            "comparison_markdown": evidence_consensus_summary["comparison_markdown"],
            "summary": evidence_consensus_summary["summary"],
            "run_manifest": evidence_consensus_summary["run_manifest"],
            "weighted_baseline": evidence_consensus_summary["weighted_baseline"],
            "all_tool_results": evidence_consensus_summary["all_tool_results"],
            "rules": str(consensus_rules_path) if consensus_rules_path.is_file() else "embedded_default",
            "rules_name": consensus_rules.get("metadata", {}).get("name", ""),
            "rules_version": consensus_rules.get("metadata", {}).get("version", ""),
            "rules_status": consensus_rules.get("metadata", {}).get("status", "PROVISIONAL_RESEARCH_ONLY"),
            "legacy_ranking_modified": evidence_consensus_summary["legacy_ranking_modified"],
        },
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
        entry_mode=entry_mode or "",
    )
    report_paths = make_dual_reports(reports, report_bundle)
    report_path = report_paths["evidence_report"]
    prov = outdir / "provenance.json"
    write_json(prov, prov_payload)
    return {
        "raw_events": str(raw_events_path),
        "raw_peptides": str(raw_peptides_path),
        "presentation_evidence": str(pres_path),
        "expression_evidence": evidence_paths["expression_evidence"],
        "rna_junction_evidence": evidence_paths["rna_junction_evidence"],
        "safety_evidence": evidence_paths["safety_evidence"],
        "peptide_safety": str(peptide_safety_path),
        "event_safety": str(event_safety_path),
        "immune_escape_summary": immune_paths["immune_escape_summary"],
        "peptide_escape_flags": immune_paths["peptide_escape_flags"],
        "appm_lite": str(appm / "appm_lite.tsv"),
        "appm_summary": str(appm / "appm_summary.tsv"),
        "appm_gene_status": str(appm / "appm_gene_status.tsv"),
        "appm_pathway_status": str(appm / "appm_pathway_status.tsv"),
        "appm_module_scores": str(appm / "appm_module_scores.tsv"),
        "appm_submodule_scores": str(appm / "appm_submodule_scores.tsv"),
        "appm_conflicts": str(appm / "appm_conflicts.tsv"),
        "appm_input_status": str(appm / "appm_input_status.tsv"),
        "appm_peptide_modifiers": str(appm / "appm_peptide_modifiers.tsv"),
        "ccf_2": str(ccf_path),
        "ccf_lite": str(ccf_lite_alias),
        "ranked_events": str(ranked_events),
        "ranked_peptides": str(ranked_peptides),
        "ranked_peptides_evidence_consensus": str(evidence_consensus_path),
        "ranked_events_evidence_consensus": evidence_consensus_summary["ranked_events"],
        "evidence_states": evidence_consensus_summary["evidence_states"],
        "evidence_conflicts": evidence_consensus_summary["evidence_conflicts"],
        "evidence_source_conflicts": comprehensive_summary["conflicts_tsv"],
        "weighted_vs_consensus_comparison": evidence_consensus_summary["comparison"],
        "ranking_compare_weighted_vs_consensus": evidence_consensus_summary["comparison"],
        "ranking_compare_weighted_vs_consensus_md": evidence_consensus_summary["comparison_markdown"],
        "evidence_consensus_summary": evidence_consensus_summary["summary"],
        "evidence_consensus_run": evidence_consensus_summary["run_manifest"],
        "ranked_peptides_weighted_baseline": evidence_consensus_summary["weighted_baseline"],
        "all_tool_results": evidence_consensus_summary["all_tool_results"],
        "comprehensive_peptide_evidence": str(comprehensive_path),
        "validation_plan": str(val_path),
        "evidence_report": str(report_path),
        "evidence_report_patient": report_paths["evidence_report_patient"],
        "evidence_report_technical": report_paths["evidence_report_technical"],
        "provenance": str(prov),
    }
