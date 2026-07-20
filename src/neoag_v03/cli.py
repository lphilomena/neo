from __future__ import annotations
import argparse
import os
from pathlib import Path
from .pipeline_v03 import run_v03
from .config import load_profile
from .adapters.pvactools_parser import parse_pvactools_outputs
from .adapters.netmhcpan import parse_netmhcpan, write_netmhcpan_evidence
from .adapters.mhcflurry import parse_mhcflurry, write_mhcflurry_evidence
from .presentation import build_presentation_evidence
from .appm_lite import build_appm_lite
from .appm_v2 import build_appm_2
from .ccf_lite import build_ccf_lite
from .ccf_v2 import build_ccf_2
from .scoring_v03 import score_v03
from .validation import make_validation_plan_v03
from .reports_dual import load_report_bundle, make_patient_report, make_technical_report
from .reports_v041 import make_report_v041
from .utils import read_tsv, write_tsv
from .tools import check_all_tools, run_tool, load_run_config, run_upstream
from .tools.postprocess import (
    lohhla_to_hla_loh_tsv,
    spechla_to_hla_loh_tsv,
    facets_to_purity_tsv,
    facets_to_cnv_tsv,
    ascat_to_purity_tsv,
    ascat_to_cnv_tsv,
)
from .tools.registry import RunContext
from .tools.runner import RUNNERS
from .benchmark_improve import DEFAULT_CEDAR, DEFAULT_INHOUSE, run_benchmark_improve
from .benchmark_system import run_system_benchmark
from .adapters.peptide_input import convert_peptide_input
from .peptide_predict import run_peptide_predict
from .sv.phase1 import build_sv_phase1_raw
from .sv.wes_adapter import build_sv_wes_phase1_5_raw
from .sv.score_pipeline import run_sv_v03_score
from .sv.peptide_builder import read_hla_alleles
from .input_router import build_raw_intermediates, resolve_entry_mode
from .evidence_layer import build_standard_evidence_layer
from .snv_call.pipeline import run_snv_wes_call, run_snv_wes_full
from .vep.annotate import (
    load_vep_annotate_config,
    resolve_vep_annotate_from_config,
    run_vep_pvacseq_annotate,
)
from .vep.extract_peptides import extract_variant_peptides_from_vcf, parse_peptide_lengths
from .peptide_safety_gate import build_peptide_safety_gate
from .immune_escape import build_immune_escape_evidence
from .hla_loh_crosscheck import write_hla_loh_crosscheck
from .comprehensive_evidence import build_comprehensive_peptide_evidence

ROOT = Path(__file__).resolve().parents[2]
def fixture(x): return ROOT/"data"/"fixtures"/x
def resource(x): return ROOT/"resources"/x

def cmd_run_demo(args):
    out = run_v03(
        outdir=args.outdir, profile_name_or_path=args.profile, sample_id=args.sample_id,
        pvac_paths=[fixture("pvacseq_aggregated.tsv"), fixture("pvacfuse_aggregated.tsv")],
        netmhcpan=fixture("netmhcpan_example.xls"), mhcflurry=fixture("mhcflurry_predictions.csv"),
        vep_appm=fixture("vep_appm.tsv"), expression=fixture("gene_expression.tsv"),
        hla_loh=fixture("hla_loh.tsv"), purity=fixture("purity.tsv"), cnv=fixture("cnv_segments.tsv"),
        normal_expression=resource("normal_expression.example.tsv"), normal_hla_ligands=resource("normal_hla_ligands.example.tsv"),
        immunogenicity_stub=True,
    )
    print("NeoAg v0.4.3 demo completed. Outputs retain .v03.tsv names for schema compatibility.")
    for k,v in out.items(): print(f"  {k}: {v}")

def cmd_run_v03(args):
    kwargs = {}
    if args.raw_events and args.raw_peptides:
        kwargs["raw_events"] = args.raw_events
        kwargs["raw_peptides"] = args.raw_peptides
    elif not args.pvac:
        raise SystemExit("Provide --pvac or both --raw-events and --raw-peptides")
    out = run_v03(
        args.outdir, args.profile, args.sample_id, args.pvac or [],
        args.netmhcpan, args.mhcflurry, args.netmhcstabpan,
        args.vep_appm, args.expression, args.hla_loh, args.purity, args.cnv,
        args.normal_expression, args.normal_hla_ligands,
        immunogenicity_stub=args.immunogenicity_stub,
        rna_junction=args.rna_junction,
        entry_mode=args.entry_mode,
        **kwargs,
    )
    print("NeoAg run completed. Outputs retain .v03.tsv names for schema compatibility.")
    for k,v in out.items(): print(f"  {k}: {v}")

def cmd_parse_pvac(args):
    parse_pvactools_outputs(args.pvac, args.sample_id, load_profile(args.profile)["_profile_name"], args.events_out, args.peptides_out)
    print("Parsed pVACtools outputs.")

def cmd_parse_netmhcpan(args):
    write_netmhcpan_evidence(args.out, parse_netmhcpan(args.input, args.sample_id))
    print("Parsed NetMHCpan output.")

def cmd_parse_mhcflurry(args):
    write_mhcflurry_evidence(args.out, parse_mhcflurry(args.input, args.sample_id))
    print("Parsed MHCflurry output.")

def cmd_build_presentation(args):
    build_presentation_evidence(
        args.raw_peptides, args.netmhcpan, args.mhcflurry,
        load_profile(args.profile), args.out, args.netmhcstabpan,
    )
    print("Built PresentationEvidence.")

def cmd_appm(args):
    build_appm_lite(args.sample_id, args.vep_tsv, args.expression, args.hla_loh, load_profile(args.profile), args.outdir, cnv_tsv=getattr(args, "cnv", None), raw_peptides=getattr(args, "raw_peptides", None))
    print("Built APPM 2.0 sidecars via compatibility entry point.")

def cmd_appm_2(args):
    out = build_appm_2(
        sample_id=args.sample_id,
        vep_tsv=args.vep_tsv,
        expression_tsv=args.expression,
        hla_loh_tsv=args.hla_loh,
        cnv_tsv=args.cnv,
        raw_peptides=args.raw_peptides,
        tumor_purity_tsv=getattr(args, "tumor_purity", None),
        profile=load_profile(args.profile),
        outdir=args.outdir,
    )
    print("Built APPM 2.0.")
    for k, v in out.items(): print(f"  {k}: {v}")

def cmd_ccf(args):
    build_ccf_lite(args.events, args.purity, args.cnv, load_profile(args.profile), args.out)
    print("Built CCF 2.0 fields via compatibility entry point.")

def cmd_ccf_2(args):
    rows = build_ccf_2(
        args.events,
        args.purity,
        args.cnv,
        load_profile(args.profile),
        args.out,
        external_clonality_tsv=getattr(args, "external_clonality", None),
        svclone_tsv=getattr(args, "svclone", None),
        sidecar_dir=getattr(args, "sidecar_dir", None),
        input_qc_out=getattr(args, "input_qc_out", None),
        conflicts_out=getattr(args, "conflicts_out", None),
        clusters_out=getattr(args, "clusters_out", None),
    )
    print(f"Built CCF 2.0: {args.out} ({len(rows)} events)")

def cmd_score(args):
    from .scoring_v03 import resolve_appm_peptide_modifiers_tsv

    appm_modifiers = resolve_appm_peptide_modifiers_tsv(
        getattr(args, "appm_peptide_modifiers", None),
        appm_summary_tsv=getattr(args, "appm_summary", None),
    )
    score_v03(
        args.raw_events, args.raw_peptides, args.presentation, args.appm_summary, args.ccf,
        args.normal_expression, args.normal_hla_ligands, load_profile(args.profile),
        args.out_events, args.out_peptides,
        peptide_safety_tsv=getattr(args, "peptide_safety", None),
        peptide_escape_flags_tsv=getattr(args, "peptide_escape_flags", None),
        appm_peptide_modifiers_tsv=appm_modifiers,
    )
    print("Scored current evidence model into v03 schema-compatible output files.")



def cmd_peptide_safety(args):
    rows, _ = build_peptide_safety_gate(
        raw_events=args.raw_events,
        raw_peptides=args.raw_peptides,
        out_peptide_safety=args.out,
        out_event_safety=args.event_out,
        profile=load_profile(args.profile),
        normal_expression=args.normal_expression,
        normal_hla_ligands=args.normal_hla_ligands,
        reference_proteome=args.reference_proteome,
        normal_junctions=args.normal_junctions,
    )
    print(f"Built peptide safety evidence: {args.out} ({len(rows)} peptides)")


def cmd_immune_escape(args):
    out = build_immune_escape_evidence(
        sample_id=args.sample_id,
        raw_peptides=args.raw_peptides,
        outdir=args.outdir,
        vep_tsv=args.vep_tsv,
        cnv_tsv=args.cnv,
        expression_tsv=args.expression,
        hla_loh_tsv=args.hla_loh,
        profile=load_profile(args.profile),
        appm_gene_status=getattr(args, "appm_gene_status", None),
        appm_pathway_status=getattr(args, "appm_pathway_status", None),
        ccf_tsv=getattr(args, "ccf", None),
        therapy_context=getattr(args, "therapy_context", None),
        ranked_peptides=getattr(args, "ranked_peptides", None),
        top_priority_threshold=getattr(args, "top_priority_threshold", "B_CAUTION"),
    )
    print("Built immune escape evidence.")
    for k, v in out.items():
        print(f"  {k}: {v}")


def cmd_report_v041(args):
    profile = load_profile(args.profile)
    events = read_tsv(args.ranked_events)
    peptides = read_tsv(args.ranked_peptides)
    appm_rows = read_tsv(args.appm_summary) if args.appm_summary else []
    appm_summary = appm_rows[0] if appm_rows else {}
    validation_rows = read_tsv(args.validation_plan) if args.validation_plan else []
    make_report_v041(
        args.out,
        profile,
        events,
        peptides,
        appm_summary=appm_summary,
        validation_rows=validation_rows,
        appm_gene_status=read_tsv(args.appm_gene_status) if args.appm_gene_status else [],
        appm_module_scores=read_tsv(args.appm_module_scores) if args.appm_module_scores else [],
        appm_submodule_scores=read_tsv(args.appm_submodule_scores) if args.appm_submodule_scores else [],
        appm_conflicts=read_tsv(args.appm_conflicts) if args.appm_conflicts else [],
        appm_peptide_modifiers=read_tsv(args.appm_peptide_modifiers) if args.appm_peptide_modifiers else [],
        immune_escape_summary=read_tsv(args.immune_escape_summary) if args.immune_escape_summary else [],
        peptide_escape_flags=read_tsv(args.peptide_escape_flags) if args.peptide_escape_flags else [],
        peptide_safety=read_tsv(args.peptide_safety) if args.peptide_safety else [],
        ccf=read_tsv(args.ccf) if args.ccf else [],
    )
    print(f"Wrote v0.4.1 evidence report to {args.out}")


def cmd_benchmark_system(args):
    out = run_system_benchmark(
        outdir=args.outdir,
        profile_name=args.profile,
        ligandome_ms=getattr(args, "ligandome_ms", None),
        mode=args.mode,
        ranked_peptides=getattr(args, "ranked_peptides", None),
        appm_summary=getattr(args, "appm_summary", None),
        appm_module_scores=getattr(args, "appm_module_scores", None),
        appm_submodule_scores=getattr(args, "appm_submodule_scores", None),
        peptide_appm_flags=getattr(args, "peptide_appm_flags", None),
        peptide_escape_flags=getattr(args, "peptide_escape_flags", None),
    )
    print("System benchmark completed.")
    for k, v in out.items():
        print(f"  {k}: {v}")


def cmd_validation_plan(args):
    peptides = read_tsv(args.ranked_peptides)
    rows = make_validation_plan_v03(
        peptides,
        peptide_catalog_tsv=getattr(args, "variant_peptides", None),
        outdir=getattr(args, "outdir", None),
    )
    from .validation import validation_plan_fieldnames
    write_tsv(args.out, rows, validation_plan_fieldnames())
    print(f"Wrote validation plan to {args.out}")

def cmd_report(args):
    profile = load_profile(args.profile)
    events = read_tsv(args.ranked_events)
    peptides = read_tsv(args.ranked_peptides)
    appm_rows = read_tsv(args.appm_summary) if args.appm_summary else []
    appm_summary = appm_rows[0] if appm_rows else {"mhc_i_integrity_score": "1.0", "mhc_ii_integrity_score": "1.0", "appm_overall_status": "UNASSESSED"}
    validation_rows = read_tsv(args.validation_plan) if args.validation_plan else []
    bundle = load_report_bundle(
        profile=profile,
        events=events,
        peptides=peptides,
        appm_summary=appm_summary,
        validation_rows=validation_rows,
        outdir=getattr(args, "outdir", None),
        provenance=_read_json_optional(getattr(args, "provenance", None)) if getattr(args, "provenance", None) else None,
        sample_id=getattr(args, "sample_id", "") or "",
    )
    audience = getattr(args, "audience", "both")
    out = Path(args.out)
    if audience in {"both", "technical"}:
        tech_out = out if audience == "technical" else out.parent / "evidence_report.technical.html"
        if audience == "both" and out.name == "evidence_report.v03.html":
            tech_out = out.parent / "evidence_report.technical.html"
        make_technical_report(tech_out, bundle)
        if audience == "both":
            out.write_text(tech_out.read_text(encoding="utf-8"), encoding="utf-8")
    if audience in {"both", "patient"}:
        patient_out = out.parent / "evidence_report.patient.html" if audience == "both" else out
        make_patient_report(patient_out, bundle)
    if audience == "both":
        print(f"Wrote patient report to {out.parent / 'evidence_report.patient.html'}")
        print(f"Wrote technical report to {out.parent / 'evidence_report.technical.html'}")
        print(f"Wrote legacy alias to {out}")
    else:
        print(f"Wrote {audience} report to {out}")


def _read_json_optional(path):
    from .reports_dual import _read_json_optional as _load
    return _load(path)


def cmd_extract_variant_peptides(args):
    normal_proteome_fasta = args.normal_proteome_fasta or os.environ.get("NEOAG_NORMAL_PROTEOME_FASTA")
    if args.filter_normal_proteome and not normal_proteome_fasta:
        raise SystemExit("--filter-normal-proteome requires --normal-proteome-fasta or NEOAG_NORMAL_PROTEOME_FASTA")
    lengths = parse_peptide_lengths(
        args.lengths,
        length_min=args.length_min,
        length_max=args.length_max,
    )
    filter_normal = bool(args.filter_normal_proteome)
    if normal_proteome_fasta and not args.annotate_normal_proteome_only:
        filter_normal = True
    hla_alleles = []
    if args.hla_alleles:
        hla_alleles = [a.strip() for a in args.hla_alleles.split(",") if a.strip()]
    result = extract_variant_peptides_from_vcf(
        args.input_vcf,
        args.output,
        lengths=lengths,
        pass_only=not args.include_filtered,
        sample_id=args.sample_id,
        exclude_multi_aa=args.exclude_multi_aa,
        single_aa_only=args.single_aa_only,
        mini_len=args.mini_len,
        normal_proteome_fasta=normal_proteome_fasta,
        filter_normal_proteome=filter_normal,
        hla_alleles=hla_alleles or None,
        netmhcpan_xls=args.netmhcpan_xls,
        mhcflurry_csv=args.mhcflurry_csv,
        annotate_netmhcpan=args.annotate_netmhcpan,
        tumor_sample_name=args.tumor_sample_name,
    )
    print("Variant short peptides extracted (sliding window, no MHC binding).")
    print(f"  Method: {result['generation_method']}")
    for key, val in result.items():
        print(f"  {key}: {val}")


def cmd_vep_annotate(args):
    if args.config:
        cfg = load_vep_annotate_config(args.config)
        kwargs = resolve_vep_annotate_from_config(cfg, root=ROOT)
    else:
        if not args.input_vcf or not args.output_vcf or not args.reference_fasta:
            raise SystemExit(
                "vep-annotate requires --input-vcf, --output-vcf, --fasta "
                "or --config conf/vep.annotate.example.toml"
            )
        kwargs = {
            "input_vcf": args.input_vcf,
            "output_vcf": args.output_vcf,
            "reference_fasta": args.reference_fasta,
            "workdir": args.workdir,
            "cache_dir": args.cache_dir,
            "plugins_dir": args.plugins_dir,
            "online": True if args.online else None,
            "fork": args.fork,
            "pick": not args.no_pick,
            "expression_custom": args.expression_custom,
            "index_vcf": not args.no_index,
        }
    result = run_vep_pvacseq_annotate(**kwargs)
    print("VEP annotation completed (pVACseq-compatible).")
    for key, val in result.items():
        print(f"  {key}: {val}")
    print("")
    print("Next: set in conf/run.*.toml:")
    print(f'  tumor_vcf = "{result["annotated_vcf"]}"')
    print(f'  variants_vcf = "{result["annotated_vcf"]}"')


def cmd_check_tools(_args):
    tools_env = ROOT / "conf" / "tools.env.sh"
    if tools_env.is_file() and not __import__("os").environ.get("NEOAG_CONDA_ENV"):
        print(f"Tip: source {tools_env}  # activates neoag-tools conda env on PATH")
        print()
    print(f"{'Tool':<14} {'Executable':<22} {'Status':<10} Message")
    print("-" * 72)
    for st in check_all_tools():
        status = "OK" if st.available else "MISSING"
        print(f"{st.name:<14} {st.executable:<22} {status:<10} {st.message}")
    print()
    print(f"Setup guide: {ROOT / 'docs' / 'TOOLS_SETUP.md'}")


def cmd_run_tool(args):
    cfg = load_run_config(args.config) if args.config else {}
    tools = cfg.get("tools", {}) if cfg else {}
    inputs = cfg.get("inputs", {}) if cfg else {}
    ctx = RunContext(
        sample_id=args.sample_id,
        outdir=Path(args.workdir),
        stub=args.stub or bool(tools.get("stub")),
        executables={k: str(v) for k, v in (tools.get("executables") or {}).items()},
        hla_alleles=list(inputs.get("hla_alleles") or []),
        raw_peptides=Path(args.raw_peptides) if args.raw_peptides else None,
        tumor_vcf=Path(args.tumor_vcf) if args.tumor_vcf else None,
        normal_vcf=Path(args.normal_vcf) if args.normal_vcf else None,
        fusion_tsv=Path(args.fusion_tsv) if args.fusion_tsv else None,
        variants_vcf=Path(args.variants_vcf) if args.variants_vcf else None,
    )
    out = run_tool(args.tool, ctx, args.output)
    print(f"Wrote {args.tool} output to {out}")


def cmd_run_upstream(args):
    outs = run_upstream(args.config, args.outdir)
    print("Upstream tools completed.")
    for k, v in sorted(outs.items()):
        print(f"  {k}: {v}")


def cmd_convert_peptide_input(args):
    summary = convert_peptide_input(args.input, args.outdir, sample_id=args.sample_id)
    print("Converted peptide input to strict peptide–HLA pair tables.")
    print(f"  Input rows: {summary.input_rows}")
    print(f"  Unique pairs: {summary.pair_rows} ({summary.unique_peptides} peptides × {summary.unique_hla} HLA)")
    print(f"  Skipped rows: {summary.skipped_rows}")
    print(f"  raw_peptides.tsv: {summary.raw_peptides_tsv}")
    print(f"  peptide_hla_pairs.tsv: {summary.pairs_tsv}")
    print(f"  hla_alleles.txt: {summary.hla_alleles_txt}")


def cmd_convert_lohhla(args):
    lohhla_to_hla_loh_tsv(args.input, args.output)
    print(f"Wrote hla_loh.tsv: {args.output}")


def cmd_convert_spechla(args):
    spechla_to_hla_loh_tsv(args.input, args.output)
    print(f"Wrote hla_loh.tsv: {args.output}")


def cmd_crosscheck_hla_loh(args):
    out = write_hla_loh_crosscheck(
        args.out,
        lohhla_hla_loh=args.lohhla_hla_loh,
        spechla_hla_loh=args.spechla_hla_loh,
        consensus_out=args.consensus_out,
        include_single_tool=not args.strict_consensus_only,
    )
    print("Built HLA LOH cross-check evidence.")
    for k, v in out.items():
        print(f"  {k}: {v}")


def cmd_convert_facets(args):
    facets_to_purity_tsv(args.purity_input, args.sample_id, args.purity_output)
    print(f"Wrote purity.tsv: {args.purity_output}")
    if args.cnv_input and args.cnv_output:
        facets_to_cnv_tsv(args.cnv_input, args.cnv_output)
        print(f"Wrote cnv_segments.tsv: {args.cnv_output}")


def cmd_convert_ascat(args):
    ascat_to_purity_tsv(args.summary_input, args.sample_id, args.purity_output)
    print(f"Wrote purity.tsv: {args.purity_output}")
    if args.segments_input and args.cnv_output:
        ascat_to_cnv_tsv(args.segments_input, args.cnv_output)
        print(f"Wrote cnv_segments.tsv: {args.cnv_output}")


def cmd_peptide_predict(args):
    outs = run_peptide_predict(
        args.input,
        args.outdir,
        sample_id=args.sample_id,
        profile_name=args.profile,
        stub=args.stub,
        skip_netmhcpan=args.skip_netmhcpan,
        skip_mhcflurry=args.skip_mhcflurry,
        skip_prime=args.skip_prime,
        skip_bigmhc_im=args.skip_bigmhc_im,
        skip_deepimmuno=args.skip_deepimmuno,
        skip_stabpan=args.skip_stabpan,
    )
    summary = outs["summary"]
    print("Peptide–HLA prediction completed.")
    print(f"  Unique pairs: {summary.pair_rows} ({summary.unique_peptides} peptides, {summary.unique_hla} HLA)")
    for key, val in outs.items():
        if key != "summary":
            print(f"  {key}: {val}")


def cmd_sv_build_raw(args):
    callers = args.callers or []
    if callers and len(callers) != len(args.sv_vcf):
        raise SystemExit("--callers length must match --sv-vcf length when provided")
    out = build_sv_phase1_raw(
        sample_id=args.sample_id,
        sv_vcfs=args.sv_vcf,
        callers=callers,
        reference_fasta=args.reference_fasta,
        gencode_gtf=args.gencode_gtf,
        hla=args.hla,
        outdir=args.outdir,
        profile_name=load_profile(args.profile)["_profile_name"],
        tumor_sample_name=args.tumor_sample_name,
        normal_sample_name=args.normal_sample_name,
        expression_tsv=args.expression,
        rna_junction_tsv=args.rna_junctions,
        normal_expression_tsv=args.normal_expression,
        normal_hla_ligands_tsv=args.normal_hla_ligands,
        merge_distance_bp=args.merge_distance_bp,
        allow_tier2=not args.tier1_only,
        capture_bed=getattr(args, "capture_bed", None),
    )
    print("SV Phase 1 raw inputs completed.")
    for k, v in out.items():
        print(f"  {k}: {v}")


def _sv_optional_path(val):
    if not val:
        return None
    p = Path(val)
    return p if p.is_file() else None


def cmd_sv_score_v03(args):
    raw_events = Path(args.raw_events) if args.raw_events else Path(args.sv_outdir) / "parsed" / "raw_events.tsv"
    raw_peptides = Path(args.raw_peptides) if args.raw_peptides else Path(args.sv_outdir) / "parsed" / "raw_peptides.tsv"
    if not raw_events.is_file() or not raw_peptides.is_file():
        raise SystemExit(f"Missing SV raw tables: {raw_events} / {raw_peptides}")
    out = run_sv_v03_score(
        outdir=args.outdir,
        profile_name_or_path=args.profile,
        sample_id=args.sample_id,
        raw_events=raw_events,
        raw_peptides=raw_peptides,
        netmhcpan=_sv_optional_path(args.netmhcpan),
        mhcflurry=_sv_optional_path(args.mhcflurry),
        vep_appm=_sv_optional_path(args.vep_appm) or ROOT / "assets" / "empty_vep.tsv",
        expression=_sv_optional_path(args.expression),
        hla_loh=_sv_optional_path(args.hla_loh) or ROOT / "assets" / "empty_hla_loh.tsv",
        purity=_sv_optional_path(args.purity) or ROOT / "assets" / "empty_purity.tsv",
        cnv=_sv_optional_path(args.cnv) or ROOT / "assets" / "empty_cnv.tsv",
        normal_expression=_sv_optional_path(args.normal_expression) or resource("normal_expression.example.tsv"),
        normal_hla_ligands=_sv_optional_path(args.normal_hla_ligands) or resource("normal_hla_ligands.example.tsv"),
        binding_stub=args.binding_stub,
        immunogenicity_stub=args.immunogenicity_stub,
        run_binding=not args.skip_binding,
        reference_proteome=getattr(args, "reference_proteome", None),
        normal_junctions=getattr(args, "normal_junctions", None),
    )
    print("SV Phase 1 scoring completed.")
    for k, v in out.items():
        print(f"  {k}: {v}")


def cmd_sv_run_full(args):
    callers = args.callers or []
    if callers and len(callers) != len(args.sv_vcf):
        raise SystemExit("--callers length must match --sv-vcf length when provided")
    adapter_out = Path(args.outdir) / "adapter"
    build_sv_phase1_raw(
        sample_id=args.sample_id,
        sv_vcfs=args.sv_vcf,
        callers=callers,
        reference_fasta=args.reference_fasta,
        gencode_gtf=args.gencode_gtf,
        hla=args.hla,
        outdir=adapter_out,
        profile_name=load_profile(args.profile)["_profile_name"],
        tumor_sample_name=args.tumor_sample_name,
        normal_sample_name=args.normal_sample_name,
        expression_tsv=args.expression,
        rna_junction_tsv=args.rna_junctions,
        normal_expression_tsv=args.normal_expression,
        normal_hla_ligands_tsv=args.normal_hla_ligands,
        merge_distance_bp=args.merge_distance_bp,
        allow_tier2=not args.tier1_only,
        capture_bed=getattr(args, "capture_bed", None),
        capture_near_bp=getattr(args, "capture_near_bp", 250),
        capture_slop_bp=getattr(args, "capture_slop_bp", 1000),
    )
    score_out = run_sv_v03_score(
        outdir=args.outdir,
        profile_name_or_path=args.profile,
        sample_id=args.sample_id,
        raw_events=adapter_out / "parsed" / "raw_events.tsv",
        raw_peptides=adapter_out / "parsed" / "raw_peptides.tsv",
        vep_appm=_sv_optional_path(args.vep_appm) or ROOT / "assets" / "empty_vep.tsv",
        expression=_sv_optional_path(args.expression),
        hla_loh=_sv_optional_path(args.hla_loh) or ROOT / "assets" / "empty_hla_loh.tsv",
        purity=_sv_optional_path(args.purity) or ROOT / "assets" / "empty_purity.tsv",
        cnv=_sv_optional_path(args.cnv) or ROOT / "assets" / "empty_cnv.tsv",
        normal_expression=_sv_optional_path(args.normal_expression) or resource("normal_expression.example.tsv"),
        normal_hla_ligands=_sv_optional_path(args.normal_hla_ligands) or resource("normal_hla_ligands.example.tsv"),
        binding_stub=args.binding_stub,
        immunogenicity_stub=args.immunogenicity_stub,
        run_binding=not args.skip_binding,
        reference_proteome=getattr(args, "reference_proteome", None),
        normal_junctions=getattr(args, "normal_junctions", None),
    )
    print("SV Phase 1 end-to-end completed.")
    print(f"  adapter: {adapter_out}")
    for k, v in score_out.items():
        print(f"  {k}: {v}")


def _sv_build_raw_common(args, *, wes: bool):
    callers = args.callers or []
    if callers and len(callers) != len(args.sv_vcf):
        raise SystemExit("--callers length must match --sv-vcf length when provided")
    build_fn = build_sv_wes_phase1_5_raw if wes else build_sv_phase1_raw
    return build_fn(
        sample_id=args.sample_id,
        sv_vcfs=args.sv_vcf,
        callers=callers,
        reference_fasta=args.reference_fasta,
        gencode_gtf=args.gencode_gtf,
        hla=args.hla,
        outdir=args.outdir,
        profile_name=load_profile(args.profile)["_profile_name"],
        tumor_sample_name=args.tumor_sample_name,
        normal_sample_name=args.normal_sample_name,
        expression_tsv=args.expression,
        rna_junction_tsv=args.rna_junctions,
        normal_expression_tsv=args.normal_expression,
        normal_hla_ligands_tsv=args.normal_hla_ligands,
        merge_distance_bp=args.merge_distance_bp,
        allow_tier2=not args.tier1_only,
        **({
            "capture_bed": getattr(args, "capture_bed", None),
            "capture_near_bp": getattr(args, "capture_near_bp", 250),
            "capture_slop_bp": getattr(args, "capture_slop_bp", 1000),
        } if wes else {}),
    )


def cmd_sv_build_raw_wes(args):
    out = _sv_build_raw_common(args, wes=True)
    print("WES SV Phase 1.5 raw inputs completed.")
    for k, v in out.items():
        print(f"  {k}: {v}")


def cmd_sv_run_full_wes(args):
    callers = args.callers or []
    if callers and len(callers) != len(args.sv_vcf):
        raise SystemExit("--callers length must match --sv-vcf length when provided")
    adapter_out = Path(args.outdir) / "adapter"
    profile = load_profile(args.profile)["_profile_name"]
    build_sv_wes_phase1_5_raw(
        sample_id=args.sample_id,
        sv_vcfs=args.sv_vcf,
        callers=callers,
        reference_fasta=args.reference_fasta,
        gencode_gtf=args.gencode_gtf,
        hla=args.hla,
        outdir=adapter_out,
        profile_name=profile,
        tumor_sample_name=args.tumor_sample_name,
        normal_sample_name=args.normal_sample_name,
        expression_tsv=args.expression,
        rna_junction_tsv=args.rna_junctions,
        normal_expression_tsv=args.normal_expression,
        normal_hla_ligands_tsv=args.normal_hla_ligands,
        merge_distance_bp=args.merge_distance_bp,
        allow_tier2=not args.tier1_only,
        capture_bed=getattr(args, "capture_bed", None),
        capture_near_bp=getattr(args, "capture_near_bp", 250),
        capture_slop_bp=getattr(args, "capture_slop_bp", 1000),
    )
    score_out = run_sv_v03_score(
        outdir=args.outdir,
        profile_name_or_path=profile,
        sample_id=args.sample_id,
        raw_events=adapter_out / "parsed" / "raw_events.tsv",
        raw_peptides=adapter_out / "parsed" / "raw_peptides.tsv",
        vep_appm=_sv_optional_path(args.vep_appm) or ROOT / "assets" / "empty_vep.tsv",
        expression=_sv_optional_path(args.expression),
        hla_loh=_sv_optional_path(args.hla_loh) or ROOT / "assets" / "empty_hla_loh.tsv",
        purity=_sv_optional_path(args.purity) or ROOT / "assets" / "empty_purity.tsv",
        cnv=_sv_optional_path(args.cnv) or ROOT / "assets" / "empty_cnv.tsv",
        normal_expression=_sv_optional_path(args.normal_expression) or resource("normal_expression.example.tsv"),
        normal_hla_ligands=_sv_optional_path(args.normal_hla_ligands) or resource("normal_hla_ligands.example.tsv"),
        binding_stub=args.binding_stub,
        immunogenicity_stub=args.immunogenicity_stub,
        run_binding=not args.skip_binding,
        reference_proteome=getattr(args, "reference_proteome", None),
        normal_junctions=getattr(args, "normal_junctions", None),
    )
    print("WES SV Phase 1.5 end-to-end completed.")
    print(f"  adapter: {adapter_out}")
    for k, v in score_out.items():
        print(f"  {k}: {v}")


def cmd_snv_call_wes(args):
    out = run_snv_wes_call(
        outdir=args.outdir,
        tumor_bam=args.tumor_bam,
        normal_bam=args.normal_bam,
        reference_fasta=args.reference_fasta,
        intervals_bed=args.intervals_bed,
        tumor_sample_name=args.tumor_sample_name,
        normal_sample_name=args.normal_sample_name,
        sample_id=args.sample_id,
        gatk=args.gatk,
        germline_resource=args.gnomad_vcf,
        panel_of_normals=args.panel_of_normals,
    )
    print("WES SNV Mutect2 calling completed.")
    for k, v in out.items():
        print(f"  {k}: {v}")


def cmd_snv_run_full_wes(args):
    hla = read_hla_alleles(args.hla)
    if not hla:
        raise SystemExit("No HLA alleles parsed from --hla")
    out = run_snv_wes_full(
        outdir=args.outdir,
        sample_id=args.sample_id,
        profile=load_profile(args.profile)["_profile_name"],
        hla_alleles=hla,
        tumor_sample_name=args.tumor_sample_name,
        normal_sample_name=None if args.tumor_only else args.normal_sample_name,
        somatic_vcf=args.somatic_vcf,
        tumor_bam=args.tumor_bam,
        normal_bam=args.normal_bam,
        reference_fasta=args.reference_fasta,
        intervals_bed=args.intervals_bed,
        upstream_stub=args.upstream_stub,
        immunogenicity_stub=args.immunogenicity_stub,
        gatk=args.gatk,
        germline_resource=args.gnomad_vcf,
        panel_of_normals=args.panel_of_normals,
        expression=_sv_optional_path(args.expression),
        cnv=_sv_optional_path(args.cnv),
        normal_expression=_sv_optional_path(args.normal_expression) or resource("normal_expression.example.tsv"),
        normal_hla_ligands=_sv_optional_path(args.normal_hla_ligands) or resource("normal_hla_ligands.example.tsv"),
        pass_only=not args.include_filtered,
    )
    print("WES SNV Phase 1 end-to-end completed.")
    for k, v in out.items():
        print(f"  {k}: {v}")


def cmd_benchmark_improve(args):
    input_path = args.input
    if args.dataset == "cedar":
        input_path = input_path or str(DEFAULT_CEDAR)
    elif args.dataset == "inhouse":
        input_path = input_path or str(DEFAULT_INHOUSE)
    else:
        input_path = input_path or str(DEFAULT_CEDAR)
    summary = run_benchmark_improve(
        input_path=input_path,
        outdir=args.outdir,
        profile_name=args.profile,
        sample_id=args.sample_id,
        stub=args.stub,
        limit=args.limit,
        skip_netmhcpan=args.skip_netmhcpan,
        skip_mhcflurry=args.skip_mhcflurry,
        skip_immunogenicity=args.skip_immunogenicity,
        skip_stabpan=args.skip_stabpan,
        reuse_tools_dir=args.reuse_tools_dir,
    )
    print("IMPROVE benchmark completed.")
    print(f"  Records: {summary['n_records']} (positive {summary['n_positive']})")
    print(f"  Report: {summary['report_md']}")
    print(f"  Metrics: {summary['metrics_tsv']}")
    if summary.get("top_auroc"):
        top = summary["top_auroc"]
        print(f"  Top predictor: {top.get('predictor')} AUROC={top.get('auroc')}")


def cmd_build_intermediates(args):
    cfg = load_run_config(args.config) if args.config else {
        "sample": {"id": args.sample_id, "profile": args.profile},
        "inputs": {
            "entry_mode": args.entry_mode,
            "pvac_files": args.pvac or [],
            "fusion_tsv": args.fusion_tsv,
            "easyfuse_tsv": getattr(args, "easyfuse_tsv", None),
            "easyfuse_pass_csv": getattr(args, "easyfuse_pass_csv", None),
            "splice_junction_tsv": args.splice_junction_tsv,
            "peptide_table": args.peptide_table,
            "raw_events": args.raw_events,
            "raw_peptides": args.raw_peptides,
            "sv_raw_events": args.sv_raw_events,
            "sv_raw_peptides": args.sv_raw_peptides,
        },
    }
    paths = build_raw_intermediates(cfg, args.outdir, root=ROOT)
    print("Built standard raw intermediates.")
    for k, v in paths.items():
        print(f"  {k}: {v}")


def cmd_build_evidence_layer(args):
    paths = build_standard_evidence_layer(
        args.outdir,
        args.profile,
        raw_events=args.raw_events,
        raw_peptides=args.raw_peptides,
        expression=args.expression,
        rna_junction=args.rna_junction,
        fusion_evidence=getattr(args, "fusion_evidence", None),
        rna_vaf=getattr(args, "rna_vaf", None),
        normal_expression=args.normal_expression,
        normal_hla_ligands=args.normal_hla_ligands,
        sample_id=args.sample_id,
    )
    print("Built standard evidence layer.")
    for k, v in paths.items():
        print(f"  {k}: {v}")


def cmd_run_full(args):
    cfg = load_run_config(args.config)
    sample = cfg.get("sample", {})
    outdir = Path(args.outdir or sample.get("outdir", "results/full"))
    upstream = run_upstream(args.config, outdir / "upstream")
    profile = sample.get("profile", "default")
    sample_id = sample.get("id", "SAMPLE001")
    entry_mode = resolve_entry_mode(cfg)
    pvac = []
    for key in ("pvacseq", "pvacfuse", "pvacsplice"):
        if key in upstream:
            pvac.append(upstream[key])
    inputs = cfg.get("inputs", {})
    for p in inputs.get("pvac_files") or []:
        path = ROOT / p if not Path(p).is_absolute() else Path(p)
        if path.is_file():
            pvac.append(str(path))

    raw_events = upstream.get("raw_events")
    raw_peptides = upstream.get("raw_peptides")
    if not (raw_events and raw_peptides):
        try:
            built = build_raw_intermediates(cfg, outdir / "intermediates", root=ROOT)
            raw_events = built["raw_events"]
            raw_peptides = built["raw_peptides"]
            if not pvac:
                pvac = []
        except ValueError:
            if not pvac:
                pvac = [str(fixture("pvacseq_aggregated.tsv")), str(fixture("pvacfuse_aggregated.tsv"))]

    tools_cfg = cfg.get("tools", {})
    run_kwargs = {
        "outdir": outdir,
        "profile_name_or_path": profile,
        "sample_id": sample_id,
        "pvac_paths": pvac if not (raw_events and raw_peptides) else [],
        "netmhcpan": upstream.get("netmhcpan"),
        "mhcflurry": upstream.get("mhcflurry"),
        "netmhcstabpan": upstream.get("netmhcstabpan"),
        "vep_appm": upstream.get("vep_appm"),
        "expression": upstream.get("expression"),
        "hla_loh": upstream.get("hla_loh"),
        "purity": upstream.get("purity"),
        "cnv": upstream.get("cnv"),
        "normal_expression": upstream.get("normal_expression"),
        "normal_hla_ligands": upstream.get("normal_hla_ligands"),
        "immunogenicity_stub": bool(tools_cfg.get("immunogenicity_stub", False)),
        "tool_executables": tools_cfg.get("executables") or {},
        "rna_junction": inputs.get("rna_junction_tsv") or inputs.get("rna_junction"),
        "entry_mode": entry_mode,
    }
    if raw_events and raw_peptides:
        run_kwargs["raw_events"] = raw_events
        run_kwargs["raw_peptides"] = raw_peptides
    final = run_v03(**run_kwargs)
    hla_alleles = list(inputs.get("hla_alleles") or [])
    variant_peptides = upstream.get("variant_peptides")
    if variant_peptides and hla_alleles:
        from .adapters.peptide_netmhcpan import annotate_variant_peptide_tsv

        prime_tsv = outdir / "tools" / "prime.tsv"
        bigmhc_tsv = outdir / "tools" / "bigmhc_im.tsv"
        iedb_tsv = outdir / "presentation" / "iedb_immunogenicity.tsv"
        ann = annotate_variant_peptide_tsv(
            variant_peptides,
            hla_alleles,
            netmhcpan_xls=upstream.get("netmhcpan"),
            mhcflurry_csv=upstream.get("mhcflurry"),
            netmhcstabpan_tsv=upstream.get("netmhcstabpan"),
            prime_tsv=prime_tsv if prime_tsv.is_file() else None,
            bigmhc_im_tsv=bigmhc_tsv if bigmhc_tsv.is_file() else None,
            iedb_immunogenicity_tsv=iedb_tsv if iedb_tsv.is_file() else None,
            output_tsv=outdir / "upstream" / "tools" / "variant_peptides.annotated.tsv",
        )
        final["variant_peptides_annotated"] = ann["output_tsv"]
        comprehensive = build_comprehensive_peptide_evidence(
            output_tsv=final["comprehensive_peptide_evidence"],
            annotated_peptides=ann["output_tsv"],
            ranked_peptides=final["ranked_peptides"],
            raw_peptides=final["raw_peptides"],
            raw_events=final["raw_events"],
            presentation_evidence=final["presentation_evidence"],
            appm_peptide_modifiers=final["appm_peptide_modifiers"],
            ccf_2=final["ccf_2"],
            expression_evidence=final["expression_evidence"],
            rna_junction_evidence=final["rna_junction_evidence"],
            peptide_safety=final["peptide_safety"],
            peptide_escape_flags=final["peptide_escape_flags"],
            validation_plan=final["validation_plan"],
        )
        final["comprehensive_peptide_evidence_rows"] = str(comprehensive["rows"])
    print("Full pipeline completed. Ranked outputs use v03 schema-compatible file names.")
    for k, v in final.items():
        print(f"  {k}: {v}")


def build_parser():
    p = argparse.ArgumentParser(prog="neoag-v03")
    sub = p.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("run-demo")
    d.add_argument("--profile", default="default"); d.add_argument("--sample-id", default="DEMO_V03"); d.add_argument("--outdir", required=True); d.set_defaults(func=cmd_run_demo)
    r = sub.add_parser("run-v03")
    r.add_argument("--profile", default="default"); r.add_argument("--sample-id", default="SAMPLE001"); r.add_argument("--outdir", required=True)
    r.add_argument("--pvac", nargs="*"); r.add_argument("--raw-events"); r.add_argument("--raw-peptides")
    r.add_argument("--netmhcpan"); r.add_argument("--mhcflurry"); r.add_argument("--netmhcstabpan"); r.add_argument("--vep-appm"); r.add_argument("--expression"); r.add_argument("--hla-loh"); r.add_argument("--purity"); r.add_argument("--cnv"); r.add_argument("--normal-expression"); r.add_argument("--normal-hla-ligands")
    r.add_argument("--rna-junction", help="Optional RNA junction support TSV")
    r.add_argument("--entry-mode", help="Entry mode label for provenance (snv_indel/fusion/splice_junction/sv/peptide_only/e2e)")
    r.add_argument("--immunogenicity-stub", action="store_true", default=False, help="Use stub immunogenicity scores (demo/testing only)")
    r.set_defaults(func=cmd_run_v03)
    pp = sub.add_parser("parse-pvac")
    pp.add_argument("--sample-id", default="SAMPLE001"); pp.add_argument("--profile", default="default"); pp.add_argument("--pvac", nargs="+", required=True); pp.add_argument("--events-out", required=True); pp.add_argument("--peptides-out", required=True); pp.set_defaults(func=cmd_parse_pvac)
    pn = sub.add_parser("parse-netmhcpan")
    pn.add_argument("--sample-id", default="SAMPLE001"); pn.add_argument("--input", required=True); pn.add_argument("--out", required=True); pn.set_defaults(func=cmd_parse_netmhcpan)
    pm = sub.add_parser("parse-mhcflurry")
    pm.add_argument("--sample-id", default="SAMPLE001"); pm.add_argument("--input", required=True); pm.add_argument("--out", required=True); pm.set_defaults(func=cmd_parse_mhcflurry)
    bp = sub.add_parser("build-presentation-evidence")
    bp.add_argument("--raw-peptides", required=True); bp.add_argument("--netmhcpan"); bp.add_argument("--mhcflurry"); bp.add_argument("--netmhcstabpan"); bp.add_argument("--profile", default="default"); bp.add_argument("--out", required=True); bp.set_defaults(func=cmd_build_presentation)
    ap = sub.add_parser("appm-lite")
    ap.add_argument("--sample-id", default="SAMPLE001"); ap.add_argument("--vep-tsv"); ap.add_argument("--expression"); ap.add_argument("--hla-loh"); ap.add_argument("--cnv"); ap.add_argument("--raw-peptides"); ap.add_argument("--profile", default="default"); ap.add_argument("--outdir", required=True); ap.set_defaults(func=cmd_appm)
    ap2 = sub.add_parser("appm-2", help="Build APPM 2.0 gene/pathway/peptide evidence sidecars")
    ap2.add_argument("--sample-id", default="SAMPLE001"); ap2.add_argument("--vep-tsv"); ap2.add_argument("--expression"); ap2.add_argument("--hla-loh"); ap2.add_argument("--cnv"); ap2.add_argument("--raw-peptides"); ap2.add_argument("--tumor-purity", help="Tumor purity TSV for APPM 2.0 input evidence status"); ap2.add_argument("--profile", default="default"); ap2.add_argument("--outdir", required=True); ap2.set_defaults(func=cmd_appm_2)
    cc = sub.add_parser("ccf-lite")
    cc.add_argument("--events", required=True); cc.add_argument("--purity"); cc.add_argument("--cnv"); cc.add_argument("--profile", default="default"); cc.add_argument("--out", required=True); cc.set_defaults(func=cmd_ccf)
    cc2 = sub.add_parser("ccf-2", help="Build copy-number/multiplicity-aware CCF 2.0 table")
    cc2.add_argument("--events", required=True); cc2.add_argument("--purity"); cc2.add_argument("--cnv"); cc2.add_argument("--external-clonality"); cc2.add_argument("--svclone"); cc2.add_argument("--sidecar-dir"); cc2.add_argument("--input-qc-out"); cc2.add_argument("--conflicts-out"); cc2.add_argument("--clusters-out"); cc2.add_argument("--profile", default="default"); cc2.add_argument("--out", required=True); cc2.set_defaults(func=cmd_ccf_2)
    sc = sub.add_parser("score-v03")
    sc.add_argument("--raw-events", required=True); sc.add_argument("--raw-peptides", required=True); sc.add_argument("--presentation", required=True); sc.add_argument("--appm-summary"); sc.add_argument("--appm-peptide-modifiers", help="APPM 2.0 peptide modifiers TSV (defaults to sibling of --appm-summary)"); sc.add_argument("--ccf"); sc.add_argument("--normal-expression"); sc.add_argument("--normal-hla-ligands"); sc.add_argument("--peptide-safety"); sc.add_argument("--peptide-escape-flags"); sc.add_argument("--profile", default="default"); sc.add_argument("--out-events", required=True); sc.add_argument("--out-peptides", required=True); sc.set_defaults(func=cmd_score)

    ps = sub.add_parser("peptide-safety", help="Build peptide_safety.tsv: reference proteome, normal ligandome, normal junction and anchor-risk safety gate")
    ps.add_argument("--raw-events", required=True)
    ps.add_argument("--raw-peptides", required=True)
    ps.add_argument("--profile", default="default")
    ps.add_argument("--presentation")
    ps.add_argument("--normal-expression")
    ps.add_argument("--normal-hla-ligands")
    ps.add_argument("--reference-proteome")
    ps.add_argument("--normal-junctions")
    ps.add_argument("--out", required=True)
    ps.add_argument("--event-out")
    ps.set_defaults(func=cmd_peptide_safety)

    ie = sub.add_parser("immune-escape", help="Build immune_escape_summary.tsv and peptide_escape_flags.tsv from HLA LOH/APM/JAK/B2M evidence")
    ie.add_argument("--sample-id", required=True)
    ie.add_argument("--raw-peptides", required=True)
    ie.add_argument("--profile", default="default")
    ie.add_argument("--vep-tsv")
    ie.add_argument("--cnv")
    ie.add_argument("--expression")
    ie.add_argument("--hla-loh")
    ie.add_argument("--appm-gene-status")
    ie.add_argument("--appm-pathway-status")
    ie.add_argument("--ccf")
    ie.add_argument("--therapy-context", choices=["vaccine", "tcr_target", "immunomonitoring", "discovery"])
    ie.add_argument("--ranked-peptides", help="Ranked peptides TSV for affected top-candidate burden counts")
    ie.add_argument("--top-priority-threshold", default="B_CAUTION", help="Top candidate cutoff for burden counts (default: B_CAUTION)")
    ie.add_argument("--outdir", required=True)
    ie.set_defaults(func=cmd_immune_escape)

    vp = sub.add_parser("validation-plan-v03")
    vp.add_argument("--ranked-peptides", required=True)
    vp.add_argument("--variant-peptides", help="Variant/fusion peptide catalog with minigene columns")
    vp.add_argument("--outdir", help="Run output dir; auto-discovers upstream/tools/*_peptides.tsv")
    vp.add_argument("--out", required=True)
    vp.set_defaults(func=cmd_validation_plan)

    rp = sub.add_parser("report-v03")
    rp.add_argument("--profile", default="default")
    rp.add_argument("--ranked-events", required=True)
    rp.add_argument("--ranked-peptides", required=True)
    rp.add_argument("--appm-summary")
    rp.add_argument("--validation-plan")
    rp.add_argument("--outdir", help="Run output dir; loads APPM/safety/escape/CCF sidecars for technical report")
    rp.add_argument("--provenance", help="provenance.v03.json for tool versions and file paths")
    rp.add_argument("--sample-id")
    rp.add_argument("--audience", choices=["both", "patient", "technical"], default="both",
                    help="Patient communication, research/technical, or both HTML reports")
    rp.add_argument("--out", required=True, help="Output path; with --audience both, also writes sibling patient/technical files")
    rp.set_defaults(func=cmd_report)

    rp41 = sub.add_parser("report-v041", help="Build v0.4.1 APPM/escape/safety/CCF evidence HTML report")
    rp41.add_argument("--profile", default="default")
    rp41.add_argument("--ranked-events", required=True)
    rp41.add_argument("--ranked-peptides", required=True)
    rp41.add_argument("--appm-summary")
    rp41.add_argument("--validation-plan")
    rp41.add_argument("--appm-gene-status")
    rp41.add_argument("--appm-module-scores")
    rp41.add_argument("--appm-submodule-scores")
    rp41.add_argument("--appm-conflicts")
    rp41.add_argument("--appm-peptide-modifiers")
    rp41.add_argument("--immune-escape-summary")
    rp41.add_argument("--peptide-escape-flags")
    rp41.add_argument("--peptide-safety")
    rp41.add_argument("--ccf")
    rp41.add_argument("--out", required=True)
    rp41.set_defaults(func=cmd_report_v041)

    bs = sub.add_parser("benchmark-system", help="Run synthetic, sensitivity, and optional ligandome/MS system benchmarks")
    bs.add_argument("--outdir", required=True)
    bs.add_argument("--profile", default="default")
    bs.add_argument("--mode", choices=["all", "synthetic", "sensitivity", "ligandome-ms", "ligandome_ms"], default="all")
    bs.add_argument("--ligandome-ms")
    bs.add_argument("--ranked-peptides")
    bs.add_argument("--appm-summary")
    bs.add_argument("--appm-module-scores")
    bs.add_argument("--appm-submodule-scores")
    bs.add_argument("--peptide-appm-flags")
    bs.add_argument("--peptide-escape-flags")
    bs.set_defaults(func=cmd_benchmark_system)

    va = sub.add_parser(
        "vep-annotate",
        help="pVACseq-compatible VEP annotation (Wildtype + Frameshift plugins)",
    )
    va.add_argument("--config", help="TOML with [vep] section (conf/vep.annotate.example.toml)")
    va.add_argument("--input-vcf", help="Unannotated somatic VCF (.vcf/.vcf.gz)")
    va.add_argument("--output-vcf", help="Annotated VCF output path")
    va.add_argument("--fasta", "--reference-fasta", dest="reference_fasta", help="GRCh38 reference FASTA")
    va.add_argument("--cache-dir", help="VEP cache dir (default: NEOAG_VEP_CACHE or ~/.vep)")
    va.add_argument("--plugins-dir", help="VEP_plugins dir (default: work/vep_plugins)")
    va.add_argument("--workdir", help="VEP working directory for logs")
    va.add_argument("--fork", type=int, default=4)
    va.add_argument("--online", action="store_true", help="Use Ensembl online VEP (no --cache)")
    va.add_argument("--no-pick", action="store_true", help="Disable VEP --pick")
    va.add_argument("--no-index", action="store_true", help="Skip tabix index on .vcf.gz output")
    va.add_argument("--expression-custom", help="Optional VEP --custom file for GX expression")
    va.set_defaults(func=cmd_vep_annotate)

    evp = sub.add_parser(
        "extract-variant-peptides",
        help="VEP VCF → variant short peptides with minigene and optional normal-proteome filter",
    )
    evp.add_argument("--input-vcf", required=True, help="VEP-annotated VCF (.vcf/.vcf.gz)")
    evp.add_argument("--output", required=True, help="Output TSV path")
    evp.add_argument("--sample-id", default="SAMPLE", help="Prefix for peptide_id")
    evp.add_argument(
        "--lengths",
        default="8,9,10,11",
        help="Peptide lengths, comma-separated (default: 8,9,10,11)",
    )
    evp.add_argument(
        "--length-min",
        type=int,
        default=None,
        help="Inclusive min peptide length (use with --length-max instead of --lengths)",
    )
    evp.add_argument(
        "--length-max",
        type=int,
        default=None,
        help="Inclusive max peptide length (use with --length-min instead of --lengths)",
    )
    evp.add_argument(
        "--mini-len",
        type=int,
        default=10,
        help="Flanking amino acids for minigene/minigene_nt (default: 10)",
    )
    evp.add_argument(
        "--normal-proteome-fasta",
        default=None,
        help="Normal/reference proteome FASTA (.fa/.fasta[.gz]); default: NEOAG_NORMAL_PROTEOME_FASTA; enables filtering",
    )
    evp.add_argument(
        "--filter-normal-proteome",
        action="store_true",
        help="Drop peptides found in --normal-proteome-fasta (default when FASTA is set)",
    )
    evp.add_argument(
        "--annotate-normal-proteome-only",
        action="store_true",
        help="With --normal-proteome-fasta: set in_normal_proteome column but do not filter",
    )
    evp.add_argument(
        "--include-filtered",
        action="store_true",
        help="Include non-PASS variants (default: PASS only)",
    )
    evp.add_argument(
        "--exclude-multi-aa",
        action="store_true",
        help="Skip multi-residue substitutions, inframe_multi, and complex consequences",
    )
    evp.add_argument(
        "--single-aa-only",
        action="store_true",
        help="Keep only single_aa missense variants (strictest filter)",
    )
    evp.add_argument(
        "--hla-alleles",
        default=None,
        help="Comma-separated sample HLA alleles (e.g. HLA-A*02:06,HLA-B*13:02)",
    )
    evp.add_argument(
        "--netmhcpan-xls",
        default=None,
        help="Optional NetMHCpan XLS/TSV to annotate MT/WT binding columns",
    )
    evp.add_argument(
        "--annotate-netmhcpan",
        action="store_true",
        help="Query IEDB NetMHCpan API for MT/WT binding (slow; prefer --netmhcpan-xls)",
    )
    evp.add_argument(
        "--mhcflurry-csv",
        default=None,
        help="Optional MHCflurry CSV to annotate MT/WT binding columns",
    )
    evp.add_argument(
        "--tumor-sample-name",
        default=None,
        help="Tumor sample column in VCF for VAF extraction (default: second sample or --sample-id)",
    )
    evp.set_defaults(func=cmd_extract_variant_peptides)

    ct = sub.add_parser("check-tools", help="Check availability of integrated bioinformatics tools")
    ct.set_defaults(func=cmd_check_tools)

    rt = sub.add_parser("run-tool", help="Run a single upstream tool (netmhcpan, mhcflurry, pvacseq, ...)")
    rt.add_argument("tool", choices=sorted(RUNNERS.keys()))
    rt.add_argument("--sample-id", default="SAMPLE001")
    rt.add_argument("--output", required=True)
    rt.add_argument("--workdir", default="work/tools")
    rt.add_argument("--config", help="Optional run TOML for stub/executable overrides")
    rt.add_argument("--stub", action="store_true")
    rt.add_argument("--raw-peptides")
    rt.add_argument("--tumor-vcf")
    rt.add_argument("--normal-vcf")
    rt.add_argument("--fusion-tsv")
    rt.add_argument("--variants-vcf")
    rt.set_defaults(func=cmd_run_tool)

    ru = sub.add_parser("run-upstream", help="Run enabled tools from conf/run.*.toml")
    ru.add_argument("--config", default=str(ROOT / "conf" / "run.stub.toml"))
    ru.add_argument("--outdir")
    ru.set_defaults(func=cmd_run_upstream)

    rf = sub.add_parser("run-full", help="Upstream tools + current scoring/reporting into v03 schema-compatible outputs")
    rf.add_argument("--config", default=str(ROOT / "conf" / "run.stub.toml"))
    rf.add_argument("--outdir")
    rf.set_defaults(func=cmd_run_full)

    bi = sub.add_parser("build-intermediates", help="Build parsed/raw_events + raw_peptides (multi-entry A–F)")
    bi.add_argument("--outdir", required=True)
    bi.add_argument("--config", help="Run TOML with [inputs.entry_mode]")
    bi.add_argument("--sample-id", default="SAMPLE001")
    bi.add_argument("--profile", default="default")
    bi.add_argument("--entry-mode", default="pvac")
    bi.add_argument("--pvac", nargs="*")
    bi.add_argument("--fusion-tsv")
    bi.add_argument("--easyfuse-tsv", help="EasyFuse fusions.csv or fusions.pass.csv")
    bi.add_argument("--easyfuse-pass-csv", help="Alias for EasyFuse fusions.pass.csv")
    bi.add_argument("--splice-junction-tsv")
    bi.add_argument("--peptide-table")
    bi.add_argument("--raw-events")
    bi.add_argument("--raw-peptides")
    bi.add_argument("--sv-raw-events")
    bi.add_argument("--sv-raw-peptides")
    bi.set_defaults(func=cmd_build_intermediates)

    be = sub.add_parser("build-evidence-layer", help="Write expression/RNA junction/safety evidence TSVs")
    be.add_argument("--outdir", required=True)
    be.add_argument("--profile", default="default")
    be.add_argument("--sample-id", default="SAMPLE001")
    be.add_argument("--raw-events")
    be.add_argument("--raw-peptides")
    be.add_argument("--expression")
    be.add_argument("--rna-junction")
    be.add_argument("--rna-vaf", help="RNA allele-count/VAF TSV keyed by event_id, gene, or chrom:pos:ref>alt")
    be.add_argument("--fusion-evidence", help="parsed/fusion_evidence.tsv from EasyFuse adapter")
    be.add_argument("--normal-expression")
    be.add_argument("--normal-hla-ligands")
    be.set_defaults(func=cmd_build_evidence_layer)

    sv = sub.add_parser("sv-build-raw", help="Build v03 schema-compatible raw_events/raw_peptides from WGS tumor-normal SV VCFs")
    sv.add_argument("--sample-id", required=True)
    sv.add_argument("--profile", default="sv_wgs_phase1")
    sv.add_argument("--sv-vcf", nargs="+", required=True, help="One or more SV VCF/VCF.GZ files from Manta/SvABA/GRIDSS/DELLY or merged calls")
    sv.add_argument("--callers", nargs="+", help="Caller names matching --sv-vcf order, e.g. Manta SvABA GRIDSS2")
    sv.add_argument("--reference-fasta", required=True)
    sv.add_argument("--gencode-gtf", required=True)
    sv.add_argument("--hla", required=True, help="HLA allele file or comma/space-separated HLA allele list")
    sv.add_argument("--outdir", required=True)
    sv.add_argument("--tumor-sample-name")
    sv.add_argument("--normal-sample-name")
    sv.add_argument("--expression")
    sv.add_argument("--rna-junctions")
    sv.add_argument("--normal-expression")
    sv.add_argument("--normal-hla-ligands")
    sv.add_argument("--merge-distance-bp", type=int, default=200)
    sv.add_argument("--tier1-only", action="store_true", help="Only export Tier1 SV events")
    sv.set_defaults(func=cmd_sv_build_raw)

    svs = sub.add_parser(
        "sv-score-v03",
        help="Run NetMHCpan/MHCflurry + score_v03 on SV Phase 1 raw tables",
    )
    svs.add_argument("--sample-id", required=True)
    svs.add_argument("--profile", default="sv_wgs_phase1")
    svs.add_argument("--outdir", required=True)
    svs.add_argument("--sv-outdir", help="Directory from sv-build-raw (uses parsed/raw_*.tsv)")
    svs.add_argument("--raw-events", help="Override path to raw_events.tsv")
    svs.add_argument("--raw-peptides", help="Override path to raw_peptides.tsv")
    svs.add_argument("--netmhcpan", help="Pre-computed NetMHCpan output (skip run)")
    svs.add_argument("--mhcflurry", help="Pre-computed MHCflurry output (skip run)")
    svs.add_argument("--vep-appm")
    svs.add_argument("--expression")
    svs.add_argument("--hla-loh")
    svs.add_argument("--purity")
    svs.add_argument("--cnv")
    svs.add_argument("--normal-expression")
    svs.add_argument("--normal-hla-ligands")
    svs.add_argument("--reference-proteome")
    svs.add_argument("--normal-junctions")
    svs.add_argument("--binding-stub", action="store_true", help="Use fixture binding outputs")
    svs.add_argument("--immunogenicity-stub", action="store_true", default=False)
    svs.add_argument("--no-immunogenicity-stub", action="store_false", dest="immunogenicity_stub")
    svs.add_argument("--skip-binding", action="store_true", help="Require --netmhcpan and --mhcflurry")
    svs.set_defaults(func=cmd_sv_score_v03)

    svf = sub.add_parser(
        "sv-run-full",
        help="SV Phase 1 adapter + NetMHCpan/MHCflurry + score_v03 in one command",
    )
    svf.add_argument("--sample-id", required=True)
    svf.add_argument("--profile", default="sv_wgs_phase1")
    svf.add_argument("--sv-vcf", nargs="+", required=True)
    svf.add_argument("--callers", nargs="+")
    svf.add_argument("--reference-fasta", required=True)
    svf.add_argument("--gencode-gtf", required=True)
    svf.add_argument("--hla", required=True)
    svf.add_argument("--outdir", required=True)
    svf.add_argument("--tumor-sample-name")
    svf.add_argument("--normal-sample-name")
    svf.add_argument("--expression")
    svf.add_argument("--rna-junctions")
    svf.add_argument("--normal-expression")
    svf.add_argument("--normal-hla-ligands")
    svf.add_argument("--reference-proteome")
    svf.add_argument("--normal-junctions")
    svf.add_argument("--vep-appm")
    svf.add_argument("--hla-loh")
    svf.add_argument("--purity")
    svf.add_argument("--cnv")
    svf.add_argument("--merge-distance-bp", type=int, default=200)
    svf.add_argument("--tier1-only", action="store_true")
    svf.add_argument("--binding-stub", action="store_true")
    svf.add_argument("--immunogenicity-stub", action="store_true", default=False)
    svf.add_argument("--no-immunogenicity-stub", action="store_false", dest="immunogenicity_stub")
    svf.add_argument("--skip-binding", action="store_true")
    svf.set_defaults(func=cmd_sv_run_full)

    svw = sub.add_parser(
        "sv-build-raw-wes",
        help="Build v03 schema-compatible raw tables from WES tumor-normal SV VCFs (Phase 1.5 tiers)",
    )
    svw.add_argument("--sample-id", required=True)
    svw.add_argument("--profile", default="sv_wes_phase1_5")
    svw.add_argument("--sv-vcf", nargs="+", required=True)
    svw.add_argument("--callers", nargs="+")
    svw.add_argument("--reference-fasta", required=True)
    svw.add_argument("--gencode-gtf", required=True)
    svw.add_argument("--hla", required=True)
    svw.add_argument("--outdir", required=True)
    svw.add_argument("--tumor-sample-name")
    svw.add_argument("--normal-sample-name")
    svw.add_argument("--expression")
    svw.add_argument("--rna-junctions")
    svw.add_argument("--normal-expression")
    svw.add_argument("--normal-hla-ligands")
    svw.add_argument("--capture-bed", help="WES capture BED; enables capture-aware Phase 1.5 status")
    svw.add_argument("--capture-near-bp", type=int, default=250)
    svw.add_argument("--capture-slop-bp", type=int, default=1000)
    svw.add_argument("--merge-distance-bp", type=int, default=200)
    svw.add_argument("--tier1-only", action="store_true", help="Only export WES_Tier1 SV events")
    svw.set_defaults(func=cmd_sv_build_raw_wes)

    svfw = sub.add_parser(
        "sv-run-full-wes",
        help="WES SV Phase 1.5 adapter + NetMHCpan/MHCflurry + score_v03 in one command",
    )
    svfw.add_argument("--sample-id", required=True)
    svfw.add_argument("--profile", default="sv_wes_phase1_5")
    svfw.add_argument("--sv-vcf", nargs="+", required=True)
    svfw.add_argument("--callers", nargs="+")
    svfw.add_argument("--reference-fasta", required=True)
    svfw.add_argument("--gencode-gtf", required=True)
    svfw.add_argument("--hla", required=True)
    svfw.add_argument("--outdir", required=True)
    svfw.add_argument("--tumor-sample-name")
    svfw.add_argument("--normal-sample-name")
    svfw.add_argument("--expression")
    svfw.add_argument("--rna-junctions")
    svfw.add_argument("--normal-expression")
    svfw.add_argument("--normal-hla-ligands")
    svfw.add_argument("--capture-bed", help="WES capture BED; enables capture-aware Phase 1.5 status")
    svfw.add_argument("--capture-near-bp", type=int, default=250)
    svfw.add_argument("--capture-slop-bp", type=int, default=1000)
    svfw.add_argument("--reference-proteome")
    svfw.add_argument("--normal-junctions")
    svfw.add_argument("--vep-appm")
    svfw.add_argument("--hla-loh")
    svfw.add_argument("--purity")
    svfw.add_argument("--cnv")
    svfw.add_argument("--merge-distance-bp", type=int, default=200)
    svfw.add_argument("--tier1-only", action="store_true")
    svfw.add_argument("--binding-stub", action="store_true")
    svfw.add_argument("--immunogenicity-stub", action="store_true", default=False)
    svfw.add_argument("--no-immunogenicity-stub", action="store_false", dest="immunogenicity_stub")
    svfw.add_argument("--skip-binding", action="store_true")
    svfw.set_defaults(func=cmd_sv_run_full_wes)

    snvc = sub.add_parser(
        "snv-call-wes",
        help="WES tumor-normal BAM → Mutect2 + FilterMutectCalls",
    )
    snvc.add_argument("--sample-id", required=True)
    snvc.add_argument("--tumor-bam", required=True)
    snvc.add_argument("--normal-bam", required=True)
    snvc.add_argument("--reference-fasta", required=True)
    snvc.add_argument("--intervals-bed", required=True, help="WES capture BED")
    snvc.add_argument("--tumor-sample-name", required=True)
    snvc.add_argument("--normal-sample-name", required=True)
    snvc.add_argument("--outdir", required=True)
    snvc.add_argument("--gatk", help="GATK executable (default: gatk on PATH)")
    snvc.add_argument("--gnomad-vcf", help="Optional gnomAD resource for FilterMutectCalls")
    snvc.add_argument("--panel-of-normals", help="Optional PoN VCF for FilterMutectCalls")
    snvc.set_defaults(func=cmd_snv_call_wes)

    snvf = sub.add_parser(
        "snv-run-full-wes",
        help="WES SNV Mutect2 (or --somatic-vcf) + pVAC stub + score_v03",
    )
    snvf.add_argument("--sample-id", required=True)
    snvf.add_argument("--profile", default="default")
    snvf.add_argument("--outdir", required=True)
    snvf.add_argument("--hla", required=True, help="HLA file or comma-separated alleles")
    snvf.add_argument("--tumor-sample-name", default="TUMOR")
    snvf.add_argument("--normal-sample-name", default="NORMAL")
    snvf.add_argument("--tumor-only", action="store_true", help="Single-sample VCF (no matched normal)")
    snvf.add_argument("--include-filtered", action="store_true", help="Include non-PASS variants in pVACseq")
    snvf.add_argument("--somatic-vcf", help="Skip calling; use existing filtered VCF")
    snvf.add_argument("--tumor-bam")
    snvf.add_argument("--normal-bam")
    snvf.add_argument("--reference-fasta")
    snvf.add_argument("--intervals-bed")
    snvf.add_argument("--gatk")
    snvf.add_argument("--gnomad-vcf")
    snvf.add_argument("--panel-of-normals")
    snvf.add_argument("--expression")
    snvf.add_argument("--cnv")
    snvf.add_argument("--normal-expression")
    snvf.add_argument("--normal-hla-ligands")
    snvf.add_argument(
        "--upstream-stub",
        action="store_true",
        default=True,
        help="Use fixture pVAC/binding outputs (default: on)",
    )
    snvf.add_argument("--no-upstream-stub", action="store_false", dest="upstream_stub")
    snvf.add_argument("--immunogenicity-stub", action="store_true", default=False)
    snvf.add_argument("--no-immunogenicity-stub", action="store_false", dest="immunogenicity_stub")
    snvf.set_defaults(func=cmd_snv_run_full_wes)

    bi = sub.add_parser(
        "benchmark-improve",
        help="Benchmark MHC presentation predictors on IMPROVE/CEDAR immunogenicity labels",
    )
    bi.add_argument(
        "--input",
        help="IMPROVE TSV path (default: data/improve/CEDAR or inhouse by --dataset)",
    )
    bi.add_argument(
        "--dataset",
        choices=["cedar", "inhouse", "custom"],
        default="cedar",
        help="Preset dataset when --input is omitted",
    )
    bi.add_argument("--outdir", required=True, help="Output directory for metrics and report")
    bi.add_argument("--profile", default="default")
    bi.add_argument("--sample-id", default="IMPROVE_CEDAR")
    bi.add_argument("--stub", action="store_true", help="Use fixture tool outputs (fast smoke test)")
    bi.add_argument("--limit", type=int, help="Max peptide–HLA pairs (for quick runs)")
    bi.add_argument("--skip-netmhcpan", action="store_true", help="Skip NetMHCpan (IEDB API is slow)")
    bi.add_argument("--skip-mhcflurry", action="store_true", help="Skip MHCflurry")
    bi.add_argument(
        "--skip-immunogenicity",
        action="store_true",
        help="Skip immunogenicity predictors (PRIME, BigMHC_IM; cached tools/ still reused)",
    )
    bi.add_argument(
        "--skip-stabpan",
        action="store_true",
        help="Skip NetMHCstabpan pMHC stability (IEDB API is slow)",
    )
    bi.add_argument(
        "--reuse-tools-dir",
        help="Reuse cached netmhcpan/mhcflurry/stabpan outputs from this tools/ directory "
        "(default: results/benchmark_cedar_v2/tools when skipping predictors)",
    )
    bi.set_defaults(func=cmd_benchmark_improve)

    cpi = sub.add_parser(
        "convert-peptide-input",
        help="Convert CSV/TSV peptide table to strict peptide–HLA pair TSVs",
    )
    cpi.add_argument("-i", "--input", required=True, help="Input CSV/TSV with peptide + HLA columns")
    cpi.add_argument("-o", "--outdir", required=True, help="Output directory")
    cpi.add_argument("--sample-id", default="SAMPLE001")
    cpi.set_defaults(func=cmd_convert_peptide_input)

    cl = sub.add_parser(
        "convert-lohhla",
        help="Convert LOHHLA HLAlossPrediction_CI output to hla_loh.tsv",
    )
    cl.add_argument("-i", "--input", required=True, help="LOHHLA *HLAlossPrediction_CI* file")
    cl.add_argument("-o", "--output", required=True, help="Output hla_loh.tsv path")
    cl.set_defaults(func=cmd_convert_lohhla)

    cs = sub.add_parser(
        "convert-spechla",
        help="Convert SpecHLA merge.hla.copy.txt to hla_loh.tsv",
    )
    cs.add_argument("-i", "--input", required=True, help="SpecHLA merge.hla.copy.txt")
    cs.add_argument("-o", "--output", required=True, help="Output hla_loh.tsv path")
    cs.set_defaults(func=cmd_convert_spechla)

    hx = sub.add_parser(
        "crosscheck-hla-loh",
        help="Cross-check normalized LOHHLA and SpecHLA hla_loh.tsv files and write a consensus table",
    )
    hx.add_argument("--lohhla-hla-loh", help="hla_loh.tsv produced by convert-lohhla")
    hx.add_argument("--spechla-hla-loh", help="hla_loh.tsv produced by convert-spechla")
    hx.add_argument("--out", required=True, help="Output hla_loh.crosscheck.tsv path")
    hx.add_argument("--consensus-out", help="Optional downstream-compatible consensus hla_loh.tsv path")
    hx.add_argument("--strict-consensus-only", action="store_true", help="Exclude single-tool LOH calls from consensus output")
    hx.set_defaults(func=cmd_crosscheck_hla_loh)

    cf = sub.add_parser(
        "convert-facets",
        help="Convert FACETS purity text and optional cncf segments to neoag TSVs",
    )
    cf.add_argument("--purity-input", required=True, help="FACETS purity.txt from runFACETS.R")
    cf.add_argument("--purity-output", required=True, help="Output purity.tsv path")
    cf.add_argument("--sample-id", default="SAMPLE001")
    cf.add_argument("--cnv-input", help="FACETS cncf.tsv from runFACETS.R")
    cf.add_argument("--cnv-output", help="Output cnv_segments.tsv path")
    cf.set_defaults(func=cmd_convert_facets)

    ca = sub.add_parser(
        "convert-ascat",
        help="Convert ASCAT summary and optional segment table to neoag TSVs",
    )
    ca.add_argument("--summary-input", required=True, help="ASCAT ascat_summary.tsv")
    ca.add_argument("--purity-output", required=True, help="Output purity.tsv path")
    ca.add_argument("--sample-id", default="SAMPLE001")
    ca.add_argument("--segments-input", help="ASCAT ascat_segments.tsv")
    ca.add_argument("--cnv-output", help="Output cnv_segments.tsv path")
    ca.set_defaults(func=cmd_convert_ascat)

    pp = sub.add_parser(
        "peptide-predict",
        help="Run binding/immunogenicity predictors from flexible peptide CSV/TSV input",
    )
    pp.add_argument("-i", "--input", required=True, help="Input CSV/TSV with peptide + HLA columns")
    pp.add_argument("-o", "--outdir", required=True, help="Output directory")
    pp.add_argument("--sample-id", default="SAMPLE001")
    pp.add_argument("--profile", default="default")
    pp.add_argument("--stub", action="store_true", help="Use stub predictors (fast smoke test)")
    pp.add_argument("--skip-netmhcpan", action="store_true")
    pp.add_argument("--skip-mhcflurry", action="store_true")
    pp.add_argument("--skip-prime", action="store_true")
    pp.add_argument("--skip-bigmhc-im", action="store_true")
    pp.add_argument("--skip-deepimmuno", action="store_true")
    pp.add_argument("--skip-stabpan", action="store_true")
    pp.set_defaults(func=cmd_peptide_predict)


    # controlled-execution Phase 0-2 execution layer
    odoc = sub.add_parser("doctor", help="controlled-execution read-only health, reference, tool and release-boundary check")
    odoc.add_argument("--project-root", default=".")
    odoc.add_argument("--outdir", required=True)
    odoc.add_argument("--tools-manifest")
    odoc.add_argument("--reference-manifest")
    odoc.add_argument("--sample-manifest")
    odoc.add_argument("--profile", default="local")
    odoc.add_argument("--run-demo", action="store_true")
    odoc.add_argument("--run-pytest", action="store_true")
    odoc.add_argument("--run-nextflow", action="store_true")
    odoc.add_argument("--mini-smoke", action="store_true")
    odoc.add_argument("--skip-release-audit", action="store_true")
    odoc.add_argument("--dry-run", action="store_true")
    def _cmd_controlled_execution_doctor(args):
        from .controlled_execution.doctor import run_doctor
        res = run_doctor(project_root=args.project_root, outdir=args.outdir, tools_manifest=args.tools_manifest, reference_manifest=args.reference_manifest, sample_manifest=args.sample_manifest, profile=args.profile, run_demo=args.run_demo, run_pytest=args.run_pytest, run_nextflow=args.run_nextflow, mini_smoke=args.mini_smoke, release_audit=not args.skip_release_audit, allow_execute=not args.dry_run)
        print(f"NeoAg Doctor status: {res.status}")
        for k, v in res.outputs.items(): print(f"  {k}: {v}")
    odoc.set_defaults(func=_cmd_controlled_execution_doctor)

    opipe = sub.add_parser("pipeline-full", help="controlled-execution manifest-driven full pipeline runner (safe dry-run by default)")
    opipe.add_argument("--sample-manifest", required=True)
    opipe.add_argument("--tools-manifest")
    opipe.add_argument("--reference-manifest")
    opipe.add_argument("--project-root", default=".")
    opipe.add_argument("--outdir", required=True)
    opipe.add_argument("--profile", default="local")
    opipe.add_argument("--execute", action="store_true")
    opipe.add_argument("--strict", action="store_true")
    opipe.add_argument("--run-demo-fixture", action="store_true")
    def _cmd_controlled_execution_pipeline(args):
        from .controlled_execution.pipeline_runner import run_pipeline_full
        run = run_pipeline_full(sample_manifest=args.sample_manifest, tools_manifest=args.tools_manifest, reference_manifest=args.reference_manifest, project_root=args.project_root, outdir=args.outdir, profile=args.profile, dry_run=not args.execute, allow_partial=not args.strict, run_demo_for_fixture=args.run_demo_fixture)
        print(f"NeoAg pipeline-full status: {run.status}")
        print(f"  run_id: {run.run_id}")
        print(f"  run_manifest: {Path(run.output_dir) / 'run_manifest.json'}")
    opipe.set_defaults(func=_cmd_controlled_execution_pipeline)

    oaudit = sub.add_parser("release-audit", help="Scan release tree for cache artifacts, private paths and patient/site hints")
    oaudit.add_argument("--root", default=".")
    oaudit.add_argument("--outdir", required=True)
    oaudit.add_argument("--scan-generated-dirs", action="store_true")
    def _cmd_controlled_execution_release_audit(args):
        from .controlled_execution.release_audit import scan_release_boundary, write_release_audit
        result = scan_release_boundary(args.root, skip_dirs=set() if args.scan_generated_dirs else None)
        outs = write_release_audit(result, args.outdir)
        print(f"Release audit status: {result['status']}")
        for k, v in outs.items(): print(f"  {k}: {v}")
    oaudit.set_defaults(func=_cmd_controlled_execution_release_audit)

    return p

def main(argv=None):
    args = build_parser().parse_args(argv)
    args.func(args)

if __name__ == "__main__":
    main()
