#!/usr/bin/env python3
"""Generate a production manifest that reuses completed upstream tool results."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def q(value) -> str:
    return json.dumps(str(value))


def require(path: str | None, label: str) -> str:
    if not path or not Path(path).exists():
        raise SystemExit(f"{label} missing: {path}")
    return str(Path(path).resolve())


def stage(lines, name, *, outputs, source="", command="", required=True, depends=None):
    lines += ["", f"[stages.{name}]", f"required = {str(required).lower()}"]
    if source: lines.append(f"source = {q(source)}")
    if depends: lines.append("depends_on = [" + ", ".join(q(value) for value in depends) + "]")
    if command: lines.append(f"command = {q(command)}")
    lines.append(f"[stages.{name}.outputs]")
    lines.extend(f"{key} = {q(value)}" for key, value in outputs.items())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--project-root", default=Path(__file__).resolve().parents[1])
    ap.add_argument("--sample-id", required=True); ap.add_argument("--outdir", required=True); ap.add_argument("--output", required=True)
    ap.add_argument("--hla-file", required=True); ap.add_argument("--somatic-vcf")
    ap.add_argument("--facets", required=True); ap.add_argument("--sequenza", required=True); ap.add_argument("--purple", required=True)
    ap.add_argument("--purity", required=True); ap.add_argument("--cnv", required=True)
    ap.add_argument("--lohhla", required=True); ap.add_argument("--spechla-loh", required=True); ap.add_argument("--hla-loh", required=True)
    ap.add_argument("--expression"); ap.add_argument("--transcript-expression")
    ap.add_argument("--easyfuse"); ap.add_argument("--star-fusion"); ap.add_argument("--arriba")
    ap.add_argument("--junctions"); ap.add_argument("--snaf"); ap.add_argument("--splicemutr"); ap.add_argument("--normal-junctions")
    ap.add_argument("--normal-expression"); ap.add_argument("--normal-hla-ligands"); ap.add_argument("--reference-proteome")
    ap.add_argument("--netchop-executable", default="netChop"); ap.add_argument("--netchop-home", default="")
    args = ap.parse_args()
    root = Path(args.project_root).resolve(); hla = require(args.hla_file, "HLA consensus")
    lines = ["# Generated from completed upstream tool results.", "[run]", f"sample_id = {q(args.sample_id)}", 'profile = "sarcoma_rna_supported_v3_provisional"', f"outdir = {q(Path(args.outdir).resolve())}", f"hla_file = {q(hla)}", "tools_stub = false", "immunogenicity_stub = false", 'presentation_predictors = ["netmhcpan", "mhcflurry", "netmhcstabpan", "netchop"]', 'required_presentation_predictors = ["netmhcpan", "mhcflurry", "netmhcstabpan", "netchop"]', 'reports = "patient,technical"', f"netchop_executable = {q(args.netchop_executable)}"]
    if args.netchop_home: lines.append(f"netchop_home = {q(args.netchop_home)}")
    lines += ["", "[run.required_tool_groups.purity_cnv]", 'tools = ["facets", "sequenza", "purple"]', "min_successful = 2", "require_all_declared = true", "", "[run.required_tool_groups.hla_loh]", 'tools = ["lohhla", "spechla"]', "min_successful = 2", "require_all_declared = true"]
    for tool in ("facets", "sequenza", "purple"):
        stage(lines, f"purity_{tool}", outputs={f"{tool}_result": require(getattr(args, tool), tool)}, required=False)
    stage(lines, "hla_loh_lohhla", outputs={"lohhla_hla_loh": require(args.lohhla, "LOHHLA")}, required=True)
    stage(lines, "hla_loh_spechla", outputs={"spechla_hla_loh": require(args.spechla_loh, "SpecHLA LOH")}, required=True)
    candidate_stages = []
    if args.somatic_vcf:
        command = f"PYTHONPATH={root / 'src'} python {root / 'scripts/run_candidate_upstream.py'} --mode snv --input {require(args.somatic_vcf, 'somatic VCF')} --hla-file {hla} --sample-id {args.sample_id} --outdir {{outdir}}/branches/snv"
        stage(lines, "snv_indel_candidates", source="SNV_INDEL", command=command, outputs={"raw_events": "{outdir}/branches/snv/parsed/raw_events.tsv", "raw_peptides": "{outdir}/branches/snv/parsed/raw_peptides.tsv"})
        candidate_stages.append("snv_indel_candidates")
    if args.easyfuse:
        easyfuse = require(args.easyfuse, "EasyFuse")
        command = f"PYTHONPATH={root / 'src'} python -m neoag.cli build-intermediates --entry-mode fusion --easyfuse-tsv {easyfuse} --sample-id {args.sample_id} --outdir {{outdir}}/branches/fusion/intermediates"
        stage(lines, "fusion_candidates", source="EasyFuse", command=command, outputs={"raw_events": "{outdir}/branches/fusion/intermediates/parsed/raw_events.tsv", "raw_peptides": "{outdir}/branches/fusion/intermediates/parsed/raw_peptides.tsv"})
        candidate_stages.append("fusion_candidates")
        review = f"python {root / 'scripts/review_rna_fusions.py'} --easyfuse {easyfuse}"
        if args.star_fusion: review += f" --star-fusion {require(args.star_fusion, 'STAR-Fusion')}"
        if args.arriba: review += f" --arriba {require(args.arriba, 'Arriba')}"
        review += " --outdir {outdir}/branches/fusion/consensus"
        stage(lines, "fusion_cross_validation", command=review, outputs={"fusion_consensus": "{outdir}/branches/fusion/consensus/fusion_consensus.tsv"}, required=True)
    if args.junctions and (args.snaf or args.splicemutr):
        command = f"PYTHONPATH={root / 'src'} python {root / 'scripts/normalize_rna_fusion_splice.py'} --sample-id {args.sample_id} --profile sarcoma_rna_supported_v3_provisional --junctions {require(args.junctions, 'junctions')}"
        if args.snaf: command += f" --snaf {require(args.snaf, 'SNAF')}"
        if args.splicemutr: command += f" --splicemutr {require(args.splicemutr, 'SpliceMutr')}"
        if args.normal_junctions: command += f" --normal-junctions {require(args.normal_junctions, 'normal junctions')}"
        command += " --outdir {outdir}/branches/splice/intermediates"
        stage(lines, "splice_candidates", source="SpliceConsensus", command=command, outputs={"raw_events": "{outdir}/branches/splice/intermediates/raw_events.tsv", "raw_peptides": "{outdir}/branches/splice/intermediates/raw_peptides.tsv"})
        candidate_stages.append("splice_candidates")
    if not candidate_stages: raise SystemExit("At least one SNV/Fusion/Splice candidate source is required")
    lines += ["", "[evidence]", f"purity = {q(require(args.purity, 'purity consensus'))}", f"cnv = {q(require(args.cnv, 'CNV consensus'))}", f"hla_loh = {q(require(args.hla_loh, 'HLA LOH consensus'))}"]
    for key in ("expression", "transcript_expression", "normal_junctions", "normal_expression", "normal_hla_ligands", "reference_proteome"):
        value = getattr(args, key)
        if value: lines.append(f"{key} = {q(require(value, key))}")
    Path(args.output).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__": raise SystemExit(main())
