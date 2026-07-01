from __future__ import annotations

import re
import shlex
from pathlib import Path

import pytest

from neoag_v03 import cli
from neoag_v03.cli import build_parser, main

ROOT = Path(__file__).resolve().parents[1]


def _parse(argv: list[str]):
    return build_parser().parse_args(argv)


NEXTFLOW_COMMANDS = [
    ["appm-2", "--sample-id", "S", "--profile", "default", "--vep-tsv", "vep.tsv", "--expression", "expr.tsv", "--hla-loh", "hla.tsv", "--cnv", "cnv.tsv", "--raw-peptides", "raw_peptides.tsv", "--tumor-purity", "purity.tsv", "--outdir", "out"],
    ["appm-lite", "--sample-id", "S", "--profile", "default", "--vep-tsv", "vep.tsv", "--expression", "expr.tsv", "--hla-loh", "hla.tsv", "--outdir", "out"],
    ["build-presentation-evidence", "--raw-peptides", "raw_peptides.tsv", "--netmhcpan", "net.tsv", "--mhcflurry", "mhc.tsv", "--profile", "default", "--out", "presentation.tsv"],
    ["ccf-2", "--events", "raw_events.tsv", "--profile", "default", "--purity", "purity.tsv", "--cnv", "cnv.tsv", "--out", "ccf_2.tsv"],
    ["ccf-lite", "--events", "raw_events.tsv", "--profile", "default", "--purity", "purity.tsv", "--cnv", "cnv.tsv", "--out", "ccf_lite.tsv"],
    ["immune-escape", "--sample-id", "S", "--raw-peptides", "raw_peptides.tsv", "--profile", "default", "--vep-tsv", "vep.tsv", "--cnv", "cnv.tsv", "--expression", "expr.tsv", "--hla-loh", "hla.tsv", "--appm-gene-status", "gene.tsv", "--appm-pathway-status", "pathway.tsv", "--ccf", "ccf_2.tsv", "--outdir", "out"],
    ["parse-mhcflurry", "--sample-id", "S", "--input", "mhc.csv", "--out", "mhc.tsv"],
    ["parse-netmhcpan", "--sample-id", "S", "--input", "net.xls", "--out", "net.tsv"],
    ["parse-pvac", "--sample-id", "S", "--profile", "default", "--pvac", "pvac.tsv", "--events-out", "raw_events.tsv", "--peptides-out", "raw_peptides.tsv"],
    ["peptide-safety", "--raw-events", "raw_events.tsv", "--raw-peptides", "raw_peptides.tsv", "--profile", "default", "--normal-expression", "normal_expr.tsv", "--normal-hla-ligands", "normal_lig.tsv", "--out", "peptide_safety.tsv", "--event-out", "event_safety.tsv"],
    ["report-v03", "--profile", "default", "--ranked-events", "ranked_events.tsv", "--ranked-peptides", "ranked_peptides.tsv", "--appm-summary", "appm_summary.tsv", "--validation-plan", "validation.tsv", "--outdir", "results", "--audience", "both", "--out", "report.html"],
    ["report-v041", "--profile", "default", "--ranked-events", "ranked_events.tsv", "--ranked-peptides", "ranked_peptides.tsv", "--appm-summary", "appm_summary.tsv", "--appm-gene-status", "gene.tsv", "--appm-module-scores", "module.tsv", "--appm-submodule-scores", "submodule.tsv", "--appm-conflicts", "conflicts.tsv", "--appm-peptide-modifiers", "mods.tsv", "--immune-escape-summary", "escape_summary.tsv", "--peptide-escape-flags", "escape_flags.tsv", "--peptide-safety", "peptide_safety.tsv", "--ccf", "ccf_2.tsv", "--validation-plan", "validation.tsv", "--out", "report.html"],
    ["run-tool", "netmhcpan", "--sample-id", "S", "--raw-peptides", "raw_peptides.tsv", "--output", "net.xls", "--workdir", "."],
    ["run-tool", "mhcflurry", "--sample-id", "S", "--raw-peptides", "raw_peptides.tsv", "--output", "mhc.csv", "--workdir", "."],
    ["run-upstream", "--config", "conf/run.example.toml", "--outdir", "upstream"],
    ["score-v03", "--raw-events", "raw_events.tsv", "--raw-peptides", "raw_peptides.tsv", "--presentation", "presentation.tsv", "--appm-summary", "appm_summary.tsv", "--ccf", "ccf_2.tsv", "--normal-expression", "normal_expr.tsv", "--normal-hla-ligands", "normal_lig.tsv", "--peptide-safety", "peptide_safety.tsv", "--peptide-escape-flags", "escape_flags.tsv", "--appm-peptide-modifiers", "mods.tsv", "--profile", "default", "--out-events", "ranked_events.tsv", "--out-peptides", "ranked_peptides.tsv"],
    ["validation-plan-v03", "--ranked-peptides", "ranked_peptides.tsv", "--outdir", "results", "--out", "validation.tsv"],
    ["sv-build-raw", "--sample-id", "S", "--profile", "default", "--sv-vcf", "sv.vcf", "--reference-fasta", "ref.fa", "--gencode-gtf", "genes.gtf", "--hla", "hla.txt", "--outdir", "out", "--merge-distance-bp", "200"],
    ["sv-build-raw-wes", "--sample-id", "S", "--profile", "default", "--sv-vcf", "sv.vcf", "--reference-fasta", "ref.fa", "--gencode-gtf", "genes.gtf", "--hla", "hla.txt", "--outdir", "out", "--capture-bed", "capture.bed"],
    ["sv-score-v03", "--sample-id", "S", "--profile", "default", "--raw-events", "raw_events.tsv", "--raw-peptides", "raw_peptides.tsv", "--outdir", "out", "--binding-stub"],
    ["sv-run-full", "--sample-id", "S", "--profile", "default", "--sv-vcf", "sv.vcf", "--reference-fasta", "ref.fa", "--gencode-gtf", "genes.gtf", "--hla", "hla.txt", "--outdir", "out", "--binding-stub"],
    ["sv-run-full-wes", "--sample-id", "S", "--profile", "default", "--sv-vcf", "sv.vcf", "--reference-fasta", "ref.fa", "--gencode-gtf", "genes.gtf", "--hla", "hla.txt", "--outdir", "out", "--capture-bed", "capture.bed", "--binding-stub"],
]


@pytest.mark.parametrize("argv", NEXTFLOW_COMMANDS)
def test_nextflow_module_cli_contract_parse(argv):
    args = _parse(argv)
    assert callable(args.func)


def _collect_doc_commands() -> list[tuple[Path, list[str]]]:
    paths = [ROOT / "README.md", ROOT / "RELEASE.md", *sorted((ROOT / "docs").glob("*.md")), *sorted((ROOT / "conf").glob("*.toml"))]
    commands: list[tuple[Path, list[str]]] = []
    for path in paths:
        if not path.exists():
            continue
        lines = path.read_text(encoding="utf-8").splitlines()
        i = 0
        while i < len(lines):
            line = lines[i]
            if "neoag-v03" not in line and "bin/neoag-v03" not in line:
                i += 1
                continue
            chunk = line.strip()
            while chunk.endswith("\\") and i + 1 < len(lines):
                i += 1
                chunk = chunk[:-1] + " " + lines[i].strip()
            for match in re.finditer(r"(?:bin/)?neoag-v03\s+([^`#]+)", chunk):
                raw = match.group(1).strip()
                raw = raw.strip("`")
                raw = re.sub(r"\s+#.*$", "", raw)
                if not raw or raw.startswith("=") or "{" in raw or "$" in raw:
                    continue
                try:
                    argv = shlex.split(raw)
                except ValueError:
                    continue
                allowed = {"run-demo", "check-tools", "run-upstream", "run-full", "build-intermediates", "build-evidence-layer", "run-v03", "appm-2", "ccf-2", "immune-escape", "sv-build-raw", "sv-build-raw-wes", "sv-run-full", "sv-run-full-wes", "sv-score-v03", "snv-run-full-wes", "snv-call-wes", "vep-annotate", "report-v041"}
                if not argv or argv[0] not in allowed:
                    continue
                if len(argv) > 1 and not any(tok.startswith("-") for tok in argv[1:]):
                    continue
                commands.append((path, argv))
            i += 1
    # Keep docs contract useful but avoid duplicate noise.
    seen = set()
    unique = []
    for path, argv in commands:
        key = tuple(argv)
        if key not in seen:
            seen.add(key)
            unique.append((path, argv))
    return unique


@pytest.mark.parametrize("path,argv", _collect_doc_commands(), ids=lambda x: str(x)[:80])
def test_documented_cli_examples_parse(path, argv):
    args = _parse(argv)
    assert callable(args.func), f"{path}: neoag-v03 {' '.join(argv)}"


def test_cli_contract_core_parameter_passthrough(monkeypatch, tmp_path):
    calls = {}
    monkeypatch.setattr(cli, "load_profile", lambda profile: {"_profile_name": profile})

    def fake_appm_2(**kwargs):
        calls["appm_2"] = kwargs
        return {"ok": "appm"}

    monkeypatch.setattr(cli, "build_appm_2", fake_appm_2)
    main(["appm-2", "--sample-id", "S", "--outdir", str(tmp_path / "appm"), "--tumor-purity", "purity.tsv"])
    assert calls["appm_2"]["tumor_purity_tsv"] == "purity.tsv"

    def fake_score_v03(*args, **kwargs):
        calls["score_v03"] = {"args": args, "kwargs": kwargs}
        return [], []

    monkeypatch.setattr(cli, "score_v03", fake_score_v03)
    import neoag_v03.scoring_v03 as scoring_v03
    monkeypatch.setattr(scoring_v03, "resolve_appm_peptide_modifiers_tsv", lambda explicit, appm_summary_tsv=None: explicit)
    main([
        "score-v03", "--raw-events", "events.tsv", "--raw-peptides", "peptides.tsv", "--presentation", "presentation.tsv",
        "--appm-summary", "appm_summary.tsv", "--ccf", "ccf_2.tsv", "--appm-peptide-modifiers", "mods.tsv",
        "--out-events", "ranked_events.tsv", "--out-peptides", "ranked_peptides.tsv",
    ])
    assert calls["score_v03"]["kwargs"]["appm_peptide_modifiers_tsv"] == "mods.tsv"

    def fake_ccf_2(*args, **kwargs):
        calls["ccf_2"] = {"args": args, "kwargs": kwargs}
        return [{"event_id": "E1"}]

    monkeypatch.setattr(cli, "build_ccf_2", fake_ccf_2)
    main([
        "ccf-2", "--events", "events.tsv", "--out", "ccf.tsv", "--external-clonality", "external.tsv", "--svclone", "svclone.tsv",
        "--sidecar-dir", "side", "--input-qc-out", "qc.tsv", "--conflicts-out", "conflicts.tsv", "--clusters-out", "clusters.tsv",
    ])
    assert calls["ccf_2"]["kwargs"] == {
        "external_clonality_tsv": "external.tsv",
        "svclone_tsv": "svclone.tsv",
        "sidecar_dir": "side",
        "input_qc_out": "qc.tsv",
        "conflicts_out": "conflicts.tsv",
        "clusters_out": "clusters.tsv",
    }

    def fake_escape(**kwargs):
        calls["immune_escape"] = kwargs
        return {"immune_escape_summary": "summary.tsv"}

    monkeypatch.setattr(cli, "build_immune_escape_evidence", fake_escape)
    main([
        "immune-escape", "--sample-id", "S", "--raw-peptides", "peptides.tsv", "--outdir", "escape",
        "--ranked-peptides", "ranked.tsv", "--top-priority-threshold", "C",
    ])
    assert calls["immune_escape"]["ranked_peptides"] == "ranked.tsv"
    assert calls["immune_escape"]["top_priority_threshold"] == "C"

    def fake_benchmark(**kwargs):
        calls["benchmark_system"] = kwargs
        return {"summary_tsv": "summary.tsv"}

    monkeypatch.setattr(cli, "run_system_benchmark", fake_benchmark)
    main([
        "benchmark-system", "--outdir", "bench", "--mode", "ligandome-ms", "--ranked-peptides", "ranked.tsv",
        "--appm-summary", "appm_summary.tsv", "--appm-module-scores", "module.tsv", "--appm-submodule-scores", "submodule.tsv",
        "--peptide-appm-flags", "mods.tsv", "--peptide-escape-flags", "escape.tsv", "--ligandome-ms", "lig.tsv",
    ])
    assert calls["benchmark_system"]["ranked_peptides"] == "ranked.tsv"
    assert calls["benchmark_system"]["peptide_escape_flags"] == "escape.tsv"
    assert calls["benchmark_system"]["ligandome_ms"] == "lig.tsv"


def test_immunogenicity_stub_defaults_to_off_for_production_commands():
    commands = [
        ["run-v03", "--sample-id", "S", "--outdir", "out", "--raw-events", "events.tsv", "--raw-peptides", "peptides.tsv"],
        ["sv-score-v03", "--sample-id", "S", "--outdir", "out", "--raw-events", "events.tsv", "--raw-peptides", "peptides.tsv"],
        ["sv-run-full", "--sample-id", "S", "--sv-vcf", "sv.vcf", "--reference-fasta", "ref.fa", "--gencode-gtf", "genes.gtf", "--hla", "hla.txt", "--outdir", "out"],
        ["sv-run-full-wes", "--sample-id", "S", "--sv-vcf", "sv.vcf", "--reference-fasta", "ref.fa", "--gencode-gtf", "genes.gtf", "--hla", "hla.txt", "--capture-bed", "capture.bed", "--outdir", "out"],
        ["snv-run-full-wes", "--sample-id", "S", "--outdir", "out", "--hla", "hla.txt", "--somatic-vcf", "somatic.vcf"],
    ]
    for argv in commands:
        assert _parse(argv).immunogenicity_stub is False
        assert _parse([*argv, "--immunogenicity-stub"]).immunogenicity_stub is True
