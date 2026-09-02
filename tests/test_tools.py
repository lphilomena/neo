import json
from pathlib import Path
import pytest

from neoag.tools import check_tool, run_tool, run_upstream, load_run_config
from neoag.tools.registry import RunContext, ROOT
from neoag.tools.runner import RUNNERS
from neoag.utils import read_tsv
from neoag.tools.prep import unique_peptide_hla_pairs, netmhcpan_allele_string
from neoag.pipeline import run
from neoag.cli import main

def test_tool_registry_covers_runners():
    assert "netmhcpan" in RUNNERS
    assert "pvacseq" in RUNNERS

def test_netmhcpan_allele_format():
    assert netmhcpan_allele_string(["HLA-A*02:01"]) == "HLA-A0201"
    assert netmhcpan_allele_string(["HLA-A02:06", "HLA-B*13:02"]) == "HLA-A0206,HLA-B1302"

def test_unique_peptides_from_fixture():
    pairs = unique_peptide_hla_pairs(ROOT / "data/fixtures/pvacseq_aggregated.tsv")
    assert ("VVVGADGVGK", "HLA-A*11:01") in pairs

def test_run_upstream_stub(tmp_path):
    cfg = ROOT / "conf/run.stub.toml"
    outs = run_upstream(cfg, tmp_path / "up")
    assert Path(outs["netmhcpan"]).exists()
    assert Path(outs["raw_peptides"]).exists()

def test_run_tool_stub_netmhcpan(tmp_path):
    pep = ROOT / "data/fixtures/pvacseq_aggregated.tsv"
    from neoag.adapters.pvactools_parser import parse_pvactools_outputs
    raw = tmp_path / "raw_peptides.tsv"
    parse_pvactools_outputs([pep], "S1", "default", None, raw)
    ctx = RunContext(sample_id="S1", outdir=tmp_path, stub=True, raw_peptides=raw, hla_alleles=["HLA-A*02:01"])
    out = run_tool("netmhcpan", ctx, tmp_path / "net.xls")
    assert out.exists()

def test_run_full_cli_stub(tmp_path):
    outdir = tmp_path / "full"
    main(["run-full", "--config", str(ROOT / "conf/run.stub.toml"), "--outdir", str(outdir)])
    assert (outdir / "scoring/ranked_peptides.tsv").exists()
    assert (outdir / "scoring/comprehensive_peptide_evidence.tsv").exists()


def test_run_full_cli_can_generate_technical_report_only(tmp_path):
    outdir = tmp_path / "technical_only"
    main([
        "run-full", "--config", str(ROOT / "conf/run.stub.toml"),
        "--outdir", str(outdir), "--reports", "technical",
    ])
    assert (outdir / "reports/evidence_report.technical.html").is_file()
    assert not (outdir / "reports/evidence_report.patient.html").exists()


def test_run_full_cli_marks_explicit_patient_report_as_pipeline_snapshot(tmp_path):
    outdir = tmp_path / "patient_snapshot"
    main([
        "run-full", "--config", str(ROOT / "conf/run.stub.toml"),
        "--outdir", str(outdir), "--reports", "patient",
    ])
    patient = outdir / "reports/evidence_report.patient.html"
    assert patient.is_file()
    assert "Pipeline 运行阶段的结果快照" in patient.read_text(encoding="utf-8")


def test_multisource_peptides_all_enter_presentation_prediction(tmp_path):
    cfg = tmp_path / "multisource.toml"
    cfg.write_text(
        f'''[sample]
id = "MULTISOURCE"
profile = "default"

[tools]
stub = true
enabled = ["netmhcpan", "mhcflurry"]

[inputs]
entry_mode = "e2e"
hla_alleles = ["HLA-A*02:01"]
pvac_files = [
  "{ROOT / 'data/fixtures/pvacseq_aggregated.tsv'}",
  "{ROOT / 'data/fixtures/pvacfuse_aggregated.tsv'}",
  "{ROOT / 'data/fixtures/pvacsplice_aggregated.tsv'}",
]
expected_peptide_sources = ["pVACseq", "pVACfuse", "pVACsplice"]
required_presentation_predictors = ["netmhcpan", "mhcflurry"]
extract_appm_from_vcf = false
''',
        encoding="utf-8",
    )

    outputs = run_upstream(cfg, tmp_path / "upstream")
    peptides = read_tsv(outputs["raw_peptides"])

    assert {row["source_tool"] for row in peptides} == {"pVACseq", "pVACfuse", "pVACsplice"}
    assert outputs["peptide_sources"] == "pVACfuse,pVACseq,pVACsplice"
    coverage = read_tsv(outputs["peptide_source_coverage"])[0]
    assert coverage["status"] == "COMPLETE"
    assert Path(outputs["netmhcpan"]).is_file()
    assert Path(outputs["mhcflurry"]).is_file()


def test_multisource_missing_fusion_reports_low_confidence(tmp_path, capsys):
    cfg = tmp_path / "missing_fusion.toml"
    cfg.write_text(
        f'''[sample]
id = "INCOMPLETE"
profile = "default"

[tools]
stub = true
enabled = []

[inputs]
entry_mode = "e2e"
pvac_files = ["{ROOT / 'data/fixtures/pvacseq_aggregated.tsv'}"]
expected_peptide_sources = ["pVACseq", "pVACfuse"]
extract_appm_from_vcf = false
''',
        encoding="utf-8",
    )

    outputs = run_upstream(cfg, tmp_path / "upstream")

    assert outputs["peptide_source_completeness"] == "LOW_CONFIDENCE"
    assert outputs["missing_peptide_sources"] == "pVACfuse"
    assert "confidence is LOW" in capsys.readouterr().out
    coverage = read_tsv(outputs["peptide_source_coverage"])[0]
    assert coverage["missing_sources"] == "pVACfuse"


def test_intermediates_prefilter_runs_before_presentation_tools(tmp_path):
    raw_events = tmp_path / "raw_events.tsv"
    raw_peptides = tmp_path / "raw_peptides.tsv"
    header = (
        "event_id\tevent_type\tsplice_alignment_qc_status\tjunction_read_qc_status\t"
        "unique_split_reads\ttotal_split_reads\ttumor_psi\tsplice_annotation_status\t"
        "crosses_junction\n"
    )
    raw_events.write_text(
        header
        + "S_PASS\tSplice\tPASS\tPASS\t5\t12\t0.20\tUNANNOTATED\ttrue\n"
        + "S_FAIL\tSplice\tPASS\tFAIL\t5\t12\t0.20\tUNANNOTATED\ttrue\n",
        encoding="utf-8",
    )
    raw_peptides.write_text(
        "peptide_id\tevent_id\tevent_type\tpeptide\thla_allele\tsplice_alignment_qc_status\t"
        "junction_read_qc_status\tunique_split_reads\ttotal_split_reads\ttumor_psi\t"
        "splice_annotation_status\tcrosses_junction\n"
        "P_PASS\tS_PASS\tSplice\tACDEFGHIK\tHLA-A*02:01\tPASS\tPASS\t5\t12\t0.20\tUNANNOTATED\ttrue\n"
        "P_FAIL\tS_FAIL\tSplice\tLMNPQRSTV\tHLA-A*02:01\tPASS\tFAIL\t5\t12\t0.20\tUNANNOTATED\ttrue\n",
        encoding="utf-8",
    )
    cfg = tmp_path / "intermediates.toml"
    cfg.write_text(
        f'''[sample]
id = "PREFILTER"
profile = "default"

[tools]
enabled = []

[inputs]
entry_mode = "intermediates"
raw_events = "{raw_events}"
raw_peptides = "{raw_peptides}"
''',
        encoding="utf-8",
    )

    outputs = run_upstream(cfg, tmp_path / "upstream")
    retained = read_tsv(outputs["raw_peptides"])
    assert [row["event_id"] for row in retained] == ["S_PASS"]
    funnel = {row["stage"]: row for row in read_tsv(outputs["splice_prefilter_funnel"])}
    assert funnel["FORMAL_JUNCTION_READ_QC"]["assessed_events"] == "2"
    assert funnel["FORMAL_JUNCTION_READ_QC"]["passed_events"] == "1"
    assert funnel["FORMAL_JUNCTION_READ_QC"]["failed_events"] == "1"


def test_upstream_prefers_purity_recommendation_over_facets(tmp_path):
    recommendation = tmp_path / "purity_recommendation.json"
    recommendation.write_text(json.dumps({
        "status": "CONCORDANT", "recommended_purity": 0.67, "range": "0.6500-0.6900",
        "n_tools": 2, "tool_values": {"FACETS": 0.65, "PURPLE": 0.69},
    }), encoding="utf-8")
    cfg = tmp_path / "purity.toml"
    cfg.write_text(
        f'''[sample]
id = "PURITY_CONSENSUS"
profile = "default"

[tools]
stub = true
enabled = ["facets"]

[inputs]
purity_recommendation = "{recommendation}"
extract_appm_from_vcf = false
''',
        encoding="utf-8",
    )

    outputs = run_upstream(cfg, tmp_path / "upstream")

    assert outputs["purity"] == str(recommendation)
    assert outputs["purity_recommendation"] == str(recommendation)
    assert outputs["facets_purity"].endswith("facets_purity.tsv")
    assert Path(outputs["facets_purity"]).is_file()

def test_fusion_tools_in_registry():
    from neoag.tools.registry import TOOL_REGISTRY

    for name in ("star_fusion", "arriba", "fusioncatcher", "easyfuse"):
        assert name in TOOL_REGISTRY
        assert TOOL_REGISTRY[name].category == "fusion"


def test_check_tools_runs():
    st = check_tool("netmhcpan")
    assert st.name == "netmhcpan"


def test_deepimmuno_external_uses_configured_python(tmp_path, monkeypatch):
    from neoag.adapters import deepimmuno as adapter
    from neoag.tools.runner import _run_deepimmuno_external

    configured_python = "/opt/neoag/envs/neoag-tools/bin/python"
    captured = {}
    monkeypatch.setenv("DEEPIMMUNO_PYTHON", configured_python)
    monkeypatch.setattr(adapter, "resolve_deepimmuno_dir", lambda _custom=None: tmp_path)

    def fake_batch(pairs, deep_dir, sample_id, python_exe=None):
        captured["python"] = python_exe
        return []

    monkeypatch.setattr(adapter, "run_deepimmuno_batch", fake_batch)
    monkeypatch.setattr(adapter, "write_deepimmuno_evidence", lambda path, rows: path.write_text("", encoding="utf-8"))
    context = RunContext(sample_id="S1", outdir=tmp_path, executables={"deepimmuno_dir": str(tmp_path)})
    output = tmp_path / "deepimmuno.tsv"

    _run_deepimmuno_external([("ACDEFGHIK", "HLA-A*02:01")], output, context)

    assert captured["python"] == configured_python
