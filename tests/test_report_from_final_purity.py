from pathlib import Path

from neoag.report_from_final import _declared_purity_results, _input_files, _purity_records
from neoag.utils import write_tsv


def test_purity_records_keep_declared_tool_without_estimate(tmp_path: Path):
    facets = tmp_path / "facets"
    purple = tmp_path / "purple"
    ascat = tmp_path / "ascat"
    final = tmp_path / "production" / "final"
    for path in (facets, purple, ascat, final):
        path.mkdir(parents=True)
    (facets / "purity.tsv").write_text("sample_id\tpurity\nevent\t0.229894006569338\n")
    (facets / "facets_omni2p5_summary.tsv").write_text("metric\tvalue\nstatus\tfinished\n")
    (purple / "sample.purple.purity.tsv").write_text("purity\tploidy\n0.31\t4.7\n")
    (purple / "sample.purple.qc").write_text("QCStatus\tPASS\n")
    (ascat / "ascat_summary.tsv").write_text("sample_id\tpurity\tploidy\nevent\tNA\tNA\n")
    manifest = {
        "stages": {
            "purity_facets": {"outputs": {"facets_result": str(facets)}},
            "purity_purple": {"outputs": {"purple_result": str(purple)}},
            "purity_ascat": {"outputs": {"ascat_result": str(ascat)}},
        }
    }

    rows, consensus = _purity_records(final, manifest, {})

    by_tool = {row["tool"]: row for row in rows}
    assert by_tool["FACETS"]["status"] == "FINISHED"
    assert by_tool["PURPLE"]["status"] == "PASS"
    assert by_tool["ASCAT"]["status"] == "NO_VALID_ESTIMATE"
    assert consensus["status"] == "MODERATE_DISCORDANCE"
    assert consensus["recommended_purity"] == "0.2699"
    assert "FACETS=0.2299" in consensus["basis"]
    assert "PURPLE=0.3100" in consensus["basis"]


def test_input_files_preserve_reference_assets_for_report_enrichment(tmp_path: Path):
    final = tmp_path / "production" / "final"
    final.mkdir(parents=True)
    generated = {
        "inputs": {
            "reference_proteome": "/assets/data/normal/proteome/gencode.fa",
            "gencode_gtf": "/assets/data/rna/gencode_v49/gencode.v49.annotation.gtf.gz",
        }
    }

    files = _input_files(final, {}, generated)

    assert files["reference_proteome"].endswith("gencode.fa")
    assert files["gencode_gtf"].endswith("gencode.v49.annotation.gtf.gz")


def test_purity_report_unions_external_consensus_and_declared_stage_results(tmp_path: Path):
    production = tmp_path / "production"
    final_dir = production / "final"
    final_dir.mkdir(parents=True)

    upstream = tmp_path / "upstream" / "purity" / "consensus"
    upstream.mkdir(parents=True)
    recommendation = upstream / "recommended_purity.tsv"
    write_tsv(recommendation, [{"purity": "0.9178"}])
    write_tsv(upstream / "purity_cnv_tool_summary.tsv", [
        {"tool": "FACETS", "status": "FOUND", "purity": "0.915562", "ploidy": "2.124256"},
        {"tool": "PURPLE", "status": "FOUND", "purity": "0.92", "ploidy": "2.12"},
    ])

    sequenza = tmp_path / "stage_outputs" / "sequenza"
    sequenza.mkdir(parents=True)
    write_tsv(sequenza / "sequenza_cellularity.tsv", [
        {"cellularity": "0.31", "ploidy": "2.5"},
    ])
    manifest = {
        "run": {"outdir": str(production)},
        "evidence": {"purity": str(recommendation)},
        "stages": {
            "purity_sequenza": {"outputs": {"sequenza_result": str(sequenza)}},
        },
    }

    tools, consensus = _purity_records(final_dir, manifest, {})

    assert [row["tool"] for row in tools] == ["FACETS", "PURPLE", "SEQUENZA"]
    assert [row["purity"] for row in tools] == ["0.915562", "0.92", "0.31"]
    assert consensus["status"] == "CONCORDANT_WITH_OUTLIER"
    assert consensus["recommended_purity"] == "0.9178"
    assert consensus["recommended_ploidy"] == "2.12213"
    assert consensus["selected_tool"] == "多工具中位数（排除离群值）"
    assert "保留全部结果" in consensus["basis"]


def test_declared_purity_stages_are_discovered_without_fixed_tool_list(tmp_path: Path):
    custom = tmp_path / "custom-tool"
    manifest = {
        "stages": {
            "purity_newcaller": {"outputs": {"newcaller_result": str(custom)}},
            "hla_typing": {"outputs": {"typing_result": str(tmp_path / "hla")}},
        },
    }

    assert _declared_purity_results(manifest) == {"NEWCALLER": custom}
