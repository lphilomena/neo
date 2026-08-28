from pathlib import Path

from neoag.report_from_final import _input_files, _purity_records


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
