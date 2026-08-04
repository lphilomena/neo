from __future__ import annotations

import os
from pathlib import Path

import pytest

from neoag.adapters.facets import parse_facets_cncf, parse_facets_purity
from neoag.adapters.lohhla import parse_lohhla_prediction
from neoag.tools.postprocess import facets_to_cnv_tsv, facets_to_purity_tsv, lohhla_to_hla_loh_tsv
from neoag.tools.upstream import load_run_config, run_upstream

ROOT = Path(__file__).resolve().parents[1]
TOOLS_ROOT = Path(os.environ.get("NEOAG_TOOLS_ROOT", ROOT))
LOHHLA_EXAMPLE = (
    TOOLS_ROOT
    / "tools/lohhla/example-file/correct-example-out/example.10.DNA.HLAlossPrediction_CI.xls"
)


@pytest.mark.skipif(not LOHHLA_EXAMPLE.is_file(), reason="LOHHLA example missing")
def test_parse_lohhla_example():
    rows = parse_lohhla_prediction(LOHHLA_EXAMPLE)
    by_allele = {r["hla_allele"]: r["loh_status"] for r in rows}
    assert "HLA-A*01:01" in by_allele
    assert "HLA-A*24:02" in by_allele
    assert by_allele["HLA-A*24:02"] == "loh"
    assert by_allele["HLA-A*01:01"] == "no"


@pytest.mark.skipif(not LOHHLA_EXAMPLE.is_file(), reason="LOHHLA example missing")
def test_lohhla_to_hla_loh_tsv(tmp_path):
    out = tmp_path / "hla_loh.tsv"
    lohhla_to_hla_loh_tsv(LOHHLA_EXAMPLE, out)
    text = out.read_text(encoding="utf-8")
    assert "HLA-A*24:02\tloh" in text
    assert "HLA-A*01:01\tno" in text


def test_parse_facets_purity_and_cncf(tmp_path):
    purity_txt = tmp_path / "purity.txt"
    purity_txt.write_text("purity\n0.72\n", encoding="utf-8")
    rows = parse_facets_purity(purity_txt, sample_id="S1")
    assert rows == [{"sample_id": "S1", "purity": "0.72"}]

    cncf = tmp_path / "cncf.tsv"
    cncf.write_text(
        "chrom\tstart\tend\ttcn.em\n"
        "1\t1000\t2000\t1.5\n"
        "2\t3000\t4000\t0.0\n",
        encoding="utf-8",
    )
    segs = parse_facets_cncf(cncf)
    assert len(segs) == 2
    assert segs[0]["chrom"] == "chr1"
    assert segs[0]["total_cn"] == "1.5000"


def test_lohhla_requires_copy_number_and_pvalue_and_preserves_raw_evidence(tmp_path):
    prediction = tmp_path / "sample.HLAlossPrediction_CI.tsv"
    prediction.write_text(
        "HLA_A_type1\tHLA_A_type2\tHLA_type1copyNum_withBAF\tHLA_type1copyNum_withBAF_lower\t"
        "HLA_type1copyNum_withBAF_upper\tHLA_type2copyNum_withBAF\tHLA_type2copyNum_withBAF_lower\t"
        "HLA_type2copyNum_withBAF_upper\tPVal\tLossAllele\tKeptAllele\tnumMisMatchSitesCov\tpropSupportiveSites\n"
        "hla_a_02_06_01\thla_a_30_01_01\t0.85\t0.75\t0.95\t0.25\t0.20\t0.30\t1e-8\t"
        "hla_a_30_01_01\thla_a_02_06_01\t142\t89.5\n"
        "hla_b_13_02_01\thla_b_48_01_01\t0.20\t0.10\t0.30\t0.70\t0.60\t0.80\t0.20\t"
        "hla_b_13_02_01\thla_b_48_01_01\t60\t80\n",
        encoding="utf-8",
    )
    rows = {row["hla_allele"]: row for row in parse_lohhla_prediction(prediction)}
    assert rows["HLA-A*30:01"]["loh_status"] == "loh"
    assert rows["HLA-A*30:01"]["lohhla_pval"] == "1e-8"
    assert rows["HLA-A*30:01"]["lohhla_copy_number_with_baf"] == "0.25"
    assert rows["HLA-B*13:02"]["loh_status"] == "no"


def test_facets_converters(tmp_path):
    purity_txt = tmp_path / "facets_purity.txt"
    purity_txt.write_text("purity\n0.65\n", encoding="utf-8")
    purity_out = tmp_path / "purity.tsv"
    facets_to_purity_tsv(purity_txt, "HCC1395", purity_out)
    assert "HCC1395" in purity_out.read_text(encoding="utf-8")

    cncf = tmp_path / "facets_cncf.tsv"
    cncf.write_text("chrom\tstart\tend\ttotal_cn\nchr17\t100\t200\t2\n", encoding="utf-8")
    cnv_out = tmp_path / "cnv_segments.tsv"
    facets_to_cnv_tsv(cncf, cnv_out)
    assert "chr17" in cnv_out.read_text(encoding="utf-8")


@pytest.mark.skipif(not LOHHLA_EXAMPLE.is_file(), reason="LOHHLA example missing")
def test_upstream_lohhla_prediction_passthrough(tmp_path):
    cfg_path = tmp_path / "run.toml"
    cfg_path.write_text(
        f"""
[sample]
id = "T1"
profile = "default"

[tools]
stub = false
enabled = []

[inputs]
lohhla_prediction = "{LOHHLA_EXAMPLE.as_posix()}"
""",
        encoding="utf-8",
    )
    outs = run_upstream(cfg_path, tmp_path / "upstream")
    assert "hla_loh" in outs
    text = Path(outs["hla_loh"]).read_text(encoding="utf-8")
    assert "HLA-A*24:02\tloh" in text


@pytest.mark.skipif(not LOHHLA_EXAMPLE.is_file(), reason="LOHHLA example missing")
def test_tierb_config_loads():
    cfg = load_run_config(ROOT / "conf/run.hcc1395.tierb.toml")
    assert cfg["inputs"]["lohhla_prediction"]
