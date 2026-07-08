
from pathlib import Path
from neoag_v03.appm_v2 import build_appm_2
from neoag_v03.utils import read_tsv


def test_appm_v042_confidence_and_submodules(tmp_path):
    vep = tmp_path / "vep.tsv"
    vep.write_text("gene\tconsequence\tevent_id\nB2M\tframeshift_variant\tEV_B2M\n", encoding="utf-8")
    cnv = tmp_path / "cnv.tsv"
    cnv.write_text("gene\tcopy_number_status\tloh_status\ttotal_cn\tmajor_cn\tminor_cn\tevent_id\nB2M\tdeep_deletion\tloh\t0\t0\t0\tEV_B2M\n", encoding="utf-8")
    expr = tmp_path / "expr.tsv"
    expr.write_text("gene\tTPM\nB2M\t0.0\nHLA-A\t12\nHLA-B\t9\nHLA-C\t8\n", encoding="utf-8")
    raw_pep = tmp_path / "raw_peptides.tsv"
    raw_pep.write_text("peptide_id\tevent_id\tpeptide\thla_allele\tmhc_class\nP1\tE1\tAAAAAAAAA\tHLA-A*02:01\tI\n", encoding="utf-8")
    paths = build_appm_2(sample_id="S1", outdir=tmp_path / "appm", vep_tsv=vep, expression_tsv=expr, cnv_tsv=cnv, raw_peptides=raw_pep)
    summary = read_tsv(paths["appm_summary"])[0]
    assert summary["appm_call_confidence"] in {"high", "medium"}
    assert "B2M" in summary["confidence_reason"] or "BIALLELIC" in summary["confidence_reason"]
    sub = read_tsv(paths["appm_submodule_scores"])
    core = next(r for r in sub if r["submodule"] == "MHC_I_CORE")
    assert float(core["score"]) <= 0.10
    assert core["appm_call_confidence"] in {"high", "medium"}
