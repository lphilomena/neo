from __future__ import annotations

from pathlib import Path

from neoag.pipeline import run
from neoag.utils import read_tsv, write_tsv

ROOT = Path(__file__).resolve().parents[1]
PVAC = ROOT / "data" / "fixtures" / "pvacseq_aggregated.tsv"
NET = ROOT / "data" / "fixtures" / "netmhcpan_example.xls"
MHC = ROOT / "data" / "fixtures" / "mhcflurry_predictions.csv"
STAB = ROOT / "data" / "fixtures" / "netmhcstabpan_example.tsv"
NORM_EXPR = ROOT / "resources" / "normal_expression.example.tsv"
NORM_LIG = ROOT / "resources" / "normal_hla_ligands.example.tsv"
PURITY = ROOT / "data" / "fixtures" / "purity.tsv"
CNV_SEG = ROOT / "data" / "fixtures" / "cnv_segments.tsv"


def _run(tmp_path, *, vep_rows=None, expr_rows=None, cnv_rows=None, hla_rows=None, label="S1"):
    vep = tmp_path / "vep.tsv"
    expr = tmp_path / "expr.tsv"
    cnv = tmp_path / "cnv.tsv"
    hla = tmp_path / "hla_loh.tsv"
    if vep_rows is not None:
        write_tsv(vep, vep_rows)
    if expr_rows is not None:
        write_tsv(expr, expr_rows)
    if cnv_rows is not None:
        write_tsv(cnv, cnv_rows)
    if hla_rows is not None:
        write_tsv(hla, hla_rows)
    return run(
        tmp_path / "run",
        "default",
        label,
        pvac_paths=[PVAC],
        netmhcpan=NET,
        mhcflurry=MHC,
        netmhcstabpan=STAB,
        vep_appm=vep if vep_rows is not None else None,
        expression=expr if expr_rows is not None else None,
        hla_loh=hla if hla_rows is not None else None,
        purity=PURITY,
        cnv=cnv if cnv_rows is not None else CNV_SEG,
        normal_expression=NORM_EXPR,
        normal_hla_ligands=NORM_LIG,
        entry_mode="cross_module_consistency",
    )


def _rows(out):
    return {
        "appm_summary": read_tsv(out["appm_summary"])[0],
        "appm_gene": {r["gene"]: r for r in read_tsv(out["appm_gene_status"])},
        "appm_mod": read_tsv(out["appm_peptide_modifiers"]),
        "escape_summary": read_tsv(out["immune_escape_summary"])[0],
        "escape_flags": read_tsv(out["peptide_escape_flags"]),
        "ranked": read_tsv(out["ranked_peptides"]),
    }


def test_b2m_biallelic_loss_is_consistent_across_appm_escape_and_scoring(tmp_path):
    out = _run(
        tmp_path,
        vep_rows=[{"gene": "B2M", "consequence": "frameshift_variant"}],
        expr_rows=[{"gene": "B2M", "TPM": "0.01"}, {"gene": "HLA-A", "TPM": "12"}],
        cnv_rows=[{"gene": "B2M", "copy_number_status": "copy_neutral_loh", "total_cn": "1", "minor_cn": "0"}],
    )
    r = _rows(out)
    assert r["appm_gene"]["B2M"]["biallelic_status"] == "BIALLELIC_LOSS"
    assert r["appm_summary"]["mhc_i_integrity_status"] == "MHC_I_DEFECTIVE"
    assert r["escape_summary"]["b2m_biallelic_loss"] == "yes"
    assert r["escape_summary"]["mhc_i_escape_status"] == "HIGH"
    assert all(flag["hla_class_i_global_status"] == "GLOBAL_MHC_I_LOSS" for flag in r["escape_flags"] if flag["mhc_class"].upper() not in {"II", "MHC-II", "CLASSII"})
    assert any(row["appm_integrity_status"] == "MHC_I_DEFECTIVE" and row["final_priority"] == "D" for row in r["ranked"])


def test_jak_ifng_defect_is_consistent_across_appm_escape_and_scoring(tmp_path):
    out = _run(
        tmp_path,
        vep_rows=[{"gene": "JAK1", "consequence": "stop_gained"}],
        expr_rows=[{"gene": "JAK1", "TPM": "0.01"}, {"gene": "HLA-A", "TPM": "15"}],
        cnv_rows=[{"gene": "JAK1", "copy_number_status": "copy_loss", "total_cn": "1", "minor_cn": "0"}],
    )
    r = _rows(out)
    assert r["appm_gene"]["JAK1"]["biallelic_status"] == "BIALLELIC_LOSS"
    assert r["appm_summary"]["ifng_response_status"] == "IFNG_RESPONSE_DEFECTIVE"
    assert r["escape_summary"]["jak1_biallelic_loss"] == "yes"
    assert r["escape_summary"]["ifng_response_status"] == "HIGH"
    assert any("IFNG" in row["appm_multiplier_reason"] or "jak_stat_defect" in row["escape_reason"] for row in r["ranked"])
    assert any(row["final_priority"] in {"B_CAUTION", "C_CAUTION", "C", "D"} for row in r["ranked"])


def test_ciita_mhc_ii_defect_affects_mhc_ii_not_mhc_i(tmp_path):
    out = _run(
        tmp_path,
        vep_rows=[{"gene": "CIITA", "consequence": "frameshift_variant"}],
        expr_rows=[{"gene": "CIITA", "TPM": "0.01"}, {"gene": "HLA-A", "TPM": "15"}],
        cnv_rows=[{"gene": "CIITA", "copy_number_status": "copy_loss", "total_cn": "1", "minor_cn": "0"}],
    )
    r = _rows(out)
    assert r["appm_gene"]["CIITA"]["biallelic_status"] == "BIALLELIC_LOSS"
    assert r["appm_summary"]["mhc_ii_integrity_status"] == "MHC_II_DEFECTIVE"
    assert r["escape_summary"]["ciita_defect"] == "yes"
    assert r["escape_summary"]["mhc_ii_escape_status"] == "MEDIUM"
    assert all(row["appm_integrity_status"] != "MHC_II_DEFECTIVE" for row in r["ranked"] if row["mhc_class"].upper() not in {"II", "MHC-II", "CLASSII"})


def test_no_appm_input_is_unassessed_and_not_penalized(tmp_path):
    out = _run(tmp_path, vep_rows=None, expr_rows=None, cnv_rows=None, hla_rows=None)
    r = _rows(out)
    assert r["appm_summary"]["appm_overall_status"] == "UNASSESSED"
    assert r["appm_summary"]["appm_evidence_completeness"] == "UNASSESSED"
    assert r["escape_summary"]["overall_immune_escape_risk"] == "INCONCLUSIVE"
    assert all(row["appm_multiplier"] == "1.0000" for row in r["ranked"])
    assert all(row["appm_review_required"] == "yes" for row in r["ranked"])


def test_class_ii_hla_loh_stays_out_of_class_i_and_requires_review(tmp_path):
    out = _run(
        tmp_path,
        hla_rows=[
            {"hla_allele": "HLA-DQA1*02:01", "loh_status": "loh"},
            {"hla_allele": "HLA-DQB1*02:02", "loh_status": "loh"},
        ],
    )
    r = _rows(out)
    summary = r["appm_summary"]
    assert summary["hla_i_loh_flag"] == "no"
    assert summary["hla_i_loh_alleles"] == ""
    assert summary["hla_ii_loh_flag"] == "yes"
    assert set(summary["hla_ii_loh_alleles"].split(",")) == {"HLA-DQA1*02:01", "HLA-DQB1*02:02"}
    assert summary["mhc_i_integrity_score"] == "1.0000"
    assert summary["mhc_i_integrity_status"] == "MHC_I_INTACT"
    assert summary["mhc_ii_integrity_score"] == "0.6500"
    assert summary["mhc_ii_integrity_status"] == "MHC_II_CAUTION"

    submodules = {row["submodule"]: row for row in read_tsv(out["appm_submodule_scores"])}
    assert submodules["MHC_I_HLA_LOH"]["score"] == "1.0000"
    assert submodules["MHC_II_CORE"]["score"] == "0.6500"

    escape = r["escape_summary"]
    assert escape["lost_hla_i_alleles"] == ""
    assert set(escape["lost_hla_ii_alleles"].split(",")) == {"HLA-DQA1*02:01", "HLA-DQB1*02:02"}
    assert escape["n_peptides_affected_by_hla_loh"] == "0"
    assert escape["n_top_peptides_affected_by_hla_loh"] == "0"
    assert escape["overall_immune_escape_risk"] == "REVIEW_REQUIRED"
    assert all(flag["escape_status"] == "ESCAPE_PASS" for flag in r["escape_flags"])
