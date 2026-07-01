from pathlib import Path
from neoag_v03.appm_lite import build_appm_lite
from neoag_v03.safety import load_normal_expression, apply_event_safety
from neoag_v03.config import load_profile

ROOT = Path(__file__).resolve().parents[1]


def test_appm_missing_expression_is_unassessed_not_low(tmp_path):
    profile = load_profile("sv_wgs_phase1")
    rows, summary = build_appm_lite("S1", None, None, None, profile, tmp_path)
    assert summary["appm_overall_status"] == "UNASSESSED"
    assert summary["mhc_i_integrity_score"] == "1.0000"
    assert all(r["expression_status"] == "unassessed" for r in rows)


def test_fusion_safety_uses_both_genes(tmp_path):
    normal = tmp_path / "normal.tsv"
    normal.write_text(
        "gene\tnormal_tissue_max_tpm\tnormal_hspc_tpm\tcritical_tissue_hit\n"
        "GENE1\t0.1\t0.1\tno\n"
        "GENE2\t12\t0.1\tyes\n",
        encoding="utf-8",
    )
    profile = load_profile("sv_wgs_phase1")
    e = {"gene": "GENE1::GENE2", "normal_tissue_max_tpm": "", "normal_hspc_tpm": "", "critical_tissue_hit": ""}
    out = apply_event_safety(e, profile, load_normal_expression(normal))
    assert out["safety_status"] == "FAIL"
    assert out["critical_tissue_hit"] == "yes"
