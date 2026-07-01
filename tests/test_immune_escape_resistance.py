from pathlib import Path

from neoag_v03.adapters.pvactools_parser import parse_pvactools_outputs
from neoag_v03.config import load_profile
from neoag_v03.immune_escape_resistance import build_immune_escape_resistance
from neoag_v03.scoring_v03 import score_v03
from neoag_v03.utils import read_tsv, write_tsv

ROOT = Path(__file__).resolve().parents[1]


def test_b2m_biallelic_rejects_mhc_i_peptides(tmp_path):
    profile = load_profile("default")
    _, peptides = parse_pvactools_outputs([ROOT / "data/fixtures/pvacseq_aggregated.tsv"], "S1", "default")
    pep_path = tmp_path / "peptides.tsv"
    write_tsv(pep_path, peptides)

    outdir = tmp_path / "immune_escape"
    _, summary, flags = build_immune_escape_resistance(
        "S1",
        pep_path,
        profile,
        outdir,
        vep_tsv=ROOT / "data/fixtures/vep_immune_escape.tsv",
        expression_tsv=ROOT / "data/fixtures/gene_expression.tsv",
        hla_loh_tsv=ROOT / "data/fixtures/hla_loh_lost.tsv",
        cnv_tsv=ROOT / "data/fixtures/cnv_segments.tsv",
    )

    assert summary["resistance_risk"] in {"HIGH", "MEDIUM"}
    assert summary["b2m_gene_status"] == "biallelic_loss"
    assert summary["jak_pathway_status"] == "biallelic_loss"
    assert "HLA-A*02:01" in summary["hla_loh_alleles"]
    assert "patient_resistant" not in summary

    events = read_tsv(outdir / "immune_escape_events.tsv")
    assert any(e["gene"] == "B2M" and e["resistance_risk"] == "HIGH" for e in events)
    assert all(e["resistance_risk"] in {"HIGH", "MEDIUM", "LOW", "INCONCLUSIVE"} for e in events)

    mhc_i_flags = [f for f in flags if str(f.get("mhc_class", "")).upper() not in {"II", "MHC-II", "CLASSII"}]
    assert mhc_i_flags
    assert all(f["escape_flag"] == "yes" for f in mhc_i_flags)
    assert all(float(f["escape_multiplier"]) <= 0.25 for f in mhc_i_flags)


def test_hla_loh_caps_restricting_allele(tmp_path):
    profile = load_profile("default")
    peptides = [{
        "peptide_id": "P1",
        "sample_id": "S1",
        "peptide": "AAAAAAA",
        "hla_allele": "HLA-A*02:01",
        "mhc_class": "I",
    }]
    pep_path = tmp_path / "peptides.tsv"
    write_tsv(pep_path, peptides)

    _, summary, flags = build_immune_escape_resistance(
        "S1",
        pep_path,
        profile,
        tmp_path / "immune_escape",
        hla_loh_tsv=ROOT / "data/fixtures/hla_loh_lost.tsv",
    )

    assert flags[0]["escape_flag"] == "yes"
    assert "hla_allele_loh" in flags[0]["escape_reason"]
    assert float(flags[0]["escape_multiplier"]) < 1.0
    assert summary["resistance_risk"] in {"HIGH", "MEDIUM", "LOW", "INCONCLUSIVE"}


def test_scoring_uses_escape_flags(tmp_path):
    profile = load_profile("default")
    e, p = parse_pvactools_outputs([ROOT / "data/fixtures/pvacseq_aggregated.tsv"], "S1", "default")
    ev_path = tmp_path / "events.tsv"
    pep_path = tmp_path / "peptides.tsv"
    write_tsv(ev_path, e)
    write_tsv(pep_path, p)

    immune_dir = tmp_path / "immune_escape"
    build_immune_escape_resistance(
        "S1",
        pep_path,
        profile,
        immune_dir,
        vep_tsv=ROOT / "data/fixtures/vep_immune_escape.tsv",
        hla_loh_tsv=ROOT / "data/fixtures/hla_loh_lost.tsv",
    )

    pres = tmp_path / "presentation.tsv"
    write_tsv(pres, [{
        "peptide_id": p[0]["peptide_id"],
        "binding_evidence_score": "0.8",
        "presentation_evidence_score": "0.8",
        "presentation_evidence_grade": "A",
    }])

    out_pep = tmp_path / "ranked_peptides.tsv"
    score_v03(
        ev_path,
        pep_path,
        pres,
        immune_dir / "immune_escape_summary.tsv",
        None,
        None,
        None,
        profile,
        tmp_path / "ranked_events.tsv",
        out_pep,
        peptide_escape_flags_tsv=immune_dir / "peptide_escape_flags.tsv",
    )
    ranked = read_tsv(out_pep)
    assert ranked
    assert "resistance_risk" in ranked[0]
    assert "escape_flag" in ranked[0]
    assert "escape_multiplier" in ranked[0]
