from __future__ import annotations

from pathlib import Path

import pytest

from neoag.sv.phase1 import build_sv_phase1_raw
from neoag.sv.score_pipeline import run_sv_score
from neoag.utils import read_tsv

ROOT = Path(__file__).resolve().parents[1]
FX = ROOT / "data" / "fixtures_sv"


def test_sv_score_stub_binding(tmp_path):
    adapter = tmp_path / "adapter"
    build_sv_phase1_raw(
        sample_id="SVMINI",
        sv_vcfs=[FX / "mini_sv.vcf"],
        callers=["GRIDSS2"],
        reference_fasta=FX / "mini_ref.fa",
        gencode_gtf=FX / "mini.gtf",
        hla=FX / "hla.txt",
        outdir=adapter,
        expression_tsv=FX / "expression.tsv",
        rna_junction_tsv=FX / "rna_junctions.tsv",
        normal_expression_tsv=FX / "normal_expression.tsv",
        normal_hla_ligands_tsv=FX / "normal_hla_ligands.tsv",
    )
    out = run_sv_score(
        outdir=tmp_path / "scored",
        profile_name_or_path="sv_wgs_phase1",
        sample_id="SVMINI",
        raw_events=adapter / "parsed" / "raw_events.tsv",
        raw_peptides=adapter / "parsed" / "raw_peptides.tsv",
        expression=FX / "expression.tsv",
        normal_expression=FX / "normal_expression.tsv",
        normal_hla_ligands=FX / "normal_hla_ligands.tsv",
        binding_stub=True,
        immunogenicity_stub=True,
    )
    assert "ccf_2" in out
    assert Path(out["ccf_2"]).exists()
    assert Path(out["ccf_lite"]).exists()
    ccf_header = Path(out["ccf_2"]).read_text().splitlines()[0]
    assert "clonality_confidence" in ccf_header
    assert "ccf_resolution" in ccf_header
    peptides = read_tsv(out["ranked_peptides"])
    events = read_tsv(out["ranked_events"])
    assert len(events) >= 1
    assert len(peptides) >= 1
    assert peptides[0].get("efficacy_score")
    assert Path(out["evidence_report"]).is_file()
