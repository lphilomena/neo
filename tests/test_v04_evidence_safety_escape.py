from pathlib import Path

from neoag.sv.phase1 import build_sv_phase1_raw
from neoag.utils import read_tsv, write_tsv
from neoag.peptide_safety_gate import build_peptide_safety_gate
from neoag.immune_escape import build_immune_escape_evidence

ROOT = Path(__file__).resolve().parents[1]
FX = ROOT / "data" / "fixtures_sv"


def test_wes_phase1_5_capture_annotation_and_cap(tmp_path):
    out = build_sv_phase1_raw(
        sample_id="SVMINI",
        sv_vcfs=[FX / "mini_sv.vcf"],
        callers=["GRIDSS2"],
        reference_fasta=FX / "mini_ref.fa",
        gencode_gtf=FX / "mini.gtf",
        hla=FX / "hla.txt",
        outdir=tmp_path,
        expression_tsv=FX / "expression.tsv",
        rna_junction_tsv=FX / "rna_junctions.tsv",
        normal_expression_tsv=FX / "normal_expression.tsv",
        normal_hla_ligands_tsv=FX / "normal_hla_ligands.tsv",
        wes_mode=True,
        capture_bed=ROOT / "data" / "fixtures_snv" / "wes_capture.bed",
    )
    sv_rows = read_tsv(out["sv_events_full"])
    raw_events = read_tsv(out["raw_events"])
    raw_peptides = read_tsv(out["raw_peptides"])
    assert sv_rows[0]["evidence_scope"] == "EXOME_CAPTURE_LIMITED"
    assert sv_rows[0]["breakend1_capture_status"] == "ON_TARGET"
    assert sv_rows[0]["wes_confidence_tier"] == "WES_Tier1"
    assert raw_events[0]["priority_cap"] == "B"
    assert raw_peptides[0]["evidence_scope"] == "EXOME_CAPTURE_LIMITED"



def test_wes_phase1_5_priority_cap_comes_from_profile(tmp_path):
    profile = tmp_path / "sv_wes_custom.toml"
    profile.write_text((ROOT / "profiles" / "sv_wes_phase1_5.toml").read_text().replace(
        'wes_tier1_priority_cap = "B"',
        'wes_tier1_priority_cap = "C"',
    ))
    out = build_sv_phase1_raw(
        sample_id="SVMINI",
        sv_vcfs=[FX / "mini_sv.vcf"],
        callers=["GRIDSS2"],
        reference_fasta=FX / "mini_ref.fa",
        gencode_gtf=FX / "mini.gtf",
        hla=FX / "hla.txt",
        outdir=tmp_path / "run",
        profile_name=str(profile),
        expression_tsv=FX / "expression.tsv",
        rna_junction_tsv=FX / "rna_junctions.tsv",
        normal_expression_tsv=FX / "normal_expression.tsv",
        normal_hla_ligands_tsv=FX / "normal_hla_ligands.tsv",
        wes_mode=True,
        capture_bed=ROOT / "data" / "fixtures_snv" / "wes_capture.bed",
    )
    raw_events = read_tsv(out["raw_events"])
    assert raw_events[0]["wes_confidence_tier"] == "WES_Tier1"
    assert raw_events[0]["priority_cap"] == "C"


def test_peptide_safety_reference_proteome_rejects_exact_match(tmp_path):
    events = tmp_path / "events.tsv"
    peptides = tmp_path / "peptides.tsv"
    ref = tmp_path / "prot.fa"
    safety = tmp_path / "peptide_safety.tsv"
    write_tsv(events, [{"event_id":"E1","sample_id":"S1","gene":"G1","event_type":"SNV","mutation_source":"SNV"}])
    write_tsv(peptides, [{"peptide_id":"P1","event_id":"E1","sample_id":"S1","gene":"G1","peptide":"AAAAAAAAA","hla_allele":"HLA-A*02:01","mhc_class":"I"}])
    ref.write_text(">prot1 gene_symbol:G1\nXXXAAAAAAAAAYYY\n")
    rows, _ = build_peptide_safety_gate(raw_events=events, raw_peptides=peptides, out_peptide_safety=safety, reference_proteome=ref)
    assert rows[0]["reference_proteome_exact_match"] == "yes"
    assert rows[0]["safety_status"] == "FAIL"
    assert "reference_proteome_exact_match" in rows[0]["safety_reason"]


def test_immune_escape_lost_hla_flags_restricting_peptide(tmp_path):
    peptides = tmp_path / "peptides.tsv"
    hla = tmp_path / "hla_loh.tsv"
    out = tmp_path / "escape"
    write_tsv(peptides, [{"peptide_id":"P1","event_id":"E1","sample_id":"S1","peptide":"AAAAAAAAA","hla_allele":"HLA-A*02:01","mhc_class":"I"}])
    write_tsv(hla, [{"hla_allele":"HLA-A*02:01","loh_status":"loh"}])
    paths = build_immune_escape_evidence(sample_id="S1", raw_peptides=peptides, outdir=out, hla_loh_tsv=hla)
    flags = read_tsv(paths["peptide_escape_flags"])
    assert flags[0]["restricting_hla_lost"] == "yes"
    assert flags[0]["escape_multiplier"] == "0.0000"
    assert flags[0]["priority_cap"] == "D"
