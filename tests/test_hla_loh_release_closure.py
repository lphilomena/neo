from pathlib import Path

from neoag.evidence_states import derive_hla_appm_state
from neoag.immune_escape import build_immune_escape_evidence
from neoag.utils import read_tsv, write_tsv


def test_single_tool_retention_cannot_promote_hla_appm():
    state = derive_hla_appm_state({
        "restricting_hla_lost": "no",
        "hla_loh_consensus_status": "SINGLE_TOOL_RETAINED",
        "appm_integrity_status": "MHC_I_INTACT",
        "appm_multiplier": "1.0",
    }, {})
    assert state["state"] == "HLA_LOH_UNASSESSED"
    assert not state["assessed"]


def test_discordant_restricting_allele_remains_conflict():
    state = derive_hla_appm_state({
        "restricting_hla_lost": "no",
        "hla_loh_consensus_status": "DISCORDANT",
        "appm_multiplier": "1.0",
    }, {})
    assert state["state"] == "CONFLICT"
    assert state["conflict"]


def test_immune_escape_projects_exact_allele_consensus(tmp_path: Path):
    peptides = tmp_path / "peptides.tsv"
    write_tsv(peptides, [
        {"peptide_id": "P1", "event_id": "E1", "peptide": "AAAAAAAAA", "hla_allele": "A*02:01", "mhc_class": "I"},
        {"peptide_id": "P2", "event_id": "E2", "peptide": "BBBBBBBBB", "hla_allele": "HLA-B*07:02", "mhc_class": "I"},
    ])
    loh = tmp_path / "hla_loh_consensus.tsv"
    write_tsv(loh, [
        {"hla_allele": "HLA-A*02:01", "consensus_status": "CONSENSUS_RETAINED", "crosscheck_status": "CONSENSUS_NO_LOH", "source_tools": "lohhla;spechla"},
        {"hla_allele": "HLA-B*07:02", "consensus_status": "SINGLE_TOOL_RETAINED", "crosscheck_status": "SINGLE_TOOL_NO_LOH", "source_tools": "spechla"},
    ])
    paths = build_immune_escape_evidence(
        sample_id="S1", raw_peptides=peptides, outdir=tmp_path / "out", hla_loh_tsv=loh,
    )
    flags = {row["peptide_id"]: row for row in read_tsv(paths["peptide_escape_flags"])}
    assert flags["P1"]["hla_loh_consensus_status"] == "CONSENSUS_RETAINED"
    assert flags["P1"]["escape_status"] == "ESCAPE_PASS"
    assert flags["P2"]["hla_loh_consensus_status"] == "SINGLE_TOOL_RETAINED"
    assert flags["P2"]["escape_status"] == "HLA_LOH_INCOMPLETE"
    assert flags["P2"]["priority_cap"] == "C_CAUTION"
