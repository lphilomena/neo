
from pathlib import Path
from neoag_v03.immune_escape import build_immune_escape_evidence
from neoag_v03.utils import read_tsv


def test_immune_escape_v042_affected_counts(tmp_path):
    raw = tmp_path / "raw_peptides.tsv"
    raw.write_text("peptide_id\tevent_id\tpeptide\thla_allele\tmhc_class\nP1\tE1\tAAAAAAAAA\tHLA-A*02:01\tI\nP2\tE2\tBBBBBBBBB\tHLA-B*07:02\tI\n", encoding="utf-8")
    ranked = tmp_path / "ranked.tsv"
    ranked.write_text("peptide_id\tfinal_priority\nP1\tA\nP2\tC\n", encoding="utf-8")
    hla = tmp_path / "hla_loh.tsv"
    hla.write_text("hla_allele\tloh_status\nHLA-A*02:01\tloh\n", encoding="utf-8")
    paths = build_immune_escape_evidence(sample_id="S1", raw_peptides=raw, hla_loh_tsv=hla, ranked_peptides=ranked, outdir=tmp_path / "ie")
    events = read_tsv(paths["immune_escape_events"])
    ev = next(r for r in events if r["mechanism"] == "HLA_ALLELE_SPECIFIC_LOSS")
    assert ev["affected_candidate_count"] == "1"
    assert ev["affected_top_candidate_count"] == "1"
    summary = read_tsv(paths["immune_escape_summary"])[0]
    assert summary["n_peptides_affected_by_hla_loh"] == "1"
