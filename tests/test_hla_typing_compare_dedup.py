from pathlib import Path

from neoag.agent_skills.hla_typing_compare import collect_typing, main


def test_one_vote_per_tool_prefers_primary_result(tmp_path: Path) -> None:
    optitype = tmp_path / "optitype"
    spechla = tmp_path / "spechla"
    consensus = tmp_path / "consensus"
    optitype.mkdir()
    spechla.mkdir()
    consensus.mkdir()

    (optitype / "sample_result.tsv").write_text(
        "A1\tA2\tB1\tB2\tC1\tC2\nA*02:07\tA*11:01\tB*40:01\tB*46:01\tC*01:02\tC*03:03\n"
    )
    (spechla / "hla.result.txt").write_text(
        "Sample\tHLA_A_1\tHLA_A_2\nS1\tA*11:01:01:01\tA*02:07:01:01\n"
    )
    (spechla / "hla.result.details.txt").write_text(
        "Sample\tHLA_A_1\tHLA_A_2\nS1\tA*11:01:01:01\tA*11:01\n"
    )
    (consensus / "hla_typing_consensus.tsv").write_text(
        "locus\tconsensus_2field\nA\tA*02:07 / A*11:01\n"
    )

    rows = collect_typing([optitype, spechla])
    a_rows = [row for row in rows if row["locus"] == "A"]
    assert len(a_rows) == 2
    assert {row["tool"] for row in a_rows} == {"OptiType", "SpecHLA"}
    spec = next(row for row in a_rows if row["tool"] == "SpecHLA")
    assert Path(spec["source_file"]).name == "hla.result.txt"

    assert main([
        "--result-dir", str(optitype),
        "--result-dir", str(spechla),
        "--sample-id", "S1",
        "--outdir", str(consensus),
    ]) == 0
    text = (consensus / "hla_typing_consensus.tsv").read_text()
    assert "2/2\t" in text
    assert "\tSINGLE_TOOL\t" in text
    assert "Unknown=" not in text
