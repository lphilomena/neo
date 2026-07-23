from pathlib import Path
import json
import pytest

from neoag.evidence_consensus import (
    build_evidence_consensus,
    load_consensus_rules,
    rank_evidence_consensus,
    rank_by_evidence_consensus,
    score_evidence_consensus,
)
from neoag.cli import build_parser, main as cli_main
from neoag.evidence_states import (
    derive_all_states,
    derive_event_authenticity,
    derive_hla_appm_state,
    derive_mutant_specificity,
    derive_presentation_consensus,
    derive_rna_support,
    derive_safety_state,
    evidence_state,
)
from neoag.pareto import nondominated_fronts
from neoag.utils import read_tsv, write_tsv


def complete_row(peptide_id="P1"):
    return {
        "peptide_id": peptide_id,
        "event_id": "E" + peptide_id[1:],
        "event_type": "SNV",
        "efficacy_score": "0.5",
        "l3_hla_presentation_score": "0.9",
        "presentation_gate_status": "PASS",
        "netmhcpan_mt_rank_el": "0.2",
        "mhcflurry_presentation_score": "0.8",
        "l3_hla_binding_score": "0.8",
        "presentation_evidence_grade": "A",
        "l3_expression_score": "0.7",
        "expression_evidence_status": "GENE_AND_TRANSCRIPT_SUPPORTED",
        "rna_alt_reads": "8",
        "rna_vaf": "0.12",
        "l3_clonality_score": "0.6",
        "ccf_status": "clonal_like",
        "l3_tumor_specificity_score": "0.9",
        "cross_platform_status": "CROSS_PLATFORM_PASS_CONCORDANT",
        "mutant_specificity_gate_status": "PASS",
        "l3_normal_tissue_safety_score": "0.8",
        "safety_evidence_completeness": "COMPLETE",
        "safety_status": "PASS",
        "l3_apm_integrity_score": "0.7",
        "appm_evidence_completeness": "HIGH",
        "restricting_hla_lost": "false",
    }


def test_evidence_state_separates_assessed_failure_from_conflict():
    assert evidence_state("FAIL") == "SUPPORTED"
    assert evidence_state("CONFLICT") == "CONFLICT"
    assert evidence_state("UNASSESSED") == "MISSING"


def test_complete_and_missing_consensus_states():
    complete = score_evidence_consensus(complete_row())
    assert complete["evidence_completeness_score"] == "1.0000"
    assert complete["evidence_grade"] == "R1"

    missing_row = complete_row("P2")
    missing_row["rna_alt_reads"] = "NA"
    missing_row["rna_vaf"] = "NA"
    missing_row["ccf_status"] = "RNA_ONLY_UNRESOLVED"
    missing = score_evidence_consensus(missing_row)
    assert "rna_support" in missing["evidence_missing_layers"]
    assert "clonality" in missing["evidence_missing_layers"]


def test_parallel_ranking_preserves_input_and_adds_deterministic_rank(tmp_path: Path):
    source = tmp_path / "ranked.tsv"
    output = tmp_path / "consensus.tsv"
    strong = complete_row("P1")
    weak = complete_row("P2")
    weak["netmhcpan_mt_rank_el"] = "1.0"
    weak["mhcflurry_presentation_score"] = "NA"
    write_tsv(source, [weak, strong])
    original = source.read_bytes()
    result = rank_by_evidence_consensus(source, output)
    rows = read_tsv(output)
    assert source.read_bytes() == original
    assert result["legacy_ranking_modified"] is False
    assert rows[0]["peptide_id"] == "P1"
    assert rows[0]["evidence_rank"] == "1"
    assert "evidence_layer_states" in rows[0]
    assert Path(result["evidence_states"]).is_file()
    assert Path(result["ranked_events"]).is_file()
    assert Path(result["comparison"]).is_file()


def test_parallel_ranking_accepts_missing_legacy_score(tmp_path: Path):
    source = tmp_path / "ranked.tsv"
    output = tmp_path / "consensus.tsv"
    row = complete_row()
    row["efficacy_score"] = "NA"
    write_tsv(source, [row])
    rank_by_evidence_consensus(source, output)
    assert read_tsv(output)[0]["evidence_rank"] == "1"


def test_hard_failure_and_priority_cap_are_auditable(tmp_path: Path):
    source = tmp_path / "comprehensive.tsv"
    failed = complete_row("P1")
    failed["reference_proteome_exact_match"] = "true"
    capped = complete_row("P2")
    capped["priority_cap"] = "C_CAUTION"
    write_tsv(source, [failed, capped])
    result = build_evidence_consensus(source, tmp_path / "consensus")
    rows = {row["peptide_id"]: row for row in read_tsv(result["ranked_peptides"])}
    assert rows["P1"]["hard_failure"] == "yes"
    assert rows["P1"]["evidence_grade"] == "R4"
    assert rows["P1"]["hard_failure_codes"] == "HARD_REFERENCE_PROTEOME_MATCH"
    assert "HARD_REFERENCE_PROTEOME_MATCH" in rows["P1"]["hard_failure_reasons"]
    assert rows["P2"]["evidence_grade_uncapped"] == "R1"
    assert rows["P2"]["evidence_grade"] == "R3"
    assert "grade_capped=R1->R3" in rows["P2"]["consensus_trace"]


def test_pareto_is_within_track_and_comparison_explains_shift(tmp_path: Path):
    source = tmp_path / "comprehensive.tsv"
    weaker = complete_row("P1")
    stronger = complete_row("P2")
    weaker["rna_alt_reads"] = "3"
    weaker["rna_vaf"] = "0.02"
    fusion = complete_row("P3")
    fusion["event_type"] = "Fusion"
    write_tsv(source, [weaker, stronger, fusion])
    result = build_evidence_consensus(source, tmp_path / "consensus")
    rows = {row["peptide_id"]: row for row in read_tsv(result["ranked_peptides"])}
    assert rows["P2"]["pareto_front"] == "1"
    assert rows["P1"]["pareto_front"] == "2"
    assert rows["P1"]["evidence_track"] == "MISSENSE"
    assert rows["P3"]["evidence_track"] == "FUSION"
    comparison = {row["peptide_id"]: row for row in read_tsv(result["comparison"])}
    assert "grade=" in comparison["P2"]["difference_reason"]


def test_pareto_uses_track_specific_dimensions_within_same_grade(tmp_path: Path):
    source = tmp_path / "comprehensive.tsv"
    missense = complete_row("P1")
    missense["mutant_specificity_gate_status"] = "MARGINAL_MT_ADVANTAGE"
    frameshift = complete_row("P2")
    frameshift.update({
        "event_type": "frameshift_variant",
        "peptide_consequence": "frameshift novel_tail",
        "contains_novel_aa": "true",
        "mutant_specificity_gate_status": "UNASSESSED",
    })
    fusion = complete_row("P3")
    fusion.update({
        "event_type": "Fusion",
        "rna_junction_reads": "12",
        "rna_frame_status": "IN_FRAME",
        "normal_junction_status": "NOT_DETECTED",
        "mutant_specificity_gate_status": "UNASSESSED",
    })
    write_tsv(source, [missense, frameshift, fusion])
    result = build_evidence_consensus(source, tmp_path / "consensus")
    rows = {row["peptide_id"]: row for row in read_tsv(result["ranked_peptides"])}
    assert rows["P1"]["evidence_track"] == "MISSENSE"
    assert rows["P2"]["evidence_track"] == "FRAMESHIFT"
    assert rows["P3"]["evidence_track"] == "FUSION"
    assert "mutant_specificity_grade" in rows["P1"]["pareto_dimensions"]
    assert "novel_tail_evidence_grade" in rows["P2"]["pareto_dimensions"]
    assert "junction_reads_grade" in rows["P3"]["pareto_dimensions"]
    assert rows["P2"]["novel_tail_evidence_grade"] == "3"
    assert rows["P3"]["junction_reads_grade"] == "3"
    assert rows["P3"]["normal_junction_safety_grade"] == "3"
    assert all(row["pareto_front"] == "1" for row in rows.values())


def test_manual_review_has_its_own_pareto_track_without_changing_biology(tmp_path: Path):
    source = tmp_path / "comprehensive.tsv"
    driver = complete_row("P1")
    driver["gene"] = "KRAS"
    ordinary = complete_row("P2")
    ordinary["gene"] = "ACLY"
    write_tsv(source, [driver, ordinary])
    result = build_evidence_consensus(source, tmp_path / "consensus")
    rows = {row["peptide_id"]: row for row in read_tsv(result["ranked_peptides"])}
    assert rows["P1"]["biological_event_track"] == "MISSENSE"
    assert rows["P1"]["evidence_track"] == "MANUAL_REVIEW"
    assert rows["P2"]["evidence_track"] == "MISSENSE"
    assert rows["P1"]["evidence_grade"] == rows["P2"]["evidence_grade"]


@pytest.mark.parametrize(("event_type", "expected_track", "expected_dimension"), [
    ("missense_variant", "MISSENSE", "mutant_specificity_grade"),
    ("frameshift_variant", "FRAMESHIFT", "novel_tail_evidence_grade"),
    ("fusion", "FUSION", "junction_authenticity_grade"),
    ("splice_junction", "SPLICE", "normal_junction_safety_grade"),
    ("structural_variant", "DNA_SV", "novel_tail_evidence_grade"),
])
def test_all_biological_tracks_select_their_own_pareto_dimensions(
    event_type: str, expected_track: str, expected_dimension: str,
):
    row = complete_row()
    row["event_type"] = event_type
    result = score_evidence_consensus(row)
    assert result["biological_event_track"] == expected_track
    assert result["evidence_track"] == expected_track
    assert expected_dimension in result["pareto_dimensions"]


def test_tie_break_is_deterministic_and_auditable(tmp_path: Path):
    first = complete_row("P_B")
    second = complete_row("P_A")
    source_a = tmp_path / "a.tsv"
    source_b = tmp_path / "b.tsv"
    write_tsv(source_a, [first, second])
    write_tsv(source_b, [second, first])
    result_a = build_evidence_consensus(source_a, tmp_path / "out_a")
    result_b = build_evidence_consensus(source_b, tmp_path / "out_b")
    rows_a = read_tsv(result_a["ranked_peptides"])
    rows_b = read_tsv(result_b["ranked_peptides"])
    assert [row["peptide_id"] for row in rows_a] == ["P_A", "P_B"]
    assert [row["peptide_id"] for row in rows_b] == ["P_A", "P_B"]
    assert [row["evidence_rank_key"] for row in rows_a] == [
        row["evidence_rank_key"] for row in rows_b
    ]
    key = rows_a[0]["evidence_rank_key"]
    assert "SAFETY_PASS" in key
    assert "RNA_CONFIRMED" in key
    assert "PRESENTATION_CONSISTENT_STRONG" in key
    assert "MT_SPECIFIC" in key
    assert "CCF_HIGH_CONFIDENCE" in key
    assert key.endswith("|P_A")


def test_tie_break_uses_netmhcpan_before_mhcflurry(tmp_path: Path):
    weaker_el = complete_row("P1")
    weaker_el["netmhcpan_mt_rank_el"] = "0.4"
    weaker_el["mhcflurry_presentation_score"] = "0.95"
    stronger_el = complete_row("P2")
    stronger_el["netmhcpan_mt_rank_el"] = "0.2"
    stronger_el["mhcflurry_presentation_score"] = "0.70"
    source = tmp_path / "input.tsv"
    write_tsv(source, [weaker_el, stronger_el])
    result = build_evidence_consensus(source, tmp_path / "out")
    rows = read_tsv(result["ranked_peptides"])
    assert rows[0]["peptide_id"] == "P2"
    assert "NETMHCPAN_EL=0.2" in rows[0]["evidence_rank_key"]
    assert "MHCFLURRY=0.7" in rows[0]["evidence_rank_key"]


def test_event_output_keeps_best_per_event_hla_and_at_most_two_representatives(tmp_path: Path):
    rows = []
    for peptide_id, hla, redundancy in (
        ("P_A", "HLA-A*02:01", "window_1"),
        ("P_B", "HLA-A*02:01", "window_2"),
        ("P_C", "HLA-B*07:02", "window_1"),
        ("P_D", "HLA-C*07:02", "window_3"),
    ):
        row = complete_row(peptide_id)
        row.update({
            "event_id": "E_SHARED", "hla_allele": hla,
            "peptide": peptide_id + "PEPTIDE", "redundancy_group": redundancy,
        })
        rows.append(row)
    source = tmp_path / "input.tsv"
    write_tsv(source, rows)
    result = build_evidence_consensus(source, tmp_path / "out")
    events = read_tsv(result["ranked_events"])
    assert len(events) == 1
    event = events[0]
    assert event["peptide_count"] == "4"
    assert event["event_hla_candidate_count"] == "3"
    assert event["representative_count"] == "2"
    assert event["representative_1_peptide_id"] == "P_A"
    assert event["representative_2_peptide_id"] == "P_D"
    assert event["representative_1_hla_allele"] == "HLA-A*02:01"
    assert event["representative_2_hla_allele"] == "HLA-C*07:02"


def test_phase_group_merges_component_events_before_representative_selection(tmp_path: Path):
    first = complete_row("P_A")
    first.update({
        "event_id": "TBR1_E1", "phase_group_id": "TBR1_PHASED_CIS",
        "hla_allele": "HLA-A*02:01", "redundancy_group": "overlap_1",
    })
    second = complete_row("P_B")
    second.update({
        "event_id": "TBR1_E2", "phase_group_id": "TBR1_PHASED_CIS",
        "hla_allele": "HLA-B*07:02", "redundancy_group": "overlap_2",
    })
    source = tmp_path / "input.tsv"
    write_tsv(source, [first, second])
    result = build_evidence_consensus(source, tmp_path / "out")
    events = read_tsv(result["ranked_events"])
    assert len(events) == 1
    assert events[0]["event_group_id"] == "PHASE:TBR1_PHASED_CIS"
    assert events[0]["member_event_ids"] == "TBR1_E1,TBR1_E2"
    assert events[0]["member_event_count"] == "2"
    assert events[0]["representative_count"] == "2"


def test_state_derivers_emit_discrete_auditable_vector():
    states = derive_all_states(complete_row(), load_consensus_rules())
    assert set(states) == {
        "event_authenticity", "rna_support", "presentation_consensus",
        "mutant_specificity", "clonality", "hla_appm", "safety",
        "evidence_completeness",
    }
    assert all(0 <= state["grade"] <= 3 for state in states.values())
    assert all(state["reason"] for state in states.values())
    assert all(state["reason_code"] for state in states.values())
    assert states["rna_support"]["reason_code"] == "RNA_ALT_READS_SUPPORTED"
    assert states["presentation_consensus"]["reason_code"] == "PRESENTATION_CORE_TOOLS_CONCORDANT"


def test_pareto_collapses_identical_vectors_and_is_deterministic():
    rows = [
        {"_pareto_id": "a", "x": 3, "y": 3},
        {"_pareto_id": "b", "x": 3, "y": 3},
        {"_pareto_id": "c", "x": 2, "y": 3},
        {"_pareto_id": "d", "x": 3, "y": 2},
        {"_pareto_id": "e", "x": 1, "y": 1},
    ]
    first = nondominated_fronts(rows, ["x", "y"])
    second = nondominated_fronts(list(reversed(rows)), ["x", "y"])
    assert first == second
    assert first["a"] == first["b"] == 1
    assert first["c"] == first["d"] == 2
    assert first["e"] == 3


def test_recommended_main_api_writes_conflicts_and_provenance(tmp_path: Path):
    source = tmp_path / "comprehensive.tsv"
    row = complete_row()
    row["cross_platform_status"] = "DISCORDANT"
    write_tsv(source, [row])
    provenance = tmp_path / "provenance.json"
    result = rank_evidence_consensus(
        source,
        tmp_path / "peptides.tsv",
        tmp_path / "events.tsv",
        tmp_path / "states.tsv",
        tmp_path / "conflicts.tsv",
        load_consensus_rules(),
        provenance,
    )
    assert result["conflicts"] == 1
    assert read_tsv(tmp_path / "conflicts.tsv")[0]["layer"] == "event_authenticity"
    assert "PROVISIONAL_RESEARCH_ONLY" in provenance.read_text()


def test_first_phase_output_contract_and_aliases(tmp_path: Path):
    comprehensive = tmp_path / "comprehensive.tsv"
    weighted = tmp_path / "ranked_peptides.tsv"
    canonical_input = complete_row()
    canonical_input.update({"peptide": "AAAAAAAAA", "hla_allele": "HLA-A*02:01"})
    write_tsv(comprehensive, [canonical_input])
    write_tsv(weighted, [{"peptide_id": "P1", "efficacy_score": "0.5"}])
    outdir = tmp_path / "scoring"
    result = build_evidence_consensus(
        comprehensive, outdir, weighted_baseline_tsv=weighted,
    )
    expected = {
        "ranked_peptides.evidence_consensus.tsv",
        "ranked_peptides.tsv",
        "ranked_events.evidence_consensus.tsv",
        "evidence_states.tsv",
        "evidence_conflicts.tsv",
        "evidence_consensus_summary.tsv",
        "ranking_compare_weighted_vs_consensus.tsv",
        "ranking_compare_weighted_vs_consensus.md",
        "evidence_consensus_run.json",
        "all_tool_results.tsv",
        "all_tool_results.manifest.json",
        "ranked_peptides.weighted_baseline.tsv",
    }
    assert expected <= {path.name for path in outdir.iterdir()}
    canonical = read_tsv(result["all_tool_results"])[0]
    assert canonical["all_tool_results_schema_version"] == "1.0"
    assert canonical["canonical_record_type"] == "PEPTIDE_HLA_EVIDENCE"
    assert len(canonical["canonical_record_id"]) == 24
    canonical_manifest = json.loads(Path(result["all_tool_results_manifest"]).read_text())
    assert canonical_manifest["canonical"] is True
    assert canonical_manifest["input"]["sha256"]
    assert canonical_manifest["output"]["sha256"]
    assert canonical_manifest["output"]["rows"] == 1
    assert canonical_manifest["validation"]["status"] == "PASS"
    assert canonical_manifest["validation"]["duplicate_record_ids"] == 0
    assert Path(result["ranked_peptides_compat"]).read_bytes() == weighted.read_bytes()
    assert Path(result["weighted_baseline"]).read_bytes() == weighted.read_bytes()
    assert Path(result["comparison_legacy"]).read_bytes() == Path(result["comparison"]).read_bytes()
    manifest = json.loads(Path(result["run_manifest"]).read_text())
    assert manifest["counts"]["peptide_hla_rows"] == 1
    assert manifest["legacy_ranking_modified"] is False
    assert manifest["outputs"]["weighted_baseline"]["sha256"]


def test_evidence_rank_cli_defaults_to_protected_parallel_mode(tmp_path: Path):
    comprehensive = tmp_path / "comprehensive.tsv"
    weighted = tmp_path / "ranked_peptides.tsv"
    provenance = tmp_path / "provenance.json"
    outdir = tmp_path / "consensus"
    write_tsv(comprehensive, [complete_row()])
    write_tsv(weighted, [{"peptide_id": "P1", "efficacy_score": "0.5"}])
    provenance.write_text("{}\n")
    cli_main([
        "evidence-rank",
        "--comprehensive-evidence", str(comprehensive),
        "--weighted-baseline", str(weighted),
        "--provenance", str(provenance),
        "--outdir", str(outdir),
        "--emit-event-ranking",
        "--compare-weighted",
        "--deterministic",
    ])
    assert (outdir / "ranked_peptides.evidence_consensus.tsv").is_file()
    assert (outdir / "ranked_events.evidence_consensus.tsv").is_file()
    assert weighted.read_text().startswith("peptide_id")
    assert "evidence_consensus" in json.loads(provenance.read_text())


def test_evidence_rank_cli_forbids_primary_replacement():
    parser = build_parser()
    defaults = parser.parse_args([
        "evidence-rank",
        "--comprehensive-evidence", "comprehensive.tsv",
        "--weighted-baseline", "ranked_peptides.tsv",
        "--outdir", "out",
    ])
    assert defaults.mode == "parallel"
    assert defaults.track == "all"
    assert defaults.emit_event_ranking is True
    assert defaults.compare_weighted is True
    assert defaults.deterministic is True
    with pytest.raises(SystemExit):
        parser.parse_args([
            "evidence-rank",
            "--comprehensive-evidence", "comprehensive.tsv",
            "--weighted-baseline", "ranked_peptides.tsv",
            "--outdir", "out",
            "--replace-primary-ranking",
        ])


def test_evidence_rank_cli_filters_track_using_normalized_state(tmp_path: Path):
    comprehensive = tmp_path / "comprehensive.tsv"
    weighted = tmp_path / "ranked_peptides.tsv"
    missense = complete_row("P1")
    fusion = complete_row("P2")
    fusion.update({"event_type": "Fusion", "rna_junction_reads": "12", "rna_frame_status": "IN_FRAME"})
    write_tsv(comprehensive, [missense, fusion])
    write_tsv(weighted, [{"peptide_id": "P1"}, {"peptide_id": "P2"}])
    outdir = tmp_path / "fusion_consensus"
    cli_main([
        "evidence-rank",
        "--comprehensive-evidence", str(comprehensive),
        "--weighted-baseline", str(weighted),
        "--outdir", str(outdir),
        "--track", "fusion",
    ])
    rows = read_tsv(outdir / "ranked_peptides.evidence_consensus.tsv")
    assert len(rows) == 1
    assert rows[0]["peptide_id"] == "P2"
    assert rows[0]["biological_event_track"] == "FUSION"


def test_source_precedence_conflicts_flow_into_final_conflict_table(tmp_path: Path):
    source = tmp_path / "comprehensive.tsv"
    row = complete_row()
    row["evidence_source_precedence_version"] = "1.0"
    row["evidence_conflict_fields"] = "ccf_estimate"
    row["evidence_conflict_details"] = json.dumps([{
        "field": "ccf_estimate", "selected_source": "ccf_2", "selected_value": "0.8",
        "other_source": "ranked_peptides", "other_value": "0.5",
        "precedence_version": "1.0", "conflict_type": "NONEMPTY_SOURCE_DISAGREEMENT",
    }])
    write_tsv(source, [row])
    result = build_evidence_consensus(source, tmp_path / "consensus")
    conflicts = read_tsv(result["evidence_conflicts"])
    assert conflicts[0]["layer"] == "source_precedence"
    assert conflicts[0]["selected_source"] == "ccf_2"
    states = read_tsv(result["evidence_states"])
    assert states[0]["manual_review_required"] == "yes"
    assert "source_precedence" in states[0]["evidence_conflict_layers"]


@pytest.mark.parametrize(("status", "expected"), [
    ("CROSS_PLATFORM_PASS_CONCORDANT", "EVENT_CONFIRMED"),
    ("ALT_PRESENT_BELOW_PASS_OR_CALLER_DIFFERENCE", "EVENT_STRONG"),
    ("COVERED_NO_ALT_SAMPLE_OR_ASSAY_DIFFERENCE", "EVENT_SAMPLE_SPECIFIC"),
    ("SOURCE_PASS_NOT_REPRODUCED_BY_PILEUP", "EVENT_CONFLICT"),
])
def test_event_authenticity_named_states(status, expected):
    row = complete_row()
    row["cross_platform_status"] = status
    assert derive_event_authenticity(row, load_consensus_rules())["state"] == expected


def test_rna_reuses_existing_status_and_gene_tpm_is_not_mutant_support():
    rules = load_consensus_rules()
    assert derive_rna_support({"event_type": "SNV", "rna_support_status": "RNA_ALT_SUPPORTED"}, rules)["state"] == "RNA_CONFIRMED"
    assert derive_rna_support({"event_type": "SNV", "rna_support_status": "RNA_ALT_NOT_DETECTED"}, rules)["state"] == "RNA_NEGATIVE"
    state = derive_rna_support({"event_type": "SNV", "gene_expression_tpm": "20"}, rules)
    assert state["state"] == "GENE_EXPRESSION_ONLY"
    assert state["grade"] == 1


def test_presentation_uses_core_groups_without_double_counting_immunogenicity_models():
    rules = load_consensus_rules()
    consistent = derive_presentation_consensus({
        "netmhcpan_mt_rank_el": "0.2", "mhcflurry_presentation_score": "0.8",
        "prime_rank": "1", "bigmhc_im_score": "0.9", "deepimmuno_score": "0.9",
    }, rules)
    discordant = derive_presentation_consensus({
        "netmhcpan_mt_rank_el": "0.2", "mhcflurry_presentation_score": "0.1",
    }, rules)
    single = derive_presentation_consensus({"netmhcpan_mt_rank_el": "0.2"}, rules)
    assert consistent["state"] == "PRESENTATION_CONSISTENT_STRONG"
    assert consistent["grade"] == 3
    assert discordant["state"] == "PRESENTATION_DISCORDANT"
    assert single["state"] == "PRESENTATION_SINGLE_TOOL"


@pytest.mark.parametrize(("specificity", "state", "hard"), [
    ("MT_SPECIFIC", "MT_SPECIFIC", False),
    ("MARGINAL_MT_ADVANTAGE", "MARGINAL_MT_ADVANTAGE", False),
    ("MT_WT_SIMILAR", "MT_WT_SIMILAR", False),
    ("WT_BETTER", "WT_BETTER", False),
    ("NON_MUTANT_SEQUENCE", "NON_MUTANT_SEQUENCE", True),
    ("UNASSESSED", "UNASSESSED", False),
])
def test_mutant_specificity_reuses_existing_states(specificity, state, hard):
    result = derive_mutant_specificity({"mutant_specificity_status": specificity}, {})
    assert result["state"] == state
    assert result["hard_fail"] is hard


def test_hla_and_safety_strong_rules_have_stable_hard_codes():
    hla = derive_hla_appm_state({"restricting_hla_lost": "true"}, {})
    safety = derive_safety_state({"reference_proteome_exact_match": "true"}, {})
    partial = derive_safety_state({
        "safety_status": "PASS", "safety_missing_layers": "normal_ligandome",
    }, {})
    event_partial = derive_safety_state({
        "safety_status": "PASS", "event_safety_status": "PASS",
        "event_safety_missing_layers": "normal_junction",
    }, {})
    event_reject = derive_safety_state({
        "safety_status": "PASS", "event_reference_proteome_exact_match": "true",
    }, {})
    assert hla["state"] == "RESTRICTING_HLA_LOST"
    assert hla["hard_code"] == "HARD_RESTRICTING_HLA_LOST"
    assert safety["hard_code"] == "HARD_REFERENCE_PROTEOME_MATCH"
    assert partial["state"] == "SAFETY_PARTIAL"
    assert event_partial["state"] == "SAFETY_PARTIAL"
    assert event_reject["hard_code"] == "HARD_REFERENCE_PROTEOME_MATCH"


def test_state_driven_caps_apply_before_pareto(tmp_path: Path):
    source = tmp_path / "comprehensive.tsv"
    row = complete_row()
    row["mhcflurry_presentation_score"] = "NA"
    write_tsv(source, [row])
    result = build_evidence_consensus(source, tmp_path / "consensus")
    ranked = read_tsv(result["ranked_peptides"])[0]
    assert ranked["evidence_grade_cap"] == "R3"
    assert "CAP_SINGLE_PRESENTATION_TOOL" in ranked["evidence_grade_cap_reasons"]
    assert ranked["evidence_grade"] == "R3"


def _rank_one(tmp_path: Path, row: dict[str, str]):
    source = tmp_path / f"{row['peptide_id']}.tsv"
    outdir = tmp_path / f"out_{row['peptide_id']}"
    write_tsv(source, [row])
    result = build_evidence_consensus(source, outdir)
    return read_tsv(result["ranked_peptides"])[0]


def test_r1_requires_all_core_evidence_and_has_first_batch_action(tmp_path: Path):
    row = _rank_one(tmp_path, complete_row("P_R1"))
    assert row["evidence_grade"] == "R1"
    assert row["consensus_action"] == "FIRST_BATCH_EXPERIMENTAL_PRIORITY"


def test_r2_allows_one_caution_without_using_legacy_letter_grade(tmp_path: Path):
    candidate = complete_row("P_R2")
    candidate["ccf_confidence"] = "low_confidence"
    row = _rank_one(tmp_path, candidate)
    assert row["evidence_grade"] == "R2"
    assert row["consensus_priority_cap"] == "R2"
    assert row["consensus_action"] == "ADVANCE_WITH_CAUTION"


def test_r3_prioritizes_evidence_completion_not_direct_elispot(tmp_path: Path):
    candidate = complete_row("P_R3")
    candidate["rna_alt_reads"] = "NA"
    candidate["rna_vaf"] = "NA"
    candidate["gene_expression_tpm"] = "15"
    row = _rank_one(tmp_path, candidate)
    assert row["evidence_grade"] == "R3"
    assert row["consensus_action"] == "EVIDENCE_COMPLETION_FIRST"
    assert "targeted RNA" in row["recommended_next_steps"]
    assert "ELISpot" not in row["recommended_next_steps"]


def test_r4_and_manual_review_do_not_upgrade_driver_event(tmp_path: Path):
    candidate = complete_row("P_R4")
    candidate["gene"] = "TP53"
    candidate["mutant_specificity_status"] = "WT_BETTER"
    candidate["mutant_specificity_gate_status"] = "WT_BETTER"
    row = _rank_one(tmp_path, candidate)
    assert row["evidence_grade"] == "R4"
    assert row["manual_review_required"] == "yes"
    assert row["consensus_action"] == "DO_NOT_ADVANCE"


def test_clear_novel_junction_can_satisfy_specificity_without_dna_ccf(tmp_path: Path):
    candidate = complete_row("P_NOVEL")
    candidate.update({
        "event_type": "Splice",
        "crosses_junction": "true",
        "contains_novel_aa": "true",
        "mutant_specificity_status": "UNASSESSED",
        "mutant_specificity_gate_status": "UNASSESSED",
        "rna_support_status": "RNA_JUNCTION_SUPPORTED",
        "rna_junction_reads": "20",
        "normal_junction_assessment_status": "NOT_DETECTED",
        "ccf_status": "RNA_ONLY_UNRESOLVED",
        "l3_clonality_score": "NA",
    })
    row = _rank_one(tmp_path, candidate)
    assert row["evidence_grade"] == "R1"


def test_rna_only_fusion_r3_recommends_orthogonal_confirmation(tmp_path: Path):
    candidate = complete_row("P_FUSION")
    candidate.update({
        "event_type": "Fusion",
        "mutation_source": "RNA_ONLY",
        "crosses_junction": "true",
        "contains_novel_aa": "true",
        "rna_support_status": "RNA_JUNCTION_SUPPORTED",
        "rna_junction_reads": "20",
        "normal_junction_assessment_status": "NOT_DETECTED",
        "ccf_status": "RNA_ONLY_UNRESOLVED",
        "l3_clonality_score": "NA",
    })
    row = _rank_one(tmp_path, candidate)
    assert row["evidence_grade"] == "R3"
    assert "CAP_RNA_ONLY_FUSION" in row["evidence_grade_cap_reasons"]
    assert "RT-PCR/Sanger" in row["recommended_next_steps"]
    assert "second fusion caller" in row["recommended_next_steps"]
