from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path


def write_tsv(path: Path, rows: list[dict[str, str]]):
    fieldnames = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def test_ranking_compare_skill(tmp_path: Path):
    rec = tmp_path / "ranked_peptides.recommendation.tsv"
    net = tmp_path / "ranked_peptides.netmhcpan42.tsv"
    rows_rec = [
        {"peptide_id": "p1", "gene": "G1", "peptide": "AAAA", "hla_allele": "HLA-A*02:01", "event_type": "SNV", "final_priority": "B", "recommended_use": "short peptide"},
        {"peptide_id": "p2", "gene": "G2", "peptide": "BBBB", "hla_allele": "HLA-A*02:01", "event_type": "Fusion", "final_priority": "C", "recommended_use": "minigene"},
        {"peptide_id": "p3", "gene": "G3", "peptide": "CCCC", "hla_allele": "HLA-B*07:02", "event_type": "InDel", "final_priority": "D", "recommended_use": "defer"},
    ]
    rows_net = [rows_rec[1], rows_rec[0], rows_rec[2]]
    write_tsv(rec, rows_rec)
    write_tsv(net, rows_net)
    out = tmp_path / "out"
    rc = subprocess.run([sys.executable, "-m", "neoag.agent_skills.ranking_compare", "--recommendation", str(rec), "--netmhcpan42", str(net), "--outdir", str(out)], text=True, capture_output=True)
    assert rc.returncode == 0, rc.stderr
    assert (out / "ranking_compare_report.md").exists()
    assert "Spearman" in (out / "ranking_compare_report.md").read_text(encoding="utf-8")


def test_ranking_compare_without_peptide_id(tmp_path: Path):
    rec = tmp_path / "ranked_peptides.recommendation.tsv"
    net = tmp_path / "ranked_peptides.netmhcpan42.tsv"
    rows_rec = [
        {"peptide_id": "", "gene": "G1", "peptide": "AAAA", "hla_allele": "HLA-A*02:01", "event_type": "SNV", "final_priority": "B", "recommended_use": "short peptide"},
        {"peptide_id": "", "gene": "G2", "peptide": "BBBB", "hla_allele": "HLA-B*07:02", "event_type": "Fusion", "final_priority": "C", "recommended_use": "minigene"},
    ]
    rows_net = [rows_rec[1], rows_rec[0]]
    write_tsv(rec, rows_rec)
    write_tsv(net, rows_net)
    out = tmp_path / "out_no_id"
    rc = subprocess.run([sys.executable, "-m", "neoag.agent_skills.ranking_compare", "--recommendation", str(rec), "--netmhcpan42", str(net), "--outdir", str(out)], text=True, capture_output=True)
    assert rc.returncode == 0, rc.stderr
    report = (out / "ranking_compare_report.md").read_text(encoding="utf-8")
    assert "Common candidate IDs" in report
    assert (out / "rank_shift.tsv").read_text(encoding="utf-8").startswith("peptide_id")


def test_generic_weighted_vs_consensus_comparison_audits_evidence(tmp_path: Path):
    left = tmp_path / "ranked_peptides.weighted_baseline.tsv"
    right = tmp_path / "ranked_peptides.evidence_consensus.tsv"
    rows = [
        {"peptide_id": "p1", "gene": "G1", "peptide": "AAAA", "hla_allele": "HLA-A*02:01", "event_type": "SNV", "final_priority": "B"},
        {"peptide_id": "p2", "gene": "G2", "peptide": "BBBB", "hla_allele": "HLA-B*07:02", "event_type": "Fusion", "final_priority": "C"},
        {"peptide_id": "p3", "gene": "G3", "peptide": "CCCC", "hla_allele": "HLA-C*07:02", "event_type": "InDel", "final_priority": "D"},
    ]
    write_tsv(left, rows)
    right_rows = [
        {**rows[1], "evidence_grade": "R4", "hard_failure": "yes", "hard_failure_codes": "HARD_NORMAL_JUNCTION", "manual_review_required": "yes", "evidence_conflict_fields": "rna_support_status", "evidence_missing_layers": "safety"},
        {**rows[0], "evidence_grade": "R1", "hard_failure": "no", "manual_review_required": "no"},
        {**rows[2], "evidence_grade": "R3", "hard_failure": "no", "manual_review_required": "yes"},
    ]
    write_tsv(right, right_rows)
    out = tmp_path / "generic_comparison"
    rc = subprocess.run([
        sys.executable, "-m", "neoag.agent_skills.ranking_compare",
        "--left", str(left), "--left-name", "weighted_baseline",
        "--right", str(right), "--right-name", "evidence_consensus",
        "--outdir", str(out),
    ], text=True, capture_output=True)
    assert rc.returncode == 0, rc.stderr
    summary = json.loads((out / "ranking_comparison_summary.json").read_text())
    assert summary["left"]["name"] == "weighted_baseline"
    assert summary["right"]["name"] == "evidence_consensus"
    assert summary["spearman_correlation"] == 0.5
    overlap = list(csv.DictReader((out / "topn_overlap.tsv").open(), delimiter="\t"))
    assert [row["top_n"] for row in overlap] == ["10", "20", "50", "100"]
    changes = {row["candidate_id"]: row for row in csv.DictReader((out / "candidate_rank_changes.tsv").open(), delimiter="\t")}
    assert changes["p2"]["change"] == "PROMOTED_IN_RIGHT"
    hard = list(csv.DictReader((out / "high_rank_hard_fail.tsv").open(), delimiter="\t"))
    assert hard[0]["candidate_id"] == "p2"
    assert (out / "top_composition.tsv").is_file()
    assert (out / "evidence_qc_summary.tsv").is_file()
    assert (out / "manual_review_candidates.tsv").is_file()
    report = (out / "ranking_compare_report.md").read_text()
    assert "Conflict and missing-evidence rates" in report
    assert "Top20 HLA coverage" in report


def test_input_qc_detects_result_files(tmp_path: Path):
    (tmp_path / "ranked_peptides.recommendation.tsv").write_text("peptide_id\nP1\n", encoding="utf-8")
    (tmp_path / "evidence_report.v041.html").write_text("<html></html>", encoding="utf-8")
    out = tmp_path / "qc"
    rc = subprocess.run([sys.executable, "-m", "neoag.agent_skills.input_qc", "--result-dir", str(tmp_path), "--outdir", str(out)], text=True, capture_output=True)
    assert rc.returncode == 0, rc.stderr
    status = json.loads((out / "input_status.json").read_text(encoding="utf-8"))
    assert status["features"]["has_ranked_peptides_recommendation"] is True
    assert status["recommended_workflow"] == "result_review"


def test_coordinator_plans_ranking_compare(tmp_path: Path):
    rec = tmp_path / "ranked_peptides.recommendation.tsv"
    net = tmp_path / "ranked_peptides.netmhcpan42.tsv"
    write_tsv(rec, [{"peptide_id": "p1", "gene": "G", "peptide": "AAAA", "hla_allele": "HLA-A*02:01", "event_type": "SNV", "final_priority": "B", "recommended_use": ""}])
    write_tsv(net, [{"peptide_id": "p1", "gene": "G", "peptide": "AAAA", "hla_allele": "HLA-A*02:01", "event_type": "SNV", "final_priority": "B", "recommended_use": ""}])
    out = tmp_path / "agent"
    rc = subprocess.run([sys.executable, "-m", "neoag.agents.coordinator", "--message", "比较 recommendation 和 NetMHCpan42 排序差异", "--file", str(rec), "--file", str(net), "--outdir", str(out), "--execute"], text=True, capture_output=True)
    assert rc.returncode == 0, rc.stderr
    state = json.loads((out / "case_state.json").read_text(encoding="utf-8"))
    assert state["intent"] == "ranking_compare"
    assert "neoag-ranking-compare" in state["planned_skills"]
    assert (out / "neoag-ranking-compare" / "ranking_compare_report.md").exists()


def test_coordinator_plans_sliding_run(tmp_path: Path):
    out = tmp_path / "agent_sliding"
    rc = subprocess.run([sys.executable, "-m", "neoag.agents.coordinator", "--message", "run sliding-window SNV/InDel VEP workflow from somatic VCF", "--outdir", str(out)], text=True, capture_output=True)
    assert rc.returncode == 0, rc.stderr
    state = json.loads((out / "case_state.json").read_text(encoding="utf-8"))
    assert state["intent"] == "sliding_run"
    assert "neoag-sliding-run" in state["planned_skills"]
    assert any(call["skill"] == "neoag-sliding-run" and call["status"] == "SKIPPED" for call in state["calls"])
