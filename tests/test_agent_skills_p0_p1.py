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
    rc = subprocess.run([sys.executable, "-m", "neoag_v03.agent_skills.ranking_compare", "--recommendation", str(rec), "--netmhcpan42", str(net), "--outdir", str(out)], text=True, capture_output=True)
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
    rc = subprocess.run([sys.executable, "-m", "neoag_v03.agent_skills.ranking_compare", "--recommendation", str(rec), "--netmhcpan42", str(net), "--outdir", str(out)], text=True, capture_output=True)
    assert rc.returncode == 0, rc.stderr
    report = (out / "ranking_compare_report.md").read_text(encoding="utf-8")
    assert "Common candidate IDs" in report
    assert (out / "rank_shift.tsv").read_text(encoding="utf-8").startswith("peptide_id")


def test_input_qc_detects_result_files(tmp_path: Path):
    (tmp_path / "ranked_peptides.recommendation.tsv").write_text("peptide_id\nP1\n", encoding="utf-8")
    (tmp_path / "evidence_report.v041.html").write_text("<html></html>", encoding="utf-8")
    out = tmp_path / "qc"
    rc = subprocess.run([sys.executable, "-m", "neoag_v03.agent_skills.input_qc", "--result-dir", str(tmp_path), "--outdir", str(out)], text=True, capture_output=True)
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
    rc = subprocess.run([sys.executable, "-m", "neoag_v03.agents.coordinator", "--message", "比较 recommendation 和 NetMHCpan42 排序差异", "--file", str(rec), "--file", str(net), "--outdir", str(out), "--execute"], text=True, capture_output=True)
    assert rc.returncode == 0, rc.stderr
    state = json.loads((out / "case_state.json").read_text(encoding="utf-8"))
    assert state["intent"] == "ranking_compare"
    assert "neoag-ranking-compare" in state["planned_skills"]
    assert (out / "neoag-ranking-compare" / "ranking_compare_report.md").exists()
