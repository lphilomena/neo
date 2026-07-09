from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


def make_rank_file(path: Path, rows: list[dict[str, str]]) -> None:
    cols = ["peptide_id", "gene", "peptide", "hla_allele", "event_type", "final_priority", "recommended_use"]
    with path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(cols)+"\n")
        for r in rows:
            fh.write("\t".join(r.get(c, "") for c in cols)+"\n")


def test_llm_coordinator_plan_ranking_compare(tmp_path: Path) -> None:
    rec = tmp_path / "ranked_peptides.recommendation.tsv"
    net = tmp_path / "ranked_peptides.netmhcpan42.tsv"
    rows = [
        {"peptide_id":"p1","gene":"A","peptide":"AAAA","hla_allele":"HLA-A*02:01","event_type":"SNV","final_priority":"B"},
        {"peptide_id":"p2","gene":"B","peptide":"BBBB","hla_allele":"HLA-A*02:01","event_type":"Fusion","final_priority":"D"},
    ]
    make_rank_file(rec, rows)
    make_rank_file(net, list(reversed(rows)))
    out = tmp_path / "out"
    cmd = [sys.executable, "-m", "neoag_v03.llm_coordinator.coordinator_agent", "--message", "比较 recommendation 和 NetMHCpan42 排序差异", "--file", str(rec), "--file", str(net), "--outdir", str(out), "--mode", "plan"]
    proc = subprocess.run(cmd, text=True, capture_output=True)
    assert proc.returncode == 0, proc.stderr
    state = json.loads((out / "case_state.json").read_text(encoding="utf-8"))
    assert state["intent"]["intent"] == "ranking_compare"
    assert "neoag-ranking-compare" in [s["skill"] for s in state["plan"]["steps"]]
    assert (out / "coordinator_plan.md").exists()


def test_llm_coordinator_execute_safe_ranking_compare(tmp_path: Path) -> None:
    rec = tmp_path / "ranked_peptides.recommendation.tsv"
    net = tmp_path / "ranked_peptides.netmhcpan42.tsv"
    rows_rec = [
        {"peptide_id":"p1","gene":"A","peptide":"AAAA","hla_allele":"HLA-A*02:01","event_type":"SNV","final_priority":"B"},
        {"peptide_id":"p2","gene":"B","peptide":"BBBB","hla_allele":"HLA-A*02:01","event_type":"Fusion","final_priority":"D"},
        {"peptide_id":"p3","gene":"C","peptide":"CCCC","hla_allele":"HLA-B*07:02","event_type":"InDel","final_priority":"C"},
    ]
    rows_net = [rows_rec[1], rows_rec[2], rows_rec[0]]
    make_rank_file(rec, rows_rec)
    make_rank_file(net, rows_net)
    out = tmp_path / "out_exec"
    cmd = [sys.executable, "-m", "neoag_v03.llm_coordinator.coordinator_agent", "--message", "比较 recommendation 和 NetMHCpan42 排序差异", "--file", str(rec), "--file", str(net), "--outdir", str(out), "--mode", "execute-safe"]
    proc = subprocess.run(cmd, text=True, capture_output=True)
    assert proc.returncode == 0, proc.stderr
    assert (out / "neoag-ranking-compare" / "ranking_compare_report.md").exists()
    state = json.loads((out / "case_state.json").read_text(encoding="utf-8"))
    statuses = {r["skill_name"]: r["status"] for r in state["skill_results"]}
    assert statuses.get("neoag-ranking-compare") == "PASS"


def test_llm_coordinator_guardrail_hpc_plan(tmp_path: Path) -> None:
    out = tmp_path / "out_hpc"
    cmd = [sys.executable, "-m", "neoag_v03.llm_coordinator.coordinator_agent", "--message", "请直接提交 Slurm 到 HPC 跑全流程", "--outdir", str(out), "--mode", "execute-safe"]
    proc = subprocess.run(cmd, text=True, capture_output=True)
    assert proc.returncode == 0, proc.stderr
    state = json.loads((out / "case_state.json").read_text(encoding="utf-8"))
    assert state["plan"]["approval_required"] is True
    assert all(r["status"] in {"PASS", "PLANNED", "APPROVAL_REQUIRED", "SKIPPED"} for r in state["skill_results"])
