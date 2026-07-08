from __future__ import annotations

import re
from pathlib import Path
from .intent_schema import INTENT_TO_SKILLS


def classify_intent(message: str) -> str:
    m = message.lower()
    if any(x in m for x in ["比较", "差异", "netmhc", "netmhcpan", "排序"]):
        return "ranking_compare"
    if any(x in m for x in ["患者", "沟通", "word", "报告", "docx", "ppt"]):
        return "patient_report_update"
    if any(x in m for x in ["appm", "免疫逃逸", "hla loh", "抗原呈递", "hlo", "lohhla"]):
        return "appm_escape_review"
    if any(x in m for x in ["ccf", "克隆", "clonality", "purity", "纯度"]):
        return "ccf_review"
    if any(x in m for x in ["工具", "安装", "部署", "reference", "参考库", "check-tools"]):
        return "tool_check"
    if any(x in m for x in ["demo", "smoke", "测试", "pytest", "nextflow"]):
        return "demo_smoke"
    if any(x in m for x in ["评分", "scoring", "score-v03", "综合评分"]):
        return "run_scoring"
    if any(x in m for x in ["实验", "验证", "elis", "minigene", "long peptide"]):
        return "experiment_design"
    if any(x in m for x in ["能不能跑", "缺", "输入", "manifest", "检查"]):
        return "input_check"
    return "input_check"


def skills_for_intent(intent: str) -> list[str]:
    return INTENT_TO_SKILLS.get(intent, [])


def find_named_files(result_dir: str | None, explicit_files: list[str] | None = None) -> dict[str, str]:
    files = [Path(f) for f in (explicit_files or [])]
    if result_dir:
        base = Path(result_dir)
        if base.exists():
            files.extend([p for p in base.rglob("*") if p.is_file()])
    out: dict[str, str] = {}
    def put(key: str, pred):
        if key in out:
            return
        for p in files:
            name = p.name.lower()
            if pred(name):
                out[key] = str(p)
                return
    put("recommendation", lambda n: "ranked_peptides.recommendation" in n)
    put("netmhcpan42", lambda n: "ranked_peptides.netmhcpan42" in n)
    put("evidence_report", lambda n: n.startswith("evidence_report") and n.endswith(".html"))
    put("hla_loh", lambda n: "hla_loh" in n or "lohhla" in n or "spechla" in n)
    put("purity_table", lambda n: "purity" in n or "facets" in n or "purple" in n)
    put("appm_gene_status", lambda n: "appm_gene_status" in n)
    put("appm_submodule_scores", lambda n: "appm_submodule_scores" in n)
    put("ranked_peptides", lambda n: n.startswith("ranked_peptides") and n.endswith(".tsv"))
    return out
