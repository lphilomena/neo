from __future__ import annotations

HIGH_RISK_KEYWORDS = ["提交", "slurm", "sbatch", "覆盖", "删除", "安装", "download", "rm -", "hpc", "正式报告"]
PROHIBITED_CLAIMS = ["一定有效", "保证有效", "临床耐药", "确定治疗方案", "已确认新抗原", "一定获益"]


def needs_human_approval(message: str, planned_skills: list[str]) -> bool:
    lower = message.lower()
    if any(k.lower() in lower for k in HIGH_RISK_KEYWORDS):
        return True
    if any(s in planned_skills for s in ["neoag-run-demo-and-smoke", "neoag-evidence-scoring"]):
        # These can be dry-run without approval, but execution should be confirmed by caller.
        return False
    return False


def sanitize_patient_language(text: str) -> str:
    out = text
    for claim in PROHIBITED_CLAIMS:
        out = out.replace(claim, "计算预测候选，需进一步验证")
    return out


def boundary_note() -> str:
    return "计算筛选结果仅用于候选优先级和实验设计参考；候选新抗原不等同于已验证新抗原，也不等同于确定治疗方案。"
