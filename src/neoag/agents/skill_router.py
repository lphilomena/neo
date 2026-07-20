from __future__ import annotations

import re
from pathlib import Path
from .intent_schema import INTENT_TO_SKILLS


def classify_intent(message: str) -> str:
    m = message.lower()
    if any(x in m for x in ["概述", "介绍项目", "项目功能", "功能介绍", "项目是做什么", "项目做什么", "project overview", "what does this project do"]):
        return "project_overview"
    if any(x in m for x in ["状态", "跑完", "进度", "日志", "status", "running", "finished"]):
        return "check_status"
    if any(x in m for x in ["报错", "错误", "失败", "traceback", "error", "debug"]):
        return "debug_error"
    if any(x in m for x in ["复制", "同步", "rsync", "scp", "移动硬盘", "copy"]):
        return "data_transfer"
    if any(x in m for x in ["提交", "推送", "发布", "release", "git", "github"]):
        return "git_release"
    if any(x in m for x in ["readme", "文档", "说明", "docs"]):
        return "update_docs"
    if any(x in m for x in ["安装", "配置", "docker", "apptainer", "容器"]):
        return "setup_tool"
    if any(x in m for x in ["facets", "purple", "sequenza", "ascat", "纯度", "ploidy", "purity", "cnv", "copy number", "拷贝数"]) and any(x in m for x in ["综合", "对比", "分析", "review", "compare", "跑", "运行", "run", "检测"]):
        return "purity_cnv_review"
    if any(x in m for x in ["optitype", "spechla", "hla-la", "hla typing", "hla分型", "hla 分型", "分型结果"]) and any(x in m for x in ["综合", "对比", "比较", "分析", "review", "compare", "跑", "运行", "run", "检测"]):
        return "hla_typing_compare"
    if any(x in m for x in ["tpm", "表达量", "表达矩阵", "rna表达", "rna 表达", "rsem", "salmon", "kallisto", "fastq生成tpm", "fastq 生成 tpm"]):
        return "rna_fastq_to_tpm"
    if any(x in m for x in ["easyfuse", "star-fusion", "star_fusion", "arriba", "fusioncatcher", "fusion", "融合"]):
        return "fusion_rna_run"
    if any(x in m for x in ["解读", "分析结果", "结果", "对比", "compare"]):
        return "inspect_results"
    workflow_terms = ["spechla", "hla-la", "optitype", "facets", "purple", "sequenza", "easyfuse", "lohhla", "netmhcpan", "bigmhc"]
    if any(x in m for x in workflow_terms) and any(x in m for x in ["运行", "跑", "run", "执行"]):
        return "workflow_run_request"
    if any(x in m for x in ["sliding", "run-full", "snv_indel", "vep", "extract-variant-peptides", "variant peptide", "variant peptides", "短肽", "滑窗", "somatic vcf"]):
        return "sliding_run"
    if any(x in m for x in ["比较", "差异", "netmhc", "netmhcpan", "排序"]):
        return "ranking_compare"
    if any(x in m for x in ["患者", "沟通", "word", "报告", "docx", "ppt"]):
        return "patient_report_update"
    if any(x in m for x in ["appm", "免疫逃逸", "hla loh", "抗原呈递", "hlo", "lohhla"]):
        return "appm_escape_review"
    if any(x in m for x in ["ccf", "克隆", "clonality"]):
        return "ccf_review"
    if any(x in m for x in ["工具", "安装", "部署", "reference", "参考库", "check-tools"]):
        return "tool_check"
    if any(x in m for x in ["demo", "smoke", "测试", "pytest", "nextflow"]):
        return "demo_smoke"
    if any(x in m for x in ["评分", "scoring", "score", "综合评分"]):
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
    put("variants_vcf", lambda n: (n.endswith(".vcf") or n.endswith(".vcf.gz")) and ("somatic" in n or "variant" in n or "mutect" in n))
    put("somatic_vcf", lambda n: (n.endswith(".vcf") or n.endswith(".vcf.gz")) and "somatic" in n)
    put("hla", lambda n: ("hla" in n) and (n.endswith(".txt") or n.endswith(".tsv") or n.endswith(".csv")))
    put("fusion", lambda n: ("fusion" in n or "arriba" in n or "easyfuse" in n or "star-fusion" in n or "fusioncatcher" in n) and (n.endswith(".tsv") or n.endswith(".csv") or n.endswith(".txt") or n.endswith(".json")))
    put("expression", lambda n: ("tpm" in n or "expression" in n or "genes.results" in n or "quant.sf" in n) and (n.endswith(".tsv") or n.endswith(".csv") or n.endswith(".txt") or n.endswith(".sf") or n.endswith(".results")))
    put("fastq1", lambda n: n.endswith(("_1.fq.gz", "_r1.fq.gz", "_1.fastq.gz", "_r1.fastq.gz", "_1.fq", "_r1.fq")))
    put("fastq2", lambda n: n.endswith(("_2.fq.gz", "_r2.fq.gz", "_2.fastq.gz", "_r2.fastq.gz", "_2.fq", "_r2.fq")))
    put("tumor_bam", lambda n: n.endswith(".bam") and any(x in n for x in ["tumor", "tumour", "237", "184", "aml"]))
    put("normal_bam", lambda n: n.endswith(".bam") and any(x in n for x in ["normal", "blood", "236", "control"]))
    return out
