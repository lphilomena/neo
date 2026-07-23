from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any


@dataclass(frozen=True)
class SkillSpec:
    name: str
    category: str
    category_label: str
    purpose: str
    description: str
    use_when: list[str]
    do_not_use_when: list[str]
    required_inputs: list[str] = field(default_factory=list)
    optional_inputs: list[str] = field(default_factory=list)
    outputs: list[str] = field(default_factory=list)
    risk_level: str = "LOW"
    approval_required: bool = False
    handler: str = ""
    external_tools: list[str] = field(default_factory=list)
    downstream_skills: list[str] = field(default_factory=list)
    boundaries: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


CATEGORY_LABELS = {
    "A": "入口适配型 Skills：把不同来源输入转换为 Project B 标准 raw_events/raw_peptides/evidence tables",
    "B": "公共证据分析型 Skills：对所有入口共用的 HLA、表达、CCF、APPM、安全和排序证据层进行标准化分析",
    "C": "审阅/报告/实验设计型 Skills：解释结果、生成报告、设计实验验证和患者沟通材料",
    "D": "工程治理/执行控制型 Skills：输入质控、环境健康检查、全流程编排、发布审计和受控执行",
}

COMMON_BOUNDARIES = [
    "Skill 不承担临床决策；不得判断患者是否适合治疗或推荐临床用药。",
    "缺失证据只能标记为 missing/unassessed，不能解释为阴性结果。",
    "高风险写入、HPC 提交、安装工具、下载参考库、删除或覆盖文件必须经过 human approval。",
    "Skill 目录不包含患者 BAM/FASTQ/VCF、大型参考库、VEP cache、NetMHCpan license、LOHHLA reference 或大型 conda env。",
]


def _spec(**kwargs: Any) -> SkillSpec:
    category = kwargs["category"]
    boundaries = list(kwargs.pop("boundaries", [])) + COMMON_BOUNDARIES
    return SkillSpec(category_label=CATEGORY_LABELS[category], boundaries=boundaries, **kwargs)


SKILL_SPECS: list[SkillSpec] = [
    # A: Entry adapters
    _spec(
        name="neoag-vcf", category="A", handler="vcf",
        purpose="VCF/SNV/InDel 入口标准化",
        description="Parse a somatic VCF or VEP-annotated VCF into Project B raw_events/raw_peptides-compatible tables.",
        use_when=["用户提供 tumor-normal somatic VCF", "需要把 SNV/InDel 输入转成 raw_events.tsv"],
        do_not_use_when=["输入是 fusion/splice/SV/peptide 表", "需要直接生成综合排序，应调用 neoag-ranking"],
        required_inputs=["vcf"], optional_inputs=["sample_id", "hla", "expression_tsv"],
        outputs=["raw_events.tsv", "raw_peptides.tsv", "vcf_parse_qc.tsv", "candidate_generation_plan.md"],
        risk_level="LOW", downstream_skills=["neoag-presentation", "neoag-ranking"],
    ),
    _spec(
        name="neoag-fusion", category="A", handler="fusion",
        purpose="Fusion caller 输出标准化",
        description="Normalize EasyFuse/Arriba/STAR-Fusion/FusionCatcher fusion outputs into raw_events/raw_peptides and fusion evidence tables.",
        use_when=["用户提供 fusions.pass.csv 或其他 fusion caller 输出", "需要生成 fusion junction 候选输入"],
        do_not_use_when=["fusion 结果尚未经过 read-through 复核时直接作为临床结论"],
        required_inputs=["fusion"], optional_inputs=["sample_id", "normal_readthrough_db"],
        outputs=["fusion_events.tsv", "raw_events.tsv", "raw_peptides.tsv", "fusion_evidence.tsv", "fusion_qc.tsv"],
        risk_level="LOW", downstream_skills=["neoag-presentation", "neoag-experiment-design"],
    ),
    _spec(
        name="neoag-splice", category="A", handler="splice",
        purpose="Splice/junction 输入标准化",
        description="Normalize RegTools/splice junction tables and optional variant VCF into splice-junction event and peptide candidate inputs.",
        use_when=["用户提供 RegTools splice junction TSV", "需要构建 splice/exon junction 候选"],
        do_not_use_when=["仅有普通 SNV/InDel VCF 时，应调用 neoag-vcf"],
        required_inputs=["junctions"], optional_inputs=["vcf", "sample_id"],
        outputs=["splice_events.tsv", "raw_events.tsv", "raw_peptides.tsv", "splice_qc.tsv"],
        risk_level="LOW", downstream_skills=["neoag-rna-evidence", "neoag-experiment-design"],
    ),
    _spec(
        name="neoag-sv-wgs", category="A", handler="sv_wgs",
        purpose="WGS DNA SV 入口标准化",
        description="Parse WGS SV VCFs (Manta/GRIDSS/SvABA/DELLY) into SV events requiring transcript/protein reconstruction.",
        use_when=["用户提供 WGS SV VCF", "需要进入 DNA SV-derived neoantigen workflow"],
        do_not_use_when=["输入来自 WES/capture-limited SV，应调用 neoag-sv-wes"],
        required_inputs=["sv_vcf"], optional_inputs=["gtf", "reference_fasta", "sample_id", "rna_junctions"],
        outputs=["sv_events.tsv", "raw_events.tsv", "raw_peptides.tsv", "sv_reconstruction_tasks.tsv"],
        risk_level="LOW", downstream_skills=["neoag-presentation", "neoag-experiment-design"],
    ),
    _spec(
        name="neoag-sv-wes", category="A", handler="sv_wes",
        purpose="WES/capture-limited SV 入口标准化",
        description="Parse WES/capture-limited SV VCFs and enforce conservative confidence caps and capture-limited flags.",
        use_when=["用户提供 WES/capture-limited SV VCF + capture BED", "需要保守解释外显子捕获 SV"],
        do_not_use_when=["输入为 WGS SV，应调用 neoag-sv-wgs"],
        required_inputs=["sv_vcf", "capture_bed"], optional_inputs=["gtf", "reference_fasta", "sample_id"],
        outputs=["sv_events.tsv", "raw_events.tsv", "raw_peptides.tsv", "sv_wes_confidence.tsv"],
        risk_level="LOW", downstream_skills=["neoag-presentation", "neoag-experiment-design"],
        boundaries=["WES SV 是 capture-limited hypothesis，默认 final priority 不应直接升至 A/B high-confidence。"],
    ),
    _spec(
        name="neoag-peptide-csv", category="A", handler="peptide_csv",
        purpose="已有 peptide-HLA 表入口标准化",
        description="Normalize an existing peptide-HLA table into raw_peptides and optional presentation_evidence tables.",
        use_when=["用户已有 peptide-HLA 候选表", "需要将外部排序/预测表接入 Project B scoring"],
        do_not_use_when=["需要从 VCF/fusion/SV 生成候选时"],
        required_inputs=["peptide_csv"], optional_inputs=["sample_id", "event_annotation"],
        outputs=["raw_peptides.tsv", "presentation_evidence.tsv", "peptide_input_qc.tsv"],
        risk_level="LOW", downstream_skills=["neoag-safety", "neoag-ranking"],
    ),
    # B: Evidence analysis
    _spec(
        name="neoag-hla-typing-loh", category="B", handler="hla_typing_loh",
        purpose="HLA typing / HLA LOH 共识与 peptide-level HLA loss flags",
        description="Normalize HLA typing and HLA LOH outputs from OptiType, SpecHLA, HLA-LA/HD and LOHHLA into consensus tables.",
        use_when=["需要标准化 HLA 分型", "需要判断 restricting HLA 是否 LOH"],
        do_not_use_when=["只需要解释已有 ranking 差异，不需要更新 HLA 状态"],
        required_inputs=["hla"], optional_inputs=["hla_loh", "ranked_peptides", "sample_id"],
        outputs=["hla_typing.normalized.tsv", "hla_typing_consensus.tsv", "hla_loh_consensus.tsv", "restricting_hla_peptide_flags.tsv", "hla_review.md"],
        risk_level="LOW", external_tools=["OptiType", "SpecHLA", "LOHHLA", "HLA-LA/HLA-HD optional"],
    ),
    _spec(
        name="neoag-presentation", category="B", handler="presentation",
        purpose="HLA binding / presentation 证据标准化",
        description="Normalize NetMHCpan/MHCflurry/PRIME/BigMHC/MixMHCpred predictions into presentation_evidence.tsv.",
        use_when=["已有 prediction 表需要标准化", "raw_peptides 需要标记 presentation 计算计划"],
        do_not_use_when=["需要解释综合推荐排序，应调用 neoag-ranking-compare"],
        required_inputs=["predictions_or_raw_peptides"], optional_inputs=["hla", "sample_id"],
        outputs=["presentation_evidence.tsv", "presentation_summary.tsv", "presentation_qc.tsv"],
        risk_level="LOW", external_tools=["NetMHCpan", "MHCflurry", "PRIME", "BigMHC", "MixMHCpred"],
    ),
    _spec(
        name="neoag-expression", category="B", handler="expression",
        purpose="gene expression/TPM 证据标准化",
        description="Normalize WTS/gene expression TPM tables into expression_evidence.tsv with expressed/low/not_detected labels.",
        use_when=["用户提供 TPM/gene expression 表", "需要把表达证据接入 ranking"],
        do_not_use_when=["需要 RNA alt reads 或 junction reads，使用 neoag-rna-evidence"],
        required_inputs=["expression_tsv"], optional_inputs=["sample_id"],
        outputs=["expression_evidence.tsv", "expression_qc.tsv"], risk_level="LOW",
    ),
    _spec(
        name="neoag-rna-evidence", category="B", handler="rna_evidence",
        purpose="RNA alt / RNA VAF / junction reads 证据标准化",
        description="Normalize RNA allele, RNA VAF and fusion/splice junction evidence into event-level RNA support tables.",
        use_when=["用户提供 RNA alt reads、RNA VAF、fusion/splice junction reads"],
        do_not_use_when=["只有 gene TPM 时，使用 neoag-expression"],
        required_inputs=["rna_tsv"], optional_inputs=["sample_id"],
        outputs=["rna_alt_evidence.tsv", "rna_junction_evidence.tsv", "rna_evidence_qc.tsv"], risk_level="LOW",
    ),
    _spec(
        name="neoag-ccf", category="B", handler="ccf",
        purpose="CCF/clonality 估计与 modifier 输出",
        description="Estimate or normalize Cancer Cell Fraction and clonality status for events, with confidence flags.",
        use_when=["需要 CCF/clonality 分层", "需要计算或解释 CCF modifier"],
        do_not_use_when=["RNA-only fusion 不应伪造 DNA CCF"],
        required_inputs=["event_table_or_ranked_peptides"], optional_inputs=["purity", "cnv_segments", "sample_id"],
        outputs=["ccf_lite.tsv", "ccf_input_qc.tsv", "ccf_modifier_summary.tsv"], risk_level="LOW",
    ),
    _spec(
        name="neoag-appm-escape", category="B", handler="appm_escape",
        purpose="APPM / HLA LOH / immune escape 证据层",
        description="Aggregate APPM gene status, HLA LOH and immune escape flags into pathway-level and peptide-level modifiers.",
        use_when=["需要 APPM、HLA LOH、B2M/TAP/JAK/NLRC5/CIITA 风险评估"],
        do_not_use_when=["缺所有 APPM 输入时，只能输出 unassessed，不得输出 intact"],
        required_inputs=["gene_status_or_appm", "peptides_optional"], optional_inputs=["hla_loh", "ranked_peptides"],
        outputs=["appm_summary.tsv", "appm_gene_status.tsv", "appm_peptide_modifiers.tsv", "immune_escape_summary.tsv", "peptide_escape_flags.tsv"], risk_level="LOW",
    ),
    _spec(
        name="neoag-safety", category="B", handler="safety",
        purpose="peptide safety gate",
        description="Check normal proteome exact matches, WT peptide flags, anchor-only risks and normal junction/ligandome overlap.",
        use_when=["需要筛查正常蛋白组匹配、WT peptide、anchor-only、normal junction 风险"],
        do_not_use_when=["不能替代湿实验安全性评估"],
        required_inputs=["raw_peptides_or_ranked_peptides"], optional_inputs=["normal_proteome", "normal_junctions", "wt_binding"],
        outputs=["peptide_safety.tsv", "event_safety.tsv", "safety_review.tsv"], risk_level="LOW",
    ),
    _spec(
        name="neoag-ranking", category="B", handler="ranking",
        purpose="生产级加权基线与证据共识并行排序",
        description="Thin SOP wrapper around the production neoag evidence-rank CLI; it does not implement ranking logic.",
        use_when=["已有 comprehensive evidence 和 weighted baseline，需要生成并行 evidence-consensus 排序"],
        do_not_use_when=["只需比较两个排序文件时，使用 neoag-ranking-compare"],
        required_inputs=["comprehensive_evidence", "weighted_baseline"], optional_inputs=["rules", "provenance", "track"],
        outputs=["all_tool_results.tsv", "ranked_peptides.weighted_baseline.tsv", "ranked_peptides.evidence_consensus.tsv", "ranked_events.evidence_consensus.tsv", "ranking_compare_weighted_vs_consensus.md"], risk_level="LOW",
        boundaries=["Skill 只调用 neoag evidence-rank；正式排序算法只由 src/neoag/evidence_consensus.py 维护。"],
    ),
    _spec(
        name="open-neo-run", category="B", handler="ranking",
        purpose="对外宏 Skill2：执行生产级并行证据排序",
        description="Public macro Skill2 alias that invokes neoag evidence-rank through the shared production CLI wrapper.",
        use_when=["外部 Agent 需要用统一 Skill2 入口生成加权基线与证据共识并行结果"],
        do_not_use_when=["尚未生成 comprehensive evidence 或 weighted baseline"],
        required_inputs=["comprehensive_evidence", "weighted_baseline"], optional_inputs=["rules", "provenance", "track"],
        outputs=["all_tool_results.tsv", "ranked_peptides.weighted_baseline.tsv", "ranked_peptides.evidence_consensus.tsv", "ranked_events.evidence_consensus.tsv", "ranking_compare_weighted_vs_consensus.md"], risk_level="LOW",
        boundaries=["Skill = SOP/包装层；不得在 Skill 中复制 R1-R4、Pareto、hard fail 或 priority-cap 算法。"],
    ),
    # C: Review/report/design
    _spec(
        name="neoag-ranking-compare", category="C", handler="ranking_compare",
        purpose="任意两个候选排序的通用比较与审计",
        description="Compare arbitrary left/right rankings using overlap, Spearman correlation, rank shifts, hard-fail audit, event/HLA composition and evidence-quality rates.",
        use_when=["需要比较 weighted baseline 与 evidence consensus", "需要比较任意两个 peptide ranking"],
        do_not_use_when=["需要生成 ranking 时，使用 open-neo-run 或 neoag-ranking"],
        required_inputs=["left", "right"], optional_inputs=["left_name", "right_name"],
        outputs=["ranking_compare_report.md", "topn_overlap.tsv", "candidate_rank_changes.tsv", "high_rank_hard_fail.tsv", "top_composition.tsv", "evidence_qc_summary.tsv", "manual_review_candidates.tsv", "ranking_comparison_summary.json"], risk_level="LOW",
    ),
    _spec(
        name="neoag-experiment-design", category="C", handler="experiment_design",
        purpose="候选实验验证设计",
        description="Assign event-deduplicated representatives to short peptide, WT control, long peptide, minigene and targeted RNA validation routes.",
        use_when=["需要第一批 10–20 个实验候选", "需要区分 short peptide / long peptide / minigene / targeted RNA"],
        do_not_use_when=["不要把设计建议写成治疗处方"],
        required_inputs=["ranked_events_or_ranked_peptides"], optional_inputs=["ranked_events", "ranked_peptides", "top_n", "therapy_context"],
        outputs=["experiment_candidates.tsv", "short_peptide_pool.tsv", "long_peptide_design.tsv", "minigene_design.tsv", "targeted_rna_validation_plan.md"], risk_level="LOW",
    ),
    _spec(
        name="neoag-patient-report", category="C", handler="patient_report",
        purpose="患者沟通版报告",
        description="Generate a patient-facing markdown/html/docx draft with research boundary statements.",
        use_when=["需要更新患者沟通版报告", "需要解释 A/B/C/D、CCF、APPM、HLA LOH"],
        do_not_use_when=["不要生成临床处方或疗效承诺"],
        required_inputs=["ranked_peptides_or_summary"], optional_inputs=["evidence_report", "ranking_compare", "appm_review", "ccf_review"],
        outputs=["patient_report.md", "patient_report.html", "patient_report.docx"], risk_level="MEDIUM",
    ),
    _spec(
        name="neoag-technical-report", category="C", handler="technical_report",
        purpose="技术审阅报告",
        description="Generate a technical report preserving tool versions, evidence completeness, fields, warnings and provenance.",
        use_when=["需要医生/生信团队审阅的技术报告"],
        do_not_use_when=["患者直接沟通应使用 neoag-patient-report"],
        required_inputs=["result_dir_or_summary"], optional_inputs=["doctor_report", "pipeline_manifest"],
        outputs=["technical_report.md", "technical_report.html"], risk_level="LOW",
    ),
    _spec(
        name="neoag-concept-explainer", category="C", handler="concept_explainer",
        purpose="术语解释与报告注释",
        description="Generate bounded explanations for APPM, CCF, HLA LOH, NetMHCpan, minigene, ELISpot and validation concepts.",
        use_when=["用户询问术语解释", "报告需要插入概念框"],
        do_not_use_when=["不要替代具体结果分析"],
        required_inputs=["concept"], optional_inputs=["audience"],
        outputs=["concept_explanation.md"], risk_level="LOW",
    ),
    # D: Governance/execution
    _spec(
        name="neoag-input-qc", category="D", handler="input_qc",
        purpose="输入状态检查与 workflow 推荐",
        description="Inspect manifests/result directories and recommend workflow while listing missing inputs.",
        use_when=["任何任务的第一步", "用户问能不能跑、缺什么输入"],
        do_not_use_when=["不能用 input-qc 的缺失信息直接做生物学阴性结论"],
        required_inputs=["manifest_or_result_dir"], outputs=["input_status.json", "input_qc_report.tsv", "missing_inputs.tsv"], risk_level="LOW",
    ),
    _spec(name="neoag-doctor", category="D", handler="doctor", purpose="只读环境健康检查", description="Run controlled-execution read-only health check for tools, references and release boundaries.", use_when=["部署/HPC/新机器验证", "用户问工具和参考库是否可用"], do_not_use_when=["不要用 doctor 安装或修改工具"], required_inputs=["project_root"], optional_inputs=["tools_manifest", "reference_manifest"], outputs=["doctor_status.json", "doctor_summary.md", "blocking_issues.tsv"], risk_level="LOW"),
    _spec(name="neoag-tool-reference-qc", category="D", handler="tool_reference_qc", purpose="工具/参考库检查", description="Check external tool entrypoints, models, caches and references; stronger than simple PATH check.", use_when=["需要检查 VEP/NetMHCpan/MHCflurry/LOHHLA/FACETS 等"], do_not_use_when=["不要将工具 missing 解释为生物学阴性"], required_inputs=["project_root"], optional_inputs=["tools_manifest", "reference_manifest"], outputs=["tool_qc_report.tsv", "reference_qc_report.tsv", "tool_smoke_report.md"], risk_level="LOW"),
    _spec(name="neoag-run-demo-and-smoke", category="D", handler="run_demo_smoke", purpose="demo/pytest/Nextflow smoke", description="Run or plan Project B run-demo, pytest and optional Nextflow smoke tests.", use_when=["release 验收", "新环境最小可运行性测试"], do_not_use_when=["不要在没有确认的生产目录覆盖结果"], required_inputs=["project_root"], outputs=["smoke_test_report.md", "demo_output_manifest.tsv"], risk_level="MEDIUM", approval_required=False),
    _spec(name="neoag-pipeline-full", category="D", handler="pipeline_full", purpose="manifest-driven full pipeline runner", description="controlled-execution explicit pipeline-full DAG planner/executor with dry-run default.", use_when=["需要从 manifest 到报告的端到端规划或执行"], do_not_use_when=["不应绕过 Doctor 或 approval 直接执行重型任务"], required_inputs=["sample_manifest"], optional_inputs=["tools_manifest", "reference_manifest"], outputs=["pipeline_plan.md", "pipeline_status.tsv", "run_manifest.json"], risk_level="HIGH", approval_required=True),
    _spec(name="neoag-release-qc", category="D", handler="release_qc", purpose="发布边界审计", description="Scan for caches, runtime artifacts, hard-coded private paths and patient/site hints.", use_when=["生成发布包前", "用户问 release 是否干净"], do_not_use_when=["不要自动删除文件，除非用户明确确认"], required_inputs=["project_root"], outputs=["release_audit_report.md", "release_boundary_findings.tsv"], risk_level="LOW"),
    _spec(name="neoag-gateway-submit", category="D", handler="gateway_submit", purpose="Gateway 受控提交", description="Submit low/medium/high-risk skill or pipeline requests to NeoAg Gateway with approval controls.", use_when=["Agent 需要通过 Gateway 受控执行任务"], do_not_use_when=["不要直接执行 shell 或绕过 approval"], required_inputs=["gateway_url", "task"], outputs=["gateway_job.json", "gateway_submission.md"], risk_level="HIGH", approval_required=True),
    _spec(
        name="neoag-remote-deploy", category="D", handler="remote_deploy",
        purpose="新机器迁移部署专用 SOP",
        description="Deploy or migrate Project B to a new machine with checksum/unpack, preflight, core install, smoke tests, local manifests, Doctor, dry-run planning and deployment report.",
        use_when=["用户要求把项目迁移到新机器", "需要让目标机器上的编程 agent 快速部署、安装、检查和验收", "需要生成 deployment_report.md 或判断 READY/PARTIAL/BLOCKED/UNSAFE"],
        do_not_use_when=["用户只是分析已有 ranked_peptides 或报告结果", "尚未完成 Doctor 和 dry-run 时不应直接运行患者级生产流程"],
        required_inputs=["project_root 或 release_tarball"],
        optional_inputs=["sha256", "deployment_tier", "tools_manifest", "reference_manifest", "sample_manifest"],
        outputs=["deployment_layout.md", "preflight_report.md", "smoke_test_report.md", "deployment_report.md", "doctor_summary.md", "local manifests"],
        risk_level="LOW_TO_HIGH_BY_STEP", approval_required=True,
        downstream_skills=["neoag-doctor", "neoag-tool-reference-qc", "neoag-run-demo-and-smoke", "neoag-pipeline-full", "neoag-release-qc"],
        boundaries=[
            "默认只做核心部署、smoke、Doctor 和 dry-run，不自动安装授权工具或下载大型参考库。",
            "不得把患者数据、本机私有路径或 license 工具打进 release。",
            "check-tools/which 不能等同于真实工具可用；需要 Doctor mini smoke 或标记 PARTIAL。",
            "HPC 提交、删除、覆盖、安装工具、下载参考库必须人工确认。",
        ],
    ),
    _spec(name="neoag-hpc-runner", category="D", handler="hpc_runner", purpose="HPC dry-run/job wrapper", description="Create Slurm/SGE/PBS job manifests and submit only after explicit approval.", use_when=["需要在 HPC 上运行 pipeline-full 或重型工具"], do_not_use_when=["默认不得直接提交，应先 dry-run"], required_inputs=["job_manifest"], outputs=["hpc_dry_run.sh", "hpc_job_manifest.json", "hpc_submission.md"], risk_level="HIGH", approval_required=True),
]

SKILLS_BY_NAME = {s.name: s for s in SKILL_SPECS}


def list_specs(category: str | None = None) -> list[SkillSpec]:
    if category:
        return [s for s in SKILL_SPECS if s.category == category]
    return list(SKILL_SPECS)


def registry_dict() -> dict[str, Any]:
    return {
        "schema_version": "neoag-skill-taxonomy-v1",
        "categories": CATEGORY_LABELS,
        "skills": [s.to_dict() for s in SKILL_SPECS],
    }
