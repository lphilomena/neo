---
name: neoag-pipeline-full
description: controlled-execution explicit pipeline-full DAG planner/executor with dry-run default.
category: D - 工程治理/执行控制型 Skills：输入质控、环境健康检查、全流程编排、发布审计和受控执行
risk_level: HIGH
approval_required: true
---

# neoag-pipeline-full

## 目标

manifest-driven full pipeline runner；上游 HLA、纯度/CNV、HLA LOH 完成后接入 VCF 排序

## 什么时候使用
- 需要从 manifest 到报告的端到端规划或执行

## 什么时候不要使用
- 不应绕过 Doctor 或 approval 直接执行重型任务

## 必需输入
- `sample_manifest`
- Step 0 `input_status.json` 必须不是 FAIL
- `recommended_hla.txt`
- `recommended_purity.tsv`
- `recommended_cnv_segments.tsv`
- `recommended_hla_loh.tsv`，或明确标记为 UNASSESSED

## 可选输入
- `tools_manifest`
- `reference_manifest`

## 输出
- `pipeline_plan.md`
- `pipeline_status.tsv`
- `run_manifest.json`

## 运行示例

```bash
neoag-skill run neoag-pipeline-full --outdir work/neoag-pipeline-full --dry-run
```

## 边界
- Skill 不承担临床决策；不得判断患者是否适合治疗或推荐临床用药。
- 缺失证据只能标记为 missing/unassessed，不能解释为阴性结果。
- 高风险写入、HPC 提交、安装工具、下载参考库、删除或覆盖文件必须经过 human approval。
- Skill 目录不包含患者 BAM/FASTQ/VCF、大型参考库、VEP cache、NetMHCpan license、LOHHLA reference 或大型 conda env。
- `pipeline-full` 不得把缺失的 HLA、纯度、CNV、RNA 或 HLA LOH 自动解释为阴性。
- 没有肿瘤 RNA 时必须采用 DNA-only profile，并将 RNA ALT/VAF/TPM 写为 `UNASSESSED`。
- 生产执行前必须检查共享工作盘剩余空间；不得在空间不足的 `/root` 中生成 WGS 中间文件。
