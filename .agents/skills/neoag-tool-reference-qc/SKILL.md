---
name: neoag-tool-reference-qc
description: Check external tool entrypoints, models, caches and references; stronger than simple PATH check.
category: D - 工程治理/执行控制型 Skills：输入质控、环境健康检查、全流程编排、发布审计和受控执行
risk_level: LOW
approval_required: false
---

# neoag-tool-reference-qc

## 目标

工具/参考库检查

## 什么时候使用
- 需要检查 VEP/NetMHCpan/MHCflurry/LOHHLA/FACETS 等

## 什么时候不要使用
- 不要将工具 missing 解释为生物学阴性

## 必需输入
- `project_root`

## 可选输入
- `tools_manifest`
- `reference_manifest`

## 输出
- `tool_qc_report.tsv`
- `reference_qc_report.tsv`
- `tool_smoke_report.md`

## 运行示例

```bash
neoag-skill run neoag-tool-reference-qc --outdir work/neoag-tool-reference-qc --dry-run
```

## 边界
- Skill 不承担临床决策；不得判断患者是否适合治疗或推荐临床用药。
- 缺失证据只能标记为 missing/unassessed，不能解释为阴性结果。
- 高风险写入、HPC 提交、安装工具、下载参考库、删除或覆盖文件必须经过 human approval。
- Skill 目录不包含患者 BAM/FASTQ/VCF、大型参考库、VEP cache、NetMHCpan license、LOHHLA reference 或大型 conda env。
