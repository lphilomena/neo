---
name: neoag-experiment-design
description: Assign event-deduplicated representatives to short peptide, WT control, long peptide, minigene and targeted RNA validation routes.
category: C - 审阅/报告/实验设计型 Skills：解释结果、生成报告、设计实验验证和患者沟通材料
risk_level: LOW
approval_required: false
---

# neoag-experiment-design

## 目标

候选实验验证设计

## 什么时候使用
- 需要第一批 10–20 个实验候选
- 需要区分 short peptide / long peptide / minigene / targeted RNA

## 什么时候不要使用
- 不要把设计建议写成治疗处方

## 必需输入
- `ranked_events_or_ranked_peptides`

## 可选输入
- `ranked_events`
- `ranked_peptides`
- `top_n`
- `therapy_context`

优先使用 `ranked_events.evidence_consensus.tsv`。该表已经按事件/HLA、
`phase_group_id` 和 `redundancy_group` 去重，并且每个事件最多提供两个代表
peptide-HLA。只有事件表不存在时才允许回退到 peptide 表，回退状态会记录为
`ranked_peptides_fallback`。

生成首批候选时先从每个事件各取代表 1；只有候选名额仍有剩余时，才按事件
排名顺序回填代表 2。这样可优先增加事件覆盖度，避免少数事件占满前 20 位。

## 输出
- `experiment_candidates.tsv`
- `short_peptide_pool.tsv`
- `long_peptide_design.tsv`
- `minigene_design.tsv`
- `targeted_rna_validation_plan.md`

## 运行示例

```bash
neoag-skill run neoag-experiment-design --outdir work/neoag-experiment-design --dry-run
```

## 边界
- Skill 不承担临床决策；不得判断患者是否适合治疗或推荐临床用药。
- 缺失证据只能标记为 missing/unassessed，不能解释为阴性结果。
- 高风险写入、HPC 提交、安装工具、下载参考库、删除或覆盖文件必须经过 human approval。
- Skill 目录不包含患者 BAM/FASTQ/VCF、大型参考库、VEP cache、NetMHCpan license、LOHHLA reference 或大型 conda env。
