---
name: neoag-ranking
description: Integrate presentation, expression, RNA evidence, CCF, APPM, escape, safety and event confidence into recommendation ranking.
category: B - 公共证据分析型 Skills：对所有入口共用的 HLA、表达、CCF、APPM、安全和排序证据层进行标准化分析
risk_level: LOW
approval_required: false
---

# neoag-ranking

## 目标

综合推荐排序生成

## 什么时候使用
- 需要生成 final_priority 和 recommendation ranked peptides

## 什么时候不要使用
- 只需比较两个排序文件时，使用 neoag-ranking-compare

## 必需输入
- `raw_peptides`
- `presentation`

## 可选输入
- `expression`
- `rna_evidence`
- `ccf`
- `appm`
- `safety`

## 输出
- `ranked_peptides.recommendation.tsv`
- `ranked_peptides.netmhcpan42.tsv`
- `ranked_events.tsv`
- `validation_plan.tsv`

## 运行示例

```bash
neoag-skill run neoag-ranking --outdir work/neoag-ranking --dry-run
```

## 边界
- Skill 不承担临床决策；不得判断患者是否适合治疗或推荐临床用药。
- 缺失证据只能标记为 missing/unassessed，不能解释为阴性结果。
- 高风险写入、HPC 提交、安装工具、下载参考库、删除或覆盖文件必须经过 human approval。
- Skill 目录不包含患者 BAM/FASTQ/VCF、大型参考库、VEP cache、NetMHCpan license、LOHHLA reference 或大型 conda env。
