---
name: neoag-safety
description: Check normal proteome exact matches, WT peptide flags, anchor-only risks and normal junction/ligandome overlap.
category: B - 公共证据分析型 Skills：对所有入口共用的 HLA、表达、CCF、APPM、安全和排序证据层进行标准化分析
risk_level: LOW
approval_required: false
---

# neoag-safety

## 目标

peptide safety gate

## 什么时候使用
- 需要筛查正常蛋白组匹配、WT peptide、anchor-only、normal junction 风险

## 什么时候不要使用
- 不能替代湿实验安全性评估

## 必需输入
- `raw_peptides_or_ranked_peptides`

## 可选输入
- `normal_proteome`
- `normal_junctions`
- `wt_binding`

## 输出
- `peptide_safety.tsv`
- `event_safety.tsv`
- `safety_review.tsv`

## 运行示例

```bash
neoag-skill run neoag-safety --outdir work/neoag-safety --dry-run
```

## 边界
- Skill 不承担临床决策；不得判断患者是否适合治疗或推荐临床用药。
- 缺失证据只能标记为 missing/unassessed，不能解释为阴性结果。
- 高风险写入、HPC 提交、安装工具、下载参考库、删除或覆盖文件必须经过 human approval。
- Skill 目录不包含患者 BAM/FASTQ/VCF、大型参考库、VEP cache、NetMHCpan license、LOHHLA reference 或大型 conda env。
