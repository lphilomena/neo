---
name: neoag-presentation
description: Normalize NetMHCpan/MHCflurry/PRIME/BigMHC/MixMHCpred predictions into presentation_evidence.tsv.
category: B - 公共证据分析型 Skills：对所有入口共用的 HLA、表达、CCF、APPM、安全和排序证据层进行标准化分析
risk_level: LOW
approval_required: false
---

# neoag-presentation

## 目标

HLA binding / presentation 证据标准化

## 什么时候使用
- 已有 prediction 表需要标准化
- raw_peptides 需要标记 presentation 计算计划

## 什么时候不要使用
- 需要解释综合推荐排序，应调用 neoag-ranking-compare

## 必需输入
- `predictions_or_raw_peptides`

## 可选输入
- `hla`
- `sample_id`

## 输出
- `presentation_evidence.tsv`
- `presentation_summary.tsv`
- `presentation_qc.tsv`

## 运行示例

```bash
neoag-skill run neoag-presentation --outdir work/neoag-presentation --dry-run
```

## 边界
- Skill 不承担临床决策；不得判断患者是否适合治疗或推荐临床用药。
- 缺失证据只能标记为 missing/unassessed，不能解释为阴性结果。
- 高风险写入、HPC 提交、安装工具、下载参考库、删除或覆盖文件必须经过 human approval。
- Skill 目录不包含患者 BAM/FASTQ/VCF、大型参考库、VEP cache、NetMHCpan license、LOHHLA reference 或大型 conda env。

## 底层工具
- NetMHCpan
- MHCflurry
- PRIME
- BigMHC
- MixMHCpred
