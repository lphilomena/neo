---
name: neoag-expression
description: Normalize WTS/gene expression TPM tables into expression_evidence.tsv with expressed/low/not_detected labels.
category: B - 公共证据分析型 Skills：对所有入口共用的 HLA、表达、CCF、APPM、安全和排序证据层进行标准化分析
risk_level: LOW
approval_required: false
---

# neoag-expression

## 目标

gene expression/TPM 证据标准化

## 什么时候使用
- 用户提供 TPM/gene expression 表
- 需要把表达证据接入 ranking

## 什么时候不要使用
- 需要 RNA alt reads 或 junction reads，使用 neoag-rna-evidence

## 必需输入
- `expression_tsv`

## 可选输入
- `sample_id`

## 输出
- `expression_evidence.tsv`
- `expression_qc.tsv`

## 运行示例

```bash
neoag-skill run neoag-expression --outdir work/neoag-expression --dry-run
```

## 边界
- Skill 不承担临床决策；不得判断患者是否适合治疗或推荐临床用药。
- 缺失证据只能标记为 missing/unassessed，不能解释为阴性结果。
- 高风险写入、HPC 提交、安装工具、下载参考库、删除或覆盖文件必须经过 human approval。
- Skill 目录不包含患者 BAM/FASTQ/VCF、大型参考库、VEP cache、NetMHCpan license、LOHHLA reference 或大型 conda env。
