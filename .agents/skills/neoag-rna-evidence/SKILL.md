---
name: neoag-rna-evidence
description: Normalize RNA allele, RNA VAF and fusion/splice junction evidence into event-level RNA support tables.
category: B - 公共证据分析型 Skills：对所有入口共用的 HLA、表达、CCF、APPM、安全和排序证据层进行标准化分析
risk_level: LOW
approval_required: false
---

# neoag-rna-evidence

## 目标

RNA alt / RNA VAF / junction reads 证据标准化

## 什么时候使用
- 用户提供 RNA alt reads、RNA VAF、fusion/splice junction reads

## 什么时候不要使用
- 只有 gene TPM 时，使用 neoag-expression

## 必需输入
- `rna_tsv`

## 可选输入
- `sample_id`

## 输出
- `rna_alt_evidence.tsv`
- `rna_junction_evidence.tsv`
- `rna_evidence_qc.tsv`

## 运行示例

```bash
neoag-skill run neoag-rna-evidence --outdir work/neoag-rna-evidence --dry-run
```

## 边界
- Skill 不承担临床决策；不得判断患者是否适合治疗或推荐临床用药。
- 缺失证据只能标记为 missing/unassessed，不能解释为阴性结果。
- 高风险写入、HPC 提交、安装工具、下载参考库、删除或覆盖文件必须经过 human approval。
- Skill 目录不包含患者 BAM/FASTQ/VCF、大型参考库、VEP cache、NetMHCpan license、LOHHLA reference 或大型 conda env。
