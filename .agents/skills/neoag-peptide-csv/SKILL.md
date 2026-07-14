---
name: neoag-peptide-csv
description: Normalize an existing peptide-HLA table into raw_peptides and optional presentation_evidence tables.
category: A - 入口适配型 Skills：把不同来源输入转换为 Project B 标准 raw_events/raw_peptides/evidence tables
risk_level: LOW
approval_required: false
---

# neoag-peptide-csv

## 目标

已有 peptide-HLA 表入口标准化

## 什么时候使用
- 用户已有 peptide-HLA 候选表
- 需要将外部排序/预测表接入 Project B scoring

## 什么时候不要使用
- 需要从 VCF/fusion/SV 生成候选时

## 必需输入
- `peptide_csv`

## 可选输入
- `sample_id`
- `event_annotation`

## 输出
- `raw_peptides.tsv`
- `presentation_evidence.tsv`
- `peptide_input_qc.tsv`

## 运行示例

```bash
neoag-skill run neoag-peptide-csv --outdir work/neoag-peptide-csv --dry-run
```

## 边界
- Skill 不承担临床决策；不得判断患者是否适合治疗或推荐临床用药。
- 缺失证据只能标记为 missing/unassessed，不能解释为阴性结果。
- 高风险写入、HPC 提交、安装工具、下载参考库、删除或覆盖文件必须经过 human approval。
- Skill 目录不包含患者 BAM/FASTQ/VCF、大型参考库、VEP cache、NetMHCpan license、LOHHLA reference 或大型 conda env。

## 推荐下游 Skill
- `neoag-safety`
- `neoag-ranking`
