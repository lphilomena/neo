---
name: neoag-splice
description: Normalize RegTools/splice junction tables and optional variant VCF into splice-junction event and peptide candidate inputs.
category: A - 入口适配型 Skills：把不同来源输入转换为 Project B 标准 raw_events/raw_peptides/evidence tables
risk_level: LOW
approval_required: false
---

# neoag-splice

## 目标

Splice/junction 输入标准化

## 什么时候使用
- 用户提供 RegTools splice junction TSV
- 需要构建 splice/exon junction 候选

## 什么时候不要使用
- 仅有普通 SNV/InDel VCF 时，应调用 neoag-vcf

## 必需输入
- `junctions`

## 可选输入
- `vcf`
- `sample_id`

## 输出
- `splice_events.tsv`
- `raw_events.tsv`
- `raw_peptides.tsv`
- `splice_qc.tsv`

## 运行示例

```bash
neoag-skill run neoag-splice --outdir work/neoag-splice --dry-run
```

## 边界
- Skill 不承担临床决策；不得判断患者是否适合治疗或推荐临床用药。
- 缺失证据只能标记为 missing/unassessed，不能解释为阴性结果。
- 高风险写入、HPC 提交、安装工具、下载参考库、删除或覆盖文件必须经过 human approval。
- Skill 目录不包含患者 BAM/FASTQ/VCF、大型参考库、VEP cache、NetMHCpan license、LOHHLA reference 或大型 conda env。

## 推荐下游 Skill
- `neoag-rna-evidence`
- `neoag-experiment-design`
