---
name: neoag-vcf
description: Parse a somatic VCF or VEP-annotated VCF into Project B raw_events/raw_peptides-compatible tables.
category: A - 入口适配型 Skills：把不同来源输入转换为 Project B 标准 raw_events/raw_peptides/evidence tables
risk_level: LOW
approval_required: false
---

# neoag-vcf

## 目标

VCF/SNV/InDel 入口标准化

## 什么时候使用
- 用户提供 tumor-normal somatic VCF
- 需要把 SNV/InDel 输入转成 raw_events.tsv

## 什么时候不要使用
- 输入是 fusion/splice/SV/peptide 表
- 需要直接生成综合排序，应调用 neoag-ranking

## 必需输入
- `vcf`

## 可选输入
- `sample_id`
- `hla`
- `expression_tsv`

## 输出
- `raw_events.tsv`
- `raw_peptides.tsv`
- `vcf_parse_qc.tsv`
- `candidate_generation_plan.md`

## 运行示例

```bash
neoag-skill run neoag-vcf --outdir work/neoag-vcf --dry-run
```

## 边界
- Skill 不承担临床决策；不得判断患者是否适合治疗或推荐临床用药。
- 缺失证据只能标记为 missing/unassessed，不能解释为阴性结果。
- 高风险写入、HPC 提交、安装工具、下载参考库、删除或覆盖文件必须经过 human approval。
- Skill 目录不包含患者 BAM/FASTQ/VCF、大型参考库、VEP cache、NetMHCpan license、LOHHLA reference 或大型 conda env。

## 推荐下游 Skill
- `neoag-presentation`
- `neoag-ranking`
