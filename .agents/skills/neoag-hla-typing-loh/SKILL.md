---
name: neoag-hla-typing-loh
description: Normalize HLA typing and HLA LOH outputs from OptiType, SpecHLA, HLA-LA/HD and LOHHLA into consensus tables.
category: B - 公共证据分析型 Skills：对所有入口共用的 HLA、表达、CCF、APPM、安全和排序证据层进行标准化分析
risk_level: LOW
approval_required: false
---

# neoag-hla-typing-loh

## 目标

HLA typing / HLA LOH 共识与 peptide-level HLA loss flags

## 什么时候使用
- 需要标准化 HLA 分型
- 需要判断 restricting HLA 是否 LOH

## 什么时候不要使用
- 只需要解释已有 ranking 差异，不需要更新 HLA 状态

## 必需输入
- `hla`

## 可选输入
- `hla_loh`
- `ranked_peptides`
- `sample_id`

## 输出
- `hla_typing.normalized.tsv`
- `hla_typing_consensus.tsv`
- `hla_loh_consensus.tsv`
- `restricting_hla_peptide_flags.tsv`
- `hla_review.md`

## 运行示例

```bash
neoag-skill run neoag-hla-typing-loh --outdir work/neoag-hla-typing-loh --dry-run
```

## 边界
- Skill 不承担临床决策；不得判断患者是否适合治疗或推荐临床用药。
- 缺失证据只能标记为 missing/unassessed，不能解释为阴性结果。
- 高风险写入、HPC 提交、安装工具、下载参考库、删除或覆盖文件必须经过 human approval。
- Skill 目录不包含患者 BAM/FASTQ/VCF、大型参考库、VEP cache、NetMHCpan license、LOHHLA reference 或大型 conda env。

## 底层工具
- OptiType
- SpecHLA
- LOHHLA
- HLA-LA/HLA-HD optional
