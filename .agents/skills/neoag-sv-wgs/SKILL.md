---
name: neoag-sv-wgs
description: Parse WGS SV VCFs (Manta/GRIDSS/SvABA/DELLY) into SV events requiring transcript/protein reconstruction.
category: A - 入口适配型 Skills：把不同来源输入转换为 Project B 标准 raw_events/raw_peptides/evidence tables
risk_level: LOW
approval_required: false
---

# neoag-sv-wgs

## 目标

WGS DNA SV 入口标准化

## 什么时候使用
- 用户提供 WGS SV VCF
- 需要进入 DNA SV-derived neoantigen workflow

## 什么时候不要使用
- 输入来自 WES/capture-limited SV，应调用 neoag-sv-wes

## 必需输入
- `sv_vcf`

## 可选输入
- `gtf`
- `reference_fasta`
- `sample_id`
- `rna_junctions`

## 输出
- `sv_events.tsv`
- `raw_events.tsv`
- `raw_peptides.tsv`
- `sv_reconstruction_tasks.tsv`

## 运行示例

```bash
neoag-skill run neoag-sv-wgs --outdir work/neoag-sv-wgs --dry-run
```

## 边界
- Skill 不承担临床决策；不得判断患者是否适合治疗或推荐临床用药。
- 缺失证据只能标记为 missing/unassessed，不能解释为阴性结果。
- 高风险写入、HPC 提交、安装工具、下载参考库、删除或覆盖文件必须经过 human approval。
- Skill 目录不包含患者 BAM/FASTQ/VCF、大型参考库、VEP cache、NetMHCpan license、LOHHLA reference 或大型 conda env。

## 推荐下游 Skill
- `neoag-presentation`
- `neoag-experiment-design`
