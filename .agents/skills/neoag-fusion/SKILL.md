---
name: neoag-fusion
description: Normalize EasyFuse/Arriba/STAR-Fusion/FusionCatcher fusion outputs into raw_events/raw_peptides and fusion evidence tables.
category: A - 入口适配型 Skills：把不同来源输入转换为 Project B 标准 raw_events/raw_peptides/evidence tables
risk_level: LOW
approval_required: false
---

# neoag-fusion

## 目标

Fusion caller 输出标准化

## 什么时候使用
- 用户提供 fusions.pass.csv 或其他 fusion caller 输出
- 需要生成 fusion junction 候选输入

## 什么时候不要使用
- fusion 结果尚未经过 read-through 复核时直接作为临床结论

## 必需输入
- `fusion`

## 可选输入
- `sample_id`
- `normal_readthrough_db`

## 输出
- `fusion_events.tsv`
- `raw_events.tsv`
- `raw_peptides.tsv`
- `fusion_evidence.tsv`
- `fusion_qc.tsv`

## 运行示例

```bash
neoag-skill run neoag-fusion --outdir work/neoag-fusion --dry-run
```

## 边界
- Skill 不承担临床决策；不得判断患者是否适合治疗或推荐临床用药。
- 缺失证据只能标记为 missing/unassessed，不能解释为阴性结果。
- 高风险写入、HPC 提交、安装工具、下载参考库、删除或覆盖文件必须经过 human approval。
- Skill 目录不包含患者 BAM/FASTQ/VCF、大型参考库、VEP cache、NetMHCpan license、LOHHLA reference 或大型 conda env。

## 推荐下游 Skill
- `neoag-presentation`
- `neoag-experiment-design`
