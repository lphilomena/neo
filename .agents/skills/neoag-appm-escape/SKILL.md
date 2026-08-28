---
name: neoag-appm-escape
description: Aggregate APPM gene status, HLA LOH and immune escape flags into pathway-level and peptide-level modifiers.
category: B - 公共证据分析型 Skills：对所有入口共用的 HLA、表达、CCF、APPM、安全和排序证据层进行标准化分析
risk_level: LOW
approval_required: false
---

# neoag-appm-escape

## 目标

APPM / HLA LOH / immune escape 证据层

## 什么时候使用
- 需要 APPM、HLA LOH、B2M/TAP/JAK/NLRC5/CIITA 风险评估

## 什么时候不要使用
- 缺所有 APPM 输入时，只能输出 unassessed，不得输出 intact

## 必需输入
- `gene_status_or_appm`
- `peptides_optional`
- `raw_events`（保留 APPM 通路变异的事件 ID、WGS/WES VAF、深度和 ALT 读段）

## 可选输入
- `hla_loh`
- `ranked_peptides`
- `ccf_2`（按 `event_id` 关联 CCF；缺失时必须标记为未评估）

## 输出
- `appm_summary.tsv`
- `appm_gene_status.tsv`
- `appm_variant_evidence.tsv`
- `appm_peptide_modifiers.tsv`
- `immune_escape_events.tsv`
- `immune_escape_summary.tsv`
- `peptide_escape_flags.tsv`

## 运行示例

```bash
neoag-skill run neoag-appm-escape --outdir work/neoag-appm-escape --dry-run
```

## 边界
- Skill 不承担临床决策；不得判断患者是否适合治疗或推荐临床用药。
- 缺失证据只能标记为 missing/unassessed，不能解释为阴性结果。
- `appm_variant_evidence.tsv` 必须保留具体变异、WGS/WES VAF、深度、ALT读段和CCF；`appm_gene_status.tsv` 必须用 `source_event_ids` 回链。
- 没有调用变异时写 `NO_CALLED_VARIANT`，不得表述为该基因所有位点均无变异或均有足够覆盖。
- 高风险写入、HPC 提交、安装工具、下载参考库、删除或覆盖文件必须经过 human approval。
- Skill 目录不包含患者 BAM/FASTQ/VCF、大型参考库、VEP cache、NetMHCpan license、LOHHLA reference 或大型 conda env。
