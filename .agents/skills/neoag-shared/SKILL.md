---
name: neoag-shared
description: NeoAg流程中presentation之后的公共段说明（APPM/CCF/安全性/免疫逃逸/打分/验证计划/报告）。供neoag-vcf、neoag-fusion、neoag-splice、neoag-sv-wgs、neoag-sv-wes、neoag-peptide-csv这6个入口skill内部引用，不建议作为用户直接触发的独立入口。
---

# 公共段：presentation之后的共用步骤

6个入口skill（VCF/Fusion/Splice/SV-WGS/SV-WES/Peptide-CSV）各自产出
`raw_events.tsv` + `raw_peptides.tsv` 之后，从这里开始完全共用同一套命令。

## 前置条件

已有：
- `<outdir>/upstream/parsed/raw_events.tsv`
- `<outdir>/upstream/parsed/raw_peptides.tsv`
- `<outdir>/presentation/`下的结合预测证据（来自各入口skill的`peptide-predict`步骤）

## 1. APPM 2.0（抗原呈递机器完整性）

```bash
neoag-v03 appm-2 \
  --sample-id <sample_id> \
  --vep-tsv <vep_appm_tsv> \
  --expression <gene_expression_tsv> \
  --hla-loh <hla_loh_tsv> \
  --cnv <cnv_segments_tsv> \
  --raw-peptides <outdir>/upstream/parsed/raw_peptides.tsv \
  --tumor-purity <purity_tsv> \
  --profile <profile> \
  --outdir <outdir>/appm
```

VEP/表达/HLA LOH/CNV/纯度任一缺失都可以省略对应参数，APPM会记录为"证据不完整"而不是报错中断。

check：`<outdir>/appm/appm_summary.tsv` 非空。

## 2. CCF 2.1（克隆性）

```bash
neoag-v03 ccf-2 \
  --events <outdir>/upstream/parsed/raw_events.tsv \
  --purity <purity_tsv> \
  --cnv <cnv_segments_tsv> \
  --profile <profile> \
  --out <outdir>/clonality/ccf_2.tsv
```

check：`<outdir>/clonality/ccf_2.tsv` 非空；若无purity/cnv，输出会带`WES_SV_CAPTURE_LIMITED_APPROX`或
`RNA_ONLY_UNRESOLVED`等方法标签，属正常降级，不算失败。

## 3. 肽段安全性门控

```bash
neoag-v03 peptide-safety \
  --raw-events <outdir>/upstream/parsed/raw_events.tsv \
  --raw-peptides <outdir>/upstream/parsed/raw_peptides.tsv \
  --profile <profile> \
  --normal-expression <normal_expression_tsv> \
  --normal-hla-ligands <normal_hla_ligands_tsv> \
  --out <outdir>/safety/peptide_safety.tsv \
  --event-out <outdir>/safety/event_safety.tsv
```

check：`<outdir>/safety/peptide_safety.tsv` 非空。这一步是计算层面的off-target风险降低，
不等于临床安全性证明，报告里要保留这条边界说明。

## 4. 免疫逃逸证据

```bash
neoag-v03 immune-escape \
  --sample-id <sample_id> \
  --raw-peptides <outdir>/upstream/parsed/raw_peptides.tsv \
  --profile <profile> \
  --vep-tsv <vep_appm_tsv> \
  --cnv <cnv_segments_tsv> \
  --expression <gene_expression_tsv> \
  --hla-loh <hla_loh_tsv> \
  --appm-gene-status <outdir>/appm/appm_gene_status.tsv \
  --appm-pathway-status <outdir>/appm/appm_pathway_status.tsv \
  --ccf <outdir>/clonality/ccf_2.tsv \
  --therapy-context discovery \
  --outdir <outdir>/immune_escape
```

`--therapy-context` 按用途选 `vaccine`/`tcr_target`/`immunomonitoring`/`discovery`，影响多重/优先级封顶策略。

check：`<outdir>/immune_escape/immune_escape_summary.tsv` 非空。这一步是机制/风险层面证据，
不是临床耐药诊断，报告里要保留这条边界说明。

以上3-4两步可以和APPM/CCF并行跑（互相之间只有弱依赖，APPM的gene/pathway status建议先跑完APPM再传给immune-escape）。

## 5. 打分排序

```bash
neoag-v03 score-v03 \
  --raw-events <outdir>/upstream/parsed/raw_events.tsv \
  --raw-peptides <outdir>/upstream/parsed/raw_peptides.tsv \
  --presentation <outdir>/presentation/presentation_evidence.tsv \
  --appm-summary <outdir>/appm/appm_summary.tsv \
  --ccf <outdir>/clonality/ccf_2.tsv \
  --normal-expression <normal_expression_tsv> \
  --normal-hla-ligands <normal_hla_ligands_tsv> \
  --peptide-safety <outdir>/safety/peptide_safety.tsv \
  --peptide-escape-flags <outdir>/immune_escape/peptide_escape_flags.tsv \
  --profile <profile> \
  --out-events <outdir>/scoring/ranked_events.v03.tsv \
  --out-peptides <outdir>/scoring/ranked_peptides.v03.tsv
```

check：两个输出文件都非空，且行数与`raw_events`/`raw_peptides`量级相符（大幅缩水通常说明上游某个证据文件是空的）。

## 6. 验证计划

```bash
neoag-v03 validation-plan-v03 \
  --ranked-peptides <outdir>/scoring/ranked_peptides.v03.tsv \
  --outdir <outdir> \
  --out <outdir>/scoring/validation_plan.v03.tsv
```

check：`validation_plan.v03.tsv` 非空。

## 7. 报告

```bash
neoag-v03 report-v03 \
  --profile <profile> \
  --ranked-events <outdir>/scoring/ranked_events.v03.tsv \
  --ranked-peptides <outdir>/scoring/ranked_peptides.v03.tsv \
  --appm-summary <outdir>/appm/appm_summary.tsv \
  --validation-plan <outdir>/scoring/validation_plan.v03.tsv \
  --outdir <outdir> \
  --sample-id <sample_id> \
  --audience both \
  --out <outdir>/reports/evidence_report.v03.html
```

`--audience both` 会同时生成 `evidence_report.patient.html`（患者沟通版）和
`evidence_report.technical.html`（研究/技术版，含dashboard和provenance）。

## 失败处理

| 阶段 | 失败时只重跑这一步 | 常见原因 |
|---|---|---|
| appm-2 | 是 | vep-tsv/expression/hla-loh路径错误或格式不对 |
| ccf-2 | 是 | purity/cnv缺失（可接受降级）或events表字段不全 |
| peptide-safety | 是 | normal_expression/normal_hla_ligands路径错误 |
| immune-escape | 是 | 依赖appm/ccf的输出，需确认这两步先跑完 |
| score-v03 | 是，但要检查所有依赖文件都非空 | 任一上游sidecar为空会导致打分结果异常缩水 |
| validation-plan-v03 | 是 | ranked_peptides为空 |
| report-v03 | 是 | outdir下sidecar路径不全，技术版报告部分卡片会缺失（不算致命错误，但要向用户说明） |

## 执行状态汇总（每个入口skill在最后都应输出）

- 已完成步骤清单 + 各步输出路径
- 中断步骤：具体哪一步、check不通过的原因
- 成功时：`ranked_peptides.v03.tsv`行数、`evidence_report.v03.html`路径、
  免疫逃逸/安全性标记的简要统计（如高风险肽段数量）
