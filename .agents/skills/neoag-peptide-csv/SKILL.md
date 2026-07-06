---
name: neoag-peptide-csv
description: 从已有的肽段-HLA表(CSV/TSV)出发，直接标准化后运行结合/免疫原性预测、APPM/CCF/安全性/免疫逃逸证据构建、打分排序与报告生成，无需变异层信息。在以下情况自动调用本skill：1. 用户提供的.csv/.tsv文件只包含肽段序列和HLA分型（没有基因组坐标/VCF等变异信息）；2. 用户明确说"这是一份肽段表，直接打分/排序"。
---

# NeoAg Peptide-only CSV（已有肽段表）新抗原分析

## 使用场景

用户已经有一份肽段-HLA表（不管来自何处），只需要标准化后接入结合预测和打分排序，
不需要从VCF/融合/SV重新产肽。这是6个入口里唯一没有独有产肽步骤的一个。

## 必需输入

- `sample_id`
- `input_csv`（含肽段序列列 + HLA分型列的CSV/TSV，具体列名不确定时先用`head`看一眼再问用户确认映射）
- `outdir`
- `profile`

## 运行前检查

```bash
neoag-v03 run-demo --entry-mode peptide_only --outdir /tmp/neoag_demo_pep --sample-id DEMO_PEP
```

本入口是纯Python标准化，不依赖任何外部生信工具，demo应该总是能跑通（除了`neoag-v03`本身需要装好）。

## 分步执行路径

### 1. 标准化

```bash
neoag-v03 build-intermediates \
  --outdir <outdir> \
  --sample-id <sample_id> \
  --profile <profile> \
  --entry-mode peptide_only \
  --peptide-table <input_csv>
```

check：

```bash
test -s <outdir>/parsed/raw_peptides.tsv
test -s <outdir>/parsed/raw_events.tsv
```

不合格 → 用`head -3 <input_csv>`确认列名/分隔符，向用户确认哪一列是肽段序列、哪一列是HLA，
必要时先用脚本重命名列后再转换。

注意：`raw_events.tsv` 会由 `raw_peptides.tsv` 自动合成一份"占位事件表"（按肽段所在基因生成一行，
VAF/CCF/克隆性等事件级字段是默认占位值，不是真实变异证据，`source`字段会标注为`peptide_input`）。
这不算失败，但在最终报告里要向用户说明：本次事件级字段（VAF、克隆性等）不是真实测序证据。

如果只需要标准化后的肽段-HLA表本身（不进入完整打分流程），也可以用更轻量的
`neoag-v03 convert-peptide-input -i <input_csv> -o <outdir>/parsed --sample-id <sample_id>`，
但它只产出`raw_peptides.tsv`/`peptide_hla_pairs.tsv`/`hla_alleles.txt`，不产出`raw_events.tsv`，
无法直接接后续打分链（`score-v03`需要`--raw-events`）。

### 2. 结合/免疫原性预测

```bash
neoag-v03 peptide-predict \
  -i <outdir>/parsed/raw_peptides.tsv \
  -o <outdir>/presentation \
  --sample-id <sample_id>
```

### 3. 公共段

参数和check标准见 `../neoag-shared/SKILL.md`。由于没有`raw_events.tsv`，
CCF/免疫逃逸中依赖事件级证据（VAF、CNV关联）的部分会标注证据不完整，不算失败。

## 关键输出

- `<outdir>/parsed/raw_peptides.tsv`
- `<outdir>/parsed/raw_events.tsv`（占位事件表，见上方注意事项）
- `<outdir>/scoring/ranked_peptides.v03.tsv`（主要参考这个）
- `<outdir>/scoring/ranked_events.v03.tsv`（事件级字段多为占位值，参考价值有限）
- `<outdir>/reports/evidence_report.v03.html`

## 参考

- 公共段完整说明：`../neoag-shared/SKILL.md`
