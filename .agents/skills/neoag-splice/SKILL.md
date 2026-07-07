---
name: neoag-splice
description: 从可变剪接事件（VCF + RegTools splice junction TSV）出发，运行剪接肽段生成、结合/免疫原性预测、APPM/CCF/安全性/免疫逃逸证据构建、打分排序与报告生成。在以下情况自动调用本skill：1. 用户提到剪接、splice junction、可变剪接新抗原、RegTools、pVACsplice；2. 用户提供了splice junction相关TSV文件。
---

# NeoAg Splice Junction（可变剪接）新抗原分析

## 使用场景

输入是变异VCF + RegTools产出的splice junction TSV，用于捕获可变剪接产生的新抗原候选。

## 必需输入

- `sample_id`
- `splice_junction_tsv`（RegTools或兼容格式的junction表）
- 可选：`pvacsplice_aggregated`（pVACsplice结果）
- `hla_alleles`
- `outdir`
- `profile`

## 运行前检查

```bash
neoag-v03 run-demo --entry-mode splice_junction --outdir /tmp/neoag_demo_splice --sample-id DEMO_SPLICE
```

用仓库自带fixture跑一遍完整链路验证代码可用，不需要真实RegTools/pVACsplice环境。

## 分步执行路径

### 1. 剪接事件标准化

```bash
neoag-v03 build-intermediates \
  --outdir <outdir> \
  --sample-id <sample_id> \
  --profile <profile> \
  --entry-mode splice_junction \
  --pvac <pvacsplice_aggregated_tsv> \
  --splice-junction-tsv <splice_junction_tsv> \
  --hla <hla_allele_1> <hla_allele_2> ...
```

`--pvac`（pVACsplice聚合表）是**真正产生肽段的来源**——`--splice-junction-tsv`只是给已有肽段做
junction支持度富集，本身不产肽，缺了`--pvac`肽段生成会静默产出0条。`--variants-vcf`是可选的富集上下文，
同样不会独立产肽。

check：

```bash
test -s <outdir>/parsed/raw_events.tsv
test -s <outdir>/parsed/raw_peptides.tsv
wc -l <outdir>/parsed/raw_peptides.tsv   # 确认>1行，不只是表头
```

不合格 → 检查junction TSV的列格式是否符合RegTools标准输出。

### 2. 结合/免疫原性预测

```bash
neoag-v03 peptide-predict \
  -i <outdir>/parsed/raw_peptides.tsv \
  -o <outdir>/presentation \
  --sample-id <sample_id>
```

### 3. 公共段

参数和check标准见 `../neoag-shared/SKILL.md`。

## 关键输出

同其余入口：`<outdir>/scoring/ranked_peptides.v03.tsv`、`<outdir>/reports/evidence_report.v03.html`等。

## 失败处理

- 标准化阶段失败：检查junction TSV列名、是否有对应的变异VCF匹配。
- 公共段失败：见 `../neoag-shared/SKILL.md`。

## 参考

- 公共段完整说明：`../neoag-shared/SKILL.md`
