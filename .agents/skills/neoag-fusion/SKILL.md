---
name: neoag-fusion
description: 从EasyFuse融合检测结果（fusions.pass.csv）出发，运行融合肽段生成、结合/免疫原性预测、APPM/CCF/安全性/免疫逃逸证据构建、打分排序与报告生成。在以下情况自动调用本skill：1. 用户提供或提到fusions.pass.csv、EasyFuse、STAR-Fusion、pVACfuse等融合检测相关文件或工具；2. 用户想跑融合基因新抗原分析。
---

# NeoAg Fusion（融合基因）新抗原分析

## 使用场景

输入是EasyFuse（或兼容格式）产出的融合事件表 `fusions.pass.csv`，可选再补一份pVACfuse聚合结果补充HLA分型表位。

## 必需输入

- `sample_id`
- `easyfuse_pass_csv`（EasyFuse `fusions.pass.csv` 路径）
- `hla_alleles`
- `outdir`
- `profile`
- 可选：`pvacfuse_aggregated`（pVACfuse结果，用于补充完整HLA分型表位）
- 可选：`normal_expression`、`normal_hla_ligands`

若用户还没跑EasyFuse（只有RNA FASTQ），要说明EasyFuse本身是外部Nextflow流程
（`bin/easyfuse-neoag`），本skill只从`fusions.pass.csv`开始接手，不负责跑EasyFuse本身。

## 运行前检查

```bash
neoag-v03 run-demo --entry-mode fusion --outdir /tmp/neoag_demo_fusion --sample-id DEMO_FUSION
```

该命令用仓库自带的`pvacfuse_aggregated.tsv`等fixture跑一遍完整链路，验证代码可用，
无需真实EasyFuse/pVACfuse环境。

若要跑真实EasyFuse上游，检查：

```bash
command -v bin/easyfuse-neoag >/dev/null || echo "MISSING: EasyFuse wrapper"
```

## 分步执行路径

### 1. 融合事件标准化

```bash
neoag-v03 build-intermediates \
  --outdir <outdir> \
  --sample-id <sample_id> \
  --profile <profile> \
  --entry-mode fusion \
  --easyfuse-pass-csv <easyfuse_pass_csv> \
  --hla <hla_allele_1> <hla_allele_2> ...
```

`--hla` 是**必需参数**，不是可选项——缺了它肽段生成会静默产出0条（只有事件表，没有肽段表，不会报错）。
若同时有pVACfuse结果，追加 `--pvac <pvacfuse_aggregated.tsv>` 补充表位信息。

check：

```bash
test -s <outdir>/parsed/raw_events.tsv
test -s <outdir>/parsed/raw_peptides.tsv
# 关键：确认raw_peptides.tsv行数>1，不只是文件存在
wc -l <outdir>/parsed/raw_peptides.tsv
```

不合格 → 检查 `fusions.pass.csv` 是否包含 `FTID`/`Fusion_Gene` 等EasyFuse标准列。

### 2. 结合/免疫原性预测

```bash
neoag-v03 peptide-predict \
  -i <outdir>/parsed/raw_peptides.tsv \
  -o <outdir>/presentation \
  --sample-id <sample_id>
```

### 3. 公共段（APPM/CCF/安全性/免疫逃逸/打分/验证计划/报告）

参数和check标准见 `../neoag-shared/SKILL.md`，直接复用，无需修改。

## 关键输出

同 `neoag-vcf`：`<outdir>/scoring/ranked_peptides.v03.tsv`、`<outdir>/reports/evidence_report.v03.html`等。

## 失败处理

- 标准化阶段失败：检查 `fusions.pass.csv` 列名/编码是否符合EasyFuse标准输出格式。
- 公共段失败：见 `../neoag-shared/SKILL.md`。
- 中断/完成都要向用户汇总已完成步骤、产出文件、失败原因（若有）。

## 参考

- 公共段完整说明：`../neoag-shared/SKILL.md`
