---
name: neoag-sv-wgs
description: 从全基因组结构变异VCF（WGS，Manta/GRIDSS/SvABA等）出发，重建融合蛋白滑窗肽段，运行结合/免疫原性预测、APPM/CCF/安全性/免疫逃逸证据构建、打分排序与报告生成。在以下情况自动调用本skill：1. 用户提供或提到SV VCF、结构变异、Manta/GRIDSS/SvABA/DELLY，且明确是全基因组(WGS)数据，没有提到capture BED/外显子捕获；2. 用户想跑WGS结构变异新抗原分析。
---

# NeoAg SV-WGS（全基因组结构变异）新抗原分析

## 使用场景

输入是全基因组SV VCF（Manta/GRIDSS/SvABA/DELLY等，可多个caller合并），无capture区间限制，
走标准SV Phase 1流程重建融合蛋白并滑窗产肽（8-11aa）。

如果用户的SV VCF来自外显子捕获(WES)测序、或明确提供了capture BED，改用 `neoag-sv-wes`，不要用本skill。

## 必需输入

- `sample_id`
- `sv_vcf`（一个或多个SV VCF路径）
- `callers`（对应的caller名，顺序需与`sv_vcf`一致，如 `Manta SvABA`）
- `reference_fasta`（GRCh38参考FASTA）
- `gencode_gtf`
- `hla`（HLA分型文件或逗号分隔列表）
- `outdir`
- 可选：`tumor_sample_name`、`normal_sample_name`、`expression`、`rna_junctions`、
  `normal_expression`、`normal_hla_ligands`

## 运行前检查

```bash
neoag-v03 run-demo --entry-mode sv_wgs --outdir /tmp/neoag_demo_svwgs --sample-id DEMO_SVWGS
```

用仓库自带SV fixture验证代码可用，不需要真实参考数据。

真实数据运行前检查参考FASTA/GTF是否存在：

```bash
test -f "<reference_fasta>" || echo "MISSING: 参考FASTA"
test -f "<gencode_gtf>" || echo "MISSING: GENCODE GTF"
```

## 分步执行路径

### 1. SV事件重建 + 滑窗产肽

```bash
neoag-v03 sv-build-raw \
  --sample-id <sample_id> \
  --profile <profile> \
  --sv-vcf <sv_vcf...> \
  --callers <callers...> \
  --reference-fasta <reference_fasta> \
  --gencode-gtf <gencode_gtf> \
  --hla <hla> \
  --outdir <outdir> \
  --tumor-sample-name <tumor_sample_name> \
  --expression <expression_tsv> \
  --rna-junctions <rna_junctions_tsv>
```

内部流程：`read_sv_inputs → 聚类 → reconstruct_cluster_protein → build_mhc1_peptides`（8-11aa滑窗）→ 展开为肽段×HLA。

check：

```bash
test -s <outdir>/parsed/raw_events.tsv
test -s <outdir>/parsed/raw_peptides.tsv
```

不合格 → 检查SV VCF是否有`SVTYPE`/`BND`等标准字段、参考FASTA/GTF版本是否匹配（都应为GRCh38）。

### 2. 结合/免疫原性预测

```bash
neoag-v03 peptide-predict \
  -i <outdir>/parsed/raw_peptides.tsv \
  -o <outdir>/presentation \
  --sample-id <sample_id>
```

### 3. 公共段

参数和check标准见 `../neoag-shared/SKILL.md`。

### 一键路径（冒烟测试用）

```bash
neoag-v03 sv-run-full \
  --sample-id <sample_id> \
  --sv-vcf <sv_vcf...> \
  --reference-fasta <reference_fasta> \
  --gencode-gtf <gencode_gtf> \
  --hla <hla> \
  --outdir <outdir> \
  --immunogenicity-stub
```

`sv-run-full` 内部同样调用分步路径里的同一批函数，不会和分步结果不一致。

## 关键输出

同其余入口：`<outdir>/scoring/ranked_peptides.v03.tsv`、`<outdir>/reports/evidence_report.v03.html`等。

## 边界声明

WGS SV重建依赖caller的断点精度和聚类合并策略，属于计算层面证据，不能替代专门的克隆结构变异分析工具，
也不是临床SV诊断，报告里要保留这条边界说明。

## 参考

- 公共段完整说明：`../neoag-shared/SKILL.md`
- WES对应版本：`../neoag-sv-wes/SKILL.md`
