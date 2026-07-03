---
name: neoag-sv-wes
description: 从外显子捕获(WES)结构变异VCF出发，用capture BED做捕获感知的SV Phase 1.5证据分层，重建融合蛋白滑窗肽段，运行结合/免疫原性预测、APPM/CCF/安全性/免疫逃逸证据构建、打分排序与报告生成。在以下情况自动调用本skill：1. 用户提供或提到SV VCF且明确是外显子捕获(WES)数据，或提到capture BED/捕获区间；2. 用户想跑WES结构变异新抗原分析。
---

# NeoAg SV-WES（外显子捕获结构变异）新抗原分析

## 使用场景

输入是WES测序得到的SV VCF，配合capture BED做"捕获感知"的Phase 1.5证据分层
（断点是否落在捕获区间内决定置信度分层和优先级封顶）。

如果没有capture BED、数据是WGS，改用 `neoag-sv-wgs`，不要用本skill。

## 必需输入

- `sample_id`
- `sv_vcf`
- `callers`
- `reference_fasta`
- `gencode_gtf`
- `hla`
- `capture_bed`（外显子捕获区间BED，本入口的强制性输入，决定WES置信度分层）
- `outdir`
- 可选：`capture_near_bp`/`capture_slop_bp`（捕获边界宽容度）、`expression`、`rna_junctions`、
  `normal_expression`、`normal_hla_ligands`、`tier1_only`（只导出WES_Tier1事件）

## 运行前检查

```bash
neoag-v03 run-demo --entry-mode sv_wes --outdir /tmp/neoag_demo_svwes --sample-id DEMO_SVWES
```

用仓库自带SV+capture fixture验证代码可用。

真实数据运行前检查capture BED是否存在（这是WES相对WGS的强约束输入，缺了这一步流程语义上就不成立）：

```bash
test -f "<capture_bed>" || echo "MISSING: capture BED（WES SV分析的必需输入）"
test -f "<reference_fasta>" || echo "MISSING: 参考FASTA"
test -f "<gencode_gtf>" || echo "MISSING: GENCODE GTF"
```

## 分步执行路径

### 1. 捕获感知SV事件重建 + 滑窗产肽

```bash
neoag-v03 sv-build-raw-wes \
  --sample-id <sample_id> \
  --profile <profile> \
  --sv-vcf <sv_vcf...> \
  --callers <callers...> \
  --reference-fasta <reference_fasta> \
  --gencode-gtf <gencode_gtf> \
  --hla <hla> \
  --outdir <outdir> \
  --capture-bed <capture_bed> \
  --capture-near-bp <capture_near_bp> \
  --capture-slop-bp <capture_slop_bp> \
  --tumor-sample-name <tumor_sample_name> \
  --expression <expression_tsv>
```

check：

```bash
test -s <outdir>/parsed/raw_events.tsv
test -s <outdir>/parsed/raw_peptides.tsv
```

不合格 → 除了通用SV检查项（VCF字段/参考版本），还要确认capture BED坐标体系（GRCh38）
与SV VCF/GTF一致，否则捕获状态判定会全部落空。

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
neoag-v03 sv-run-full-wes \
  --sample-id <sample_id> \
  --sv-vcf <sv_vcf...> \
  --reference-fasta <reference_fasta> \
  --gencode-gtf <gencode_gtf> \
  --hla <hla> \
  --capture-bed <capture_bed> \
  --outdir <outdir> \
  --immunogenicity-stub
```

## 关键输出

同其余入口，另外`raw_events.tsv`会带WES置信度分层字段（`WES_Tier1`/`WES_Tier2`/`WES_Tier3`/`WES_UNINTERPRETABLE`），
优先级封顶按profile里的`[wes_confidence_caps]`配置（默认 Tier1→B, Tier2→B_CAUTION, Tier3→C, UNINTERPRETABLE→D）。

## 边界声明

WES SV Phase 1.5 是"捕获限定"的证据层，不等价于WGS级别的SV发现完整性，报告里要明确这条边界。

## 参考

- 公共段完整说明：`../neoag-shared/SKILL.md`
- WGS对应版本：`../neoag-sv-wgs/SKILL.md`
