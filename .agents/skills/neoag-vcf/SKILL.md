---
name: neoag-vcf
description: 从体细胞VCF（SNV/InDel）出发，运行滑窗产肽、结合/免疫原性预测、APPM/CCF/安全性/免疫逃逸证据构建、打分排序与报告生成。在以下情况自动调用本skill：1. 用户提供或提到.vcf/.vcf.gz体细胞变异文件，或提到Mutect2/VEP/pVACseq；2. 用户想跑SNV/InDel新抗原分析、滑动窗口产肽、或NetMHCpan/MHCflurry打分排序。
---
# NeoAg VCF（SNV/InDel）新抗原分析

## 使用场景

输入是体细胞VCF（Mutect2等call出来的SNV/InDel），需要一路跑到候选新抗原肽段排序表和报告。

优先使用分步执行路径（可中途检查、可局部重跑）；只有做冒烟测试或用户明确要"一条命令跑完"时，
才使用 `run-full` 一键路径。

## 必需输入

执行前确认或询问：

- `sample_id`
- `variants_vcf`（体细胞VCF路径）
- `tumor_sample_name`（VCF里的肿瘤样本列名）
- `hla_alleles`（如 `HLA-A*02:01,HLA-B*07:02,HLA-C*07:02`）
- `outdir`
- `profile`（不指定则用 `default`）
- VCF是否已有VEP的 `CSQ` 注释
- 可选：`normal_expression`、`normal_hla_ligands`、`normal_proteome_fasta`（安全过滤用，）

不要臆造真实患者路径或HLA分型，缺失就直接问用户。

每一个分步执行前检查参数，有任何参数模糊或缺失，直接问用户，不要自己推测参数。

## 运行前检查

本入口只需要检查这条链路用到的工具，用带 `--entry-mode` 的demo先验证代码链路是否可用：

```bash
neoag-v03 run-demo --entry-mode snv_indel --outdir /tmp/neoag_demo_snv --sample-id DEMO_SNV
```

该命令会先打印这条链路需要的工具清单，再用仓库自带fixture跑一遍完整链路（VEP/NetMHCpan/MHCflurry
用fixture桩数据代替，不需要真实安装即可验证代码本身是否工作）。

若要跑真实VCF且VCF没有CSQ注释，额外需要本地VEP，此时检查这几个环境变量：

```bash
test -n "$NEOAG_VEP_BIN" && test -x "$NEOAG_VEP_BIN" || echo "MISSING: VEP可执行文件"
test -d "$NEOAG_VEP_CACHE" || echo "MISSING: VEP cache目录"
test -f "$NEOAG_VEP_PLUGINS/Wildtype.pm" && test -f "$NEOAG_VEP_PLUGINS/Frameshift.pm" || echo "MISSING: VEP Wildtype/Frameshift插件"
test -f "$NEOAG_REFERENCE_FASTA" || echo "MISSING: 参考FASTA"
```

若要跑本地NetMHCpan，检查：

```bash
command -v netMHCpan >/dev/null || echo "MISSING: NetMHCpan（可用 --stub 跳过，或用户已有预计算结果）"
```

环境变量或工具检查不存在时，向用户说明缺失的工具/数据，询问是否需要进行工具安装或配置。

## 分步执行路径

### 1. VEP注释（仅当VCF无CSQ时需要）

```bash
neoag-v03 vep-annotate \
  --input-vcf <variants_vcf> \
  --output-vcf <outdir>/upstream/tools/<sample_id>.vep.annotated.vcf.gz \
  --sample-id <sample_id> \
  --fasta "$NEOAG_REFERENCE_FASTA" \
  --cache-dir "$NEOAG_VEP_CACHE" \
  --plugins-dir "$NEOAG_VEP_PLUGINS" \
  --fork 4
```

check：

```bash
test -s <outdir>/upstream/tools/<sample_id>.vep.annotated.vcf.gz
```

不合格（文件为空/命令失败）→ 停止，只重跑这一步；检查 `NEOAG_VEP_BIN`/`NEOAG_VEP_CACHE`/`NEOAG_VEP_PLUGINS`/`NEOAG_REFERENCE_FASTA` 并向用户确认后再试。

### 2. 滑窗产肽

```bash
neoag-v03 extract-variant-peptides \
  --input-vcf <annotated_vcf或已有CSQ的原始vcf> \
  --output <outdir>/upstream/tools/variant_peptides.tsv \
  --sample-id <sample_id> \
  --lengths 8,9,10,11 \
  --mini-len 27 \
  --hla-alleles <hla_csv> \
  --tumor-sample-name <tumor_sample_name> \
  --normal-proteome-fasta <normal_proteome_fasta> \
  --filter-normal-proteome
```

normal_proteome-fasta为安全过滤，真实数据下必须开启，向用户确认可参考的normal proteome文件。

check：

```bash
test -s <outdir>/upstream/tools/variant_peptides.tsv
```

不合格 → 只重跑这一步，检查CSQ注释是否存在、`tumor_sample_name`是否与VCF列名一致、HLA格式是否正确。

### 3. 

### 4. 结合/免疫原性预测（公共段第一步）

```bash
neoag-v03 peptide-predict \
  -i <outdir>/upstream/tools/variant_peptides.tsv \
  -o <outdir>/presentation \
  --sample-id <sample_id> 

```

check：`presentation/` 目录下应有 `netmhcpan_evidence.tsv`/`mhcflurry_evidence.tsv` 等文件非空。

### 5. 打分排序（公共段第二步，见 `neoag-shared/SKILL.md` 完整说明）

依次调用 `appm-2`/`ccf-2`/`peptide-safety`/`immune-escape`（可并行）→ `score-v03` →
`validation-plan-v03` → `report-v03`。参数和check标准见 `../neoag-shared/SKILL.md`，
这里不重复。

### 一键路径（冒烟测试/简单重跑用）

```toml
# conf/run.<sample_id>.private.toml
[sample]
id = "<sample_id>"
profile = "default"

[tools]
stub = false
enabled = ["netmhcpan", "mhcflurry"]
immunogenicity_stub = false

[inputs]
entry_mode = "snv_indel"
variant_peptide_extraction = true
variants_vcf = "<variants_vcf>"
tumor_sample_name = "<tumor_sample_name>"
hla_alleles = ["HLA-A*02:01", "HLA-B*07:02", "HLA-C*07:02"]
normal_expression = "resources/normal_expression.example.tsv"
normal_hla_ligands = "resources/normal_hla_ligands.example.tsv"
```

```bash
neoag-v03 run-full --config conf/run.<sample_id>.private.toml --outdir results/<sample_id>
```

`run-full` 内部调用的是和分步路径完全相同的底层函数（`appm_lite.py`/`ccf_v2.py`/`peptide_safety_gate.py`/
`immune_escape.py`/`validation.py`），不会产生和分步路径不一致的结果。

## 关键输出

- `<outdir>/upstream/parsed/raw_events.tsv`
- `<outdir>/upstream/parsed/raw_peptides.tsv`
- `<outdir>/scoring/ranked_events.v03.tsv`
- `<outdir>/scoring/ranked_peptides.v03.tsv`
- `<outdir>/scoring/validation_plan.v03.tsv`
- `<outdir>/reports/evidence_report.v03.html`（`.patient.html`/`.technical.html` 双受众版本）

## 失败处理与报告

- VEP阶段失败：只重跑第1步，不动后面步骤。
- 产肽阶段失败：只重跑第2步，检查CSQ/样本名/HLA格式。
- 打分/公共段失败：见 `../neoag-shared/SKILL.md` 的失败处理表。
- 无论成功/中断，最后都向用户汇总：已完成到哪一步、产出了哪些文件、若中断则给出具体不合格原因和修复建议。

## 参考

- 公共段（presentation之后的全部步骤）：`../neoag-shared/SKILL.md`
