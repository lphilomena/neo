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
- 可选：`normal_expression`、`normal_hla_ligands`、`normal_proteome_fasta`（安全过滤用）

不要臆造真实患者路径或HLA分型，缺失就直接问用户。

每一个分步执行前检查参数，有任何参数模糊或缺失，直接问用户，不要自己推测参数。

## 运行前检查

参考 `../pipeline-get/SKILL.md`检查neoag-v03是否已有虚拟环境可用。

在所有命令前需要激活工具的path：

```bash
source conf/tools.env.sh
```

本入口只需要检查这条链路用到的工具，用带 `--entry-mode` 的demo先验证流程链路是否可用：

```bash
neoag-v03 run-demo --entry-mode snv_indel --outdir /tmp/neoag_demo_snv --sample-id DEMO_SNV
```

该命令会先打印这条链路需要的工具清单，并按真实需要的工具跑一个demo。

若用户提供的VCF没有CSQ注释，额外需要本地VEP, 注意版本为vep105，此时检查这几个环境变量：

```bash
test -n "$NEOAG_VEP_BIN" && test -x "$NEOAG_VEP_BIN" || echo "MISSING: VEP可执行文件"
test -d "$NEOAG_VEP_CACHE" || echo "MISSING: VEP cache目录"
test -f "$NEOAG_VEP_PLUGINS/Wildtype.pm" && test -f "$NEOAG_VEP_PLUGINS/Frameshift.pm" || echo "MISSING: VEP Wildtype/Frameshift插件"
test -f "$NEOAG_REFERENCE_FASTA" || echo "MISSING: 参考FASTA"
```

环境变量不存在或错误时，让用户提供准确的环境变量，并将用户的环境变量写入 `conf/tools.env.local.sh`复用。

结合/免疫原性预测需要netmhcpan, mhcflurry, bigmhc-im, prime，netmhcstabpan和deepimmuno，其中netmhcpan, mhcflurry, bigmhc-im, prime必选，netmhcstabpan和deepimmuno询问用户。

确定需要的工具后，检查工具是否已安装可用。

工具检查不存在时，向用户说明缺失的工具，询问是否需要进行工具安装或配置。

工具安装严格参考 `../pipeline-get/reference/INSTALL_AND_DATA.md`，不要自行写命令安装，没有安装命令的工具警告用户。

## 分步执行路径

### 1. VEP注释（仅当VCF无CSQ时需要）

```bash
neoag-v03 vep-annotate \
  --input-vcf <variants_vcf> \
  --output-vcf <outdir>/upstream/tools/<sample_id>.vep.annotated.vcf.gz \
  --fasta $NEOAG_REFERENCE_FASTA \
  --cache-dir $NEOAG_VEP_CACHE \
  --plugins-dir $NEOAG_VEP_PLUGINS \
  --fork 4
```

check：

```bash
test -s <outdir>/upstream/tools/<sample_id>.vep.annotated.vcf.gz
```

不合格（文件为空/命令失败）→ 停止，只重跑这一步；检查 `NEOAG_VEP_BIN`/`NEOAG_VEP_CACHE`/`NEOAG_VEP_PLUGINS`/`NEOAG_REFERENCE_FASTA` 并向用户确认后再试。

### 2. 滑窗产肽 + 标准化

```bash
neoag-v03 snv-build-raw \
  --variants-vcf <annotated_vcf或已有CSQ的原始vcf> \
  --outdir <outdir>/upstream \
  --sample-id <sample_id> \
  --hla <hla_allele_1> <hla_allele_2> ... \
  --tumor-sample-name <tumor_sample_name> \
  --lengths 8,9,10,11 \
  --mini-len 27 \
  --normal-proteome-fasta <normal_proteome_fasta> \
  --filter-normal-proteome
```

`--hla` 是必需参数。normal_proteome_fasta为安全过滤，真实数据下必须开启，向用户确认可参考的normal proteome文件。

**必须用这条命令**，`snv-build-raw`内部会先做滑窗产肽，再自动转换成 `peptide-predict`能消费的标准 `raw_peptides.tsv`（含 `peptide`列）。

check：

```bash
test -s <outdir>/upstream/tools/variant_peptides.tsv
test -s <outdir>/upstream/parsed/raw_peptides.tsv
```

不合格 → 只重跑这一步，检查CSQ注释是否存在、`tumor_sample_name`是否与VCF列名一致、HLA格式是否正确。

### 3. 结合/免疫原性预测（公共段第一步）

```bash
neoag-v03 peptide-predict \
  -i <outdir>/upstream/parsed/raw_peptides.tsv \
  -o <outdir>/presentation \
  --sample-id <sample_id>
```

对应工具中netmhcpan, mhcflurry, bigmhc-im, prime必须运行，不可skip。

| 工具          | 输入文件         | 输出文件                                                   |
| ------------- | ---------------- | ---------------------------------------------------------- |
| netmhcpan     | raw_peptides.tsv | tools/netmhcpan.xls → presentation/netmhcpan_evidence.tsv |
| mhcflurry     | raw_peptides.tsv | tools/mhcflurry.csv → mhcflurry_evidence.tsv              |
| bigmhc-im     | raw_peptides.tsv | bigmhc_im_evidence.tsv                                     |
| prime         | raw_peptides.tsv | presentation/prime_evidence.tsv                            |
| netmhcstabpan | raw_peptides.tsv | tools/netmhcstabpan.tsv                                    |
| deepimmuno    | raw_peptides.tsv | deepimmuno_evidence.tsv                                    |

check：各输出文件非空。

### 4. 打分排序（公共段第二步，见 `neoag-shared/SKILL.md` 完整说明）

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
