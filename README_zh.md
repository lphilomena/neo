# NeoAg Event Pipeline v0.4.3 Online Release

NeoAg Event Pipeline 是研究型肿瘤新抗原候选优先级分析流程。它将 SNV/InDel、fusion、splice、结构变异以及 peptide-only 候选统一转换为标准事件表和 peptide-HLA 表，并叠加 HLA 呈递、APPM、CCF、safety、immune escape、validation plan 和 report 证据层。

本包是轻量级 online release，包含源代码、CLI 入口、Nextflow workflow、测试、fixture、profile、安装脚本和文档。不包含大型参考库、授权工具、conda 环境、缓存 work 目录、真实患者数据或生产结果。

重要边界：本流程输出的是 computational triage 和实验验证规划结果，不构成临床诊断、临床耐药判定或已验证治疗建议。

## 功能概览

本流程可以：

- 将 pVACtools-like SNV/fusion/splice 输出解析为 `raw_events.tsv` 和 `raw_peptides.tsv`。
- 从 VEP 注释 VCF 生成 sliding-window variant peptides；当 VCF 缺少 CSQ 注释时，可在配置完整时自动先运行 VEP 注释。
- 使用 NetMHCpan、MHCflurry 以及可选稳定性/免疫原性工具构建 MHC 呈递证据。
- 构建 APPM 2.0 证据，包括输入完整性、冲突、肽段修饰因子和 immune-context 注释。
- 基于 purity、CNV 和 VAF context 估计 CCF/clonality。
- 基于 normal expression、normal ligandome、normal junction、matched-normal 和 reference proteome context 构建 peptide safety 证据。
- 基于 HLA LOH、APPM、CCF、B2M/JAK/APM context 及相关证据表构建 immune escape 证据。
- 为 frameshift、splice、exon-junction、fusion 和 SV 候选生成 long peptide / minigene 验证设计。
- 生成患者沟通版和科研技术版 HTML 报告。
- 通过 CLI 或内置 Nextflow wrapper 运行 fixture workflow。

`.v03.tsv` 后缀是 schema 兼容标签，不代表软件版本。当前 release 版本是 v0.4.3，但仍写出 v03-compatible 表格，便于旧下游脚本继续读取相同文件名。

## 快速开始

从项目根目录运行：

```bash
python -m pip install -e .
neoag-v03 run-demo --outdir work/demo_v043 --sample-id DEMO001
```

重要 demo 输出包括：

- `work/demo_v043/scoring/ranked_peptides.v03.tsv`
- `work/demo_v043/scoring/ranked_events.v03.tsv`
- `work/demo_v043/scoring/validation_plan.v03.tsv`
- `work/demo_v043/reports/evidence_report.v03.html`
- `work/demo_v043/reports/evidence_report.patient.html`
- `work/demo_v043/reports/evidence_report.technical.html`
- `work/demo_v043/appm/appm_summary.tsv`
- `work/demo_v043/appm/appm_peptide_modifiers.tsv`
- `work/demo_v043/clonality/ccf_lite.tsv`
- `work/demo_v043/safety/peptide_safety.tsv`
- `work/demo_v043/immune_escape/peptide_escape_flags.tsv`

运行测试：

```bash
python -m pip install -e '.[test]'
pytest -q
```

默认测试命令会刻意跳过 integration、benchmark 和 external-tool 测试。

## 常用运行命令

### 准备环境

仅用于 fixture 开发：

```bash
python -m pip install -e '.[test]'
pytest -q
neoag-v03 run-demo --outdir work/demo_v043 --sample-id DEMO001
```

需要外部工具的运行：

```bash
bash scripts/setup_tools_env.sh
source conf/tools.env.sh
python -m pip install -e '.[test]'
neoag-v03 check-tools
```

较小的开发/测试环境：

```bash
NEOAG_TOOLS_LITE=1 bash scripts/setup_tools_env.sh
source conf/tools.env.sh
python -m pip install -e '.[test]'
pytest -q
```

### 从已有 pVAC-like 表运行

当你已有 pVACseq/pVACfuse/pVACsplice-like aggregated tables 时使用：

```bash
neoag-v03 run-v03 \
  --outdir results/sample \
  --sample-id SAMPLE001 \
  --profile default \
  --pvac data/fixtures/pvacseq_aggregated.tsv \
  --immunogenicity-stub
```

### 从预生成 raw intermediates 运行

当 `parsed/raw_events.tsv` 和 `parsed/raw_peptides.tsv` 已经存在时使用：

```bash
neoag-v03 run-v03 \
  --outdir results/sample \
  --sample-id SAMPLE001 \
  --profile default \
  --raw-events results/sample/parsed/raw_events.tsv \
  --raw-peptides results/sample/parsed/raw_peptides.tsv \
  --netmhcpan results/sample/presentation/netmhcpan.xls \
  --mhcflurry results/sample/presentation/mhcflurry.csv \
  --expression results/sample/parsed/expression.tsv \
  --hla-loh results/sample/tools/hla_loh.tsv \
  --purity results/sample/tools/purity.tsv \
  --cnv results/sample/tools/cnv_segments.tsv
```

### Sliding-window Variant Peptides 到排序结果

当你有 somatic SNV/InDel VCF，并希望通过 sliding window 生成 mutant peptides、预测 peptide-HLA 呈递并生成 event/peptide 排序结果时使用。

如果 VCF 已含 VEP `CSQ` 注释，pipeline 会直接使用。如果缺少 `CSQ`，在 VEP、cache、reference FASTA 和 plugins 配置完整时，`run-full` 会先自动运行 VEP 注释。

```bash
cat > conf/run.sliding.private.toml <<'TOML'
[sample]
id = "SAMPLE001"
profile = "default"

[tools]
stub = false
enabled = ["netmhcpan", "mhcflurry"]
immunogenicity_stub = false

[inputs]
entry_mode = "snv_indel"
variant_peptide_extraction = true
variants_vcf = "/path/to/sample.somatic.pass.vcf.gz"
tumor_sample_name = "TUMOR"
hla_alleles = ["HLA-A*02:01", "HLA-B*07:02", "HLA-C*07:02"]
extract_appm_from_vcf = false
normal_expression = "resources/normal_expression.example.tsv"
normal_hla_ligands = "resources/normal_hla_ligands.example.tsv"
TOML

neoag-v03 run-full \
  --config conf/run.sliding.private.toml \
  --outdir results/SAMPLE001_sliding
```

关键输出：

- `results/SAMPLE001_sliding/upstream/tools/variant_peptides.tsv`
- `results/SAMPLE001_sliding/upstream/tools/variant_peptides.annotated.tsv`
- `results/SAMPLE001_sliding/upstream/parsed/raw_events.tsv`
- `results/SAMPLE001_sliding/upstream/parsed/raw_peptides.tsv`
- `results/SAMPLE001_sliding/scoring/ranked_events.v03.tsv`
- `results/SAMPLE001_sliding/scoring/ranked_peptides.v03.tsv`
- `results/SAMPLE001_sliding/scoring/validation_plan.v03.tsv`
- `results/SAMPLE001_sliding/reports/evidence_report.v03.html`
- `results/SAMPLE001_sliding/reports/evidence_report.patient.html`
- `results/SAMPLE001_sliding/reports/evidence_report.technical.html`

手动调试 variant peptide extraction：

```bash
neoag-v03 extract-variant-peptides \
  --input-vcf /path/to/sample.vep.annotated.vcf.gz \
  --output results/SAMPLE001_sliding/upstream/tools/variant_peptides.tsv \
  --sample-id SAMPLE001 \
  --lengths 8,9,10,11 \
  --mini-len 27 \
  --hla-alleles HLA-A*02:01,HLA-B*07:02,HLA-C*07:02 \
  --tumor-sample-name TUMOR \
  --normal-proteome-fasta /path/to/Homo_sapiens.GRCh38.pep.all.fa \
  --filter-normal-proteome
```

没有授权 predictor 的 smoke test 可以在 TOML 中设置 `stub = true`，或在直接调用 `run-v03` 时添加 `--immunogenicity-stub`。生产排序应使用真实 NetMHCpan/MHCflurry 输出以及真实 normal-expression/normal-ligand evidence，不应使用 fixture resources。

### 构建标准证据 sidecars

```bash
neoag-v03 build-evidence-layer \
  --outdir results/sample \
  --profile default \
  --sample-id SAMPLE001 \
  --raw-events results/sample/parsed/raw_events.tsv \
  --raw-peptides results/sample/parsed/raw_peptides.tsv \
  --expression results/sample/parsed/gene_expression.tsv \
  --rna-vaf results/sample/parsed/rna_vaf.tsv \
  --rna-junction results/sample/parsed/rna_junctions.tsv \
  --fusion-evidence results/sample/parsed/fusion_evidence.tsv \
  --normal-expression resources/normal_expression.example.tsv \
  --normal-hla-ligands resources/normal_hla_ligands.example.tsv
```

### HLA LOH 转换与交叉检查

```bash
neoag-v03 convert-lohhla \
  -i results/sample/tools/LOHHLA.HLAlossPrediction_CI.xls \
  -o results/sample/tools/lohhla.hla_loh.tsv

neoag-v03 convert-spechla \
  -i results/sample/tools/merge.hla.copy.txt \
  -o results/sample/tools/spechla.hla_loh.tsv

neoag-v03 crosscheck-hla-loh \
  --lohhla-hla-loh results/sample/tools/lohhla.hla_loh.tsv \
  --spechla-hla-loh results/sample/tools/spechla.hla_loh.tsv \
  --out results/sample/tools/hla_loh.crosscheck.tsv \
  --consensus-out results/sample/tools/hla_loh.consensus.tsv
```

### 生成报告

生成默认综合报告，同时输出患者沟通版和科研技术版报告：

```bash
neoag-v03 report-v03 \
  --profile default \
  --ranked-events results/sample/scoring/ranked_events.v03.tsv \
  --ranked-peptides results/sample/scoring/ranked_peptides.v03.tsv \
  --appm-summary results/sample/appm/appm_summary.tsv \
  --validation-plan results/sample/scoring/validation_plan.v03.tsv \
  --outdir results/sample \
  --audience both \
  --out results/sample/reports/evidence_report.v03.html
```

### Nextflow Fixture 运行

建议使用项目 wrapper，不要直接调用 `nextflow`。该 wrapper 会优先使用当前 checkout 的 `bin/neoag-v03`，设置项目路径，并避免将 Nextflow metadata 写入 root-owned 位置。

```bash
export NXF_HOME=/path/to/writable/nextflow_cache
bin/neoag-nextflow -version
bin/neoag-nextflow run workflows/main.nf \
  -w /tmp/neoag_nf_work \
  --pvac_files data/fixtures/pvacseq_aggregated.tsv \
  --outdir results/demo_nf \
  --sample_id NF_DEMO
```

查看某个命令的完整参数：

```bash
neoag-v03 <command> --help
```

## 配置文件

真实部署路径应保存在 local/private 文件中。建议先复制模板，再编辑站点相关路径。

| 文件 | 用途 | 真实数据需要编辑？ | 是否提交/打包？ |
| --- | --- | --- | --- |
| `conf/tools.env.sh` | 主环境入口。设置项目路径、conda env 名称、工具根目录、VEP cache fallback 和 wrapper `PATH`。 | 通常不直接改；用 local override 覆盖。 | 是 |
| `conf/tools.env.local.example.sh` | 私有站点路径模板，例如患者数据根目录、共享参考库、授权工具安装路径和 cache 目录。 | 复制为 `conf/tools.env.local.sh` 后编辑。 | 示例文件提交；复制后的 local 文件不提交 |
| `conf/site.config.example` | 站点/集群/Nextflow executor 模板。 | 复制为 `conf/site.config` 后编辑。 | 示例文件提交；复制后的 local 文件不提交 |
| `conf/run.private.example.toml` | 真实样本私有运行配置模板。 | 复制为 private TOML 后编辑。 | 示例文件提交；复制后的 local 文件不提交 |
| `conf/run.snv_wes.example.toml` | WES SNV workflow 示例，包含 Mutect2/annotation 输入。 | 使用真实 BAM/VCF 前复制并编辑。 | 示例文件提交 |
| `conf/run.stub.toml` | 轻量 stub/demo upstream 配置。 | 生产不使用；适合 smoke test。 | 是 |
| `conf/*.example.toml` | SV、fusion、splice、peptide-only 或 site mode 的 workflow-specific 示例。 | 复制并编辑。 | 示例文件提交 |

典型设置流程：

```bash
cp conf/tools.env.local.example.sh conf/tools.env.local.sh
cp conf/run.private.example.toml conf/run.sample.private.toml
# 编辑两个文件，填入站点路径和样本输入。
source conf/tools.env.sh
neoag-v03 check-tools
neoag-v03 run-full --config conf/run.sample.private.toml --outdir results/sample
```

`conf/tools.env.sh` 或 local override 中常见的重要变量：

| 变量 | 含义 |
| --- | --- |
| `NEOAG_PROJECT_ROOT` | 项目 checkout 根目录。 |
| `NEOAG_TOOLS_ROOT` | 外部工具、wrapper 和本地 artifact bundle 根目录。 |
| `NEOAG_CONDA_BASE` | Miniforge/Mambaforge 安装路径。 |
| `NEOAG_CONDA_ENV` | 主 Python CLI 环境。 |
| `NEOAG_VEP_ENV` | VEP conda 环境名称/路径。 |
| `NEOAG_VEP_BIN` | VEP 可执行文件或 wrapper 路径。 |
| `NEOAG_VEP_CACHE` | VEP offline cache 根目录，预期包含 `homo_sapiens/105_GRCh38` 或等价目录。 |
| `NEOAG_VEP_PLUGINS` | VEP plugin 目录，包含 `Wildtype.pm` 和 `Frameshift.pm`。 |
| `NEOAG_REFERENCE_FASTA` | VEP/GATK/SV peptide workflow 使用的 GRCh38 FASTA。 |
| `NEOAG_NORMAL_PROTEOME_FASTA` | peptide safety filtering 使用的 normal/reference proteome FASTA。 |
| `NETMHCPAN_HOME` / `NEOAG_NETMHCPAN_BIN` | NetMHCpan 安装和可执行文件路径。 |
| `NEOAG_NETMHCPAN_TMPDIR` | NetMHCpan 使用的短临时目录。 |
| `NETMHCSTABPAN_HOME` | NetMHCstabpan 安装路径。 |
| `LOHHLA_HOME`, `POLYSOLVER_HOME`, `NOVOALIGN_LICENSE_FILE` | LOHHLA 及依赖路径/license。 |
| `FACETS_HOME`, `NEOAG_DBSNP_VCF` | FACETS 脚本和 `snp-pileup` 用 dbSNP/common SNP VCF。 |
| `NEOAG_ASCAT_ENV`, `ASCAT_HOME` | ASCAT conda env 和 wrapper 路径。 |
| `NXF_HOME` | 可写 Nextflow cache；干净 online/offline workflow 运行需要设置。 |

不要提交或打包以下 local/private 文件：

- `conf/tools.env.local.sh`
- `conf/site.config`
- `conf/private/*`
- `conf/*.private.toml`
- 包含患者标识、临床数据绝对路径、集群凭据或授权工具路径的文件

## 安装、工具与数据

### 基础系统环境

推荐基线：

| 组件 | 推荐 | 说明 |
| --- | --- | --- |
| OS | Linux x86_64 | 内置 wrapper 和多数 upstream bioinformatics 工具假设 Linux。 |
| Shell | Bash | 脚本使用 Bash。 |
| CPU/RAM | 真实数据建议 8+ CPU、32+ GB RAM | Fixture demo 很小；患者级 WES/WGS/fusion 需要更多资源。 |
| Disk | 真实数据部署建议 200 GB+ | VEP cache、hg38 reference、fusion reference 和 Nextflow work dir 较大。 |
| Network | 首次安装/下载需要 | Offline 安装应预先准备 tarball、conda package、reference 和 Nextflow cache。 |
| Java | Java 11 或更新 | Nextflow 和部分工具需要。 |
| Conda/Mamba | 推荐 Miniforge/Mambaforge | 工具环境定义在 `conda/` 下。 |

Ubuntu/Debian 建议安装：

```bash
sudo apt-get update
sudo apt-get install -y \
  bash coreutils curl wget git tar gzip unzip bzip2 xz-utils \
  ca-certificates build-essential openjdk-17-jre-headless rsync file
```

如果迁移后的压缩包丢失执行权限：

```bash
find bin -maxdepth 1 -type f -exec chmod +x {} \;
find scripts -maxdepth 1 -type f -name '*.sh' -exec chmod +x {} \;
```

### 外部工具安装表

Fixture demo 不需要多数外部工具；真实数据模式按 workflow 需要安装。只安装所选 workflow 需要的工具即可。

| 工具 | 用途 | 是否必需 | 安装/下载命令 | 关键变量 | 验证 |
| --- | --- | --- | --- | --- | --- |
| pVACtools (`pvacseq`, `pvacfuse`, `pvacsplice`) | upstream SNV/fusion/splice 候选生成 | 使用 pVAC upstream 模式时需要 | `bash scripts/setup_tools_env.sh` | `NEOAG_PVAC_DOCKER`, `NEOAG_PVAC_WORKDIR` | `neoag-v03 check-tools` |
| NetMHCpan 4.2 | Binding/presentation prediction | 本地真实 NetMHCpan 运行需要，除非使用 fallback/stub | `bash scripts/install_netmhcpan.sh /path/to/netMHCpan-4.2*.tar.gz` | `NETMHCPAN_HOME`, `NETMHCpan`, `NEOAG_NETMHCPAN_BIN`, `NEOAG_NETMHCPAN_BACKEND` | `neoag-v03 check-tools` |
| MHCflurry | Binding/presentation prediction | NetMHCpan 的可选替代/补充 | `bash scripts/setup_tools_env.sh`；需要时再运行 `mhcflurry-downloads fetch` | `NEOAG_CONDA_ENV`, `NEOAG_FORCE_CPU` | `neoag-v03 check-tools` |
| NetMHCstabpan | pMHC stability evidence | 可选 | `bash scripts/install_netmhcstabpan.sh --iedb` 或授权 tarball 安装 | `NETMHCSTABPAN_HOME` | `neoag-v03 check-tools` |
| PRIME / MixMHCpred / BigMHC | Immunogenicity evidence | 可选 | `bash scripts/install_immunogenicity_tools.sh` | `PRIME_HOME`, `MIXMHCPRED_HOME`, `BIGMHC_DIR`, `NEOAG_PRIME_JOBS` | `neoag-v03 check-tools` |
| DeepImmuno | 可选 immunogenicity evidence | 可选 | `bash scripts/install_deepimmuno.sh` | `DEEPIMMUNO_DIR` | `neoag-v03 check-tools` |
| VEP | VCF annotation 和 peptide extraction | `vep-annotate` / `run-full` 自动注释需要 | `bash scripts/install_vep.sh`；cache 用 `bash scripts/install_vep_cache.sh` | `NEOAG_VEP_ENV`, `NEOAG_VEP_BIN`, `NEOAG_VEP_CACHE`, `NEOAG_VEP_PLUGINS`, `NEOAG_VEP_ONLINE` | `neoag-v03 check-tools` |
| GATK4 / Mutect2 | WES/WGS SNV calling | 从 BAM 开始时需要 | `bash scripts/install_gatk.sh` | `NEOAG_GATK_ENV` | `neoag-v03 check-tools` |
| LOHHLA | HLA LOH evidence | 可选；immune-escape evidence 推荐 | `bash scripts/install_lohhla.sh`；Polysolver/Novoalign 需另行配置 | `LOHHLA_HOME`, `POLYSOLVER_HOME`, `NOVOALIGN_LICENSE_FILE` | `neoag-v03 check-tools` |
| SpecHLA | HLA copy/LOH conversion | 可选 | 外部安装，并将输出提供给 `convert-spechla` | 站点自定义 | `neoag-v03 convert-spechla --help` |
| FACETS | Purity/CNV/LOH evidence | 可选；CCF/escape 推荐 | `bash scripts/install_facets.sh` | `FACETS_HOME`, `NEOAG_DBSNP_VCF` | `neoag-v03 check-tools` |
| ASCAT | CNV/LOH evidence | 可选 | `bash scripts/install_ascat_pyclone.sh` | `NEOAG_ASCAT_ENV`, `ASCAT_HOME` | `neoag-v03 check-tools` |
| PyClone-VI | Clonality context | 可选 | `bash scripts/install_ascat_pyclone.sh` | `NEOAG_PYCLONE_ENV`, `NEOAG_PYCLONE_BIN` | `neoag-v03 check-tools` |
| STAR-Fusion / FusionCatcher / Arriba / EasyFuse | Fusion discovery | 可选；对应 fusion workflow 需要 | 外部安装/挂载；仅当已有 Nextflow conda cache 时 seed EasyFuse env | `NEOAG_FUSION_ENV`, `NEOAG_STAR_FUSION_HOME`, `NEOAG_CTAT_LIB_DIR`, `NEOAG_EASYFUSE_HOME`, `NEOAG_EASYFUSE_REF` | `neoag-v03 check-tools` |
| Manta / GRIDSS / SvABA / Sniffles2 | SV discovery | 可选 upstream SV caller | 外部安装或站点 conda/module 提供 | `NEOAG_SV_ENV`, `NEOAG_MANTA_ENV` | `neoag-v03 check-tools` |
| PURPLE / AMBER / COBALT | Purity、ploidy、CNV、LOH evidence | 可选 | 见 `docs/TOOLS_SETUP.md` 和本地 wrapper | `HMFTOOLS_HOME`、站点参考库 | 工具 wrapper `--help` |
| DASH | HLA LOH / allele-specific deletion evidence | 可选 | 见 `docs/TOOLS_SETUP.md`；模型可能需另行提供 | DASH env/model path | 工具 wrapper |

NetMHCpan、NetMHCstabpan、LOHHLA、Novoalign/Polysolver 等授权工具可能需要学术或机构许可。不要将其二进制文件重新分发到 online release 中。

### NetMHCpan 说明

从 DTU 授权 Linux tarball 安装：

```bash
mkdir -p vendor
cp /path/to/netMHCpan-4.2c.Linux.tar.gz vendor/
export NEOAG_CONDA_BASE="$(conda info --base)"
bash scripts/install_netmhcpan.sh vendor/netMHCpan-4.2c.Linux.tar.gz
source conf/tools.env.sh
neoag-v03 check-tools
netMHCpan -h | head
```

如果 NetMHCpan 已安装但 wrapper 损坏：

```bash
bash scripts/install_netmhcpan.sh --repair
```

在较老 host glibc 上运行 patched NetMHCpan binary 时，`neoag-tools` 环境必须保留 `sysroot_linux-64` 和 `patchelf`。lite 环境文件已经包含这些依赖。

### VEP 说明

安装 VEP 并配置 cache/plugins：

```bash
bash scripts/install_vep.sh
source conf/tools.env.sh
neoag-v03 check-tools
```

安装或指定 VEP cache：

```bash
# 在线安装 cache，可能较慢且占用空间大。
bash scripts/install_vep_cache.sh

# 或使用已有 cache。
export NEOAG_VEP_CACHE=/path/to/vep_cache
export NEOAG_VEP_CACHE_VERSION=105
source conf/tools.env.sh
```

`NEOAG_VEP_PLUGINS` 应指向包含 `Wildtype.pm` 和 `Frameshift.pm` 的目录。普通 VEP 可以不使用这些 plugins，但本项目会用它们生成 pVACseq-compatible WT/frameshift 信息，并让 peptide extraction 更完整。

### 参考数据表

大型数据应放在 `NEOAG_TOOLS_ROOT`、`NEOAG_SHARED_REF_DIR` 或其他站点管理的 reference 区域，而不是源代码 checkout 内。

| 数据/参考库 | 用途 | 下载/设置命令 | 期望变量/路径 | 验证 |
| --- | --- | --- | --- | --- |
| VEP cache, GRCh38 release 105 | Offline VEP annotation | `bash scripts/install_vep_cache.sh` | `NEOAG_VEP_CACHE=/path/to/data/vep`，包含 `homo_sapiens/105_GRCh38` 或等价目录 | `test -d "$NEOAG_VEP_CACHE/homo_sapiens"` |
| VEP plugins | WT 和 frameshift plugin annotation | 由 VEP/pVAC tooling 安装，或从站点 bundle 复制 | `NEOAG_VEP_PLUGINS=/path/to/work/vep_plugins` | `test -f "$NEOAG_VEP_PLUGINS/Wildtype.pm"` |
| GRCh38 FASTA 和索引 | VEP peptide extraction、GATK、SV peptide building | `bash scripts/download_ref_hg38.sh /path/to/ref/hg38` | `NEOAG_REFERENCE_FASTA=/path/to/Homo_sapiens_assembly38.fasta` | `test -f "$NEOAG_REFERENCE_FASTA"` |
| dbSNP/common SNP VCF | FACETS `snp-pileup` 和部分 CNV workflow | 站点 reference bundle 或 hg38 bundle 下载 | `NEOAG_DBSNP_VCF=/path/to/dbsnp_chr.vcf.gz` | `test -f "$NEOAG_DBSNP_VCF"` |
| gnomAD AF VCF 和 PoN | GATK Mutect2 filtering | `bash scripts/download_ref_hg38.sh /path/to/ref/hg38` 或站点 bundle | 所选 run config 中的路径 | `test -f /path/to/af-only-gnomad.hg38.vcf.gz` |
| Ensembl protein FASTA | Peptide safety normal/reference proteome screen | 手动下载 Ensembl GRCh38 peptide FASTA 或使用站点 bundle | `NEOAG_NORMAL_PROTEOME_FASTA=/path/to/Homo_sapiens.GRCh38.pep.all.fa` | `test -f "$NEOAG_NORMAL_PROTEOME_FASTA"` |
| Normal expression table | Peptide safety evidence | 站点生成 TSV 或 fixture example | CLI 参数或 run config 路径 | 检查 TSV header |
| Normal HLA ligand table | Peptide safety evidence | 站点生成 TSV 或 fixture example | CLI 参数或 run config 路径 | 检查 TSV header |
| CTAT genome lib | STAR-Fusion | 按 STAR-Fusion/CTAT 文档下载，或挂载站点 bundle | `CTAT_GENOME_LIB`, `NEOAG_CTAT_LIB_DIR`, `NEOAG_SHARED_REF_DIR` | `test -d "$CTAT_GENOME_LIB"` |
| EasyFuse reference | EasyFuse workflow | 按 EasyFuse 文档下载，或挂载站点 bundle | `NEOAG_EASYFUSE_REF`, `NEOAG_SHARED_REF_DIR` | `test -d "$NEOAG_EASYFUSE_REF"` |
| GTF annotation | SV/fusion peptide generation | 使用与 reference FASTA 匹配的 GENCODE/Ensembl GTF | CLI 参数 `--gencode-gtf` | `test -f /path/to/genes.gtf` |
| Capture BED | WES SV Phase 1.5 | 使用 panel/exome capture BED | CLI 参数 `--capture-bed` | `test -f /path/to/capture.bed` |
| HLA allele file | Peptide prediction 和 SV workflow | 站点 HLA typing 输出，转换为每行一个 allele | CLI 参数 `--hla` | `head /path/to/hla.txt` |

推荐外部 bundle 布局：

```text
/path/to/neoag_artifact_bundle/
  tools/
    netMHCpan/
    netMHCstabpan/
    DeepImmuno/
    prime/
    mixMHCpred_install/
  data/
    ref/hg38/
    vep/
    ref/ctat/
  work/
    vep_plugins/
```

然后配置：

```bash
export NEOAG_TOOLS_ROOT=/path/to/neoag_artifact_bundle
export NEOAG_SHARED_REF_DIR=/path/to/shared_refs
source conf/tools.env.sh
neoag-v03 check-tools
```

本地路径清单详见 `docs/PROJECT_DATA_PATHS.md`。

## Workflow 依赖矩阵

| Workflow / command | 最小输入 | 工具 | 参考库/数据 |
| --- | --- | --- | --- |
| Fixture demo: `neoag-v03 run-demo --outdir work/demo_v043 --sample-id DEMO001` | 内置 fixture | Python package 之外不需要 | 内置 fixtures/resources |
| Parsed pVAC results: `neoag-v03 run-v03 --outdir results/sample --sample-id SAMPLE001 --pvac data/fixtures/pvacseq_aggregated.tsv --immunogenicity-stub` | pVAC-like TSVs | 若输入已存在则不需外部工具 | 可选 normal expression/ligand tables |
| Raw intermediates: `neoag-v03 run-v03 --raw-events ... --raw-peptides ...` | `raw_events.tsv`, `raw_peptides.tsv` | 如提供则使用 NetMHCpan/MHCflurry 输出；可选 evidence tools | 可选 expression、LOH、purity、CNV、normal evidence |
| Full upstream run: `neoag-v03 run-full --config conf/run.sample.private.toml --outdir results/sample` | Run config | 取决于 enabled tools | 取决于 enabled tools |
| Binding prediction only: `peptide-predict` | Peptide/HLA table | 按选择使用 NetMHCpan、MHCflurry、PRIME/BigMHC/DeepImmuno | HLA alleles；predictor model data |
| VEP annotation: `vep-annotate` | VCF | VEP | VEP cache、reference FASTA、plugins |
| Variant peptide extraction: `extract-variant-peptides` | VEP-annotated VCF | Python；可选 VEP pre-step | Reference FASTA、可选 normal proteome |
| WES SNV calling: `snv-call-wes` | Tumor/normal BAM | GATK4 | GRCh38 FASTA、gnomAD AF VCF、PoN、intervals |
| WES SNV full: `snv-run-full-wes` | Somatic VCF 或 BAMs | BAM 模式需要 GATK；如启用则需 pVAC/binding tools | GRCh38 FASTA、HLA、可选 normal evidence |
| SV WGS raw build: `sv-build-raw` | SV VCF、FASTA、GTF、HLA | Python | Reference FASTA、GTF、HLA file |
| SV WES raw build: `sv-build-raw-wes` | SV VCF、FASTA、GTF、HLA、capture BED | Python | Reference FASTA、GTF、capture BED、HLA file |
| SV score: `sv-score-v03` | Raw events/peptides | 除非 `--binding-stub`，否则需要 NetMHCpan/MHCflurry | HLA alleles、可选 evidence tables |
| Long-read SV wrapper | FASTQ/BAM 或 Sniffles2 VCF | 按选择使用 minimap2/samtools/Sniffles2 | Reference FASTA、GTF、HLA |
| Fusion discovery | FASTQ/BAM 或 caller outputs | 按选择使用 STAR-Fusion、FusionCatcher、Arriba、EasyFuse | CTAT/EasyFuse/fusion caller references |
| Immune escape evidence: `immune-escape` | Raw peptides、APPM/CCF/LOH evidence | 可选 LOHHLA/FACETS upstream | HLA LOH、CNV、VEP/APM/JAK/B2M evidence |
| Nextflow fixture | 内置 pVAC fixture | Java/Nextflow runtime | 内置 fixtures；可写 `NXF_HOME` |

## 测试

默认 pytest 只运行快速 unit tests：

```bash
pytest -q
```

显式运行更大的测试分组：

```bash
pytest -q --run-integration
pytest -q --run-benchmark
pytest -q --run-external
pytest -q --run-all
```

也支持 marker 形式：

```bash
pytest -q -m unit
pytest -q -m integration --run-integration
pytest -q -m benchmark --run-benchmark
pytest -q -m external --run-external
```

这样可以避免 lightweight release 用户在普通 `pytest` 下误运行耗时 Nextflow、benchmark 或 external-tool 测试。

## 安装验收命令

安装后从项目根目录运行。

### 基础包验收

```bash
source conf/tools.env.sh
python -m pip install -e '.[test]'
pytest -q
neoag-v03 run-demo --outdir work/demo_v043 --sample-id DEMO001
```

### 工具可见性验收

```bash
source conf/tools.env.sh
neoag-v03 check-tools
bash scripts/check_tools_env.sh
```

如果所选 workflow 不需要某些工具，`check-tools` 可能会报告可选工具 missing。生产运行时，所选 workflow 需要的每个工具都应为 `OK`。

### Nextflow 验收

```bash
export NXF_HOME=/path/to/writable/nextflow_cache
bin/neoag-nextflow -version
bin/neoag-nextflow run workflows/main.nf \
  -w /tmp/neoag_nf_work \
  --pvac_files data/fixtures/pvacseq_aggregated.tsv \
  --outdir results/demo_nf \
  --sample_id NF_DEMO
```

预期输出包括：

- `results/demo_nf/scoring/ranked_peptides.v03.tsv`
- `results/demo_nf/scoring/ranked_events.v03.tsv`
- `results/demo_nf/reports/evidence_report.v03.html`
- `results/demo_nf/provenance/workflow_provenance.yml`

### Reference 文件验收

```bash
test -f "$NEOAG_REFERENCE_FASTA"
test -d "$NEOAG_VEP_CACHE/homo_sapiens"
test -f "$NEOAG_VEP_PLUGINS/Wildtype.pm"
test -f "$NEOAG_NORMAL_PROTEOME_FASTA"
```

只运行与你所选 workflow 和已配置路径相关的检查。

## 常见错误与处理

| 现象 | 可能原因 | 处理 |
| --- | --- | --- |
| `neoag-v03: command not found` | 包未安装，或项目 `bin/` 不在 `PATH` | 运行 `source conf/tools.env.sh`，再运行 `python -m pip install -e '.[test]'`。 |
| `No module named neoag_v03` | 缺少 `PYTHONPATH` 或 editable install | 运行 `python -m pip install -e .`，或用 `PYTHONPATH=src python -m neoag_v03.cli ...`。 |
| `pytest: command not found` | 未安装 test extra | 运行 `python -m pip install -e '.[test]'`。 |
| `bin/neoag-nextflow: Permission denied` | archive/migration 后丢失执行位 | `find bin -maxdepth 1 -type f -exec chmod +x {} \;`。 |
| `conda not found` | 未安装或未初始化 Miniforge/Mambaforge | 安装 Miniforge 并打开新 shell，或 source 其 `etc/profile.d/conda.sh`。 |
| `mamba: unrecognized arguments -n ...` | 安装了 pip package `mamba`，不是 conda-forge mamba | 不要用 `pip install mamba`；使用脚本默认 conda 模式或安装真正的 conda/mamba。 |
| `CXXABI_1.3.15 not found` | MHCflurry/scipy 加载了旧系统 `libstdc++` | `conda install -n neoag-tools -c conda-forge 'libstdcxx-ng>=13'`，并确保 conda `lib` 在 `LD_LIBRARY_PATH` 最前。 |
| `mhcflurry-downloads fetch failed` | 网络或模型下载问题 | 激活 env 后重跑 `mhcflurry-downloads fetch`；offline 部署需预先准备 model data。 |
| `NetMHCpan MISSING` | 授权 tarball 未安装或 `NETMHCPAN_HOME` 错误 | 用 `bash scripts/install_netmhcpan.sh /path/to/tar.gz` 安装，然后 `source conf/tools.env.sh`。 |
| NetMHCpan wrapper 指向旧 conda 路径 | wrapper 用旧 conda base 创建 | 设置 `NEOAG_CONDA_BASE="$(conda info --base)"`，然后运行 `bash scripts/install_netmhcpan.sh --repair`。 |
| NetMHCpan 报缺少 `ld-linux-x86-64.so.2` | `neoag-tools` 中 conda sysroot 被移除 | 在 `neoag-tools` 中安装 `sysroot_linux-64` 和 `patchelf`；lite env 文件已保留这些依赖。 |
| `VEP cache not found` | Offline cache 缺失或 `NEOAG_VEP_CACHE` 错误 | 运行 `bash scripts/install_vep_cache.sh`，或在 `conf/tools.env.local.sh` 中设置 `NEOAG_VEP_CACHE`。 |
| `vep MISSING` 但 VEP 已安装 | VEP env/wrapper 未配置到路径 | 运行 `bash scripts/install_vep.sh`，source `conf/tools.env.sh`，并检查 `NEOAG_VEP_BIN`。 |
| VEP 运行时 `Can't locate DBI.pm` | 其他 conda env 的 Perl 环境污染 VEP | 使用 `bin/vep-neoag`，该 wrapper 会清理冲突 Perl 环境变量。 |
| `No CSQ annotations` | 输入 VCF 未经过 VEP 注释 | 配置 VEP 后使用 `run-full` 自动注释，或先运行 `neoag-v03 vep-annotate`。 |
| `.nextflow/history.lock (Permission denied)` | root-owned `.nextflow` metadata | 设置 `export NXF_HOME=/path/to/writable/cache` 并运行 `bin/neoag-nextflow`。 |
| `Downloading nextflow dependencies` 卡住 | 首次启动没有 cache 或网络受阻 | 预先填充 `NXF_HOME`，使用共享 cache，或允许网络完成下载。 |
| `Java not found` 或 Java 版本不支持 | Java 缺失或过旧 | 安装 OpenJDK 11+；用 `java -version` 验证。 |
| `work/`、`results/` 或 `tools/` 下 `Permission denied` | 目录归其他用户/root 所有 | 使用用户可写 out/work 目录，或请管理员修正 ownership。 |
| `GATK reference dictionary missing` | FASTA index/dict 缺失 | 运行 `bash scripts/download_ref_hg38.sh /path/to/ref/hg38`，或用 samtools/picard 创建 `.fai`/`.dict`。 |
| 真实数据 workflow 仍在使用 fixture paths | private run config 未编辑 | 复制 example config 到 private local config，并在生产前更新所有路径。 |
| 可选工具 missing 但 demo 可运行 | fixture demo 不需要该工具 | 仅在所选 workflow 需要时安装；见依赖矩阵。 |

## Release Boundary

不要提交或打包：

- `.git`, `.venv`, `.nextflow`, `.pytest_cache`
- `tools/`, `results/`, `work/`, `dist/`, `conda_packs/`
- `conf/tools.env.local.sh`
- `conf/site.config`
- `conf/private/*`
- `conf/*.private.toml`
- 真实患者数据或样本标识符
- 授权工具二进制文件
- `data/ref`、`data/vep` 等大型参考库

准备 online release 前，运行 `scripts/check_release_boundary.sh`。

## 更多文档

- `docs/TOOLS_SETUP.md`：外部工具安装详解。
- `docs/PROJECT_DATA_PATHS.md`：项目参考库/数据路径清单。
- `docs/INSTALL_AND_DATA.md`：安装和数据准备说明。
- `docs/V043_CCF21.md`：CCF 2.1 说明。
- `docs/V042_P1_APPM_EXPLAINABILITY.md`：APPM explainability 说明。
- `docs/V04_EVIDENCE_SAFETY_ESCAPE.md`：safety 和 immune-escape evidence 说明。
- `RELEASE.md`：release boundary 和测试总结。

## 解释边界

本流程用于研究分层和验证规划。排序候选需要结合 assay validation、疾病背景、HLA typing、tumor purity、expression/protein support、safety evidence、immune-escape context 以及适当的临床或科研治理流程进行复核。
