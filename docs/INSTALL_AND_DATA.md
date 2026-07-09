# Installation, Tool, and Data Setup

部署检查清单：基础环境 → Python/conda → Nextflow → **按入口需要装什么工具** → 参考数据 → 验收命令 → 常见错误。

在线发布包本身很轻量：只有源码、workflow、测试、示例配置、小fixture。不含licensed工具、大参考数据、
真实患者数据、`results/`/`work/`/`conda_packs/`/Nextflow缓存。

> **最快摸底方式**：不确定自己需要装哪些工具，先装好Python包，直接对着你要用的入口跑一次
> `run-demo --entry-mode`（6选1：`snv_indel`/`fusion`/`splice_junction`/`sv_wgs`/`sv_wes`/`peptide_only`）——
> 它会先打印这个入口需要的全部工具的OK/MISSING状态，
> 再用bundled fixture跑一遍全流程冒烟测试。这比读完整张工具表快得多，本文档第4节的表格和它读的是同一份
> `conf/tools.toml`。

## 1. Base System Environment

| Component | Recommended | Notes |
| --- | --- | --- |
| OS | Linux x86_64 | 项目里的wrapper脚本和大多数上游生信工具假定Linux |
| Shell | Bash | 脚本都是Bash写的 |
| CPU/RAM | 8+ CPU, 32+ GB RAM（真实数据） | fixture demo很小；患者级WES/WGS/fusion需要更多 |
| Disk | 200 GB+（真实数据部署） | VEP cache、hg38参考、fusion参考、Nextflow work目录都很大 |
| Network | 首次安装/下载需要 | 离线部署需预先准备好tarball/conda包/参考数据/Nextflow缓存 |
| Java | Java 11+ | Nextflow和部分工具需要 |
| Conda/Mamba | 推荐Miniforge/Mambaforge | 工具环境定义在`conda/`下 |

```bash
sudo apt-get update
sudo apt-get install -y bash coreutils curl wget git tar gzip unzip bzip2 xz-utils ca-certificates build-essential openjdk-17-jre-headless
```

没有`sudo`权限的话，找站点管理员把这些工具装到`PATH`上。

## 2. Python And Conda Environment

### 2.1 建主环境

```bash
bash scripts/setup_tools_env.sh
source conf/tools.env.sh
python -m pip install -e '.[test]'
neoag-v03 check-tools
```

只想跑测试、不需要重量级工具栈的轻量开发环境：

```bash
NEOAG_TOOLS_LITE=1 bash scripts/setup_tools_env.sh
source conf/tools.env.sh
python -m pip install -e '.[test]'
pytest -q
```

默认`pytest -q`是发布安全的快速检查；`pytest -q --run-all`是维护者级别的检查（含integration/benchmark/external-tool/长workflow），
只有在工具、参考数据、Nextflow缓存都装好之后再跑。

| 变量 | 作用 | 典型值 |
| --- | --- | --- |
| `NEOAG_PROJECT_ROOT` | 源码checkout路径 | 由`conf/tools.env.sh`自动设置 |
| `NEOAG_TOOLS_ROOT` | 外部工具/参考数据根目录 | `/path/to/neoag_artifact_bundle` |
| `NEOAG_CONDA_BASE` | Miniforge/Mambaforge根目录 | `/opt/miniforge3` |
| `NEOAG_CONDA_ENV` | 主环境名 | `neoag-tools` |
| `NEOAG_FORCE_CPU` | 禁用GPU探测 | CPU-only节点上设为`1` |

站点专属路径不要写进仓库，复制一份本地override（已在`.gitignore`里，不会被提交）：

```bash
cp conf/tools.env.local.example.sh conf/tools.env.local.sh
# 然后编辑 conf/tools.env.local.sh
```

### 2.2 验证CLI装好了

```bash
python -m pip install -e '.[test]'
neoag-v03 run-demo --entry-mode snv_indel --outdir work/demo_snv --sample-id DEMO001
```

预期产出：`scoring/ranked_peptides.v03.tsv`、`scoring/ranked_events.v03.tsv`、`scoring/validation_plan.v03.tsv`、
`reports/evidence_report.v03.html`、`provenance.v03.json`。其余5个入口的验证命令见
[`USAGE_GUIDE.md`](USAGE_GUIDE.md)。

## 3. Nextflow, Java, And Cache Setup

用项目自带的wrapper，不要直接调`nextflow`——它会优先用当前checkout的`bin/neoag-v03`、设置`PYTHONPATH=src`，
并把Nextflow元数据存到仓库外面。

```bash
export NXF_HOME=/path/to/writable/nextflow_cache
bin/neoag-nextflow -version
bin/neoag-nextflow run workflows/main.nf -w /tmp/neoag_nf_work \
  --pvac_files data/fixtures/pvacseq_aggregated.tsv --outdir results/demo_nf --sample_id NF_DEMO
```

- **在线模式**：包里自带轻量`bin/nextflow`launcher，首次运行可能需要联网下载运行时依赖到`NXF_HOME`。
- **离线模式**：预先准备好Java、Nextflow运行时缓存、conda/容器资源和全部外部参考数据，把`NXF_HOME`指向这份预填好的可写缓存。
- 共享集群上把`NXF_HOME`指到共享可写缓存，避免重复下载。
- 不要从root拥有的`.nextflow`目录启动。如果看到`.nextflow/history.lock (Permission denied)`，
  设置`NXF_HOME`到可写路径并用`bin/neoag-nextflow`。

## 4. 按入口需要装什么工具

整个流程分三层：**入口专属上游工具**（每个入口不一样）→ **结合/免疫原性预测**（6个入口共用）→
**样本上下文证据**（6个入口共用，全部可选）。下面按这个结构列，而不是甩一张大表让你自己筛——
想知道某个具体工具属于哪层/哪个入口，看`conf/tools.toml`里对应的`[tool]`段和`entries`字段，
这份文档和那份TOML是同一个信息源，没有第二套。

### 4.1 结合/免疫原性预测（presentation，6个入口全部共用）

| 工具 | 必需？ | 安装 | 关键变量 | 验证 |
| --- | --- | --- | --- | --- |
| NetMHCpan | 真实结合预测必需，除非用`--stub`/IEDB API回退 | `bash scripts/install_netmhcpan.sh /path/to/netMHCpan-4.2*.tar.gz` | `NEOAG_NETMHCPAN_BIN`, `NEOAG_NETMHCPAN_BACKEND` | `neoag-v03 check-tools` |
| MHCflurry | 可选，NetMHCpan的替代/补充 | `bash scripts/setup_tools_env.sh`；需要时`mhcflurry-downloads fetch` | `NEOAG_CONDA_ENV` | `neoag-v03 check-tools` |
| NetMHCstabpan | 可选，稳定性证据 | `bash scripts/install_netmhcstabpan.sh --iedb`（或本地tarball） | `NETMHCSTABPAN_HOME` | `neoag-v03 check-tools` |
| PRIME / MixMHCpred / BigMHC | 可选，免疫原性证据 | `bash scripts/install_immunogenicity_tools.sh` | `PRIME_HOME`, `MIXMHCPRED_HOME`, `BIGMHC_DIR` | `neoag-v03 check-tools` |
| DeepImmuno | 可选，免疫原性证据 | `bash scripts/install_deepimmuno.sh` | `DEEPIMMUNO_DIR` | `neoag-v03 check-tools` |

一个都不装也能跑：`peptide-predict --stub` 或 `run-demo`（demo默认用bundled预测结果fixture代替）。

### 4.2 样本上下文证据（sample context，6个入口全部共用，全部可选）

| 工具 | 用途 | 安装 | 关键变量 | 验证 |
| --- | --- | --- | --- | --- |
| LOHHLA + Polysolver | HLA杂合性丢失证据（免疫逃逸层用） | 按LOHHLA license/文档安装，本地env设路径 | `LOHHLA_HOME`, `POLYSOLVER_HOME`, `NOVOALIGN_LICENSE_FILE` | `neoag-v03 check-tools` |
| SpecHLA | HLA拷贝数/LOH转换，和LOHHLA交叉验证 | 外部安装 | Site-specific | `neoag-v03 convert-spechla ...` / `crosscheck-hla-loh ...` |
| FACETS | 纯度/CNV/LOH证据 | `bash scripts/install_ascat_pyclone.sh`或站点R安装 | `FACETS_HOME`, `NEOAG_DBSNP_VCF` | `neoag-v03 check-tools` |
| ASCAT 2.5.2 | CNV/LOH证据（legacy基线） | `bash scripts/install_ascat_pyclone.sh` | `NEOAG_ASCAT_ENV`, `ASCAT_HOME` | `neoag-v03 check-tools` |
| ASCAT 3.2.0 | CNV/LOH交叉验证 | `conda env create -f conda/env.neoag-ascat-v3.yml` | `NEOAG_ASCAT_V3_ENV`, `NEOAG_ASCAT_V3_BIN` | `bin/ascat-v3 --check` |
| PyClone-VI | 外部克隆性反卷积（CCF 2.1可选输入） | `bash scripts/install_ascat_pyclone.sh` | `NEOAG_PYCLONE_ENV`, `NEOAG_PYCLONE_BIN` | `neoag-v03 check-tools` |

不装这些，APPM/CCF/免疫逃逸证据层会自动降级为"证据不完整"，不会报错中断。

### 4.3 各入口专属工具

| 入口 | 专属工具 | 必需？ | 安装 | 关键变量 |
| --- | --- | --- | --- | --- |
| **A｜SNV/InDel** | VEP | 走滑窗产肽路径（`snv-build-raw`）时必需，除非VCF已有CSQ注释 | `bash scripts/install_vep.sh`；cache用`bash scripts/install_vep_cache.sh` | `NEOAG_VEP_BIN`, `NEOAG_VEP_CACHE`, `NEOAG_VEP_PLUGINS` |
| | GATK4/Mutect2 | 仅当从BAM起（`snv-call-wes`）时必需 | `bash scripts/install_gatk.sh` | `NEOAG_GATK_ENV` |
| | pVACtools（pvacseq） | 可选，走pVACseq产肽路径时用（P1路径，二选一，见下方说明） | `bash scripts/setup_tools_env.sh`或Docker `scripts/pull_docker_tools.sh` | `NEOAG_PVAC_DOCKER`, `NEOAG_PVAC_WORKDIR` |
| **B｜Fusion** | EasyFuse | 可选，跑EasyFuse本身（Nextflow流程）时需要；如果已有`fusions.pass.csv`就不需要装 | 外部安装；EasyFuse conda envs用`bash scripts/seed_easyfuse_conda_envs.sh` | `NEOAG_EASYFUSE_HOME`, `NEOAG_EASYFUSE_REF` |
| | STAR-Fusion / FusionCatcher | 可选，EasyFuse内部调用 | 外部安装/挂载 | `NEOAG_STAR_FUSION_HOME`, `NEOAG_FUSIONCATCHER_BIN` |
| | pVACfuse | 可选，补充HLA分型表位 | 同pVACtools | 同上 |
| **C｜Splice junction** | VEP | 同Entry A（给变异层做CSQ注释） | 同上 | 同上 |
| | pVACsplice | **真正产生肽段的来源**，必需（`--splice-junction-tsv`只做富集，不产肽） | 同pVACtools | 同上 |
| **D1｜SV-WGS** | Manta / GRIDSS / SvABA | 可选，跑SV caller本身时需要；如果已有SV VCF就不需要装 | 外部安装或站点conda/module | `NEOAG_SV_ENV`, `NEOAG_MANTA_ENV` |
| **D2｜SV-WES** | 同D1 + capture BED | 同上，另加capture BED文件（不是"工具"，是数据） | 同上 | 同上 |
| **E｜Peptide-only CSV** | 无 | 纯Python标准化，不依赖任何外部工具 | — | — |

licensed工具（NetMHCpan、NetMHCstabpan、LOHHLA、Novoalign/Polysolver组件）可能需要学术/机构授权，不要把它们的二进制打进在线发布包。

### 4.4 一次性核对全部工具

```bash
source conf/tools.env.sh
neoag-v03 check-tools          # 全量核对，覆盖上面4.1-4.3全部工具，适合部署/CI
neoag-v03 run-demo --entry-mode snv_indel \
  --outdir /tmp/demo --sample-id DEMO   # 只核对某一个入口需要的工具，适合日常开发（其余5个入口同理，把entry-mode换掉）
```

## 5. Reference Data Download Table

大数据放在`NEOAG_TOOLS_ROOT`或站点管理的参考数据区，不要放进源码checkout里。

| 数据/参考 | 用途 | 下载/准备命令 | 变量/路径 | 验证 |
| --- | --- | --- | --- | --- |
| VEP cache（GRCh38 release 105） | 离线VEP注释 | `bash scripts/install_vep_cache.sh` | `NEOAG_VEP_CACHE` | `test -d "$NEOAG_VEP_CACHE/homo_sapiens"` |
| GRCh38 FASTA+索引 | VEP产肽、GATK、SV产肽 | `bash scripts/download_ref_hg38.sh /path/to/ref/hg38` | `NEOAG_REFERENCE_FASTA` | `test -f "$NEOAG_REFERENCE_FASTA"` |
| dbSNP/常见SNP VCF | FACETS `snp-pileup`等CNV流程 | 随hg38 bundle一起或站点提供 | `NEOAG_DBSNP_VCF` | `test -f "$NEOAG_DBSNP_VCF"` |
| gnomAD AF VCF + PoN | GATK Mutect2过滤 | `bash scripts/download_ref_hg38.sh /path/to/ref/hg38` | 写在run config里 | `test -f /path/to/af-only-gnomad.hg38.vcf.gz` |
| Ensembl protein FASTA | 肽段安全性正常/参考蛋白质组过滤 | 手动下载或站点bundle | `NEOAG_NORMAL_PROTEOME_FASTA` | `test -f "$NEOAG_NORMAL_PROTEOME_FASTA"` |
| Normal expression表 | 肽段安全性证据 | 站点生成，或用`resources/normal_expression.example.tsv`跑fixture | CLI参数/run config路径 | 表头匹配预期schema |
| Normal HLA ligand表 | 肽段安全性证据 | 站点生成，或用`resources/normal_hla_ligands.example.tsv`跑fixture | CLI参数/run config路径 | 表头匹配预期schema |
| RNA allele-count/RNA VAF TSV | RNA变异支持证据 | `scripts/rna_allele_counts_pysam.py`或站点RNA genotyper | `build-evidence-layer --rna-vaf` | 含`event_id`/`rna_alt_reads`/`rna_vaf`等列 |
| CTAT genome lib | STAR-Fusion | 按STAR-Fusion/CTAT文档下载或挂载站点bundle | `NEOAG_CTAT_LIB_DIR` | `test -d "$CTAT_GENOME_LIB"` |
| EasyFuse参考 | EasyFuse流程 | 按EasyFuse文档下载或挂载站点bundle | `NEOAG_EASYFUSE_REF` | `test -d "$NEOAG_EASYFUSE_REF"` |
| GTF注释 | SV/融合产肽 | 用和参考FASTA匹配的GENCODE/Ensembl GTF | CLI参数`--gencode-gtf` | `test -f /path/to/genes.gtf` |
| Capture BED | WES SV Phase 1.5 | 用panel/exome捕获BED | CLI参数`--capture-bed` | `test -f /path/to/capture.bed` |
| HLA分型文件 | 结合预测、SV流程 | 站点HLA分型结果，每行一个allele | CLI参数`--hla` | `head /path/to/hla.txt` |

推荐的外部bundle目录结构：

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

```bash
export NEOAG_TOOLS_ROOT=/path/to/neoag_artifact_bundle
export NEOAG_SHARED_REF_DIR=/path/to/shared_refs
source conf/tools.env.sh
neoag-v03 check-tools
```

## 6. 常用命令一览

需要更完整的"6个入口每个入口的分步命令+中间文件+参数说明"，见[`USAGE_GUIDE.md`](USAGE_GUIDE.md)；
这里只列日常最常用的几条：

| 命令 | 用途 | 依赖 |
| --- | --- | --- |
| `neoag-v03 run-demo --entry-mode snv_indel --outdir <dir> --sample-id <id>` | fixture冒烟测试+按入口工具核对（`--entry-mode`可换成其余5个入口） | 无（bundled fixture） |
| `neoag-v03 check-tools` | 全部工具的可用性总览 | 无 |
| `neoag-v03 run-upstream --config conf/run.stub.toml --outdir results/upstream` | 只跑上游生信工具 | 取决于config里enable的工具 |
| `neoag-v03 run-full --config <toml> --outdir results/<sample>` | 上游+打分一键跑完 | 同上 |
| `neoag-v03 peptide-predict -i <raw_peptides.tsv> -o <outdir>/presentation` | 只跑结合/免疫原性预测 | 见4.1（缺了就用`--stub`） |
| `bin/neoag-nextflow run workflows/main.nf -w /tmp/nf_work --pvac_files data/fixtures/pvacseq_aggregated.tsv --outdir results/demo_nf --sample_id NF_DEMO` | Nextflow方式跑fixture | Java/Nextflow运行时 |

`neoag-v03 <command> --help`查看完整参数。

## 7. Installation Acceptance Commands

安装完之后从项目根目录跑这些。

### 7.1 基础包验收

```bash
source conf/tools.env.sh
python -m pip install -e '.[test]'
pytest -q
neoag-v03 run-demo --entry-mode snv_indel --outdir work/demo_snv --sample-id DEMO001
# RNA VAF / junction证据验收（用你自己的raw表）
neoag-v03 build-evidence-layer --outdir results/sample --profile default \
  --raw-events results/sample/parsed/raw_events.tsv \
  --raw-peptides results/sample/parsed/raw_peptides.tsv \
  --rna-vaf results/sample/parsed/rna_vaf.tsv \
  --rna-junction results/sample/parsed/rna_junctions.tsv
# HLA LOH交叉验证（转换完LOHHLA和SpecHLA输出之后）
neoag-v03 crosscheck-hla-loh \
  --lohhla-hla-loh results/sample/tools/lohhla.hla_loh.tsv \
  --spechla-hla-loh results/sample/tools/spechla.hla_loh.tsv \
  --out results/sample/tools/hla_loh.crosscheck.tsv \
  --consensus-out results/sample/tools/hla_loh.consensus.tsv
```

### 7.2 工具可见性验收

```bash
source conf/tools.env.sh
neoag-v03 check-tools              # 全量核对（部署/CI用）
bash scripts/check_tools_env.sh
neoag-v03 run-demo --entry-mode snv_indel --outdir /tmp/demo --sample-id DEMO   # 按入口核对（日常开发用）
```

`check-tools`对当前workflow不需要的可选工具报MISSING是正常的；生产运行前，确认**你选定的workflow**需要的每个工具都是`OK`。

### 7.3 Nextflow验收

```bash
export NXF_HOME=/path/to/writable/nextflow_cache
bin/neoag-nextflow -version
bin/neoag-nextflow run workflows/main.nf -w /tmp/neoag_nf_work \
  --pvac_files data/fixtures/pvacseq_aggregated.tsv --outdir results/demo_nf --sample_id NF_DEMO
```

预期产出：`results/demo_nf/scoring/ranked_peptides.v03.tsv`、`ranked_events.v03.tsv`、
`reports/evidence_report.v041.html`、`provenance/workflow_provenance.yml`。

### 7.4 参考文件验收

```bash
test -f "$NEOAG_REFERENCE_FASTA"
test -d "$NEOAG_VEP_CACHE/homo_sapiens"
test -f "$NEOAG_NORMAL_PROTEOME_FASTA"
```

只跑和你选定workflow、已配置路径相关的检查项。

## 8. Common Errors And Fixes

| 症状 | 可能原因 | 修复 |
| --- | --- | --- |
| `neoag-v03: command not found` | 包没装，或项目`bin/`不在`PATH`上 | 先`source conf/tools.env.sh`，再`python -m pip install -e '.[test]'` |
| `No module named neoag_v03` | `PYTHONPATH`或editable install缺失 | `python -m pip install -e .`，或`PYTHONPATH=src python -m neoag_v03.cli ...` |
| `pytest: command not found` | 测试依赖没装 | `python -m pip install -e '.[test]'` |
| `conda not found` | Miniforge/Mambaforge没装或没初始化 | 安装Miniforge后开新shell，或`source`它的`etc/profile.d/conda.sh` |
| `mhcflurry-downloads fetch failed` | 网络/模型下载问题 | 激活环境重跑`mhcflurry-downloads fetch`；离线部署需预先准备模型数据 |
| `NetMHCpan MISSING` | licensed tarball没装，或`NETMHCPAN_HOME`路径不对 | `bash scripts/install_netmhcpan.sh /path/to/tar.gz`，然后`source conf/tools.env.sh` |
| `VEP cache not found` | 离线cache缺失，或`NEOAG_VEP_CACHE`路径不对 | `bash scripts/install_vep_cache.sh`，或在`conf/tools.env.local.sh`里设`NEOAG_VEP_CACHE` |
| `Cannot detect peptide column`（跑`peptide-predict`时） | 喂给它的是`extract-variant-peptides`的原始输出（`variant_peptides.tsv`），列名不兼容 | 改用`snv-build-raw`——它会自动转换成`peptide-predict`能识别的`raw_peptides.tsv`，见[`USAGE_GUIDE.md`](USAGE_GUIDE.md) |
| fusion/splice入口`raw_peptides.tsv`只有表头，0条数据 | `build-intermediates`缺了`--hla`（fusion）或`--pvac`（splice，真正的产肽来源） | 补上对应参数，见`USAGE_GUIDE.md`里Entry B/C的表格 |
| `.nextflow/history.lock (Permission denied)` | `.nextflow`元数据被root占用 | `export NXF_HOME=/path/to/writable/cache`，用`bin/neoag-nextflow` |
| `Downloading nextflow dependencies`卡住 | 首次运行没缓存，或网络受限 | 预填`NXF_HOME`、用共享缓存，或放行网络直到下载完成 |
| `Java not found`或版本太旧 | Java缺失/太老 | 装OpenJDK 11+，`java -version`验证 |
| `work/`/`results/`/`tools/`下`Permission denied` | 目录属于别的用户/root | 换一个当前用户可写的输出/work目录，或找管理员改权限 |
| `GATK reference dictionary missing` | FASTA索引/dict缺失 | `bash scripts/download_ref_hg38.sh /path/to/ref/hg38`，或用samtools/picard建`.fai`/`.dict` |
| 真实数据workflow却在用fixture路径 | 私有run config没改 | 复制一份示例配置到私有本地配置，改完全部路径再跑生产 |
| 某个可选工具显示missing但demo能跑 | fixture demo本来就不需要这个工具 | 只有你选定的workflow真的需要它时才装，见第4节 |

## 9. Release Boundary Reminder

不要提交或打包：

- `conf/tools.env.local.sh`
- `conf/site.config`
- `conf/private/*`
- 真实患者数据、样本标识符、患者专属脚本、站点本地绝对路径
- licensed工具二进制
- 大参考数据
- `tools/`、`results/`、`work/`、`dist/`、`conda_packs/`、`.nextflow*`

准备在线发布前跑`scripts/check_release_boundary.sh`。完整的打包边界见[`RELEASE_BOUNDARY.md`](RELEASE_BOUNDARY.md)。
