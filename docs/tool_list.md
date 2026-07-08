# 项目工具清单

## 变异注释

| 工具 | 可执行文件 | Conda 环境 | 说明 | Docker 镜像 | pip/conda 安装 |
|------|-----------|-----------|------|-------------|---------------|
| VEP 115 | `vep` | `neoag-vep` | Ensembl Variant Effect Predictor，变异功能注释 | ✅ `ensemblorg/ensembl-vep` (官方) | ✅ `conda install -c bioconda ensembl-vep` |

## 新抗原预测

| 工具 | 可执行文件 | Conda 环境 | 说明 | Docker 镜像 | pip/conda 安装 |
|------|-----------|-----------|------|-------------|---------------|
| pVACseq | `pvacseq` | `neoag-tools` | SNV/InDel 新抗原预测 | ✅ `griffithlab/pvactools` (官方) | ✅ `pip install pvactools` |
| pVACfuse | `pvacfuse` | `neoag-tools` | 融合新抗原预测 | ✅ `griffithlab/pvactools` (官方) | ✅ `pip install pvactools` |
| pVACsplice | `pvacsplice` | `neoag-tools` | 剪接位点新抗原预测 | ✅ `griffithlab/pvactools` (官方) | ✅ `pip install pvactools` |

> pVACseq/pVACfuse/pVACsplice 均为 pVACtools 套件的一部分，统一通过 `pip install pvactools` 安装。
> 另有 Biocontainers: `quay.io/biocontainers/pvacseq`

## MHC 结合与呈递预测

| 工具 | 可执行文件 | 环境/路径 | 说明 | Docker 镜像 | pip/conda 安装 |
|------|-----------|----------|------|-------------|---------------|
| NetMHCpan 4.2 | `netMHCpan` | `tools/netMHCpan/` | MHC-I 结合和 EL 预测（DTU 许可） | ❌ 无官方镜像 (许可限制) | ❌ 无；需从 DTU 下载源码包手动安装 |
| MHCflurry | `mhcflurry-predict` | `neoag-tools` | MHC-I 呈递预测 | ❌ 无官方镜像 | ✅ `pip install mhcflurry` / `conda install -c bioconda mhcflurry` |
| NetMHCstabpan | `netMHCstabpan` | `tools/netMHCstabpan/` | pMHC 稳定性预测 | ❌ 无官方镜像 (许可限制) | ❌ 无；需从 DTU 下载源码包手动安装 |

## 免疫原性预测

| 工具 | 可执行文件 | 环境/路径 | 说明 | Docker 镜像 | pip/conda 安装 |
|------|-----------|----------|------|-------------|---------------|
| PRIME | `PRIME` | `tools/prime/` | MHC-I 免疫原性预测 | ❌ 无 | ❌ 无；需从 [GfellerLab/PRIME](https://github.com/GfellerLab/PRIME) 下载 ZIP 并编译 C++ 源码 |
| BigMHC_IM | `bigmhc_predict` | `tools/bigmhc/` | 新表位免疫原性预测 | ❌ 无 | ✅ `pip install git+https://github.com/karchinlab/bigmhc.git` (非 PyPI) |
| DeepImmuno-CNN | `deepimmuno-cnn.py` | `tools/DeepImmuno/` | 9/10-mer MHC-I 免疫原性预测 | ❌ 无 | ✅ `pip install git+https://github.com/griffithlab/deepimmuno.git#egg=deepimmuno` (非 PyPI) |

## HLA 分型与 LOH

| 工具 | 可执行文件 | Conda 环境 | 说明 | Docker 镜像 | pip/conda 安装 |
|------|-----------|-----------|------|-------------|---------------|
| LOHHLA | `LOHHLA` | `tools/lohhla/` | HLA 等位基因特异性 LOH | ✅ `quay.io/biocontainers/lohhla` (Biocontainers) | ✅ `conda install -c bioconda lohhla` |
| SpecHLA | — | — | HLA 分型与拷贝数（仅解析输出） | ❌ 无 | ✅ `conda install -c bioconda spechla` |
| OptiType | — | `neoag-optitype` | HLA 分型 | ✅ `fred2/optitype` (官方) | ✅ `conda install -c bioconda optitype` |

> OptiType 另有 Biocontainers: `quay.io/biocontainers/optitype`

## 肿瘤纯度与拷贝数

| 工具 | 可执行文件 | Conda 环境 | 说明 | Docker 镜像 | pip/conda 安装 |
|------|-----------|-----------|------|-------------|---------------|
| FACETS | `runFACETS.R` | `tools/facets/` | 肿瘤纯度与拷贝数 | ⚠️ `uclahs-cds/docker-cnv_facets` (社区) | ✅ `conda install -c bioconda r-facets` |
| ASCAT | `ascat.R` | `neoag-ascat` | 等位基因特异性拷贝数与纯度 | ✅ `quay.io/biocontainers/ascat` (Biocontainers) | ✅ `conda install -c bioconda ascat` |
| ASCAT v3 | `ascat-v3` | `neoag-ascat-v3` | ASCAT v3 版本 | ✅ `quay.io/biocontainers/ascat` (Biocontainers) | ✅ `conda install -c bioconda ascat` |

## 克隆性分析

| 工具 | 可执行文件 | Conda 环境 | 说明 | Docker 镜像 | pip/conda 安装 |
|------|-----------|-----------|------|-------------|---------------|
| PyClone-VI | `pyclone-vi` | `neoag-pyclone` | 克隆簇 CCF 推断 | ❌ 无 | ✅ `conda install -c bioconda pyclone-vi` / `pip install git+https://github.com/Roth-Lab/pyclone-vi.git` |

## 变异检测

| 工具 | 可执行文件 | Conda 环境 | 说明 | Docker 镜像 | pip/conda 安装 |
|------|-----------|-----------|------|-------------|---------------|
| GATK4 | `gatk` | `neoag-gatk` | Mutect2 体细胞变异检测 + FilterMutectCalls | ✅ `broadinstitute/gatk` (官方) | ✅ `conda install -c bioconda gatk4` |

> GATK4 官方 Docker 镜像内已包含 samtools、bcftools、tabix、Python3 及常用 Python 包

## 融合检测

| 工具 | 可执行文件 | 路径 | 说明 | Docker 镜像 | pip/conda 安装 |
|------|-----------|------|------|-------------|---------------|
| STAR-Fusion | `star-fusion-neoag` | `tools/STAR-Fusion/` | 嵌合转录本检测 | ✅ `trinityctat/ctatfusion` (官方) | ✅ `conda install -c bioconda star-fusion` |
| FusionCatcher | `fusioncatcher-neoag` | `tools/fusioncatcher/` | 体细胞融合检测 | ✅ `quay.io/biocontainers/fusioncatcher` (Biocontainers) | ✅ `conda install -c bioconda fusioncatcher` |
| Arriba | `arriba` | `neoag-fusion` | RNA-seq 融合检测 | ✅ `uhrigs/arriba` (官方) | ✅ `conda install -c bioconda arriba` |
| EasyFuse | `easyfuse-neoag` | `tools/EasyFuse/` | 融合 meta-caller（整合以上三个 + ML） | ✅ `tronbioinformatics/easyfuse:1.3.4` (官方) | ✅ `pip install pyeasyfuse` / `conda install -c bioconda pyeasyfuse` |

> Arriba 另有 Biocontainers: `quay.io/biocontainers/arriba`
> STAR-Fusion 另有 Biocontainers: `quay.io/biocontainers/star-fusion`

## 结构变异检测

| 工具 | 环境 | 说明 | Docker 镜像 | pip/conda 安装 |
|------|------|------|-------------|---------------|
| Manta | `neoag-manta` | SV 检测 | ✅ `quay.io/biocontainers/manta` (Biocontainers) | ✅ `conda install -c bioconda manta` |
| SvABA | — | SV 检测 | ✅ `quay.io/biocontainers/svaba` (Biocontainers) | ✅ `conda install -c bioconda svaba` |
| GRIDSS2 | — | SV 检测 | ✅ `gridss/gridss` (官方) | ✅ `conda install -c bioconda gridss` |
| DELLY | — | SV 检测 | ✅ `dellytools/delly` (官方) | ✅ `conda install -c bioconda delly` |

> GRIDSS2 另有 Biocontainers: `quay.io/biocontainers/gridss`
> DELLY 另有社区镜像: `getwilds/delly`

## VEP 插件（自定义）

| 插件 | 路径 | 说明 |
|------|------|------|
| Wildtype | `work/vep_plugins/Wildtype.pm` | 输出突变位点附近的野生型蛋白序列 |
| Frameshift | `work/vep_plugins/Frameshift.pm` | 输出移码变异的替代蛋白序列 |

## 基础依赖

| 组件 | 环境 | 说明 | Docker 镜像 | pip/conda 安装 |
|------|------|------|-------------|---------------|
| IEDB API | 网络 | NetMHCpan/NetMHCstabpan 远程回退 | — | — |
| tabix | `neoag-tools` | VCF 索引 | ✅ `quay.io/biocontainers/tabix` (Biocontainers) | ✅ `conda install -c bioconda tabix` |
| samtools | `neoag-tools` | BAM/VCF 操作 | ✅ `quay.io/biocontainers/samtools` (Biocontainers) | ✅ `conda install -c bioconda samtools` |
| Python 3.11 | `neoag-tools` | 管线主运行环境 | — | — |
| Perl 5.32 | `neoag-vep` | VEP 运行环境（含 perl-dbi） | — | — |
| TensorFlow 2.21 | `neoag-tools` | MHCflurry 模型推理 | — | ✅ `pip install tensorflow` |
| R | `neoag-r` | FACETS/LOHHLA 等 R 脚本运行 | — | ✅ `conda install -c conda-forge r-base` |

## Conda 环境一览

| 环境 | 用途 |
|------|------|
| `neoag-tools` | 主环境：pVACtools、MHCflurry、Python 3.11 |
| `neoag-vep` | VEP + Perl 5.32 + perl-dbi |
| `neoag-gatk` | GATK4 Mutect2 |
| `neoag-fusion` | 融合检测（Arriba 等） |
| `neoag-ascat` | ASCAT 拷贝数 |
| `neoag-ascat-v3` | ASCAT v3 |
| `neoag-pyclone` | PyClone-VI |
| `neoag-optitype` | OptiType HLA 分型 |
| `neoag-r` | R 语言环境 |
| `neoag-bm` | Benchmark 评估 |

---

## 安装方式汇总

### Docker 镜像可用性

| 状态 | 工具 |
|------|------|
| ✅ 有官方 Docker 镜像 | VEP 115, pVACtools, OptiType, GATK4, STAR-Fusion, Arriba, EasyFuse, GRIDSS2, DELLY |
| ✅ 有 Biocontainers 镜像 | LOHHLA, ASCAT, ASCAT v3, FusionCatcher, Manta, SvABA, samtools, tabix |
| ⚠️ 仅有社区镜像 | FACETS (`uclahs-cds/docker-cnv_facets`), NetMHCpan 4.2 (`ghcr.io/macromnex/netmhcpan_mcp`) |
| ❌ 无 Docker 镜像 | NetMHCstabpan, PRIME, BigMHC_IM, DeepImmuno-CNN, SpecHLA, PyClone-VI, MHCflurry |

全部镜像拉取完成，共 17 个镜像，约 48 GB：

  官方 Docker 镜像（8 个）

  ┌─────────────────────────────┬────────────────┬─────────┐
  │            镜像             │      Tag       │  大小   │
  ├─────────────────────────────┼────────────────┼─────────┤
  │ ensemblorg/ensembl-vep      │ release_115.2  │ 766 MB  │
  ├─────────────────────────────┼────────────────┼─────────┤
  │ griffithlab/pvactools       │ 6.1.1          │ 19.5 GB │
  ├─────────────────────────────┼────────────────┼─────────┤
  │ fred2/optitype              │ release-v1.3.1 │ 1.87 GB │
  ├─────────────────────────────┼────────────────┼─────────┤
  │ broadinstitute/gatk         │ 4.6.2.0        │ 5.75 GB │
  ├─────────────────────────────┼────────────────┼─────────┤
  │ uhrigs/arriba               │ 2.5.1          │ 1.14 GB │
  ├─────────────────────────────┼────────────────┼─────────┤
  │ tronbioinformatics/easyfuse │ 1.3.7          │ 9.56 GB │
  ├─────────────────────────────┼────────────────┼─────────┤
  │ gridss/gridss               │ 2.13.2         │ 5.22 GB │
  ├─────────────────────────────┼────────────────┼─────────┤
  │ dellytools/delly            │ v2.3.0         │ 20 MB   │
  └─────────────────────────────┴────────────────┴─────────┘

  Biocontainers 镜像（9 个）

  ┌─────────────────────────────────────┬──────────────────────┬─────────┐
  │                镜像                 │         Tag          │  大小   │
  ├─────────────────────────────────────┼──────────────────────┼─────────┤
  │ quay.io/biocontainers/lohhla        │ 20171108--hdfd78af_3 │ 1.51 GB │
  ├─────────────────────────────────────┼──────────────────────┼─────────┤
  │ quay.io/biocontainers/ascat         │ 2.5.2--r40hdfd78af_3 │ 753 MB  │
  ├─────────────────────────────────────┼──────────────────────┼─────────┤
  │ quay.io/biocontainers/ascat         │ 3.2.0--r44hdfd78af_1 │ 1.16 GB │
  ├─────────────────────────────────────┼──────────────────────┼─────────┤
  │ quay.io/biocontainers/star-fusion   │ 1.15.1--hdfd78af_1   │ 1.60 GB │
  ├─────────────────────────────────────┼──────────────────────┼─────────┤
  │ quay.io/biocontainers/fusioncatcher │ 1.33b--hdfd78af_0    │ 1.13 GB │
  ├─────────────────────────────────────┼──────────────────────┼─────────┤
  │ quay.io/biocontainers/manta         │ 1.6.0--h9ee0642_3    │ 178 MB  │
  ├─────────────────────────────────────┼──────────────────────┼─────────┤
  │ quay.io/biocontainers/svaba         │ 1.2.0--h69ac913_1    │ 73.5 MB │
  ├─────────────────────────────────────┼──────────────────────┼─────────┤
  │ quay.io/biocontainers/samtools      │ 1.23.1--ha83d96e_0   │ 80.3 MB │
  ├─────────────────────────────────────┼──────────────────────┼─────────┤
  │ quay.io/biocontainers/tabix         │ 1.11--hdfd78af_0     │ 94.3 MB │
  └─────────────────────────────────────┴──────────────────────┴─────────┘



### pip 安装可用性

| 状态 | 工具 |
|------|------|
| ✅ PyPI 官方发布 | pVACtools (`pip install pvactools`), MHCflurry (`pip install mhcflurry`), EasyFuse (`pip install pyeasyfuse`) |
| ✅ GitHub pip 安装 | PyClone-VI, BigMHC_IM, DeepImmuno-CNN |
| ❌ 无 pip 安装 | VEP 115, NetMHCpan 4.2, NetMHCstabpan, PRIME, LOHHLA, SpecHLA, OptiType, FACETS, ASCAT, ASCAT v3, GATK4, STAR-Fusion, FusionCatcher, Arriba, Manta, SvABA, GRIDSS2, DELLY |

### conda 安装可用性

| 状态 | 工具 |
|------|------|
| ✅ Bioconda/Conda-forge | VEP 115, MHCflurry, LOHHLA, SpecHLA, OptiType, FACETS, ASCAT, ASCAT v3, PyClone-VI, GATK4, STAR-Fusion, FusionCatcher, Arriba, EasyFuse, Manta, SvABA, GRIDSS2, DELLY, samtools, tabix |
| ❌ 无 conda 安装 | NetMHCpan 4.2, NetMHCstabpan, PRIME, BigMHC_IM, DeepImmuno-CNN, pVACtools (仅 pip) |

### 许可限制（需单独申请许可证）

| 工具 | 许可方式 |
|------|---------|
| NetMHCpan 4.2 | DTU Health Tech 学术许可，需从 [services.healthtech.dtu.dk](https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/) 申请下载 |
| NetMHCstabpan | DTU Health Tech 学术许可，需从 [services.healthtech.dtu.dk](https://services.healthtech.dtu.dk/services/NetMHCstabpan-1.0/) 申请下载 |

### 需手动编译/安装的工具

| 工具 | 安装方式 |
|------|---------|
| PRIME | 从 [GfellerLab/PRIME](https://github.com/GfellerLab/PRIME) 下载 ZIP，`g++ -O3 PRIME.cc -o PRIME.x` 编译 |
| DeepImmuno-CNN | `pip install git+https://github.com/griffithlab/deepimmuno.git#egg=deepimmuno`，依赖 Python 3.6~3.7 / TensorFlow<2.16 |
| BigMHC_IM | `pip install git+https://github.com/karchinlab/bigmhc.git`，依赖 PyTorch |
