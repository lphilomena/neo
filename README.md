# NeoAg Event Pipeline v0.4.3 Online Release

NeoAg Event Pipeline 是研究型肿瘤新抗原候选优先级分析流程。它将 SNV/InDel、fusion、splice、SV 或 peptide-only 候选统一转换为事件表和 peptide-HLA 表，并叠加 HLA 呈递、APPM、CCF、safety、immune escape、validation plan 和 evidence report。

> **边界说明**：本流程输出的是 computational triage / 候选优先级结果，不是临床诊断、临床耐药判定或已验证治疗方案。

---

## 1. 本版迁移测试结论

本 README 根据 169 机器迁移测试结果更新。测试中已经确认：

- `python -m pip install -e .` 后，`neoag-v03 run-demo --outdir work/demo_v043 --sample-id DEMO001` 可成功运行。
- `pytest -q` 可完成轻量测试；测试记录中为 `175 passed, 95 skipped`。
- `bin/neoag-nextflow run workflows/main.nf --pvac_files data/fixtures/pvacseq_aggregated.tsv --outdir results/demo_nf` 在给 `bin/` 下文件补执行权限后可成功运行。
- 基础系统包、`setup_tools_env.sh`、GATK、VEP、DeepImmuno、NetMHCstabpan、NetMHCpan 在修复路径/环境问题后可安装或识别。
- 测试中暴露的问题包括：`bin/` 文件权限、MHCflurry 模型下载时 C++ runtime、NetMHCpan 脚本硬编码 conda 路径、ASCAT/PyClone 脚本误用 pip 版 mamba、VEP env 未加入 PATH、PRIME/BigMHC/MixMHCpred 依赖和变量缺失、LOHHLA 缺安装说明、FACETS/ASCAT wrapper 缺失。

本版已将上述问题对应的安装说明写入 README 和 `docs/TOOLS_SETUP.md`，并修正相关安装脚本的默认行为。

---

## 2. 推荐安装顺序

请从项目根目录执行。推荐使用 conda/miniconda，不建议将所有工具安装在系统 Python。

### 2.1 系统依赖

Ubuntu/Debian：

```bash
sudo apt-get update
sudo apt-get install -y \
  bash coreutils curl wget git tar gzip unzip bzip2 xz-utils \
  ca-certificates build-essential openjdk-17-jre-headless rsync file
```

RHEL/CentOS 请安装等价组件：`bash coreutils curl wget git tar gzip unzip bzip2 xz ca-certificates gcc gcc-c++ make java-17-openjdk-headless rsync file`。

### 2.2 解压后先修复可执行权限

部分迁移环境会丢失 `bin/` 执行权限。若出现 `Permission denied`，执行：

```bash
find bin -maxdepth 1 -type f -exec chmod +x {} \;
find scripts -maxdepth 1 -type f -name '*.sh' -exec chmod +x {} \;
```

### 2.3 安装 Python 包并跑 demo

```bash
python -m pip install -e .
neoag-v03 run-demo --outdir work/demo_v043 --sample-id DEMO001
```

关键输出：

```text
work/demo_v043/scoring/ranked_peptides.v03.tsv
work/demo_v043/scoring/ranked_events.v03.tsv
work/demo_v043/reports/evidence_report.v03.html
work/demo_v043/reports/evidence_report.patient.html
work/demo_v043/reports/evidence_report.technical.html
work/demo_v043/appm/appm_summary.tsv
work/demo_v043/clonality/ccf_lite.tsv
work/demo_v043/safety/peptide_safety.tsv
work/demo_v043/immune_escape/peptide_escape_flags.tsv
```

`.v03.tsv` 是 schema 兼容文件名，不代表软件版本。当前软件版本为 v0.4.3。

### 2.4 运行轻量测试

```bash
python -m pip install -e '.[test]'
pytest -q
```

迁移测试记录中结果为：

```text
175 passed, 95 skipped
```

`skipped` 多为外部工具、benchmark 或真实数据测试。真实生产运行前仍需要逐个安装并检查外部工具。

### 2.5 Nextflow fixture 测试

```bash
export NXF_HOME=${PWD}/.nextflow_cache
bin/neoag-nextflow run workflows/main.nf \
  --pvac_files data/fixtures/pvacseq_aggregated.tsv \
  --outdir results/demo_nf \
  --sample_id NF_DEMO
```

如果出现 `bin/neoag-nextflow: Permission denied`，执行第 2.2 节的 `chmod` 命令后重试。

---

## 3. 外部工具安装总览

先安装基础环境：

```bash
bash scripts/setup_tools_env.sh
source conf/tools.env.sh
neoag-v03 check-tools
```

`setup_tools_env.sh` 主要安装/检查 pVACtools 和 MHCflurry。真实样本还需要根据入口模式安装其他工具。

### 工具与用途

| 工具 | 用途 | 是否必需 | 安装命令/说明 |
|---|---|---:|---|
| pVACtools | 兼容 pVACseq/pVACfuse/pVACsplice 表格 | 可选 | `bash scripts/setup_tools_env.sh` |
| MHCflurry | HLA-I binding/presentation | 推荐 | `bash scripts/setup_tools_env.sh` 后运行 `mhcflurry-downloads fetch` |
| NetMHCpan 4.2 | 主要 HLA binding/EL 预测 | 生产推荐 | 需 DTU 学术许可包；见 4.2 |
| NetMHCstabpan | pMHC 稳定性 | 可选 | `bash scripts/install_netmhcstabpan.sh --iedb` |
| VEP | VCF 注释、蛋白序列提取 | SNV/InDel 推荐 | `bash scripts/install_vep.sh`；cache 见 4.3 |
| GATK | Mutect2/FilterMutectCalls | WES/WGS upstream 可选 | `bash scripts/install_gatk.sh` |
| DeepImmuno | 免疫原性辅助 | 可选 | `bash scripts/install_deepimmuno.sh` |
| PRIME / MixMHCpred / BigMHC | 免疫原性辅助 | 可选 | `bash scripts/install_immunogenicity_tools.sh` |
| FACETS | purity/CNV/LOH | CCF/APPM 推荐 | `bash scripts/install_facets.sh` |
| ASCAT / PyClone-VI | CNV/克隆性辅助 | 可选 | `bash scripts/install_ascat_pyclone.sh` |
| LOHHLA | HLA LOH | vaccine/TCR 场景推荐 | `bash scripts/install_lohhla.sh` 后配置 Polysolver/Novoalign |
| STAR-Fusion / FusionCatcher / EasyFuse | fusion detection | Fusion 流程需要 | 测试包含 wrapper；真实运行需参考库 |
| Arriba | fusion cross-check | 可选 | 可用 conda 或官方安装；测试中未自动安装 |

---

## 4. 关键工具安装细节与测试中发现的修复点

### 4.1 MHCflurry 模型下载和 C++ runtime

测试中 `mhcflurry-downloads fetch` 曾因 `libstdc++.so.6: CXXABI_1.3.15 not found` 失败。处理方式：

```bash
conda activate neoag-tools
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
conda install -n neoag-tools -c conda-forge -y 'libstdcxx-ng>=13'
mhcflurry-downloads fetch
```

若下载仍失败，可稍后手动重试。`neoag-v03 check-tools` 只能确认 `mhcflurry-predict` 是否存在，不能完全证明模型文件已下载。建议额外运行：

```bash
mhcflurry-predict --help | head
```

### 4.2 NetMHCpan 4.2

NetMHCpan 需要单独从 DTU 获取 Linux tarball。流程：

```bash
# 1. 从 DTU 邮件链接下载，例如 netMHCpan-4.2c.Linux.tar.gz
mkdir -p vendor
cp /path/to/netMHCpan-4.2c.Linux.tar.gz vendor/

# 2. 安装
export NEOAG_CONDA_BASE="$(conda info --base)"
bash scripts/install_netmhcpan.sh vendor/netMHCpan-4.2c.Linux.tar.gz

# 3. 检查
source conf/tools.env.sh
neoag-v03 check-tools
netMHCpan -h | head
```

迁移测试中旧脚本曾硬编码 `/home/na/miniforge3`，导致 wrapper 找不到 conda sysroot。本版脚本默认使用 `conda info --base` 或 `NEOAG_CONDA_BASE`，不要再手动改成固定路径。

如已安装但 wrapper 失效：

```bash
bash scripts/install_netmhcpan.sh --repair
```

### 4.3 VEP 和 cache

安装 VEP：

```bash
bash scripts/install_vep.sh
source conf/tools.env.sh
neoag-v03 check-tools
```

安装 cache 可能很慢，可二选一：

```bash
# 在线下载 cache，耗时和空间较大
bash scripts/install_vep_cache.sh

# 或使用已有 cache
export NEOAG_VEP_CACHE=/path/to/vep_cache
export NEOAG_VEP_CACHE_VERSION=105
source conf/tools.env.sh
```

测试中 VEP 已安装但 `check-tools` 显示 missing，原因是 VEP env 未加入 PATH。本版 `install_vep.sh` 会把 `NEOAG_VEP_BIN` 和 VEP env `bin` 写入 `conf/tools.env.sh`。

### 4.4 GATK

```bash
bash scripts/install_gatk.sh
source conf/tools.env.sh
neoag-v03 check-tools
gatk --help | head
```

迁移测试中 GATK 安装成功。

### 4.5 ASCAT / PyClone-VI

测试中旧脚本调用 `mamba`，用户执行 `pip install mamba` 后得到的是 Python 测试框架 `mamba`，不是 conda-forge 的 mamba solver，导致 `mamba env create -n ...` 参数错误。

本版脚本默认使用 `conda`，只有在 `NEOAG_USE_MAMBA=1` 且真正的 mamba 可用时才使用 mamba。

```bash
bash scripts/install_ascat_pyclone.sh
source conf/tools.env.sh
neoag-v03 check-tools
```

安装后会创建：

```text
bin/ascat.R
bin/pyclone
```

`ascat.R` 是 wrapper；真实 ASCAT 分析仍建议使用项目内或自定义 R 脚本调用 ASCAT 包。

### 4.6 FACETS

测试中 FACETS 没有安装入口，`check-tools` 显示 `runFACETS.R MISSING`。本版新增：

```bash
bash scripts/install_facets.sh
source conf/tools.env.sh
neoag-v03 check-tools
```

说明：FACETS 真实样本通常还需要 `snp-pileup`、dbSNP/common SNP VCF、purity/ploidy fitting 参数和样本 BAM；安装 wrapper 只证明 R 包和入口可用。

### 4.7 LOHHLA

测试中 LOHHLA 缺少安装命令。本版新增基本安装入口：

```bash
bash scripts/install_lohhla.sh
source conf/tools.env.sh
neoag-v03 check-tools
```

生产运行还必须配置：

```bash
export POLYSOLVER_HOME=/path/to/polysolver
export NOVOALIGN_LICENSE_FILE=/path/to/novoalign.lic
```

LOHHLA 还需要 HLA calls、HLA FASTA、tumor/normal BAM、purity/ploidy 等输入。`check-tools` OK 只代表脚本入口存在，不代表所有运行资源齐全。

### 4.8 PRIME / MixMHCpred / BigMHC

```bash
bash scripts/install_immunogenicity_tools.sh
source conf/tools.env.sh
neoag-v03 check-tools
```

迁移测试中旧脚本的问题包括：

- `NEOAG_PRIME_BIN`、`MIXMHCPRED_BIN`、`BIGMHC_DIR` 未赋值；
- MixMHCpred 缺 `numpy`；
- BigMHC 缺 `torch`、`pandas`、`psutil`；
- PRIME 编译成 `PRIME.x.bin`，而 PRIME wrapper 实际调用 `PRIME.x`。

本版脚本已改为：

- 编译 `tools/prime/lib/PRIME.x`；
- 自动安装 `numpy pandas psutil torch`；
- 创建 `bin/bigmhc_predict` wrapper；
- 将 PRIME/MixMHCpred/BigMHC 路径写入 `conf/tools.env.sh`。

如果 BigMHC clone 因网络中断失败，直接重跑同一命令即可；BigMHC 仓库较大，慢速网络下可能需要多次尝试或预先离线拷贝。

### 4.9 NetMHCstabpan

测试中 IEDB shim 路径可安装：

```bash
bash scripts/install_netmhcstabpan.sh --iedb
source conf/tools.env.sh
neoag-v03 check-tools
```

### 4.10 EasyFuse / Arriba / STAR-Fusion / FusionCatcher

测试中 `seed_easyfuse_conda_envs.sh` 因默认读取 `work/.nextflow_conda` 失败。该脚本适合“已有 Nextflow conda cache”场景，不是全新安装脚本。

建议：

```bash
# 若已有 EasyFuse/STAR-Fusion/FusionCatcher wrappers，直接检查
source conf/tools.env.sh
neoag-v03 check-tools

# 若无 work/.nextflow_conda，不要把 seed_easyfuse_conda_envs.sh 当作全新安装命令。
# 先通过 Nextflow fixture 或站点工具包准备 conda env / reference，再 seed。
```

Arriba 在测试中仍为 missing，建议用 conda 或官方二进制单独安装，并把 `arriba` 放入 PATH。

---

## 5. 推荐安装检查顺序

最小开发/fixture：

```bash
python -m pip install -e '.[test]'
pytest -q
neoag-v03 run-demo --outdir work/demo_v043 --sample-id DEMO001
```

常用真实样本 scoring：

```bash
bash scripts/setup_tools_env.sh
bash scripts/install_vep.sh
bash scripts/install_gatk.sh
bash scripts/install_netmhcstabpan.sh --iedb
bash scripts/install_deepimmuno.sh
source conf/tools.env.sh
neoag-v03 check-tools
```

完整增强证据层：

```bash
bash scripts/install_netmhcpan.sh vendor/netMHCpan-4.2c.Linux.tar.gz
bash scripts/install_facets.sh
bash scripts/install_ascat_pyclone.sh
bash scripts/install_lohhla.sh
bash scripts/install_immunogenicity_tools.sh
source conf/tools.env.sh
neoag-v03 check-tools
```

期望 `check-tools` 至少显示：

```text
pvacseq OK
pvacfuse OK
pvacsplice OK
netmhcpan OK   # 若安装了 DTU tarball
mhcflurry OK
netmhcstabpan OK
vep OK
gatk OK
pyclone OK
```

可选工具可能仍为 `MISSING`，取决于是否需要对应流程：

```text
lohhla      HLA LOH，需要时安装和配置
facets      purity/CNV，需要时安装
ascat       备用 CNV/LOH，需要时安装
prime       免疫原性辅助，需要时安装
bigmhc_im   免疫原性辅助，需要时安装
arriba      fusion cross-check，需要时安装
```

---

## 6. 常见错误与处理

| 错误 | 原因 | 处理 |
|---|---|---|
| `bin/neoag-nextflow: Permission denied` | bin 权限丢失 | `find bin -maxdepth 1 -type f -exec chmod +x {} \;` |
| `CXXABI_1.3.15 not found` | MHCflurry/Scipy 加载系统 libstdc++ | `conda install -n neoag-tools -c conda-forge 'libstdcxx-ng>=13'` 并设置 `LD_LIBRARY_PATH=$CONDA_PREFIX/lib` |
| `mamba: unrecognized arguments -n ...` | 安装了 pip 版 mamba | 不要 `pip install mamba`；使用本版脚本默认 conda，或安装 conda-forge mamba |
| `vep MISSING` 但 VEP 已安装 | VEP env 未进 PATH | 运行本版 `scripts/install_vep.sh` 后 `source conf/tools.env.sh` |
| `NetMHCpan wrapper 找 /home/na/miniforge3` | 旧脚本硬编码 conda base | 本版脚本使用 `NEOAG_CONDA_BASE`/`conda info --base`；重新运行 install 或 `--repair` |
| `PRIME.x Syntax error` | PRIME 编译目标不对 | 本版脚本编译 `lib/PRIME.x` |
| BigMHC 缺 `torch/pandas/psutil` | Python 依赖缺失 | 本版脚本自动安装；也可手动 `python -m pip install torch pandas psutil` |
| `LOHHLA MISSING` | 未安装/未配置 | `bash scripts/install_lohhla.sh`，并配置 Polysolver/Novoalign |
| `FACETS MISSING` | 未安装/无 wrapper | `bash scripts/install_facets.sh` |

---

## 7. 生产运行注意事项

1. 不要在 production 模式使用 stub 结果；设置 `NEOAG_STRICT_MODE=1`。
2. NetMHCpan、LOHHLA、FACETS、VEP cache 等大型工具/参考库不随轻量包分发。
3. 患者 BAM/FASTQ/VCF 不应写入公开 release 包。
4. 每次正式运行前保存：`neoag-v03 check-tools` 输出、工具版本、reference manifest、运行配置和 provenance。
5. 如果某个工具在 `check-tools` 中显示 `MISSING`，对应证据层应标记为 missing/unassessed，而不是解释为阴性结果。

---

## 8. 常用命令

### 从 pVAC-like 表运行

```bash
neoag-v03 run-v03 \
  --outdir results/sample \
  --sample-id SAMPLE001 \
  --profile default \
  --pvac data/fixtures/pvacseq_aggregated.tsv \
  --immunogenicity-stub
```

### 从 raw intermediates 运行

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

### 生成 report

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

---

## 9. 更多文档

- `docs/TOOLS_SETUP.md`：外部工具安装详解。
- `docs/V043_CCF21.md`：CCF 2.1 说明。
- `docs/V042_P1_APPM_EXPLAINABILITY.md`：APPM explainability 说明。
- `docs/V04_EVIDENCE_SAFETY_ESCAPE.md`：Safety / Immune Escape 说明。
- `RELEASE.md`：发布边界和测试说明。
