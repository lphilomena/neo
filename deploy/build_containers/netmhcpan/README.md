# NetMHCpan 4.2c 容器镜像

## 镜像定位

**纯净工具镜像** — 只包含 NetMHCpan 4.2c 本身，不包含 `neoag-v03` CLI。

```
neoag-netmhcpan:4.2c-ubuntu22.04
    ├── /opt/netMHCpan/netMHCpan      ← bash wrapper（替换了原始 tcsh 脚本）
    ├── /opt/netMHCpan/Linux_x86_64/  ← ELF 二进制 + bin/
    └── /opt/netMHCpan/data/          ← MHC 等位基因数据
```

## 构建

```bash
# 1. 准备 tarball（需要 DTU 学术许可）
cp /path/to/netMHCpan-4.2c.Linux.tar.gz vendor/
# 2. 构建
bash deploy/build_containers/netmhcpan/build.sh
# 3. 验证
bash deploy/build_containers/netmhcpan/test.sh
```

## 直接使用（独立于管线）

```bash
# FASTA 输入 + 指定等位基因
docker run --rm -v $(pwd):/data neoag-netmhcpan:4.2c-ubuntu22.04 \
  netMHCpan -f /data/peptides.fa -a HLA-A02:01

# PMHC 格式（肽段 + MHC 在同一文件）
docker run --rm -v $(pwd):/data neoag-netmhcpan:4.2c-ubuntu22.04 \
  netMHCpan -pmhc -BA -f /data/input.pmhc
```

## 在 Nextflow 管线中使用

### 当前状态

管线通过 `neoag-v03 run-tool netmhcpan` 调用工具（见 `modules/run_binding_predictors/main.nf`）。由于本镜像是纯净工具镜像（不含 `neoag-v03`），**不能直接替换为管线默认容器**。

现有两种运行模式：

| 模式 | 工具提供方式 | neoag-v03 来源 |
|------|------------|---------------|
| **conda**（当前可用） | conda 环境安装到宿主机 | conda 环境 |
| **docker**（待完成） | 各工具独立镜像 | 需要额外提供（见下方） |

### 未来集成方案

要在 Docker 模式下使用本镜像，需要解决 `neoag-v03` 的来源问题：

**方案 A：bind-mount（开发/测试）**
```bash
NEOAG_RUNNER_MODE=docker \
neoag-nextflow run workflows/main_all.nf \
  -profile docker \
  --with-process-container RUN_BINDING_PREDICTORS=neoag-netmhcpan:4.2c-ubuntu22.04
```
需确保 `neoag-v03` 在 Nextflow 进程可访问的路径（如宿主 PATH 经 bind mount）。

**方案 B：组合镜像（生产）**
```dockerfile
# 将多个工具镜像 + neoag-v03 打包成一个镜像
FROM neoag-netmhcpan:4.2c-ubuntu22.04 AS netmhcpan
FROM neoag-mhcflurry:xxx AS mhcflurry
FROM python:3.12-slim
COPY --from=netmhcpan /opt/netMHCpan /opt/netMHCpan
COPY --from=mhcflurry /opt/mhcflurry /opt/mhcflurry
RUN pip install neoag-v03
```
构建后推送到 `ghcr.io/org/neoag-tools:latest`（即 `nextflow.config` 中引用的镜像）。

## 镜像标签

| 标签 | 用途 |
|------|------|
| `neoag-netmhcpan:4.2c-ubuntu22.04` | 本镜像（纯 netMHCpan） |
| `ghcr.io/org/neoag-tools:latest` | 组合镜像（管线 Docker 模式目标，尚未构建） |

## Apptainer (HPC)

```bash
apptainer build netmhcpan-4.2c.sif netmhcpan-4.2c.apptainer.def
apptainer run netmhcpan-4.2c.sif netMHCpan -f input.fa -a HLA-A02:01
```

## 许可

NetMHCpan 4.2c 需要 DTU 学术许可。本镜像不包含软件本体 — tarball 需从 DTU 下载，构建时通过 `vendor/` 目录注入。
