# 部署功能修复记录

日期:2026-07-06

## 背景

在 `test/run_demo.sh` 中运行 `neoag-v03 run-full` 完整管线，逐一解决工具发现、运行时依赖和数据缺失问题，最终管线成功完成。

---

## 修复内容

### 1. VEP 未找到

**现象**: `vep (vep) not available: not found on PATH`

**原因**: `run_demo.sh` 未 source `conf/tools.env.sh`，`NEOAG_VEP_BIN` 未设置，`vep` 不在 PATH 中。

**修复**: `run_demo.sh` 中添加 `source /home/na/project/neo/conf/tools.env.sh`。

### 2. VEP Perl DBI.pm 缺失

**现象**: `Can't locate DBI.pm in @INC`

**原因**: `tools.env.sh` 中 `neoag-tools/bin` 排在 `neoag-vep/bin` 之前，`#!/usr/bin/env perl` 解析到无 `perl-dbi` 包的 neoag-tools Perl。

**修复**: `run_demo.sh` 中 source 后重新 prepend `neoag-vep/bin` 到 PATH：
```bash
export PATH="/home/na/miniforge3/envs/neoag-vep/bin:$PATH"
```

### 3. VEP `--online` 参数不存在

**影响文件**: `src/neoag_v03/vep/annotate.py:74`, `src/neoag_v03/tools/runner.py:739`

**原因**: VEP 115 不支持 `--online` 标志。在线模式需使用 `--database`（或不传 `--cache`/`--offline`）。

**修复**: 将 `--online` 改为 `--database --species homo_sapiens`，同时支持 `NEOAG_VEP_ONLINE` 环境变量：
```python
if online is True or (online is None and os.environ.get("NEOAG_VEP_ONLINE", "").lower() in {"1", "true", "yes"}):
    cmd.extend(["--database", "--species", "homo_sapiens"])
```

### 4. VEP 缓存路径

**影响文件**: `conf/tools.env.sh`

**原因**: `NEOAG_VEP_CACHE` 默认指向 `${NEOAG_TOOLS_ROOT}/data/vep`（不存在），quarantine 回退路径也不存在。

**修复**: 添加共享存储回退路径：
```bash
if [[ ! -d "${NEOAG_VEP_CACHE}/homo_sapiens" && -d "/mnt/zjl-bgi-zzb/peixunban/gl/data/reference/homo_sapiens" ]]; then
  export NEOAG_VEP_CACHE="/mnt/zjl-bgi-zzb/peixunban/gl/data/reference"
```

缓存版本为 105（`NEOAG_VEP_CACHE_VERSION=105`），结构为 `homo_sapiens/105_GRCh38/`。

### 5. VEP 自定义插件

**影响文件**: `work/vep_plugins/Wildtype.pm`, `work/vep_plugins/Frameshift.pm`（新文件）

**原因**: 项目管线依赖 `WildtypeProtein` 和 `FrameshiftSequence` 两个 CSQ 字段，需通过 VEP 插件生成。这两个插件不在 Ensembl 公共 VEP_plugins 仓库中。

**修复**: 在 `work/vep_plugins/` 中创建两个 VEP 115 兼容插件：
- `Wildtype.pm` — 从翻译序列提取突变位点附近的野生型蛋白窗口（`WildtypeProtein` 字段）
- `Frameshift.pm` — 为移码变异生成替代蛋白序列（`FrameshiftSequence` 字段）

`NEOAG_VEP_PLUGINS` 自动指向 `${NEOAG_TOOLS_ROOT}/work/vep_plugins`。

### 6. 无插件时的肽段提取回退

**影响文件**: `src/neoag_v03/vep/extract_peptides.py`

**原因**: `pick_csq_transcript` 要求 `WildtypeProtein` 或 `FrameshiftSequence` 非空才通过过滤；`build_mutant_protein` 在 `wt_protein` 为空时直接返回空。当 VEP 插件不可用时，所有变异被跳过。

**修复**:

- `pick_csq_transcript`: 允许 `Amino_acids` 包含 "/" 的变异通过（错义/框内 indel）
- `build_mutant_protein`: 当 `wt_protein` 为空但 `Amino_acids` 有置换信息时，构建合成蛋白骨架用于 k-mer 滑窗肽段生成

### 7. NetMHCpan 安装

**原因**: `netMHCpan` 可执行文件不存在，管线在 binding 预测阶段失败。

**修复**: 使用安装包 `netMHCpan-4.2c.Linux.tar.gz`，通过 `scripts/install_netmhcpan.sh` 安装到 `tools/netMHCpan/`。依赖 `neoag-tools` conda 环境中的 `patchelf` + `sysroot_linux-64`。

### 8. MHCflurry 模型下载与路径配置

**原因**: MHCflurry 首次运行需要下载 `models_class1_presentation` 模型（~135MB 压缩包，198MB 解压后）。

**修复**:
- 下载模型：`mhcflurry-downloads fetch models_class1_presentation`
- 移动模型到项目共享目录：`/home/na/project/tool_data/mhcflurry/`
- 在 `conf/tools.env.sh` 中配置环境变量避免重复下载：
  ```bash
  export MHCFLURRY_DOWNLOADS_DIR="${NEOAG_TOOLS_ROOT}/../tool_data/mhcflurry"
  ```

---

## 修改文件清单

| 文件 | 变更类型 |
|------|----------|
| `conf/tools.env.sh` | 修改 — VEP 缓存路径、MHCflurry 模型路径、NetMHCpan 自动追加 |
| `src/neoag_v03/vep/annotate.py` | 修改 — VEP `--online` → `--database` |
| `src/neoag_v03/tools/runner.py` | 修改 — 同上 |
| `src/neoag_v03/vep/extract_peptides.py` | 修改 — 无插件时肽段提取回退逻辑 |
| `work/vep_plugins/Wildtype.pm` | 新增 — VEP 野生型蛋白插件 |
| `work/vep_plugins/Frameshift.pm` | 新增 — VEP 移码蛋白插件 |
| `test/run_demo.sh` | 修改 — source tools.env.sh + PATH 重排序 |

## 新增工具/数据依赖

| 依赖 | 路径 | 说明 |
|------|------|------|
| NetMHCpan 4.2c | `tools/netMHCpan/` | DTU 学术许可，patchelf 适配 conda sysroot |
| MHCflurry 模型 | `../tool_data/mhcflurry/` | `MHCFLURRY_DOWNLOADS_DIR` 环境变量 |

## 管线验证结果

`neoag-v03 run-full` 完整运行成功，产出包括：
- 变异肽段提取（5930 变异 → 3852 肽段-HLA 对）
- NetMHCpan + MHCflurry binding 预测
- APPM 评分、CCF 克隆性分析、免疫逃逸分析
- 排序评分、验证计划、患者/技术 HTML 报告


# nextflow功能修复和优化

日期：0710

 ---
  优化总结

  根据 docs/workflow.md 报告的改进点，对 workflows/main_full.nf 及相关文件进行了以下优化：

  新增文件

  ┌───────────────────────┬─────────────────────────────┐
  │         文件          │            说明             │
  ├───────────────────────┼─────────────────────────────┤
  │ nextflow.config       │ 根级配置，自动加载，无需 -c │
  ├───────────────────────┼─────────────────────────────┤
  │ conf/main_full.config │ main_full.nf 配套配置       │
  └───────────────────────┴─────────────────────────────┘

  修改文件

  ┌──────────────────────────────┬───────────────────────────────────────────────────────────────────────────────────────────────────────┐
  │             文件             │                                               改进内容                                                │
  ├──────────────────────────────┼───────────────────────────────────────────────────────────────────────────────────────────────────────┤
  │ workflows/main_full.nf       │ 重写：输入验证、文档头、strict_mode、proper channel 模式、file() 包装                                 │
  ├──────────────────────────────┼───────────────────────────────────────────────────────────────────────────────────────────────────────┤
  │ workflows/neoag_v03_rc.nf    │ 弃用说明、CCF_2.out.ccf_lite 通道修复（3处）                                                          │
  ├──────────────────────────────┼───────────────────────────────────────────────────────────────────────────────────────────────────────┤
  │ modules/run_upstream/main.nf │ label 'large'、container 指令、time '24h'、输出 glob 修复                                             │
  ├──────────────────────────────┼───────────────────────────────────────────────────────────────────────────────────────────────────────┤
  │ 6 个关键模块                 │ 添加 label 和 container 指令（score_v041, ccf_2, appm_2, immune_escape, peptide_safety, report_v041） │
  ├──────────────────────────────┼───────────────────────────────────────────────────────────────────────────────────────────────────────┤
  │ 2 个旧版模块                 │ @Deprecated 标记（ccf_lite, appm_lite）                                                               │
  └──────────────────────────────┴───────────────────────────────────────────────────────────────────────────────────────────────────────┘

  对应报告改进点

  ┌────────────────────┬─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
  │     报告改进点     │                                                        实现                                                         │
  ├────────────────────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
  │ 1. 仅 local 执行器 │ ✅ nextflow.config 中定义 slurm、sge、awsbatch profiles，模块添加 label（small/medium/large/huge）用于 HPC 资源调度 │
  ├────────────────────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
  │ 2. 容器化缺失      │ ✅ 添加 docker、singularity、conda profiles；所有关键模块添加 container 指令                                        │
  ├────────────────────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
  │ 3. 主配置文件缺失  │ ✅ 创建根 nextflow.config（自动加载）+ conf/main_full.config                                                        │
  ├────────────────────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
  │ 4. CI 测试         │ ⬚ 超出范围（需要 GitHub Actions）                                                                                   │
  ├────────────────────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
  │ 5. 模块版本碎片化  │ ✅ ccf_lite/appm_lite 添加弃用标记；neoag_v03_rc.nf 统一使用 CCF_2.out.ccf_lite                                     │
  └────────────────────┴─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘

  验证结果

  [SUCCESS] completed=13 failed=0 cached=0

  完整流水线 13 个进程全部成功，产出文件：
  - ranked_events.v03.tsv / ranked_peptides.v03.tsv
  - evidence_report.v041.html
  - timeline.html / report.html / trace.txt（流水线执行报告）

