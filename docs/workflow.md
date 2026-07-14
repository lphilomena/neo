# NeoAg Nextflow 工作流

## 概述

项目使用 [Nextflow](https://www.nextflow.io/) (DSL2, v26.04.3) 作为工作流编排引擎，管理从原始测序数据（BAM/VCF）到新抗原（neoantigen）预测评分的完整计算流水线。

**启动器**：`bin/neoag-nextflow`，封装了 `bin/nextflow`，自动设置 `PYTHONPATH=src`、`PATH` 和 `NXF_HOME`。

```bash
# 设置 Nextflow 缓存目录（避免权限问题）
export NXF_HOME=/path/to/writable/nextflow_cache

# 运行 Demo
bin/neoag-nextflow run workflows/main.nf -c conf/demo.config

# 或直接使用 bin/nextflow
bin/nextflow run workflows/main.nf -c conf/demo.config

```

## 目录与文件总览

```
  项目共有 36 个 Nextflow 相关文件（1525 行代码），分布在四个层级：

  bin/nextflow              ← Nextflow 官方启动器 (v26.04.3, Apache 2.0 协议)
  bin/neoag-nextflow        ← 项目封装器：设置 PYTHONPATH/PATH/NXF_HOME

  conf/                     ← 6 个运行配置文件
  workflows/                ← 9 个工作流定义文件
  modules/                  ← 26 个模块定义文件

```

## 目录结构

```
bin/nextflow              ← Nextflow 官方启动器 (v26.04.3)
bin/neoag-nextflow        ← 项目封装器：设置环境变量 + 可写 NXF_HOME

conf/                     ← 运行配置文件（6 个）
  demo.config               - 快速演示（leukemia profile + 内置 fixture）
  local.config              - 本地执行器基础配置
  tools.config              - 上游工具路径配置
  sv_demo.config            - SV WGS fixture 演示
  sv_wes_demo.config        - SV WES fixture 演示
  snv_wes_demo.config       - SNV WES fixture 演示

workflows/                ← 工作流定义文件（10 个）
  main.nf                   - 核心评分流程（从 pVAC TSV 输入）
  main_all.nf               - 完整端到端（BAM→HLA Typing→Mutect2→Upstream→评分）
  main_fromVCF.nf           - VCF 入口端到端（已知 VCF+HLA→Upstream→评分→报告）
  neoag_v03_rc.nf           - 核心评分链子流程（被 include）
  snv_phase1_wes.nf         - WES SNV Phase1（BAM → Mutect2 → 评分）
  snv_phase1_wes_fixture.nf - WES SNV fixture（已有 VCF，跳过变异检出）
  sv_phase1_wgs.nf          - WGS SV Phase1（BAM → Manta/SvABA/GRIDSS → 评分）
  sv_phase1_fixture.nf      - SV fixture（已有 VCF，跳过 BAM 调用器）
  sv_phase1_5_wes.nf        - WES SV Phase1.5（捕获感知分级）
  sv_score_v03.nf           - SV 评分链子流程（被 include）

modules/                  ← 模块定义文件（26 个）
  parse_pvac/               - 解析 pVACtools 输出 → raw_events + raw_peptides
  parse_netmhcpan/          - 解析 NetMHCpan 输出
  parse_mhcflurry/          - 解析 MHCflurry 输出
  run_upstream/             - 运行上游工具（pVAC/NetMHCpan/MHCflurry/VEP）
  run_binding_predictors/   - 运行 NetMHCpan + MHCflurry 结合预测
  build_presentation/       - 构建 HLA 呈递证据
  appm_2/                   - 构建 APPM 2.0（抗原呈递加工机制评估）
  appm_lite/                - APPM lite（兼容入口）
  ccf_2/                    - 构建 CCF 2.0（癌细胞分数/克隆性）
  ccf_lite/                 - CCF lite（兼容入口）
  peptide_safety/           - 肽段安全性门控（参考蛋白组/正常配体组）
  immune_escape/            - 免疫逃逸评估（HLA LOH/APM/JAK/B2M）
  score_v03/                - v0.3 评分（基础版）
  score_v041/               - v0.4.1 评分（完整版，含 safety + escape）
  report_v03/               - v0.3 HTML 报告（患者+技术双受众）
  report_v041/              - v0.4.1 HTML 报告（含 APPM/escape/safety/CCF）
  validation_plan/          - 验证计划生成
  workflow_provenance/      - 版本追溯聚合
  snv_write_run_config/     - 生成 SNV WES 的 run config TOML
  sv_build_raw/             - SV VCF → raw_events + raw_peptides
  sv_manta/                 - 运行 Manta SV 调用器
  sv_svaba/                 - 运行 SvABA SV 调用器
  sv_gridss/                - 运行 GRIDSS SV 调用器
  sv_normalize_merge/       - 合并多个 SV 调用器的 VCF
  gatk_mutect2/             - 运行 GATK Mutect2 SNV 检出
  gatk_filter_mutect_calls/ - 过滤 Mutect2 原始 VCF
```

## 核心流程

### 1. SNV/InDel 评分流程（`main.nf`）

从 pVACtools TSV 文件直接启动，使用 fixture 数据快速验证：

```bash
bin/neoag-nextflow run workflows/main.nf -c conf/demo.config
```

```
pVAC TSV 文件
    │
    ▼
PARSE_PVAC ───→ raw_events.tsv + raw_peptides.tsv
    │
    ├──→ PARSE_NETMHCPAN ──→ netmhcpan_evidence.tsv
    ├──→ PARSE_MHCFLURRY  ──→ mhcflurry_evidence.tsv
    │         │
    │         ▼
    ├──→ BUILD_PRESENTATION ──→ presentation_evidence.tsv
    │
    ├──→ APPM_2 ──→ appm_summary.tsv + gene/pathway/modifier 侧文件
    ├──→ CCF_2  ──→ ccf_2.tsv（克隆性评估）
    ├──→ PEPTIDE_SAFETY ──→ peptide_safety.tsv
    ├──→ IMMUNE_ESCAPE ──→ immune_escape_summary.tsv + peptide_escape_flags.tsv
    │
    ▼
SCORE_V041 ──→ ranked_events.v03.tsv + ranked_peptides.v03.tsv
    │
    ├──→ VALIDATION_PLAN ──→ validation_plan.v03.tsv
    └──→ REPORT_V041 ──→ evidence_report.v041.html
```

### 2. WES SNV 完整流程（`snv_phase1_wes.nf`）

从 tumor/normal BAM 开始的全自动流程：

```bash
bin/neoag-nextflow run workflows/snv_phase1_wes.nf \
  --sample_id P001 \
  --tumor_bam tumor.bam --normal_bam normal.bam \
  --reference_fasta GRCh38.fa --intervals_bed capture.bed \
  --tumor_sample_name TUMOR --normal_sample_name NORMAL \
  -c conf/snv_wes_demo.config
```

```
tumor/normal BAM
    │
    ▼
GATK_MUTECT2 ──→ raw VCF
    │
    ▼
GATK_FILTER_MUTECT_CALLS ──→ filtered VCF
    │
    ▼
SNV_WRITE_RUN_CONFIG ──→ run.snv_wes.generated.toml
    │
    ▼
RUN_UPSTREAM ──→ pVAC/NetMHCpan/MHCflurry/VEP 输出
    │
    ▼
NEOAG_V03_RC（评分链，同上）
```

**Fixture 模式**（`snv_phase1_wes_fixture.nf`）：跳过 Mutect2 调用，直接从已有 somatic VCF 启动。

### 3. SV WGS 流程（`sv_phase1_wgs.nf`）

三个 SV 调用器并行运行 → 合并 → 评分：

```bash
bin/neoag-nextflow run workflows/sv_phase1_wgs.nf \
  --sample_id P001 \
  --tumor_bam tumor.bam --normal_bam normal.bam \
  --reference_fasta GRCh38.fa --gencode_gtf gencode.gtf \
  --hla 'HLA-A*02:01,HLA-B*07:02' \
  -c conf/sv_demo.config
```

```
tumor/normal BAM
    │
    ├──→ SV_MANTA  ──→ somaticSV.vcf.gz    ┐
    ├──→ SV_SVABA  ──→ somatic.sv.vcf      ├─ 并行
    └──→ SV_GRIDSS ──→ gridss.vcf.gz       ┘
    │
    ▼
SV_NORMALIZE_MERGE ──→ sv_inputs.list
    │
    ▼
SV_BUILD_RAW ──→ raw_events.tsv + raw_peptides.tsv
    │
    ▼
NEOAG_SV_SCORE_V03
  ├──→ RUN_BINDING_PREDICTORS（NetMHCpan + MHCflurry）
  ├──→ BUILD_PRESENTATION
  ├──→ APPM_2 / CCF_2 / PEPTIDE_SAFETY / IMMUNE_ESCAPE
  ├──→ SCORE_V041
  └──→ REPORT_V041
```

### 4. SV WES Phase 1.5 流程（`sv_phase1_5_wes.nf`）

针对外显子组捕获的捕获感知分析，使用 `sv_wes_phase1_5` profile：

```bash
bin/neoag-nextflow run workflows/sv_phase1_5_wes.nf -c conf/sv_wes_demo.config
```

**与 WGS 版本的关键差异**：
- 使用 `sv-build-raw-wes`（捕获感知）替代 `sv-build-raw`
- 引入 WES Tier 分级：WES_Tier1 / WES_Tier2 / WES_Tier3
- 支持 `--capture-bed` 参数进行捕获区间过滤
- 对 RNA 连接 reads 证据有更严格的要求

## 配置文件

### conf/*.config

| 文件 | 用途 | 默认 profile | 适用工作流 |
|---|---|---|---|
| `demo.config` | 快速演示 | `leukemia` | `main.nf` |
| `local.config` | 本地执行器（最小配置） | — | 任意 |
| `tools.config` | 上游工具路径 | — | `main_fromVCF.nf` |
| `sv_demo.config` | SV WGS fixture 演示 | `sv_wgs_phase1` | `sv_phase1_wgs.nf`, `sv_phase1_fixture.nf` |
| `sv_wes_demo.config` | SV WES fixture 演示 | `sv_wes_phase1_5` | `sv_phase1_5_wes.nf` |
| `snv_wes_demo.config` | SNV WES fixture 演示 | `default` | `snv_phase1_wes.nf`, `snv_phase1_wes_fixture.nf` |

所有 demo/fixture config 均使用 `executor = 'local'` 和 `errorStrategy = 'terminate'`。

### profiles/*.toml

Profile 文件控制评分参数（权重、阈值、策略），与工作流正交：

| Profile | 癌种/场景 | 关键特点 |
|---|---|---|
| `default.toml` | 通用实体瘤 | 全突变来源，三源免疫原性 |
| `leukemia.toml` | 白血病 | 去 CNV/SV，收紧 HSPC 阈值 |
| `sarcoma.toml` | 肉瘤 | Fusion 权重最高 |
| `sv_wgs_phase1.toml` | SV WGS | 仅 SV 子类型，关闭免疫原性 |
| `sv_wes_phase1_5.toml` | SV WES | 捕获感知分级，RNA 支持加权 |
| `immunogenicity_extended.toml` | 扩展免疫 | 三工具免疫原性（含 DeepImmuno） |

Profile 通过 `--profile_name` 参数从 workflow 层透传到 Python 模块层，在模块内部调用 `load_profile()` 加载。

## 模块设计模式

所有模块遵循统一的 Nextflow DSL2 模式：

```nextflow
// modules/<name>/main.nf
process MODULE_NAME {
    tag "$sample_id"
    publishDir params.outdir, mode: 'copy'

    input:
    val sample_id
    path some_input_file
    val profile_name

    output:
    path "output.tsv", emit: output
    path "versions.yml", emit: versions

    script:
    """
    neoag-v03 <subcommand> \
      --sample-id '${sample_id}' \
      --profile '${profile_name}' \
      --input '${some_input_file}' \
      --out output.tsv
    """
}
```

**关键设计原则**：
1. **薄封装**：每个 module 是 `neoag-v03` 单个子命令的 Nextflow 封装，Python 业务逻辑不在 `.nf` 文件中
2. **profile 透传**：`profile_name` 作为 val 类型在所有模块间传递，不做解析
3. **版本追溯**：每个 module 输出 `versions.yml`，由 `workflow_provenance` 模块汇总
4. **空数据 fallback**：可选输入文件默认指向 `assets/empty_*.tsv`，缺失数据不会阻塞流程

## 与 profiles/ 的配合关系

```
用户指定 --profile_name leukemia
         │
         ▼
Nextflow workflow（编排层）
    params.profile_name = 'leukemia'
    → 透传给所有 module 的 script 块
         │
         ▼
neoag-v03 <subcommand> --profile 'leukemia'
         │
         ▼
Python load_profile("leukemia")
    → profiles/leukemia.toml
    → 返回评分权重/阈值/策略 dict
```

Nextflow 不解析 TOML，只负责参数透传。Profile 切换不影响流程结构。

## 运行模式

### Stub 模式

用于快速验证流程正确性，使用预置 fixture 数据替代真实工具输出：

| 参数 | 作用 |
|---|---|
| `--binding_stub` | 使用 fixture 的 NetMHCpan/MHCflurry 输出 |
| `--upstream_stub` | 使用 fixture 的 pVAC 输出 |
| `--immunogenicity_stub` | 使用占位免疫原性评分 |

```bash
# SV fixture 快速测试（使用预构建 VCF + stub 结合预测）
bin/neoag-nextflow run workflows/sv_phase1_fixture.nf -c conf/sv_demo.config

# SNV fixture 快速测试
bin/neoag-nextflow run workflows/snv_phase1_wes_fixture.nf -c conf/snv_wes_demo.config
```

### 严格模式

设置 `--strict_mode true` 可禁止 stub 模式，适用于生产环境：

```bash
bin/neoag-nextflow run workflows/sv_phase1_fixture.nf -c conf/sv_demo.config --strict_mode true
```

## 空数据 Fallback 机制

可选输入文件默认指向 `assets/empty_*.tsv`（仅含表头、无数据行的占位文件）：

```nextflow
// 示例：workflows/main.nf
params.vep_appm = params.vep_appm ?: "${projectDir}/../assets/empty_vep.tsv"
params.hla_loh  = params.hla_loh  ?: "${projectDir}/../assets/empty_hla_loh.tsv"
params.purity   = params.purity   ?: "${projectDir}/../assets/empty_purity.tsv"
```

这样即使不提供某些可选数据，下游也能正常读取文件（格式正确），只是无匹配结果。

## 当前局限性

| 方面 | 状态 | 说明 |
|---|---|---|
| **执行器** | 仅 local | 所有 config 使用 `executor = 'local'`，不支持 HPC（SLURM/SGE）或云平台 |
| **容器化** | 缺失 | 无 Docker/Singularity 镜像定义，工具依赖系统 PATH |
| **主配置文件** | 缺失 | 根目录无 `nextflow.config`，需每次通过 `-c` 指定 |
| **CI 测试** | 缺失 | 无 GitHub Actions 或其他 CI 验证流水线可运行 |
| **模块版本** | 碎片化 | 同时存在 `ccf_lite`/`ccf_2`、`appm_lite`/`appm_2` 等新旧版本 |

## 相关文档

- [SV Phase 1 WGS](SV_PHASE1_WGS.md) — SV WGS 工作流详细说明
- [SV Phase 1.5 WES](SV_PHASE1_5_WES.md) — SV WES 捕获感知工作流
- [SNV Phase 1 WES](SNV_PHASE1_WES.md) — SNV WES 工作流详细说明
- [安装与数据准备](INSTALL_AND_DATA.md) — NXF_HOME 配置与运行示例
- [线上发布说明](ONLINE_RELEASE.md) — Nextflow 缓存与离线部署
- [Tools Setup](TOOLS_SETUP.md) — 上游工具安装
