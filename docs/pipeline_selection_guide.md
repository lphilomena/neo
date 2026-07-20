# Pipeline Selection Guide — NeoAg Nextflow Workflows

> 本文档帮助用户和 AI 助手根据输入数据情况，从 5 个核心 Nextflow 工作流中选择正确的一个。

## 总览

```
                        ┌──────────────────────────┐
                        │   你有什么输入数据？       │
                        └────────────┬─────────────┘
                                     │
              ┌──────────────────────┼──────────────────────┐
              │                      │                      │
              ▼                      ▼                      ▼
      ┌──────────────┐     ┌──────────────┐       ┌──────────────┐
      │ 只有 TOML    │     │ TOML +       │       │ tumor BAM    │
      │ 配置文件     │     │ BAM/FASTQ    │       │ + normal BAM │
      └──────┬───────┘     └──────┬───────┘       └──────┬───────┘
             │                    │                      │
             ▼                    ▼              ┌───────┴───────┐
    main_fromVCF.nf    main_fromVCF_nohla.nf    │ HLA 来源？     │
                                               └───────┬───────┘
                                           ┌───────────┼───────────┐
                                           │           │           │
                                           ▼           ▼           ▼
                                      OptiType     OptiType     手动提供
                                      无 QC        有 QC
                                           │           │           │
                                           ▼           ▼           ▼
                                     main_all   main_all_qc  main_all_nohla
```

## 工作流详情

### 1. `main_fromVCF` — TOML 配置入口（已有 VCF + HLA）

**适用场景**：已经完成了变异检出（VCF 已生成），HLA 分型结果已知，都写在 TOML 配置文件中。只需运行上游工具 + 评分。

```
输入: TOML 配置 (含 VCF 路径 + HLA alleles + 样本信息)
输出: 新抗原排序报告
不执行: HLA 分型、变异检出
```

| 参数 | 必需 | 说明 |
|------|------|------|
| `--run_config` | ✅ | TOML 配置文件路径 |
| `--sample_id` | ✅ | 样本标识 |
| `--outdir` | | 输出目录 (默认: results/full) |

### 2. `main_fromVCF_nohla` — TOML 配置 + 自动 HLA 分型

**适用场景**：已有 TOML 配置和 VCF，但需要从测序数据中自动做 HLA 分型。OptiType 会自动将 HLA 结果注入 TOML。

```
输入: TOML 配置 + BAM/FASTQ (用于 HLA 分型)
输出: 新抗原排序报告
执行: OptiType HLA 分型
不执行: 变异检出
```

| 参数 | 必需 | 说明 |
|------|------|------|
| `--run_config` | ✅ | TOML 配置文件路径 |
| `--input_bam` | ✅* | 用于 HLA 分型的 BAM |
| `--input_fq1/2` | ✅* | 或 FASTQ pair（替代 BAM） |
| `--sample_id` | ✅ | 样本标识 |
| `--skip_hla_typing` | | 跳过 HLA 分型（使用 TOML 中已有的） |

### 3. `main_all` — BAM 端到端（自动 HLA + 变异检出 + 评分）

**适用场景**：只有 tumor/normal BAM 文件，需要从头开始：HLA 分型 → 变异检出 → 上游工具 → 新抗原评分。

```
输入: tumor BAM + normal BAM + 参考基因组
输出: 新抗原排序报告 + HLA 分型结果
执行: OptiType + Mutect2 + 上游工具 + 评分
```

| 参数 | 必需 | 说明 |
|------|------|------|
| `--normal_bam` | ✅ | 正常/血液 BAM |
| `--tumor_bam` | ✅ | 肿瘤 BAM |
| `--sample_id` | ✅ | 样本标识 |
| `--reference_fasta` | ✅ | 参考基因组 FASTA |
| `--tumor_sample_name` | | BAM 中肿瘤样本名 (默认: TUMOR) |
| `--normal_sample_name` | | BAM 中正常样本名 (默认: NORMAL) |

### 4. `main_all_qc` — BAM 端到端 + 质量控制

**适用场景**：与 `main_all` 相同，但额外需要肿瘤纯度/倍性（FACETS + PURPLE）和 HLA LOH（LOHHLA + SpecHLA）的 QC 分析。

> 这是**最全面**的工作流，推荐用于生产环境。

```
输入: 同 main_all + (可选) dbSNP VCF
输出: 新抗原排序报告 + HLA 分型 + QC 报告 (LOH + 纯度)
执行: OptiType + Mutect2 + (并行: LOH_CHECK + PURITY_CHECK) + 上游 + 评分
```

| 参数 | 必需 | 说明 |
|------|------|------|
| (同 main_all) | ✅ | |
| `--dbsnp_vcf` | 推荐 | dbSNP VCF，用于 FACETS snp-pileup |
| `--skip_qc` | | 跳过 QC，仅运行主流程 |

### 5. `main_all_nohla` — BAM 端到端（手动 HLA）

**适用场景**：有 BAM 文件，且 HLA 分型结果已经通过其他方式获得（或来自之前的运行），无需 OptiType。

```
输入: tumor BAM + normal BAM + HLA alleles + 参考基因组
输出: 新抗原排序报告
执行: Mutect2 + 上游工具 + 评分
不执行: HLA 分型
```

| 参数 | 必需 | 说明 |
|------|------|------|
| (同 main_all) | ✅ | |
| `--hla_alleles` | ✅ | 逗号分隔，如 `"HLA-A*02:01,HLA-B*07:02,HLA-C*07:02"` |

## 快速决策表

| 你的情况 | 使用工作流 | 需要额外提供 |
|---------|-----------|-------------|
| 已有 TOML + VCF + HLA | `main_fromVCF` | 无 |
| 已有 TOML + VCF，需要自动 HLA | `main_fromVCF_nohla` | BAM 或 FASTQ |
| 只有 BAM，需要全部自动，无 QC | `main_all` | 参考基因组 |
| 只有 BAM，需要全部自动 + QC | `main_all_qc` | 参考基因组 + (推荐) dbSNP VCF |
| 只有 BAM，已知 HLA | `main_all_nohla` | 参考基因组 + HLA alleles |

## 命令速查

### 使用统一启动脚本 (推荐)

```bash
# TOML 入口
bash scripts/run_pipeline.sh --workflow main_fromVCF \
  --run_config conf/run.mycase.toml --sample_id SAMPLE001

# TOML + 自动 HLA
bash scripts/run_pipeline.sh --workflow main_fromVCF_nohla \
  --run_config conf/run.mycase.toml --input_bam sample.bam --sample_id SAMPLE001

# BAM 端到端
bash scripts/run_pipeline.sh --workflow main_all \
  --normal_bam normal.bam --tumor_bam tumor.bam \
  --sample_id SAMPLE001 --reference_fasta ref.fa

# BAM 端到端 + QC
bash scripts/run_pipeline.sh --workflow main_all_qc \
  --normal_bam normal.bam --tumor_bam tumor.bam \
  --sample_id SAMPLE001 --reference_fasta ref.fa \
  --dbsnp_vcf dbsnp.vcf.gz

# BAM 端到端 + 手动 HLA
bash scripts/run_pipeline.sh --workflow main_all_nohla \
  --normal_bam normal.bam --tumor_bam tumor.bam \
  --hla_alleles "HLA-A*02:01,HLA-B*07:02,HLA-C*07:02" \
  --sample_id SAMPLE001 --reference_fasta ref.fa
```

### 使用 Nextflow 直接调用

```bash
# 设置运行模式（Docker 或 Conda）
export NEOAG_RUNNER_MODE=docker

# 所有工作流共享的通用格式
bin/neoag-nextflow run workflows/<workflow>.nf \
  --sample_id SAMPLE001 \
  <...workflow-specific params...> \
  -c conf/main_full.config \
  -profile docker
```

### 其他常用选项

```bash
--dry-run              # 只打印命令，不执行
--resume               # 续跑上次中断的任务
--outdir /path/to/out  # 指定输出目录
--profile_name leukemia # 使用特定癌种评分 profile
--strict_mode          # 生产模式（禁止 stub）
```

## 环境变量

| 变量 | 用途 | 默认值 |
|------|------|--------|
| `NEOAG_RUNNER_MODE` | 工具运行方式: `conda` 或 `docker` | `conda` |
| `NEOAG_REFERENCE_FASTA` | 参考基因组路径 | (需设置) |
| `NEOAG_DBSNP_VCF` | dbSNP VCF 路径 (FACETS 用) | (可选) |
| `NEOAG_PROFILE` | Nextflow 集群 profile | (可选) |
| `NXF_HOME` | Nextflow 元数据目录 | `work/.nextflow_home` |

## 输出目录结构

不同工作流的输出结构：

```
<outdir>/
├── scored/                    # 评分: 排序后的新抗原事件/肽段
├── report/                    # HTML 报告
├── upstream_*/                # 上游工具输出 (pVAC/NetMHCpan/MHCflurry/VEP)
├── hla_typing/                # HLA 分型结果 (仅 main_all / main_all_qc)
├── gatk_mutect2/              # 原始变异检出 (BAM 入口流程)
├── gatk_filter_mutect_calls/  # 过滤后 VCF (BAM 入口流程)
├── qc/                        # QC 报告 (仅 main_all_qc)
│   ├── lohhla/
│   ├── spechla/
│   ├── facets/
│   ├── purple/
│   └── quality_control.tsv
└── pipeline_info/             # Nextflow 执行报告
    ├── timeline.html
    ├── report.html
    └── trace.txt
```

## 相关文档

- [Workflow 技术文档](workflow.md) — Nextflow 工作流详细说明
- [输入架构](INPUT_ARCHITECTURE.md) — 多入口输入架构
- [安装与数据准备](INSTALL_AND_DATA.md) — 环境配置
- [Tools Setup](TOOLS_SETUP.md) — 外部工具安装
