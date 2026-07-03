---
name: pipeline-get
description: 检查并安装nextflow, 并获取可用的 pipeline 列表。在以下情况自动调用本 skill：1. 当用户提到介绍这个工具/仓库/pipeline，或者"帮我安装这个工具"、"这个工具怎么用"、"这个仓库里有什么"、"有什么功能"、"怎么分析"等意图时；2. 当用户提到或提供了.csv/.tsv/数据文件。
---

# NeoAg Pipeline 总览 / 环境自检 / 入口分诊

本 skill 不跑任何分析流程，只做三件事：环境自检（Python包+可选Nextflow）、
列出可用的子 pipeline（skill），以及根据用户提供的输入文件判断该走哪个子 skill。

## 第一步：环境自检

```bash
cd <repo_root>
python -m pip install -e '.[test]' -q
neoag-v03 --help >/dev/null && echo "OK: neoag-v03 CLI 可用"
```

只有用户明确要用 Nextflow 方式运行（而不是直接调用 `neoag-v03` CLI）时，才检查/安装 Nextflow：

```bash
source conf/tools.env.sh
bin/neoag-nextflow -version
```

`bin/neoag-nextflow` 会自动把 `NXF_HOME` 指到 `work/.nextflow_home`（避免污染仓库根目录），
第一次运行可能需要下载 Nextflow 运行时依赖，需要网络。

不要在这一步安装 VEP/NetMHCpan/EasyFuse 等重量级生信工具——那些工具的检查和安装
分散在各个子 skill 自己的运行前检查里，只有真正跑到那个入口时才需要。

## 第二步：告诉用户/AI有哪些子 pipeline 可用

本仓库把新抗原分析流程按输入文件类型拆成 6 个独立入口 skill，每个入口从各自的原始输入
一路跑到 `ranked_peptides.v03.tsv` / `evidence_report.v03.html`：

| Skill 名称 | 适用输入 | 典型触发场景 |
|---|---|---|
| `neoag-vcf` | 体细胞VCF（SNV/InDel），可选已有 pVACseq 结果 | 用户给了 `.vcf`/`.vcf.gz`，或提到"变异""VCF""Mutect2" |
| `neoag-fusion` | EasyFuse `fusions.pass.csv`，可选 pVACfuse 结果 | 用户提到"融合基因""fusion""EasyFuse" |
| `neoag-splice` | 变异VCF + RegTools splice junction TSV | 用户提到"剪接""splice junction""可变剪接" |
| `neoag-sv-wgs` | 全基因组SV VCF（Manta/GRIDSS/SvABA等）+ GTF + FASTA | 用户提到"结构变异""SV""WGS" |
| `neoag-sv-wes` | 同上 + capture BED（外显子捕获区间） | 用户提到"WES的SV""外显子捕获""capture" |
| `neoag-peptide-csv` | 已有的肽段-HLA表（CSV/TSV，两列起） | 用户直接上传/给了一份肽段表，没有变异层信息 |

## 第三步：入口判定规则

按以下优先级判断用户应该走哪个入口，判断不了就直接问：

1. 用户上传/提到的文件后缀是 `.csv`/`.tsv`，且看起来只有肽段序列+HLA两类信息（没有基因组坐标）
   → `neoag-peptide-csv`。
2. 文件是 `.vcf`/`.vcf.gz` 且看起来是SV记录（`SVTYPE`、`BND`、`INV`等字段）
   → 询问是WGS还是WES（是否有对应capture BED）→ `neoag-sv-wgs` 或 `neoag-sv-wes`。
3. 文件名/内容像 EasyFuse 的 `fusions.pass.csv`（含 `FTID`/`Fusion_Gene`等列）
   → `neoag-fusion`。
4. 文件是 RegTools/junction相关TSV，或用户明确提到splice/junction
   → `neoag-splice`。
5. 普通体细胞VCF（SNV/InDel，无SV特征字段）
   → `neoag-vcf`（默认/最常见入口）。
6. 无法从文件判断，或用户只是想了解项目功能、还没有文件
   → 展示上表，让用户选择，或直接问"你现在手上是什么类型的输入文件？"

判定后，直接切到对应子 skill 执行；不要在本 skill 里跑具体的产肽/打分命令。

## 边界声明

- 本仓库是研究级/验证-优先级的计算流程，不做临床治疗建议。候选排序结果需要湿实验验证、
  HLA分型、纯度、表达/蛋白支持和临床治理流程配合。
- 每个子 skill 都会在自己的运行前检查阶段告诉用户"这个入口具体需要哪些外部工具"，
  本 skill 不预先假设/安装任何重量级生信工具。
