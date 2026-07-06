# NeoAg Event Pipeline v0.4.3

NeoAg Event Pipeline 是一个面向研究的新抗原优先级排序流程。它把SNV/InDel、融合、剪接、结构变异、
以及已有肽段表这几类候选，统一转换为标准化的事件表和肽段-HLA表，再叠加结合呈递、APPM、CCF、安全性、
免疫逃逸、验证计划、报告等证据层。

**边界声明**：本流程产出的是计算层面的三级分诊和验证规划结果，不做临床诊断、临床耐药判定，
也不给出经过验证的治疗建议。

## 能做什么

- 解析类pVACtools的SNV/融合/剪接输出，或直接从VEP注释的VCF生成滑动窗口变异肽段
  （VCF无CSQ注释时可自动跑VEP）。
- 对NetMHCpan、MHCflurry及可选的稳定性/免疫原性工具做结合呈递证据打分。
- 构建APPM 2.0（抗原呈递机器完整性）、CCF 2.1（克隆性）、肽段安全性、免疫逃逸证据层。
- 产出排序后的事件/肽段表、验证计划，以及患者版+技术版双受众HTML报告。
- 支持直接用`neoag-v03` CLI运行，也支持通过内置Nextflow封装运行。

排序输出里的`.v03.tsv`后缀是schema兼容性标签，不代表软件版本号。

## 6类输入入口

整个流程按6类独立输入入口组织，每个入口从各自原始输入出发，最后汇入同一套打分/报告公共段。
完整的命令级教程（模块命令组合、中间文件、最终输出、参数、所需环境）见
**[`docs/USAGE_GUIDE.md`](docs/USAGE_GUIDE.md)**。

| 入口 | 输入 | 对应AI skill（`.agent/skills/`） |
|---|---|---|
| SNV/InDel | 体细胞VCF（可选pVACseq） | `neoag-vcf` |
| Fusion | EasyFuse `fusions.pass.csv`（可选pVACfuse） | `neoag-fusion` |
| Splice junction | VCF + RegTools junction TSV | `neoag-splice` |
| SV（WGS） | 结构变异VCF + GTF + FASTA | `neoag-sv-wgs` |
| SV（WES） | 同上 + capture BED | `neoag-sv-wes` |
| Peptide-only | 肽段+HLA的CSV/TSV | `neoag-peptide-csv` |

如果你在用支持`.agent/skills/`的AI编程/agent工具，从`pipeline-get`开始——它会检查环境、
列出上面6个入口，并根据你的输入文件路由到对应入口。每个入口skill和`docs/USAGE_GUIDE.md`
描述的是同一套命令，两边任一处改动都要同步更新另一处。

## 快速开始

```bash
python -m pip install -e '.[test]'
neoag-v03 run-demo --entry-mode snv_indel --outdir work/demo_snv --sample-id DEMO001
```

`run-demo --entry-mode {snv_indel,fusion,splice_junction,sv_wgs,sv_wes,peptide_only}`
会针对单个入口跑一遍完整的fixture冒烟测试，并只打印该入口相关的工具清单——不需要在尝试流程前
装齐所有可选工具。

## 测试

```bash
pytest -q                       # 默认：快速单元测试 + 发布安全性契约检查
pytest -q --run-integration
pytest -q --run-benchmark
pytest -q --run-external
pytest -q --run-all             # 全部；需要外部工具/网络/Nextflow缓存
```

## 文档索引

| 文档 | 内容 |
|---|---|
| [`docs/USAGE_GUIDE.md`](docs/USAGE_GUIDE.md) | 完整教程：各入口命令链 + 各模块参数/输出/环境参考 |
| [`docs/INSTALL_AND_DATA.md`](docs/INSTALL_AND_DATA.md) | 基础环境、Python/conda、Nextflow、外部工具安装表、参考数据、验收检查 |
| [`docs/CHANGELOG.md`](docs/CHANGELOG.md) | 版本历史（v0.4 → v0.4.3，及skills化重构） |
| [`docs/RELEASE_BOUNDARY.md`](docs/RELEASE_BOUNDARY.md) | 轻量在线发布包里包含什么/排除什么 |
| [`docs/SITE_CONFIG_BOUNDARY.md`](docs/SITE_CONFIG_BOUNDARY.md) | 哪些配置文件应提交、哪些应保留在本地 |

## 许可与引用

见 [`LICENSE`](LICENSE)、[`NOTICE`](NOTICE)、[`CITATION.cff`](CITATION.cff)。
