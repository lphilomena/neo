# NeoAg Event Pipeline 使用教程

本文档是对最初重构需求的正式交付：

1. **第一部分**：6类输入各自需要运行的流程示例（模块命令组合），列出每步的中间文件和check标准。
2. **第二部分**：每个模块的参数说明、输出文件说明、所需环境说明。

配套的AI可路由skill定义在 `.agents/skills/`（`pipeline-get` 做入口判定，6个 `neoag-*` skill 对应下面6类入口，
`neoag-shared` 是被它们共用的公共段说明）。本文档面向人工阅读，skill文档面向AI agent执行，两者内容一致，
不需要重复维护——skill文档改动时请同步更新本文档对应章节，反之亦然。

所有命令假设已执行：

```bash
cd <repo_root>
python -m pip install -e '.[test]' -q
```

需要真实生信工具（VEP/NetMHCpan/EasyFuse等）时，先 `source conf/tools.env.sh`。工具安装细节见
[`INSTALL_AND_DATA.md`](INSTALL_AND_DATA.md)。

---

# 第一部分：6类输入的流程组合命令

本仓库把新抗原分析拆成6个独立入口。每个入口从各自的原始输入一路跑到
`ranked_peptides.v03.tsv` / `evidence_report.v03.html`，从 presentation 之后的步骤（第4步起）
全部6个入口完全共用。

每类入口的表格结构一致：**步骤 → 命令 → 中间文件 → check标准**，最后列出该入口的最终输出。
参数细节/输出字段/环境变量见第二部分。

## 6类入口对照速查表

| Entry              | 对应skill             | 独有起点命令                                         | 是否有专用一键路径               | 独有强制输入                                              |
| ------------------ | --------------------- | ---------------------------------------------------- | -------------------------------- | --------------------------------------------------------- |
| A SNV/InDel        | `neoag-vcf`         | `vep-annotate`→`snv-build-raw`                  | `run-full`（通用一键，读TOML） | `variants_vcf`、`tumor_sample_name`                   |
| B Fusion           | `neoag-fusion`      | `build-intermediates --entry-mode fusion`          | 无专用一键命令，用 `run-full`  | `easyfuse_pass_csv`                                     |
| C Splice junction  | `neoag-splice`      | `build-intermediates --entry-mode splice_junction` | 无专用一键命令，用 `run-full`  | `splice_junction_tsv`                                   |
| D1 SV-WGS          | `neoag-sv-wgs`      | `sv-build-raw`                                     | `sv-run-full`                  | `sv_vcf`、`reference_fasta`、`gencode_gtf`、`hla` |
| D2 SV-WES          | `neoag-sv-wes`      | `sv-build-raw-wes`                                 | `sv-run-full-wes`              | 同上 +`capture_bed`                                     |
| E Peptide-only CSV | `neoag-peptide-csv` | `build-intermediates --entry-mode peptide_only`    | 无专用一键命令，用 `run-full`  | `peptide_table`                                         |

> 原先规划中的"Entry F / e2e 全自动"入口已从用户可见流程中去掉：`entry_mode="e2e"`仍存在于
> 底层路由代码里（供多来源混合场景内部使用），但不作为独立入口对外暴露，本教程和 `.agents/skills/`都不会引导到它。

## Entry A｜SNV/InDel（体细胞VCF） → skill: `neoag-vcf`

| 步骤 | 命令                                                                                                                                                                                                                                                                                                                                                                                                  | 中间文件                                                                                                                      | check                     |
| ---- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------- | ------------------------- |
| 0    | `neoag-v03 run-demo --entry-mode snv_indel --outdir /tmp/demo --sample-id DEMO`                                                                                                                                                                                                                                                                                                                     | 全套fixture跑通                                                                                                               | 命令退出码0               |
| 1    | `neoag-v03 vep-annotate --input-vcf <vcf> --output-vcf <outdir>/upstream/tools/<id>.vep.annotated.vcf.gz --fasta $NEOAG_REFERENCE_FASTA --cache-dir $NEOAG_VEP_CACHE --plugins-dir $NEOAG_VEP_PLUGINS`（仅当VCF无CSQ时）                                                                                                                                                                            | `<outdir>/upstream/tools/<id>.vep.annotated.vcf.gz`                                                                         | 文件非空                  |
| 2    | `neoag-v03 snv-build-raw --variants-vcf <annotated_vcf> --outdir <outdir>/upstream --sample-id <id> --hla <HLA-A*02:01> <HLA-B*07:02> --tumor-sample-name <name>`（**必须用这条命令**，不要直接用 `extract-variant-peptides`的输出接下一步——它产出的 `variant_peptides.tsv`列名和 `peptide-predict`不兼容，会报 `Cannot detect peptide column`；`snv-build-raw`内部做了这层转换） | `upstream/tools/variant_peptides.tsv`、`upstream/parsed/raw_events.tsv`、`upstream/parsed/raw_peptides.tsv`             | `raw_peptides.tsv`非空  |
| 3    | `neoag-v03 peptide-predict -i <outdir>/upstream/parsed/raw_peptides.tsv -o <outdir>/presentation --sample-id <sample_id>`                                                                                                                                                                                                                                                                           | `presentation/presentation_evidence.tsv` 等                                                                                 | 非空                      |
| 4    | `appm-2` / `ccf-2` / `peptide-safety` / `immune-escape`（可并行）                                                                                                                                                                                                                                                                                                                             | `appm/appm_summary.tsv`、`clonality/ccf_2.tsv`、`safety/peptide_safety.tsv`、`immune_escape/peptide_escape_flags.tsv` | 各自非空                  |
| 5    | `score-v03`                                                                                                                                                                                                                                                                                                                                                                                         | `scoring/ranked_events.v03.tsv`、`scoring/ranked_peptides.v03.tsv`                                                        | 非空，行数与raw表量级相符 |
| 6    | `validation-plan-v03`                                                                                                                                                                                                                                                                                                                                                                               | `scoring/validation_plan.v03.tsv`                                                                                           | 非空                      |
| 7    | `report-v03 --audience both`                                                                                                                                                                                                                                                                                                                                                                        | `reports/evidence_report.{v03,patient,technical}.html`                                                                      | 三个html都生成            |

`snv-build-raw`是这轮新增的独立子命令（之前这一步的转换逻辑只能通过 `run-upstream`/`run-full`一键命令间接触发，
没有办法单独调用、单独check），参数上兼容 `extract-variant-peptides`的大部分选项
（`--lengths`/`--mini-len`/`--normal-proteome-fasta`/`--filter-normal-proteome`等）。
如果只是想预览滑窗产肽的原始结果（不需要接入打分链），仍然可以用 `extract-variant-peptides`单独跑，
但它的输出**不能**直接喂给 `peptide-predict`。

**一键替代路径**：`neoag-v03 run-full --config conf/run.<id>.private.toml --outdir results/<id>`
（内部调用与上表完全相同的底层函数，结果一致，见第二部分"双轨说明"）

**最终输出**：`scoring/ranked_events.v03.tsv`、`scoring/ranked_peptides.v03.tsv`、
`scoring/validation_plan.v03.tsv`、`reports/evidence_report.v03.html`（+patient/technical）、`provenance.v03.json`

## Entry B｜Fusion（EasyFuse融合） → skill: `neoag-fusion`

| 步骤 | 命令                                                                                                                                                                                                                                                                   | 中间文件                                                                               | check                                                       |
| ---- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------- | ----------------------------------------------------------- |
| 0    | `neoag-v03 run-demo --entry-mode fusion --outdir /tmp/demo --sample-id DEMO`                                                                                                                                                                                         | 全套fixture跑通                                                                        | 退出码0                                                     |
| 1    | `neoag-v03 build-intermediates --outdir <outdir> --entry-mode fusion --easyfuse-pass-csv <fusions.pass.csv> --hla <HLA-A*02:01> <HLA-B*07:02> --pvac <pvacfuse_aggregated.tsv>`（`--hla` **必需**，否则肽段生成为0条；`--pvac`可选，无pVACfuse结果时省略） | `parsed/raw_events.tsv`、`parsed/raw_peptides.tsv`、`parsed/fusion_evidence.tsv` | 三者非空，尤其确认 `raw_peptides.tsv`行数>1（不只是表头） |
| 2    | `neoag-v03 peptide-predict -i <outdir>/parsed/raw_peptides.tsv -o <outdir>/presentation`                                                                                                                                                                             | `presentation/presentation_evidence.tsv`                                             | 非空                                                        |
| 3-6  | 同Entry A的步骤4-7（公共段）                                                                                                                                                                                                                                           | 同上                                                                                   | 同上                                                        |

**最终输出**：同Entry A结构。

## Entry C｜Splice Junction（可变剪接） → skill: `neoag-splice`

真正产生肽段的是**pVACsplice聚合表**（`--pvac`），`--splice-junction-tsv`（+可选 `--variants-vcf`）
只是给已有肽段做junction支持度富集，**不会**独立从VCF生成肽段——这点和README早期版本、
以及我们上一版文档里的说法不一致，已核实并更正。

| 步骤 | 命令                                                                                                                                                                                                                                                                                                                                                        | 中间文件                                               | check                                     |
| ---- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------ | ----------------------------------------- |
| 0    | `neoag-v03 run-demo --entry-mode splice_junction --outdir /tmp/demo --sample-id DEMO`                                                                                                                                                                                                                                                                     | 全套fixture跑通                                        | 退出码0                                   |
| 1    | `neoag-v03 build-intermediates --outdir <outdir> --entry-mode splice_junction --pvac <pvacsplice_aggregated.tsv> --splice-junction-tsv <regtools_junctions.tsv> --hla <HLA-A*02:01> <HLA-B*07:02>`（`--pvac`是**真正产生肽段**的来源，`--splice-junction-tsv`只做junction支持度富集，不是产肽必需项；`--variants-vcf`同样只是可选富集上下文） | `parsed/raw_events.tsv`、`parsed/raw_peptides.tsv` | 非空，尤其确认 `raw_peptides.tsv`行数>1 |
| 2    | `neoag-v03 peptide-predict -i <outdir>/parsed/raw_peptides.tsv -o <outdir>/presentation`                                                                                                                                                                                                                                                                  | `presentation/presentation_evidence.tsv`             | 非空                                      |
| 3-6  | 同Entry A的步骤4-7（公共段）                                                                                                                                                                                                                                                                                                                                | 同上                                                   | 同上                                      |

**最终输出**：同Entry A结构。

## Entry D1｜SV-WGS（全基因组结构变异） → skill: `neoag-sv-wgs`

| 步骤 | 命令                                                                                                                                                                | 中间文件                                                                                                                                                                                                               | check                       |
| ---- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------- |
| 0    | `neoag-v03 run-demo --entry-mode sv_wgs --outdir /tmp/demo --sample-id DEMO`                                                                                      | 全套fixture跑通                                                                                                                                                                                                        | 退出码0                     |
| 1    | `neoag-v03 sv-build-raw --sample-id <id> --sv-vcf <sv.vcf...> --callers <Manta ...> --reference-fasta <ref.fa> --gencode-gtf <gtf> --hla <hla> --outdir <outdir>` | `parsed/raw_events.tsv`、`parsed/raw_peptides.tsv`、`sv/sv_events.full.tsv`、`sv/sv_protein_reconstruction.tsv`、`sv/sv_mutant_proteins.fa`、`sv/sv_event_to_peptide.tsv`、`sv/sv_validation_design.tsv` | raw_events/raw_peptides非空 |
| 2    | `neoag-v03 peptide-predict -i <outdir>/parsed/raw_peptides.tsv -o <outdir>/presentation`                                                                          | `presentation/presentation_evidence.tsv`                                                                                                                                                                             | 非空                        |
| 3-6  | 同Entry A的步骤4-7（公共段）                                                                                                                                        | 同上                                                                                                                                                                                                                   | 同上                        |

**一键替代路径**：`neoag-v03 sv-run-full --sample-id <id> --sv-vcf <sv.vcf> --reference-fasta <ref.fa> --gencode-gtf <gtf> --hla <hla> --outdir <outdir>`

**最终输出**：同Entry A结构，另加 `provenance.sv_phase1.json`。

**限制**（来自SV Phase 1原始文档，仍然成立）：

- BND/融合重建是启发式的CDS前缀/后缀重建，不是完整的重排图重建。
- 复杂重排图重建未实现。
- RNA junction证据是可选的，但纯DNA证据的候选应保守解读。
- `sv-build-raw`本身不做MHC结合预测，需要后续 `peptide-predict`。

## Entry D2｜SV-WES（外显子捕获结构变异） → skill: `neoag-sv-wes`

| 步骤 | 命令                                                                                                                                                                       | 中间文件                                   | check                       |
| ---- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------ | --------------------------- |
| 0    | `neoag-v03 run-demo --entry-mode sv_wes --outdir /tmp/demo --sample-id DEMO`                                                                                             | 全套fixture跑通                            | 退出码0                     |
| 1    | `neoag-v03 sv-build-raw-wes --sample-id <id> --sv-vcf <sv.vcf> --reference-fasta <ref.fa> --gencode-gtf <gtf> --hla <hla> --capture-bed <capture.bed> --outdir <outdir>` | 同SV-WGS，另加WES置信度分层字段            | raw_events/raw_peptides非空 |
| 2    | `neoag-v03 peptide-predict -i <outdir>/parsed/raw_peptides.tsv -o <outdir>/presentation`                                                                                 | `presentation/presentation_evidence.tsv` | 非空                        |
| 3-6  | 同Entry A的步骤4-7（公共段）                                                                                                                                               | 同上                                       | 同上                        |

**一键替代路径**：`neoag-v03 sv-run-full-wes --sample-id <id> --sv-vcf <sv.vcf> --reference-fasta <ref.fa> --gencode-gtf <gtf> --hla <hla> --capture-bed <capture.bed> --outdir <outdir>`

**最终输出**：同SV-WGS，`raw_events.tsv`含 `wes_confidence_tier`/`priority_cap`字段。

**WES置信度分层规则**：

| Tier                    | 判定标准                                       |
| ----------------------- | ---------------------------------------------- |
| `WES_Tier1`           | RNA junction reads ≥ 3，或达到WGS Tier1置信度 |
| `WES_Tier2`           | RNA junction reads ≥ 1，或达到WGS Tier2置信度 |
| `WES_Tier3`           | 仅当 `--tier1-only`未启用时保留              |
| `WES_UNINTERPRETABLE` | 均不满足                                       |

优先级封顶（`priority_cap`）默认映射：`Tier1→B`，`Tier2→B_CAUTION`，`Tier3→C`，`UNINTERPRETABLE→D`；
可在 `profiles/*.toml` 的 `[wes_confidence_caps]` 里覆盖。

**边界声明**：WES SV Phase 1.5 是"捕获限定"的证据层（`EXOME_CAPTURE_LIMITED`），断点外的区域可能漏检，
不等价于WGS级别的SV发现完整性。

## Entry E｜Peptide-only CSV（已有肽段表） → skill: `neoag-peptide-csv`

| 步骤 | 命令                                                                                                         | 中间文件                                                             | check    |
| ---- | ------------------------------------------------------------------------------------------------------------ | -------------------------------------------------------------------- | -------- |
| 0    | `neoag-v03 run-demo --entry-mode peptide_only --outdir /tmp/demo --sample-id DEMO`                         | 全套fixture跑通                                                      | 退出码0  |
| 1    | `neoag-v03 build-intermediates --outdir <outdir> --entry-mode peptide_only --peptide-table <peptides.csv>` | `parsed/raw_peptides.tsv`、`parsed/raw_events.tsv`（占位事件表） | 两者非空 |
| 2    | `neoag-v03 peptide-predict -i <outdir>/parsed/raw_peptides.tsv -o <outdir>/presentation`                   | `presentation/presentation_evidence.tsv`                           | 非空     |
| 3-6  | 同Entry A的步骤4-7（公共段）                                                                                 | 同上                                                                 | 同上     |

**最终输出**：同Entry A结构，但 `ranked_events.v03.tsv`的事件级字段（VAF/CCF等）是自动合成的占位值
（`source=peptide_input`），不是真实测序证据，报告中需注明。

如果只需要标准化后的肽段-HLA表本身、不进入完整打分链，也可以用更轻量的
`neoag-v03 convert-peptide-input -i <input_csv> -o <outdir>/parsed --sample-id <id>`，
但它只产出 `raw_peptides.tsv`/`peptide_hla_pairs.tsv`/`hla_alleles.txt`，不产出 `raw_events.tsv`，
无法直接接后续打分链（`score-v03`需要 `--raw-events`）。

---

# 第二部分：模块参数 / 输出 / 环境说明

## 双轨说明（一键命令 vs 分步命令）

`run-full`/`run-upstream`/`run-v03`（以及 `sv-run-full`/`sv-run-full-wes`）是"一键"封装，内部直接
`import`并调用与分步子命令背后**同一批Python函数**（`build_appm_lite`/`build_ccf_2`/`build_peptide_safety_gate`/
`build_immune_escape_evidence`/`make_validation_plan_v03`等，定义在 `pipeline_v03.py`），不存在两套逻辑，
分步跑和一键跑的结果不会不一致。分步命令的优势在于可以在任意一步之间检查中间文件、单独重跑失败的一步。

## 环境管理：`conf/tools.toml`（按工具分段）

之前环境变量散落在 `conf/tools.env.sh`里三十多个 `NEOAG_*`变量，一次性 `source`全部，
很难看出"这个入口到底需要哪几个"。现在改为两层：

- **`conf/tools.env.sh`**：仍然负责真正的环境解析（conda环境探测、PATH拼接、GPU/CPU切换、
  quarantine回退链）——这些是有真实条件分支的shell逻辑，不适合硬套成声明式配置，保留不动。
- **`conf/tools.toml`**：新增的结构化描述文件，每个工具一段（`[vep]`/`[netmhcpan]`/`[mhcflurry]`/
  `[easyfuse]`/`[gatk]`/`[manta]`/`[lohhla]`/`[facets]`……），声明：
  - `check`：怎么判定可用（`bin`查可执行文件、`dir`查目录、`env`查环境变量已设置）
  - 对应读取哪个 `NEOAG_*`环境变量
  - `entries`：这个工具被哪些 `entry_mode`用到
  - `optional`：缺失时是否可以降级（true）还是必须有（false，目前只有VEP在snv_indel/splice_junction标为必需）

`src/neoag_v03/tools_config.py` 是配套的加载器（`load_tools_toml`/`check_entry_tools`），
`run-demo --entry-mode`/`check-tools`都读它来做检查，不再各自维护一份工具清单。

## 环境检查：`run-demo --entry-mode` 与 `check-tools` 的分工

```bash
neoag-v03 run-demo --entry-mode {snv_indel,fusion,splice_junction,sv_wgs,sv_wes,peptide_only} --outdir <dir> --sample-id <id>
```

现在做两件事，而不只是打印工具名单：

1. **真实检测**：读 `conf/tools.toml`里该入口对应的工具段，实际检查二进制/目录/环境变量是否存在
   （不再是之前那种"打印一个固定字符串，不做任何检测"的假检查），打印每个工具 OK/MISSING。
2. **REAL/STUB双模式执行**：不再是"无论如何都用预处理好的pVAC聚合表fixture走捷径"。现在：
   - **产肽阶段**：全部6个入口都从各自**真实格式**的输入出发（体细胞VCF、EasyFuse
     `fusions.pass.csv`、pVACsplice聚合表+RegTools junction TSV、SV VCF、肽段CSV），
     和第一部分教程里教的命令组合是**同一段代码**，不再是走另一条"仅用于demo"的捷径。
   - **VEP阶段**（仅snv_indel）：VEP可用则真的跑 `vep-annotate`；不可用则退化为使用仓库自带的
     一份**已带CSQ注释**的fixture VCF继续走真实的 `extract-variant-peptides`（该步骤本身是纯Python，
     不需要VEP，缺的只是"谁来产生CSQ"这一步）。
   - **呈递预测阶段**（presentation，所有入口共用）：NetMHCpan/MHCflurry任一真实可用，就调用真实
     预测（`run_netmhcpan`/`run_mhcflurry`，和 `peptide-predict`背后同一批函数）；都不可用时才回退到
     仓库自带的预计算预测结果fixture。
   - 执行完会打印一行 `verification[<阶段>]: REAL/STUB/PARTIAL`，明确告诉你这次demo到底验证了什么，
     不会出现"demo通过≠教程里的命令能用"的情况。

```bash
neoag-v03 check-tools
```

仍然保留，给出**全部**已注册工具（`src/neoag_v03/tools/registry.py`里定义的完整工具集，覆盖面比
`conf/tools.toml`目前声明的更广，包括pVACtools/PRIME/BigMHC_IM/DeepImmuno等）的可用性总览，
适合部署/CI阶段一次性核对整套环境。**注意**：`check-tools`和 `run-demo --entry-mode`目前是两套独立的
工具清单（前者读 `tools/registry.py`，后者读 `conf/tools.toml`），尚未合并成一份——这是已知的后续优化项，
不是当前的既定行为，日常按入口开发调试建议优先用 `run-demo --entry-mode`。

## Entry独有模块

### `vep-annotate`

pVACseq兼容的VEP注释（含Wildtype/Frameshift插件）。

| 参数                                | 说明              |
| ----------------------------------- | ----------------- |
| `--input-vcf`                     | 未注释的体细胞VCF |
| `--output-vcf`                    | 注释后VCF输出路径 |
| `--fasta` / `--reference-fasta` | GRCh38参考FASTA   |
| `--cache-dir`                     | VEP cache目录     |
| `--plugins-dir`                   | VEP_plugins目录   |
| `--fork`                          | 并行度            |

环境：`NEOAG_VEP_BIN`、`NEOAG_VEP_CACHE`、`NEOAG_VEP_PLUGINS`、`NEOAG_REFERENCE_FASTA`（均在 `conf/tools.env.sh`里定义）。
输出：注释后的 `.vcf.gz`，含 `CSQ`/`WildtypeProtein`/`Frameshift`字段。

### `extract-variant-peptides`

VEP CSQ滑动窗口全枚举产肽（不经pVACseq）。

| 参数                                                       | 说明                     |
| ---------------------------------------------------------- | ------------------------ |
| `--input-vcf`                                            | VEP注释过的VCF           |
| `--output`                                               | 输出TSV路径              |
| `--sample-id`                                            | peptide_id前缀           |
| `--lengths` / `--length-min` / `--length-max`        | 肽段长度范围（默认8-11） |
| `--mini-len`                                             | minigene长度             |
| `--normal-proteome-fasta` + `--filter-normal-proteome` | normal proteome安全过滤  |
| `--hla-alleles`                                          | HLA分型                  |
| `--tumor-sample-name`                                    | VCF肿瘤样本列名          |

环境：无外部工具依赖，纯Python。输出：`variant_peptides.tsv`（滑窗肽段全集）。

> **注意**：这个命令的输出 `variant_peptides.tsv`用的是 `mutant_peptide`列名（变异专属schema），
> **不能**直接喂给 `peptide-predict`（后者内部靠列名别名识别 `peptide`列，认不出 `mutant_peptide`，
> 会报 `Cannot detect peptide column`）。要接入打分链，请用下面的 `snv-build-raw`，
> 不要把这个命令的输出直接传给 `peptide-predict`。

### `snv-build-raw`（Entry A的产肽→标准化，一步到位）

| 参数                                                                                             | 说明                                                                                        |
| ------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------- |
| `--variants-vcf`                                                                               | VEP注释过的VCF（必须已有CSQ）                                                               |
| `--outdir`                                                                                     | 写入 `parsed/raw_events.tsv`、`parsed/raw_peptides.tsv`、`tools/variant_peptides.tsv` |
| `--hla`                                                                                        | 一个或多个HLA分型，**必需**                                                           |
| `--tumor-sample-name` / `--rna-sample-name`                                                  | VCF样本列名                                                                                 |
| `--lengths` / `--length-min` / `--length-max` / `--mini-len`                             | 同 `extract-variant-peptides`                                                             |
| `--normal-proteome-fasta` + `--filter-normal-proteome` / `--annotate-normal-proteome-only` | normal proteome安全过滤                                                                     |
| `--exclude-multi-aa` / `--single-aa-only`                                                    | 同 `extract-variant-peptides`                                                             |
| `--easyfuse-pass-csv` / `--easyfuse-tsv`                                                     | 可选，snv_indel+fusion合并场景一并产出                                                      |

环境：无外部工具依赖，纯Python。内部先做滑窗产肽（同 `extract-variant-peptides`），
再经 `variant_peptide_adapter`转换成标准 `raw_peptides.tsv`（含 `peptide`列，`peptide-predict`能直接消费）。
这是本轮新增的命令——之前这一步转换逻辑只能通过 `run-upstream`/`run-full`一键命令间接触发，
无法单独调用、单独check中间结果，是"分步执行"设计里的一个真实缺口，现在补上。

### `build-intermediates`（Entry B/C/E共用的标准化入口）

| 参数                                         | 说明                                                                                                                |
| -------------------------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| `--entry-mode`                             | `fusion` / `splice_junction` / `peptide_only` 等                                                              |
| `--easyfuse-pass-csv` / `--easyfuse-tsv` | Fusion入口用                                                                                                        |
| `--splice-junction-tsv`                    | Splice入口用                                                                                                        |
| `--peptide-table`                          | Peptide-only入口用                                                                                                  |
| `--hla`                                    | **Fusion/Splice入口必需**，一个或多个HLA分型（如 `HLA-A*02:01 HLA-B*07:02`）；fusion入口不传则肽段生成为0条 |
| `--pvac`                                   | **Splice入口的真实产肽来源**（pVACsplice聚合表）；fusion入口可选，补充pVACfuse表位                            |
| `--variants-vcf`                           | 可选，仅给splice入口的junction支持度富集提供上下文，**不会**独立从VCF生成肽段（真正产肽靠 `--pvac`）        |
| `--splice-junction-tsv`                    | Splice入口用，给已有肽段（来自 `--pvac`）做junction支持度富集，本身不产肽                                         |
| `--raw-events` / `--raw-peptides`        | 已有raw表时直接透传                                                                                                 |

环境：纯Python，无外部工具依赖。输出：`parsed/raw_events.tsv` + `parsed/raw_peptides.tsv`（+入口专属sidecar，如 `fusion_evidence.tsv`）。

> **验证提示**：`build-intermediates`跑完之后，务必检查 `raw_peptides.tsv`的**行数**而不只是"文件是否存在"——
> fusion入口缺 `--hla`、splice入口缺 `--pvac`时会**静默**只产出事件表、肽段表为空（不会报错，容易被忽略）。
> `run-demo --entry-mode fusion`/`splice_junction`现在走的是和本文档一致的真实命令组合
> （不再是先前版本文档里的捷径），demo跑通即代表这条命令组合本身没问题。

### `sv-build-raw` / `sv-build-raw-wes`

| 参数                                                     | 说明                          |
| -------------------------------------------------------- | ----------------------------- |
| `--sv-vcf` (可多个) + `--callers`                    | SV VCF及对应caller名          |
| `--reference-fasta` / `--gencode-gtf`                | GRCh38参考+注释               |
| `--hla`                                                | HLA分型文件或列表             |
| `--expression` / `--rna-junctions`                   | 可选RNA证据                   |
| `--capture-bed`（仅wes版）                             | WES捕获区间，启用捕获感知分层 |
| `--capture-near-bp` / `--capture-slop-bp`（仅wes版） | 捕获边界宽容度                |
| `--tier1-only`                                         | 只导出Tier1事件               |

内部流程：`read_sv_inputs → 聚类 → reconstruct_cluster_protein → build_mhc1_peptides`（8-11aa滑窗）→ 展开为肽段×HLA。
环境：无外部工具依赖（上游SV caller如Manta/GRIDSS/SvABA是独立跑的，本命令只处理其VCF输出）。
输出：`parsed/raw_events.tsv`/`raw_peptides.tsv` + `sv/sv_events.full.tsv`/`sv_protein_reconstruction.tsv`/`sv_mutant_proteins.fa`/`sv_event_to_peptide.tsv`/`sv_validation_design.tsv`。

### `convert-peptide-input`

轻量CSV/TSV→标准肽段表转换器（不产出 `raw_events.tsv`，见Entry E小节说明）。

| 参数            | 说明                      |
| --------------- | ------------------------- |
| `-i/--input`  | 含肽段序列+HLA列的CSV/TSV |
| `-o/--outdir` | 输出目录                  |
| `--sample-id` | 样本ID                    |

输出：`raw_peptides.tsv`、`peptide_hla_pairs.tsv`、`hla_alleles.txt`。

## 公共段模块（presentation之后，6个入口共用）

### `peptide-predict`

| 参数                                                                                                                               | 说明                                   |
| ---------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------- |
| `-i/--input`                                                                                                                     | `raw_peptides.tsv`                   |
| `-o/--outdir`                                                                                                                    | 输出到 `presentation/`               |
| `--stub`                                                                                                                         | 用桩预测器快速冒烟测试，不需要真实工具 |
| `--skip-netmhcpan` / `--skip-mhcflurry` / `--skip-prime` / `--skip-bigmhc-im` / `--skip-deepimmuno` / `--skip-stabpan` | 按需跳过某个预测器                     |

对应工具中netmhcpan, mhcflurry, bigmhc-im, prime必须运行。

环境：NetMHCpan（`NEOAG_NETMHCPAN_BIN`）、MHCflurry（conda环境）、NetMHCstabpan、PRIME/BigMHC_IM/DeepImmuno，均可选，
缺失时用 `--stub`或对应 `--skip-*`降级。输出：`presentation/presentation_evidence.tsv`等。

### `appm-2`（抗原呈递机器完整性）

| 参数               | 说明                    |
| ------------------ | ----------------------- |
| `--vep-tsv`      | VEP衍生的APPM相关变异表 |
| `--expression`   | 基因表达                |
| `--hla-loh`      | HLA杂合性丢失表         |
| `--cnv`          | 拷贝数分段              |
| `--raw-peptides` | 肽段表                  |
| `--tumor-purity` | 肿瘤纯度                |

任一证据缺失都可以省略，APPM会记录为"证据不完整"而不是报错。输出：`appm_summary.tsv`、`appm_gene_status.tsv`、
`appm_pathway_status.tsv`、`appm_module_scores.tsv`、`appm_peptide_modifiers.tsv`等。

### `ccf-2`（克隆性）

| 参数                                     | 说明                                                     |
| ---------------------------------------- | -------------------------------------------------------- |
| `--events`                             | `raw_events.tsv`                                       |
| `--purity` / `--cnv`                 | 纯度/拷贝数（可选）                                      |
| `--external-clonality` / `--svclone` | 可选外部克隆性工具（PyClone-VI/PhylogicNDT/SVclone）输入 |

输出：`ccf_2.tsv`（+`ccf_input_qc.tsv`/`ccf_conflicts.tsv`/`ccf_cluster.tsv`侧车）。无purity/cnv时方法标签会降级为
`WES_SV_CAPTURE_LIMITED_APPROX`或 `RNA_ONLY_UNRESOLVED`，属正常降级不算失败。

### `peptide-safety`（肽段安全性门控）

| 参数                                               | 说明                  |
| -------------------------------------------------- | --------------------- |
| `--raw-events` / `--raw-peptides`              | 必需                  |
| `--normal-expression` / `--normal-hla-ligands` | 正常组织表达/配体证据 |
| `--reference-proteome`                           | 参考蛋白质组比对      |
| `--normal-junctions`                             | 正常剪接junction比对  |

输出：`peptide_safety.tsv`、`event_safety.tsv`。**边界**：计算层面降低off-target风险，不是临床安全性证明。

### `immune-escape`（免疫逃逸证据）

| 参数                                                         | 说明                                                                              |
| ------------------------------------------------------------ | --------------------------------------------------------------------------------- |
| `--raw-peptides`                                           | 必需                                                                              |
| `--vep-tsv` / `--cnv` / `--expression` / `--hla-loh` | 各类证据（可选）                                                                  |
| `--appm-gene-status` / `--appm-pathway-status`           | 建议先跑完 `appm-2`再传入，避免APPM重复判定                                     |
| `--ccf`                                                    | 克隆性证据                                                                        |
| `--therapy-context`                                        | `vaccine`/`tcr_target`/`immunomonitoring`/`discovery`，影响优先级封顶策略 |

输出：`immune_escape_summary.tsv`、`peptide_escape_flags.tsv`。**边界**：机制/风险层面证据，不是临床耐药诊断。

### `score-v03`（打分排序）

| 参数                                                       | 说明                |
| ---------------------------------------------------------- | ------------------- |
| `--raw-events` / `--raw-peptides` / `--presentation` | 必需                |
| `--appm-summary` / `--appm-peptide-modifiers`          | APPM证据            |
| `--ccf`                                                  | 克隆性证据          |
| `--normal-expression` / `--normal-hla-ligands`         | 安全性证据          |
| `--peptide-safety` / `--peptide-escape-flags`          | 安全性/免疫逃逸证据 |

输出：`ranked_events.v03.tsv`、`ranked_peptides.v03.tsv`。任一上游sidecar为空会导致打分结果异常缩水，
是最常见的"最终结果行数偏少"的排查起点。

### `validation-plan-v03`

| 参数                   | 说明                                       |
| ---------------------- | ------------------------------------------ |
| `--ranked-peptides`  | 必需                                       |
| `--variant-peptides` | 可选，minigene列信息来源                   |
| `--outdir`           | 自动发现 `upstream/tools/*_peptides.tsv` |

输出：`validation_plan.v03.tsv`。

### `report-v03` / `report-v041`

| 参数                                        | 说明                                             |
| ------------------------------------------- | ------------------------------------------------ |
| `--ranked-events` / `--ranked-peptides` | 必需                                             |
| `--appm-summary` / `--validation-plan`  | 补充证据卡片                                     |
| `--outdir`                                | 自动加载APPM/安全性/逃逸/CCF侧车（技术版报告用） |
| `--provenance`                            | `provenance.v03.json`，工具版本追溯            |
| `--audience`                              | `both`/`patient`/`technical`               |

输出：`evidence_report.v03.html`（`--audience both`时另加 `.patient.html`/`.technical.html`）。
`report-v041`是v0.4.1时期新增的APPM/逃逸/安全性/CCF证据卡片专版报告，字段更细，可按需替代或补充 `report-v03`。

## 常见失败排查

| 症状                                                       | 排查方向                                                                                                                           |
| ---------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------- |
| `raw_events.tsv`/`raw_peptides.tsv`为空                | 检查对应Entry的独有输入文件列名/编码是否符合期望格式（见第一部分各Entry表格）                                                      |
| VEP步骤报错                                                | 检查 `NEOAG_VEP_BIN`/`NEOAG_VEP_CACHE`/`NEOAG_VEP_PLUGINS`/`NEOAG_REFERENCE_FASTA`是否指向真实路径                         |
| `ranked_peptides.v03.tsv`行数远小于 `raw_peptides.tsv` | 依次检查 `presentation_evidence.tsv`/`appm_summary.tsv`/`peptide_safety.tsv`是否为空——任一上游证据缺失都会导致打分结果缩水 |
| SV捕获状态全部落空（WES）                                  | 确认 `capture_bed`坐标体系（GRCh38）与SV VCF/GTF一致                                                                             |
| `run-demo`某个entry-mode跑不通                           | 先确认 `pip install -e .`已重新执行（尤其是刚pull过代码之后）                                                                    |
