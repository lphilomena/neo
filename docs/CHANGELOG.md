# CHANGELOG

本文件合并了原先分散在 `docs/V04_EVIDENCE_SAFETY_ESCAPE.md`、`docs/V041_APPM_CCF_IMMUNE_ESCAPE.md`、
`docs/V042_P1_APPM_EXPLAINABILITY.md`、`docs/V043_CCF21.md` 四份文档里的内容，按版本倒序排列。
四份原始文档已删除，避免重复维护。

## v0.4.3 — CCF 2.1（Quality / Multiplicity Confidence / External Clonality）

- 新增CCF输入QC侧车 `ccf_input_qc.tsv`。
- 新增multiplicity候选与置信度字段、克隆性置信度字段、粗粒度CCF区间概率摘要。
- 新增事件类型感知的CCF方法标签，包括 `WES_SV_CAPTURE_LIMITED_APPROX`、`RNA_ONLY_UNRESOLVED`。
- `ccf-2` CLI新增 `--external-clonality`（PyClone-VI/PhylogicNDT类输出）、`--svclone`（SVclone类输出）。
- 新增CCF冲突侧车 `ccf_conflicts.tsv`、聚类侧车 `ccf_cluster.tsv`。
- `build_ccf_2` 向后兼容旧调用方式，默认同时写出全部新侧车。

**边界**：CCF 2.1是面向计算三级分诊的克隆性证据层，不能替代专门的克隆结构解卷积工具，不产出临床克隆性诊断。

## v0.4.2 P1 — APPM Explainability

- 新增APPM调用置信度计算、MHC-I/MHC-II/IFNG子模块打分（`appm_submodule_scores.tsv`）。
- 新增 `reports_v041.py`：APPM证据卡片、MHC-I子模块卡片、Top APPM driver缺陷卡片、
  结合APPM/免疫逃逸/安全性/CCF的肽段机制卡片。
- 新增 `report-v041` 命令（完整参数示例见`docs/USAGE_GUIDE.md`）。
- `immune-escape` 新增 `--ranked-peptides`、`--top-priority-threshold`，输出受影响候选负担与top候选计数。
- `benchmark-system` 新增APPM相关输入选项及输出：`appm_ms_stratified_validation.tsv`、
  `appm_multiplier_delta.tsv`、`hla_ligand_detection_by_appm.tsv`。
- `appm_summary.tsv`/`appm_module_scores.tsv`/`appm_pathway_status.tsv`/`appm_peptide_modifiers.tsv`
  均新增APPM置信度字段；`immune_escape_events.tsv`新增逃逸事件CCF上下文和受影响肽段计数。

## v0.4.1 — APPM 2.0 + CCF 2.0 + Immune Escape 2.0

- 新增 `appm_v2.py`：APPM基因/通路状态、肽段级APPM标志、由损伤性变异+CN/LOH/表达证据推导的双等位基因状态逻辑。
- 新增 `ccf_v2.py`：拷贝数感知的CCF估计（含突变多重性枚举）、CCF min/max/best及置信度/方法字段。
- `immune_escape.py` 新增APPM 2.0/CCF 2.0侧车消费能力，新增治疗上下文策略
  （`vaccine`/`tcr_target`/`immunomonitoring`/`discovery`）。
- 新增CLI：`appm-2`、`ccf-2`；`immune-escape`扩展APPM/CCF/治疗上下文选项。
- `appm-lite`/`ccf-lite`改为由APPM 2.0/CCF 2.0驱动实现，同时保留原有返回结构和输出文件名，
  向后兼容`score_v03`所需的旧字段。

## v0.4 — Evidence / Safety / Escape Layer

- 新增 **捕获感知WES SV Phase 1.5**（`sv/wes_capture.py`）：capture BED解析、扩展BED sidecar、
  breakend捕获状态、WES置信度分层、优先级封顶。CLI新增 `--capture-bed`/`--capture-near-bp`/`--capture-slop-bp`。
- 新增**肽段安全性门控**（`peptide_safety_gate.py`）：参考蛋白质组精确匹配、正常配体组匹配、
  正常剪接junction匹配、仅锚点风险判定。输出 `safety/peptide_safety.tsv`、`safety/event_safety.tsv`。
- 新增**免疫逃逸/HLA LOH层**（`immune_escape.py`）：肽段级"呈递HLA丢失"标志、B2M/JAK/APM/CIITA风险摘要。
  输出 `immune_escape/immune_escape_events.tsv`、`immune_escape/immune_escape_summary.tsv`、
  `immune_escape/peptide_escape_flags.tsv`。
- `score-v03` 新增 `--peptide-safety`、`--peptide-escape-flags`；`run-full`/`sv-run-full`系列会自动构建并消费这些证据。

**边界**：
- WES SV Phase 1.5 是捕获限定的证据层，不等价于WGS级别的SV发现完整性。
- 肽段安全门控降低计算层面的off-target风险，不证明临床安全性。
- 免疫逃逸输出是机制/风险证据，不是临床耐药诊断。

## Skills 化重构（本轮）

- 新增 `.agents/skills/`：`pipeline-get`（入口路由/环境自检）+ 6个入口skill
  （`neoag-vcf`/`neoag-fusion`/`neoag-splice`/`neoag-sv-wgs`/`neoag-sv-wes`/`neoag-peptide-csv`）+
  `neoag-shared`（presentation之后的公共段说明）。
- `run-demo` 新增 `--entry-mode {snv_indel,fusion,splice_junction,sv_wgs,sv_wes,peptide_only}`，
  按入口给出范围内的工具检查+fixture冒烟测试，取代"要么全装要么不装"的单一全局检查方式
  （`check-tools`仍保留作为部署/CI阶段的全量核对手段）。
- 修复 `run-demo --entry-mode fusion/splice_junction/peptide_only` 的`shutil.SameFileError`崩溃
  （`pipeline_v03.py`中raw_events/raw_peptides源路径与目标路径相同时的自拷贝问题）。
- `docs/` 精简：合并/删除多份重叠文档（本CHANGELOG、`docs/USAGE_GUIDE.md`、`docs/RELEASE_BOUNDARY.md`即为整合结果）。
