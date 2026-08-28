# NeoAg 候选来源链置信度 C1–C4

## 1. 目的与边界

本层回答：**候选肽能否可靠回溯到所声称的 SNV、InDel、Fusion 或 Splice 事件，并形成事件→转录本/ORF→肽段→HLA 结果的有效来源链。**

它不回答候选是否一定能自然呈递、被 T 细胞识别或产生临床获益。

- `C1–C4`：候选来源链置信度。
- `R1–R4`：综合 RNA、呈递、CCF、HLA/APPM、安全性等后的实验推荐等级。
- `C1 != R1`；C1 候选仍可因 HLA 丢失、安全性或呈递失败而成为 R4。

## 2. 统一分级

| 等级 | 定义 | 默认处理 |
|---|---|---|
| C1 | 核心来源链全部成立，并有独立/跨模态正交确认 | 可进入完整 R1–R4 评估 |
| C2 | 核心计算和原始 reads 证据链完整、无重大冲突，但无独立正交确认 | 研究 profile 可进入完整评估；转化 profile 默认最高 R2 |
| C3 | 来源链合理但至少一个适用核心环节未评估、低功效或未闭合 | integrated profile 最高 R3，先补证据 |
| C4 | 事件、ORF、肽段回链或样本/坐标被反证或无效 | Hard fail，通常 R4 |

确定性判定：

```text
任一 fatal NEGATIVE/CONFLICT → C4
全部适用核心项 SUPPORTED + orthogonal SUPPORTED → C1
全部适用核心项 SUPPORTED → C2
其余 → C3
```

## 3. 状态语义

每个要求必须区分：

- `NOT_APPLICABLE`：对该事件类型或当前场景不适用。
- `UNASSESSED`：适用，但没有完成评估。
- `INDETERMINATE_LOW_POWER`：已尝试评估，但覆盖、reads 或质量不足，无法可靠判断。
- `SUPPORTED`：评估后获得支持。
- `NEGATIVE`：检测能力充分，但没有支持；是否构成 C4 由该要求的 fatal 属性决定。
- `CONFLICT`：不同来源相互矛盾，需要人工审阅；关键冲突可触发 C4。

禁止把 `NOT_APPLICABLE`、`UNASSESSED` 或低功效未检出改写为生物学阴性。

## 4. 四类事件专属核心要求

### 4.1 SNV

核心项：

1. 肿瘤 DNA 精确 SNV 和 read-level QC。
2. 匹配正常阴性，或明确进入独立胚系路线。
3. RNA ALT 或直接蛋白/肽证据。
4. 转录本、codon、蛋白改变可解析。
5. 存在邻近变异时完成 read-backed phasing；否则为 `NOT_APPLICABLE`。
6. 候选肽包含突变氨基酸。
7. MT、WT、HLA 结果可精确回链。

典型 C4：原始 DNA reads 不支持、mapping artifact、候选肽不含突变、phasing 证明序列不可能存在、坐标或 REF/ALT 错误、样本混淆。

### 4.2 InDel

核心项：

1. left-normalized 最简表示。
2. 局部重比对、重复区、微同源和 homopolymer QC。
3. 肿瘤 DNA ALT 和匹配正常证据。
4. RNA InDel 或直接蛋白/肽证据。
5. in-frame/frameshift ORF 正确重建。
6. frameshift 的 novel tail、终止位置和 NMD 注释；非 frameshift 为 `NOT_APPLICABLE`。
7. 邻近变异 phasing。
8. 肽段包含 InDel 或新尾部序列。
9. 肽段-HLA 结果可回溯到正确 ORF。

典型 C4：局部重比对否定事件、归一化后原事件/肽段身份失效、reading frame 错误、候选肽不在所称蛋白中。

### 4.3 Fusion

核心项：

1. partner、精确 breakpoint、5'/3'方向和链信息。
2. unique split/junction reads。
3. spanning pairs、多独立 reads、caller 共识或其他结构支持。
4. 排除 read-through、正常 junction、高同源/重复区伪影。
5. reading frame 和融合 ORF。
6. 肽段跨 breakpoint 或包含融合新序列。
7. HLA 结果可回溯到融合 ORF。

独立 RT-PCR/Sanger、长读长、DNA-SV 或独立 RNA 文库可将完整 C2 提升为 C1。同一 RNA BAM 的多个 caller 仅是计算一致性，不等同正交确认。

### 4.4 Splice

核心项：

1. 精确 build/chromosome/intron/junction/strand。
2. 精确 unique split reads。
3. anchor/overhang、MAPQ、重复区/伪基因风险 QC。
4. SE/A3SS/A5SS/RI/cryptic exon 等结构可解析。
5. reference path 与 alternative path 关系。
6. transcript hypothesis、frame 和可翻译 ORF。
7. 肽段跨异常 junction 或包含异常剪接新序列。
8. HLA 结果回溯到对应 ORF/junction。
9. 正常 junction 背景强制评估。

典型 C4：无法定位精确 junction、用同基因另一 junction reads 回填、mapping/坐标伪影、ORF 无效或肽段不能由事件产生。

## 5. 正交确认

可计为独立或跨模态确认的例子：

- WES/WGS 对同一 SNV/InDel 一致，并满足各自 QC。
- DNA 事件与 RNA ALT/直接蛋白证据一致。
- RT-PCR + Sanger 确认 Fusion/Splice 断点。
- 长读长、独立文库、DNA-SV 或靶向深测确认。

不自动视为正交确认：

- 两个 caller 读取同一 BAM。
- 同一结果表的重复解析。
- gene TPM 对突变或 junction 的间接支持。

## 6. 与 Evidence Consensus 的两种接入模式

### 兼容模式（推荐先使用）

配置：`sarcoma_evidence_consensus_v2_1_source_chain.toml`

- 输出 C1–C4 和逐项审计。
- 不改变现有 R1–R4、Pareto front 或 rank。
- 用于历史病例回顾和规则验证。

### 集成模式（研究性）

配置：`sarcoma_evidence_consensus_v3_source_chain.toml`

- C4 → hard fail/R4。
- C3 → 最高 R3。
- translational profile 下 C2 → 默认最高 R2。
- C1–C4 grade 进入赛道内 Pareto 和 tie-break。

集成模式仍为 `PROVISIONAL_RESEARCH_ONLY`，需要历史病例和实验反馈验证后才能替代兼容模式。

## 7. 主要输出

```text
source_chain_confidence.tsv
source_chain_requirements.long.tsv
ranked_peptides.evidence_consensus.tsv
ranked_events.evidence_consensus.tsv
```

关键字段：

```text
source_chain_track
source_chain_confidence_tier
source_chain_confidence_label
source_chain_orthogonal_status
source_chain_orthogonal_sources
source_chain_hard_failure_codes
source_chain_reason_codes
source_chain_missing_requirements
source_chain_low_power_requirements
source_chain_negative_requirements
source_chain_conflict_requirements
source_chain_requirement_statuses
source_chain_rule_version
```

## 8. CLI

独立评估：

```bash
neoag source-chain \
  --input scoring/all_tool_results.tsv \
  --output scoring/source_chain_confidence.tsv \
  --requirements-out scoring/source_chain_requirements.long.tsv \
  --rules configs/ranking/sarcoma_evidence_consensus_v2_1_source_chain.toml
```

兼容模式双排序：

```bash
neoag evidence-rank \
  --comprehensive-evidence scoring/all_tool_results.tsv \
  --weighted-baseline scoring/ranked_peptides.weighted_baseline.tsv \
  --rules configs/ranking/sarcoma_evidence_consensus_v2_1_source_chain.toml \
  --outdir scoring/evidence_consensus_v2_1
```

集成模式：

```bash
neoag evidence-rank \
  --comprehensive-evidence scoring/all_tool_results.tsv \
  --weighted-baseline scoring/ranked_peptides.weighted_baseline.tsv \
  --rules configs/ranking/sarcoma_evidence_consensus_v3_source_chain.toml \
  --outdir scoring/evidence_consensus_v3_alpha
```

报告字段对齐审计：

```bash
neoag report-dimension-audit \
  --input scoring/evidence_consensus_v2_1/ranked_peptides.evidence_consensus.tsv \
  --output scoring/report_dimension_audit.tsv
```

## 9. v0.5.2 source-chain rule corrections

- SNV now audits tumor/normal depth, ALT support, base quality, mapping quality,
  strand bias, FFPE artifacts, low-complexity and paralogous-region risk. Phasing
  is `NOT_APPLICABLE` unless a proximal variant can alter the candidate sequence.
- InDel representation changes are not failures by themselves. A hard failure is
  emitted only when canonical biological-event identity changes or the reconstructed
  ORF/peptide is invalid.
- RNA-only fusion is not automatically capped at R3 when a complete C1/C2
  breakpoint-to-ORF-to-peptide chain is present. Incomplete C3 chains remain capped.
- Fusion split reads and spanning pairs are evaluated jointly with independent
  start sites, duplicate handling, caller support and junction sequence uniqueness.
- Splice assessment includes exact-junction reads, anchor/overhang, MAPQ,
  PSI/junction usage, normal-isoform context, translation direction and NMD risk.
- Plain `NOT_DETECTED` is not normal-negative evidence. Normal-background support
  requires adequate coverage; insufficient coverage is
  `INDETERMINATE_LOW_POWER`, not `NEGATIVE`.
- Requirement outputs include applicability, status, value, reason code, source
  fields and conflict flag. Reports display source-chain C tier separately from
  the final R recommendation tier.
