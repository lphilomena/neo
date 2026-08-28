# Splice P0 修复说明（v0.5.2-p0）

本版本把此前仅存在于预过滤或人工审查中的关键限制，落实到正式
Splice Provenance Layer、证据链和最终 R 分级中。目标是阻止低质量
junction、错误事件借证、未解析的新氨基酸边界及局部肽产物被升级。

## 1. 正式 junction read QC

新增输出 `splice_junction_read_qc.tsv`。默认 `complete` 策略要求每条
候选来源 junction 同时满足：

- 精确且有链的 canonical junction；
- unique split reads >= 3；
- unique fragment starts >= 2；
- max overhang >= 10；
- mapping quality >= 20；
- multimapping fraction <= 0.20；
- tumor PSI >= 0.05；
- caller 未明确标记 FAIL/REJECT/BLACKLIST/ARTIFACT。

缺失必需字段得到 `INCOMPLETE`，阈值失败得到 `FAIL`。只有显式 `PASS`
可以进入正式 RNA junction 支持集合。所有阈值均可通过 CLI 修改。

建议主 junction TSV 提供以下扩展列：

```text
unique_split_reads
multi_split_reads
unique_fragment_starts
max_overhang
median_mapq
multimapping_fraction
tumor_psi
caller_filter
```

历史数据可显式使用 `--junction-qc-policy reads_only`，但该策略仅检查
有链精确坐标、unique split reads 和显式 caller failure，属于兼容模式，
不建议用于候选自动升级。

## 2. peptide-origin 级 junction 支持

`splice_peptide_origins.tsv` 新增 `required_junction_ids`。共识层要求该肽
来源所需的全部 junction 都有 QC PASS 和精确 RNA 证据。事件内其他
junction 的 reads 不再能够借给当前肽来源。

旧记录没有该字段时，按以下顺序保守解析：

1. origin 的 `required_junction_ids`；
2. origin 的 `junction_ids`；
3. event 的 `alternative_junction_ids`。

不会退回“事件中任意 junction 有支持即可”的旧规则。

## 3. pVACsplice 残基级新颖性

pVACsplice 适配器不再无条件写入 `crosses_junction=true` 和
`contains_novel_aa=true`。只有报告提供可落在表位内部的显式突变/新残基
位置时，才写入 `contains_novel_aa=true`；只有提供有效的肽内 junction
边界时，才写入 `crosses_junction=true`。否则均保留为 `UNASSESSED`，并
触发 R3 上限。

当前支持的残基位置列包括 `Pos`、`Mutation Position`、
`Mutation Position(s)`、`Novel AA Positions`；junction 边界列包括
`Junction AA Position`、`Junction Position` 和
`Junction Offset in Peptide`。超出表位长度的位置不会被猜测转换。

## 4. 正常背景状态统一

证据链和共识层共同识别：

- `DETECTED_MATCHED_NORMAL`
- `DETECTED_CRITICAL_TISSUE`
- `DETECTED_BROAD_NORMAL`
- `LOW_LEVEL_NONCRITICAL_NORMAL`
- 历史兼容状态 `NORMAL_DETECTED`

关键组织检出在两个层面均产生硬性负面结论；覆盖不足的未检出仍保持
未完成，k4neo 阴性不等同于 locus coverage 阴性。

## 5. fallback 事件不再伪称 novel

未被正式事件图连接的单个 junction 使用事件类型
`JUNCTION_ONLY_UNCLASSIFIED`，junction role 为 `UNCLASSIFIED_EDGE`，
且 `alternative_junction_ids` 为空。只有经过注释和正常背景判定的来源
适配器才能正式声明 novel junction。

## 6. partial/epitope-only 产物上限

以下 ORF 状态仍可保留供研究和 HLA 预测，但最终等级最高为 R3：

- `PARTIAL_TRANSLATED_SEGMENT`
- `VALID_PEPTIDE_PRODUCT_ONLY`
- `VALID_EPITOPE_PRODUCT_ONLY`

其共识中会记录 `CAP_PARTIAL_OR_EPITOPE_ONLY_R3`。双计算生成器产生相同
肽仍可记录为 peptide-level reconstruction consensus，但不再等价于完整
ORF 共识。

## 运行示例

```bash
PYTHONPATH=src python -m neoag.splice.cli build \
  --sample-id SAMPLE01 \
  --outdir results/SAMPLE01/splice_provenance \
  --junctions tumor.regtools.enriched.tsv \
  --junction-coordinate-system regtools_annotated \
  --junction-qc-policy complete \
  --min-junction-unique-reads 3 \
  --min-junction-unique-fragment-starts 2 \
  --min-junction-overhang 10 \
  --min-junction-mapping-quality 20 \
  --max-junction-multimapping-fraction 0.20 \
  --min-junction-tumor-psi 0.05 \
  --tool-version RegTools=VERSION \
  --strict
```

## 验证

P0 专项回归测试：

```bash
PYTHONPATH=src pytest -q \
  tests/test_splice_p0_fixes.py \
  tests/test_splice_v050_layer.py \
  tests/test_splice_v051_three_chains.py
```

这些测试覆盖低 reads、缺失 QC、完整 QC PASS、事件内错误借证、
pVACsplice 未评估/显式残基位置、正常背景状态一致性、fallback 命名及
partial/epitope-only R3 上限。
