# NeoAg Splice P0 修复验证报告

## 交付版本

- Splice schema：`0.5.2-p0`
- 基础项目：NeoAg Event Pipeline v0.5.1 / project version 0.5.2
- 修复目录：`neo-na0707_splice_p0_fixed`

## 已实现范围

1. 正式 junction read QC 表及可配置阈值；
2. pVACsplice 残基级新颖性和 junction 边界解析，禁止默认写 true；
3. peptide-origin 的 `required_junction_ids` 及逐来源全量支持要求；
4. evidence chain 与 consensus 共用正常背景检出状态函数；
5. fallback 事件更名为 `JUNCTION_ONLY_UNCLASSIFIED`；
6. partial、peptide-only、epitope-only 产物强制 R3 上限。

## 专项验证

执行：

```bash
PYTHONPATH=/tmp/neoag_p0_testdeps:src python -m pytest -q \
  tests/test_splice_p0_fixes.py \
  tests/test_splice_v050_layer.py \
  tests/test_splice_v051_three_chains.py
```

结果：`44 passed`。

另外完成：

- `python -m compileall -q src/neoag/splice`：通过；
- `python -m json.tool resources/splice_provenance_v052_p0.schema.json`：通过；
- `bash -n scripts/run_splice_provenance_v051.sh`：通过；
- CLI help 中已暴露全部 junction QC 参数。

## 全仓快速测试

执行 `python -m pytest -q` 的结果为：

- 697 passed
- 113 skipped
- 4 failed

4个失败项不位于本次修改的 Splice P0 模块：

- `tests/test_production_runner.py` 两个状态断言；
- `tests/test_reports_parallel_consensus.py` 一个固定英文文案断言；
- `tests/test_splice_v051_tool_wrappers.py` 的 EasyQuant wrapper 伪 BAM/外部检查。

因此，本报告只声明 P0 专项和直接相关回归测试通过，不声明全仓测试全绿。

## 关键行为变化

- 默认 `junction_qc_policy=complete`。旧的仅含坐标和 read count 的 TSV 会
  产生 `INCOMPLETE`，不会贡献正式 RNA 支持；
- 如必须回放历史数据，需要显式指定 `--junction-qc-policy reads_only`；
- 没有 residue position/junction offset 的 pVACsplice 表位会保留，但新颖性
  为 `UNASSESSED` 并受 R3 上限；
- 历史 partial/peptide-only 候选的最终 tier 可能由 R2 降至 R3；
- 未分类单 junction 事件 ID 会因 event type 变化而改变，不能与旧版
  `NOVEL_JUNCTION` fallback ID 混用。
