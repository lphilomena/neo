# open-neo-run 与 open-neo-review 的关系

`open-neo-run` 和 `open-neo-review` 是前后衔接、职责分离的两个宏 Skill。

```mermaid
flowchart LR
    A["原始数据或已有结果"] --> B["open-neo-run<br/>计算与证据整合"]
    B --> C["Evidence-consensus<br/>R1-R4 正式排序"]
    C --> D["open-neo-review<br/>审阅与实验设计"]
    D --> E["患者报告、技术报告<br/>实验候选清单"]
```

## open-neo-run

`open-neo-run` 负责运行分析并形成正式结果：

- 识别 VCF、BAM、FASTQ、表达量及已有工具结果。
- 根据输入数据、可用工具和参考资源自动编排可运行模块。
- 运行或复用 HLA 分型、纯度/CNV、HLA LOH、RNA、Fusion、Splice 和呈递预测等结果。
- 建立综合工具与证据表 `all_tool_results.tsv`。
- 生成 weighted baseline 和 Evidence-consensus 两套并行排序。
- 保存运行清单、证据冲突、验证计划及可追溯信息。

主要输出包括：

- `all_tool_results.tsv`
- `ranked_peptides.weighted_baseline.tsv`
- `ranked_peptides.evidence_consensus.tsv`
- `ranked_events.evidence_consensus.tsv`
- `validation_plan.tsv`
- `run_manifest.json`

`open-neo-run` 决定候选的正式证据等级和 Pipeline 排名。

为避免产生两份患者报告，`open-neo-run` 默认只生成 Pipeline 技术报告：

- `reports/evidence_report.technical.html`

底层 `neoag-v03 run-full` 为兼容既有调用，仍支持通过 `--reports patient,technical` 显式生成患者报告；该患者报告会标记为 Pipeline 结果快照，不作为最终患者报告。

## open-neo-review

`open-neo-review` 负责读取、检查和解释 `open-neo-run` 的结果，不重新运行生信 Pipeline：

- 检查必需结果是否完整且属于同一次运行。
- 以事件级 Evidence-consensus 排序作为主要审阅入口。
- 比较 weighted baseline 与 Evidence-consensus 的差异。
- 按事件、单倍型和重复窗口去重，避免同一事件占据过多实验名额。
- 在不修改 Pipeline 等级的前提下生成审阅和实验设计字段。
- 生成第一批研究验证集合，以及短肽、长肽、minigene 和 targeted RNA 验证计划。
- 调用项目正式患者报告生成器，输出唯一的 `reports/patient_report.html`。
- 生成技术报告和一页式研究摘要等辅助结果。

新增的审阅字段主要包括：

- `review_status`
- `review_reason`
- `experiment_priority`
- `recommended_validation`

`open-neo-review` 不会修改以下正式 Pipeline 字段：

- `pipeline_r_grade`
- `pipeline_event_rank`

最终对外患者报告统一为：

- `open-neo-review` 输出的 `reports/patient_report.html`

## 推荐执行顺序

```bash
open-neo run \
  --sample-manifest configs/sample.yaml \
  --mode execute \
  --approved \
  --outdir results/CASE001

open-neo review \
  --result-dir results/CASE001 \
  --top-n 12 \
  --reports patient,technical,onepage \
  --outdir reviews/CASE001
```

第一步建立综合证据和正式排序，第二步完成结果审阅、实验设计和报告生成。

## 前置条件与失败边界

`open-neo-review` 至少需要以下文件：

- `run_manifest.json`
- `ranked_events.evidence_consensus.tsv`
- `ranked_peptides.evidence_consensus.tsv`
- `ranked_peptides.weighted_baseline.tsv`
- `all_tool_results.tsv`
- `validation_plan.tsv`

如果缺少 `ranked_events.evidence_consensus.tsv`，`open-neo-review` 返回 `NEEDS_RANKING`，不会静默改用 weighted baseline 的 Top 候选。

## 一句话总结

- `open-neo-run` 回答：候选是什么、证据是否完整、正式排序是什么。
- `open-neo-review` 回答：为什么升降级、先验证哪些事件、采用什么实验，以及如何向患者和研究人员解释。

两者均属于研究性分析流程。计算候选不等于已经确认的新抗原，也不等于确定治疗方案。
