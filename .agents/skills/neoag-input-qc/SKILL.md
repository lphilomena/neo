---
name: neoag-input-qc
description: Inspect manifests/result directories and recommend workflow while listing missing inputs.
category: D - 工程治理/执行控制型 Skills：输入质控、环境健康检查、全流程编排、发布审计和受控执行
risk_level: LOW
approval_required: false
---

# neoag-input-qc

## 目标

输入状态检查与 workflow 推荐

## 什么时候使用
- 任何任务的第一步
- 用户问能不能跑、缺什么输入

## 什么时候不要使用
- 不能用 input-qc 的缺失信息直接做生物学阴性结论

## 必需输入
- `manifest_or_result_dir`

## 可选输入
- `无`

## 输出
- `input_status.json`
- `input_qc_report.tsv`
- `missing_inputs.tsv`

## 运行示例

```bash
neoag-skill run neoag-input-qc --outdir work/neoag-input-qc --dry-run
```

## 边界
- Skill 不承担临床决策；不得判断患者是否适合治疗或推荐临床用药。
- 缺失证据只能标记为 missing/unassessed，不能解释为阴性结果。
- 高风险写入、HPC 提交、安装工具、下载参考库、删除或覆盖文件必须经过 human approval。
- Skill 目录不包含患者 BAM/FASTQ/VCF、大型参考库、VEP cache、NetMHCpan license、LOHHLA reference 或大型 conda env。
