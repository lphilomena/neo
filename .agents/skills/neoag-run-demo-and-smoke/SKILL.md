---
name: neoag-run-demo-and-smoke
description: Run or plan Project B run-demo, pytest and optional Nextflow smoke tests.
category: D - 工程治理/执行控制型 Skills：输入质控、环境健康检查、全流程编排、发布审计和受控执行
risk_level: MEDIUM
approval_required: false
---

# neoag-run-demo-and-smoke

## 目标

demo/pytest/Nextflow smoke

## 什么时候使用
- release 验收
- 新环境最小可运行性测试

## 什么时候不要使用
- 不要在没有确认的生产目录覆盖结果

## 必需输入
- `project_root`

## 可选输入
- `无`

## 输出
- `smoke_test_report.md`
- `demo_output_manifest.tsv`

## 运行示例

```bash
neoag-skill run neoag-run-demo-and-smoke --outdir work/neoag-run-demo-and-smoke --dry-run
```

## 边界
- Skill 不承担临床决策；不得判断患者是否适合治疗或推荐临床用药。
- 缺失证据只能标记为 missing/unassessed，不能解释为阴性结果。
- 高风险写入、HPC 提交、安装工具、下载参考库、删除或覆盖文件必须经过 human approval。
- Skill 目录不包含患者 BAM/FASTQ/VCF、大型参考库、VEP cache、NetMHCpan license、LOHHLA reference 或大型 conda env。
