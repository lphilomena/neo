---
name: neoag-gateway-submit
description: Submit low/medium/high-risk skill or pipeline requests to NeoAg Gateway with approval controls.
category: D - 工程治理/执行控制型 Skills：输入质控、环境健康检查、全流程编排、发布审计和受控执行
risk_level: HIGH
approval_required: true
---

# neoag-gateway-submit

## 目标

Gateway 受控提交

## 什么时候使用
- Agent 需要通过 Gateway 受控执行任务

## 什么时候不要使用
- 不要直接执行 shell 或绕过 approval

## 必需输入
- `gateway_url`
- `task`

## 可选输入
- `无`

## 输出
- `gateway_job.json`
- `gateway_submission.md`

## 运行示例

```bash
neoag-skill run neoag-gateway-submit --outdir work/neoag-gateway-submit --dry-run
```

## 边界
- Skill 不承担临床决策；不得判断患者是否适合治疗或推荐临床用药。
- 缺失证据只能标记为 missing/unassessed，不能解释为阴性结果。
- 高风险写入、HPC 提交、安装工具、下载参考库、删除或覆盖文件必须经过 human approval。
- Skill 目录不包含患者 BAM/FASTQ/VCF、大型参考库、VEP cache、NetMHCpan license、LOHHLA reference 或大型 conda env。
