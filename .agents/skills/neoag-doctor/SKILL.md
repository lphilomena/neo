---
name: neoag-doctor
description: Run controlled-execution read-only health check for tools, references and release boundaries.
category: D - 工程治理/执行控制型 Skills：输入质控、环境健康检查、全流程编排、发布审计和受控执行
risk_level: LOW
approval_required: false
---

# neoag-doctor

## 目标

只读环境健康检查

## 什么时候使用
- 部署/HPC/新机器验证
- 用户问工具和参考库是否可用

## 什么时候不要使用
- 不要用 doctor 安装或修改工具

## 必需输入
- `project_root`

## 可选输入
- `tools_manifest`
- `reference_manifest`

## 输出
- `doctor_status.json`
- `doctor_summary.md`
- `blocking_issues.tsv`

## 运行示例

```bash
neoag-skill run neoag-doctor --outdir work/neoag-doctor --dry-run
```

## 边界
- Skill 不承担临床决策；不得判断患者是否适合治疗或推荐临床用药。
- 缺失证据只能标记为 missing/unassessed，不能解释为阴性结果。
- 高风险写入、HPC 提交、安装工具、下载参考库、删除或覆盖文件必须经过 human approval。
- Skill 目录不包含患者 BAM/FASTQ/VCF、大型参考库、VEP cache、NetMHCpan license、LOHHLA reference 或大型 conda env。
