---
name: neoag-hpc-runner
description: Create Slurm/SGE/PBS job manifests and submit only after explicit approval.
category: D - 工程治理/执行控制型 Skills：输入质控、环境健康检查、全流程编排、发布审计和受控执行
risk_level: HIGH
approval_required: true
---

# neoag-hpc-runner

## 目标

HPC dry-run/job wrapper

## 什么时候使用
- 需要在 HPC 上运行 pipeline-full 或重型工具

## 什么时候不要使用
- 默认不得直接提交，应先 dry-run

## 必需输入
- `job_manifest`

## 可选输入
- `无`

## 输出
- `hpc_dry_run.sh`
- `hpc_job_manifest.json`
- `hpc_submission.md`

## 运行示例

```bash
neoag-skill run neoag-hpc-runner --outdir work/neoag-hpc-runner --dry-run
```

## 边界
- Skill 不承担临床决策；不得判断患者是否适合治疗或推荐临床用药。
- 缺失证据只能标记为 missing/unassessed，不能解释为阴性结果。
- 高风险写入、HPC 提交、安装工具、下载参考库、删除或覆盖文件必须经过 human approval。
- Skill 目录不包含患者 BAM/FASTQ/VCF、大型参考库、VEP cache、NetMHCpan license、LOHHLA reference 或大型 conda env。
