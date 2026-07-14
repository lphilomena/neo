---
name: neoag-concept-explainer
description: Generate bounded explanations for APPM, CCF, HLA LOH, NetMHCpan, minigene, ELISpot and validation concepts.
category: C - 审阅/报告/实验设计型 Skills：解释结果、生成报告、设计实验验证和患者沟通材料
risk_level: LOW
approval_required: false
---

# neoag-concept-explainer

## 目标

术语解释与报告注释

## 什么时候使用
- 用户询问术语解释
- 报告需要插入概念框

## 什么时候不要使用
- 不要替代具体结果分析

## 必需输入
- `concept`

## 可选输入
- `audience`

## 输出
- `concept_explanation.md`

## 运行示例

```bash
neoag-skill run neoag-concept-explainer --outdir work/neoag-concept-explainer --dry-run
```

## 边界
- Skill 不承担临床决策；不得判断患者是否适合治疗或推荐临床用药。
- 缺失证据只能标记为 missing/unassessed，不能解释为阴性结果。
- 高风险写入、HPC 提交、安装工具、下载参考库、删除或覆盖文件必须经过 human approval。
- Skill 目录不包含患者 BAM/FASTQ/VCF、大型参考库、VEP cache、NetMHCpan license、LOHHLA reference 或大型 conda env。
