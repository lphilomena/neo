---
name: neoag-ranking-compare
description: Compute overlap, rank shifts and interpretation between NetMHCpan42 and recommendation ranking tables.
category: C - 审阅/报告/实验设计型 Skills：解释结果、生成报告、设计实验验证和患者沟通材料
risk_level: LOW
approval_required: false
---

# neoag-ranking-compare

## 目标

NetMHCpan42 与综合推荐排序比较

## 什么时候使用
- 用户问两个排序文件差异
- 需要解释 NetMHCpan42 top fusion 为什么被综合降权

## 什么时候不要使用
- 需要生成 ranking 时，使用 neoag-ranking

## 必需输入
- `recommendation`
- `netmhcpan42`

## 可选输入
- `无`

## 输出
- `ranking_compare_report.md`
- `topn_overlap.tsv`
- `rank_shift.tsv`

## 运行示例

```bash
neoag-skill run neoag-ranking-compare --outdir work/neoag-ranking-compare --dry-run
```

## 边界
- Skill 不承担临床决策；不得判断患者是否适合治疗或推荐临床用药。
- 缺失证据只能标记为 missing/unassessed，不能解释为阴性结果。
- 高风险写入、HPC 提交、安装工具、下载参考库、删除或覆盖文件必须经过 human approval。
- Skill 目录不包含患者 BAM/FASTQ/VCF、大型参考库、VEP cache、NetMHCpan license、LOHHLA reference 或大型 conda env。
