---
name: neoag-patient-report
description: Generate a patient-facing markdown/html/docx draft with research boundary statements.
category: C - 审阅/报告/实验设计型 Skills：解释结果、生成报告、设计实验验证和患者沟通材料
risk_level: MEDIUM
approval_required: false
---

# neoag-patient-report

## 目标

患者沟通版报告

## 什么时候使用
- 需要更新患者沟通版报告
- 需要解释 A/B/C/D、CCF、APPM、HLA LOH

## 什么时候不要使用
- 不要生成临床处方或疗效承诺

## 必需输入
- `ranked_peptides_or_summary`

## 可选输入
- `evidence_report`
- `ranking_compare`
- `appm_review`
- `ccf_review`

## 输出
- `patient_report.md`
- `patient_report.html`
- `patient_report.docx`

## 运行示例

```bash
neoag-skill run neoag-patient-report --outdir work/neoag-patient-report --dry-run
```

## 边界
- Skill 不承担临床决策；不得判断患者是否适合治疗或推荐临床用药。
- 缺失证据只能标记为 missing/unassessed，不能解释为阴性结果。
- 高风险写入、HPC 提交、安装工具、下载参考库、删除或覆盖文件必须经过 human approval。
- Skill 目录不包含患者 BAM/FASTQ/VCF、大型参考库、VEP cache、NetMHCpan license、LOHHLA reference 或大型 conda env。
