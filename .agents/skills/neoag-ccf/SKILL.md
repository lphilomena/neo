---
name: neoag-ccf
description: Estimate or normalize Cancer Cell Fraction and clonality status for events, with confidence flags.
category: B - 公共证据分析型 Skills：对所有入口共用的 HLA、表达、CCF、APPM、安全和排序证据层进行标准化分析
risk_level: LOW
approval_required: false
---

# neoag-ccf

## 目标

CCF/clonality 估计与 modifier 输出

## 什么时候使用
- 需要 CCF/clonality 分层
- 需要计算或解释 CCF modifier

## 什么时候不要使用
- RNA-only fusion 不应伪造 DNA CCF

## 必需输入
- `event_table_or_ranked_peptides`

## 可选输入
- `purity`
- `cnv_segments`
- `sample_id`

## 输出
- `ccf_lite.tsv`
- `ccf_input_qc.tsv`
- `ccf_modifier_summary.tsv`

## 运行示例

```bash
neoag-skill run neoag-ccf --outdir work/neoag-ccf --dry-run
```

## 边界
- Skill 不承担临床决策；不得判断患者是否适合治疗或推荐临床用药。
- 缺失证据只能标记为 missing/unassessed，不能解释为阴性结果。
- 高风险写入、HPC 提交、安装工具、下载参考库、删除或覆盖文件必须经过 human approval。
- Skill 目录不包含患者 BAM/FASTQ/VCF、大型参考库、VEP cache、NetMHCpan license、LOHHLA reference 或大型 conda env。
