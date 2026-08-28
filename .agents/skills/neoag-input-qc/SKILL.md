---
name: neoag-input-qc
description: Inventory case inputs and perform real GRCh38 BAM/BAI, VCF sample, sequencing QC, and tumor-normal identity checks before workflow selection.
category: D - 工程治理/执行控制型 Skills：输入质控、环境健康检查、全流程编排、发布审计和受控执行
risk_level: LOW
approval_required: false
---

# neoag-input-qc

## 目标

输入状态、样本配对、测序 QC、肿瘤-正常身份检查与 workflow 推荐

## 什么时候使用
- 任何任务的第一步
- 用户问能不能跑、缺什么输入

## 什么时候不要使用
- 不能用 input-qc 的缺失信息直接做生物学阴性结论

## 必需输入
- `input_dir`

## 可选输入
- `tumor_bam`
- `normal_bam`
- `somatic_vcf`
- `reference_fasta`
- `bam_matcher_loci`

## 输出
- `input_status.json`
- `input_inventory.tsv`
- `sample_pairing.tsv`
- `bam_qc.tsv`
- `vcf_qc.tsv`
- `input_qc_review.md`

## 运行示例

```bash
python scripts/run_case_input_qc.py \
  --input-dir /path/to/case \
  --tumor-bam /path/to/tumor.bam \
  --normal-bam /path/to/normal.bam \
  --somatic-vcf /path/to/somatic.pass.vcf.gz \
  --reference "$NEOAG_REFERENCE_FASTA" \
  --bam-matcher-loci "$BAM_MATCHER_LOCI" \
  --outdir results/input_qc
```

默认执行 `samtools quickcheck/header/flagstat/coverage` 和 BAM-matcher。只有显式
使用 `--quick` 时才跳过全量指标，跳过的深度、比对率、重复率或污染证据必须写为
`UNASSESSED`。文件名推断不能证明同一患者，生产运行必须以指纹结果为准。

## 边界
- Skill 不承担临床决策；不得判断患者是否适合治疗或推荐临床用药。
- 缺失证据只能标记为 missing/unassessed，不能解释为阴性结果。
- 高风险写入、HPC 提交、安装工具、下载参考库、删除或覆盖文件必须经过 human approval。
- Skill 目录不包含患者 BAM/FASTQ/VCF、大型参考库、VEP cache、NetMHCpan license、LOHHLA reference 或大型 conda env。
