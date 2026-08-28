---
name: neoag-hla-typing-loh
description: Normalize HLA typing and HLA LOH outputs from OptiType, SpecHLA, HLA-LA/HD and LOHHLA into consensus tables.
category: B - 公共证据分析型 Skills：对所有入口共用的 HLA、表达、CCF、APPM、安全和排序证据层进行标准化分析
risk_level: LOW
approval_required: false
---

# neoag-hla-typing-loh

## 目标

HLA typing / HLA LOH 共识与 peptide-level HLA loss flags

## 什么时候使用
- 需要标准化 HLA 分型
- 需要判断 restricting HLA 是否 LOH

## 什么时候不要使用
- 只需要解释已有 ranking 差异，不需要更新 HLA 状态

## 必需输入
- `hla`

## 可选输入
- `hla_loh`
- SpecHLA LOH 执行需要 SpecHLA 自身生成的 `hla.result.txt`、`HLA_*_freq.txt`、肿瘤纯度和倍性；外部分型文件不能代替频率文件
- `ranked_peptides`
- `sample_id`

## 输出
- `hla_typing.normalized.tsv`
- `hla_typing.standardized.tsv`
- `hla_typing_consensus.tsv`
- `recommended_hla.txt`
- `hla_loh_consensus.tsv`
- `recommended_hla_loh.tsv`
- `restricting_hla_peptide_flags.tsv`
- `hla_review.md`

## 运行示例

```bash
neoag-skill run neoag-hla-typing-loh --outdir work/neoag-hla-typing-loh --dry-run
```

## 边界
- Skill 不承担临床决策；不得判断患者是否适合治疗或推荐临床用药。
- 缺失证据只能标记为 missing/unassessed，不能解释为阴性结果。
- 高风险写入、HPC 提交、安装工具、下载参考库、删除或覆盖文件必须经过 human approval。
- Skill 目录不包含患者 BAM/FASTQ/VCF、大型参考库、VEP cache、NetMHCpan license、LOHHLA reference 或大型 conda env。

## 底层工具
- OptiType
- SpecHLA
- LOHHLA
- HLA-LA（DNA BAM 图谱分型复核；已安装时可直接运行）
- HLA-HD optional

## HLA-LA DNA 分型执行

优先使用正常样本的、GRCh38 比对且已建立索引的短读长 BAM。HLA-LA 需要完整的
`PRG_MHC_GRCh38_withIMGT` 图谱；仅有运行时容器不能视为可用安装。

```bash
bash scripts/verify_hla_la_container.sh

bash scripts/run_hla_la_sample.sh \
  --bam NORMAL.GRCh38.bam \
  --sample-id SAMPLE_NORMAL \
  --graph "$HLALA_GRAPH" \
  --threads 8 \
  --outdir work/hla_la/SAMPLE_NORMAL
```

成功运行必须同时得到非空的 `R1_bestguess_G.txt`（或兼容版本的
`R1_bestguess.txt`）、`run_metadata.tsv` 和 `.complete`。只有进程启动、目录存在或
图谱存在均不能视为分型成功。肿瘤 BAM 可用于交叉检查，但不能替代正常样本的
胚系 HLA 分型主结果。

## OptiType DNA 分型执行

```bash
bash scripts/run_optitype_sample.sh \
  --bam NORMAL.GRCh38.bam \
  --sample-id SAMPLE_NORMAL \
  --threads 4 \
  --outdir work/hla_typing/optitype
```

## SpecHLA DNA 分型执行

```bash
bash scripts/run_spechla_sample.sh \
  --bam NORMAL.GRCh38.bam \
  --sample-id SAMPLE_NORMAL \
  --threads 5 \
  --outdir work/hla_typing/spechla
```

三工具结束后运行 `neoag-hla-typing-compare`，必须同时输出 2-field 与
高分辨率一致性、`hla_typing.standardized.tsv`、`hla_typing_consensus.tsv` 和
`recommended_hla.txt`。OptiType 只参与 HLA-A/B/C 共识。

## SpecHLA LOH 执行

```bash
bash scripts/run_spechla_loh.sh \
  --sample-id SAMPLE \
  --typing-dir work/spechla/SAMPLE \
  --purity 0.5 \
  --ploidy 2.0 \
  --outdir work/spechla_loh/SAMPLE
```

同型位点或杂合 SNP 数低于阈值的位点必须标记为 `unassessed`，不得解释为未发生 LOH。

LOHHLA 与 SpecHLA 归一化后使用 `scripts/build_hla_loh_consensus.py`。最终表只合并
HLA-A/B/C，并明确保留 `LOST`、`RETAINED`、`UNASSESSED` 和 `CONFLICT`；供排序
使用的 `recommended_hla_loh.tsv` 只包含推荐 LOST 等位基因。
