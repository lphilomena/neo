---
name: neoag-scoring
description: Run neoag_v03's peptide ranking/scoring pipeline on a variant_peptides_annotated.tsv (extract-variant-peptides wide annotation table, already carrying NetMHCpan/MHCflurry/PRIME/BigMHC/IEDB/NetMHCstabPan columns) plus a matching raw_events.tsv, and interpret the resulting ranked_peptides.tsv for the user. Use this whenever the user asks to rank, score, prioritize, or sort neoantigen/variant peptide candidates, mentions "variant_peptides_annotated.tsv" or "raw_events.tsv" by name, or asks "which of these peptides should I validate first" / "trust this ranking" / "why is this peptide ranked here". Always run the four-step workflow below in order (explain input fields -> explain weights from the live TOML -> run the scorer -> interpret the output) rather than jumping straight to running the scorer.
---

# neoag 排序与结果解读

## 总览

这个skill覆盖四步，**必须按顺序执行，不要跳过前两步直接跑排序**：

1. 解释 `variant_peptides_annotated.tsv` 里的字段属于哪个评价维度，让用户能自己判断要不要信任后续排序
2. 解释排序用的权重和阈值——**这些数字必须现场从用户实际使用的 profile TOML 文件里读出来再讲给用户**，禁止在对话里凭记忆背诵一份可能已经过时的数字
3. 运行 `scripts/score_from_annotated_table.py` 生成 `ranked_peptides.tsv` / `ranked_events.tsv`
4. 打开生成的排序表，评估这次排序结果本身有多可信，给出下一步建议

---

## 第一步：解释输入文件——字段对应哪个评价维度

`variant_peptides_annotated.tsv` 是 `extract-variant-peptides`（已附加NetMHCpan/MHCflurry等预测）产出的宽表，每一行是一个"肽段×HLA"组合。在跑任何排序之前，先用 `view`/`head` 看一眼这个文件表头，然后按下表把列名讲给用户，**并明确告诉用户：如果不信任后面的综合排序，完全可以直接在这张表里按这些列自己筛选**（比如"只看 netmhcpan_mt_rank_el < 0.5 且 mhcflurry_mt_presentation_score > 0.8 的行"）。

| 评价维度 | 对应字段 | 怎么看 |
|---|---|---|
| **HLA结合** | `netmhcpan_mt_ic50`、`netmhcpan_mt_rank_ba`、`mhcflurry_mt_affinity`、`mhcflurry_mt_affinity_percentile` | ic50/rank/percentile 都是**越小越强**；ic50单位nM，<500通常算能结合，<50算强结合 |
| **HLA呈递** | `netmhcpan_mt_rank_el`、`mhcflurry_mt_processing_score`、`mhcflurry_mt_presentation_score`、`netmhcstabpan_score`/`netmhcstabpan_rank` | EL rank越小越好（呈递配体排名）；processing/presentation score越接近1越好；stabpan_score是复合物半衰期(小时)，越大越稳定 |
| **免疫原性** | `prime_score`、`bigmhc_im_score`、`iedb_immunogenicity_score` | 都是T细胞识别概率预测，0-1，越大越好；三个来自不同模型，不一定同时有值 |
| **安全性/自身相似性** | `wildtype_peptide`、`netmhcpan_wt_*`、`mhcflurry_wt_*`、`in_normal_proteome` | 对比突变肽和野生型肽的结合力差异；`in_normal_proteome=yes` 说明这条肽序列本身在正常蛋白组里能找到，有诱发自身免疫的风险，要重点关注 |
| **表达证据** | `gene_tpm`、`transcript_tpm`、`event_expression`、`expression_status`、`rna_vaf`、`rna_alt_reads`、`rna_depth`、`rna_allele_evidence_status` | TPM是基因/转录本表达量；rna_alt_reads是RNA层面实际测到携带这个突变的read数——DNA测到不代表真的转录出来了，这几列是独立佐证 |
| **变异检出可信度** | `vaf`、`tumor_depth`、`tumor_alt_count` | DNA测序支持这个变异的证据强度，不是"致病性"，是"这个变异是不是真的存在" |
| **肽段结构/来源** | `peptide_length`、`mutation_position_in_peptide`、`multi_aa_flag`、`crosses_junction`、`contains_novel_aa`、`consequence` | 描述这条肽是怎么产生的（错义/移码/融合/剪接），不直接决定分数高低，但影响后面走哪种验证方案 |

---

## 第二步：解释排序规则——权重必须现场读TOML，不要背数字

neoag的打分是"L3十个维度加权平均 → 综合分" 再 "× 安全性/APPM/呈递门控/克隆性/免疫逃逸 五个独立乘子" 得到每条肽的 `efficacy_score`。**具体权重和阈值因 profile（疾病场景配置）而异，且这些配置文件本身可能被用户改过**，所以必须现场读取以下文件并把读到的数字直接讲给用户，不要用训练记忆里的数字：

```bash
# 先确认这次运行用的是哪个profile（跑第三步时会指定），然后：
cat <neo仓库>/profiles/<profile名>.toml
```

读取后按下面这个话术模板讲给用户（把 `{...}` 替换成实际读到的数字）：

> 这次用的是 `{profile名}.toml`。十个维度里权重最高的是 `{从[l3_weights]区块读出权重最大的两三项及其数值}`，加起来占了综合分的大头；`normal_tissue_safety`/`apm_integrity` 这两项权重是 `{读到的数值}`（如果是0，说明安全性和呈递通路完整性不参与加权平均，而是作为下面的独立乘子直接生效，避免重复计算）。
>
> 综合分算出来后，还会依次乘上这几个独立乘子：安全性乘子（PASS/CAUTION/FAIL对应 `{从代码/toml读到的值}`）、呈递门控乘子（结合等级、EL rank、稳定性任一项不达标就打 `{gates.failure_multiplier}` 折）、克隆性乘子（`{ccf相关阈值}`）、免疫逃逸乘子（如果这条肽依赖的HLA丢失或抗原呈递通路基因缺陷，最低可以打到 `{immune_escape.multipliers里最低的值}` 折甚至清零）。这几个乘子里**任何一个接近0，都会让这条肽不管前面综合分多高，最终排序都会掉到最后**。

如果用户问"为什么某条肽排名靠后"，同样要去读它的详细分数列（第四步会生成），逐项对照这份权重表解释，不要泛泛而谈。

---

## 第三步：运行排序

先检查前置条件：

```bash
# 1. 确认两个文件的关联键对得上——宽表的variant_key必须能在raw_events.tsv的event_id里找到
cut -f1 raw_events.tsv | tail -n +2 | sort -u > /tmp/eids.txt
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="variant_key") vk=i} NR>1{print $vk}' variant_peptides_annotated.tsv | sort -u > /tmp/vks.txt
comm -23 /tmp/vks.txt /tmp/eids.txt | wc -l   # 这个数字如果很大，说明两个文件对不上，见下面"数据不匹配"处理
```

然后运行打分脚本：

```bash
python3 scripts/score_from_annotated_table.py \
    --annotated variant_peptides_annotated.tsv \
    --raw-events raw_events.tsv \
    --profile <profile名，不指定则用default> \
    --neoag-root <neo仓库根目录路径> \
    --outdir <输出目录>
```

脚本会打印类似：
```
Scored 1655 peptide x HLA rows across 45 events.
WARNING: skipped N rows whose variant_key had no matching event_id in --raw-events.
Wrote: <outdir>/ranked_peptides.tsv
Wrote: <outdir>/ranked_events.tsv
```

**数据不匹配处理**：如果 `skipped` 的数量占比很高（比如前置检查里`comm`命令数字接近variant_key总数），如实告诉用户"这两个文件对不上，raw_events.tsv 里的事件数明显少于宽表里的变异数，可能是 raw_events.tsv 经过了预筛选，或者两个文件不是同一批次生成的"，不要悄悄跳过不提。

这个脚本做的事（可以照实讲给用户，不用回避细节）：对宽表每一行，直接复用 `netmhcpan_mt_rank_ba`/`netmhcpan_mt_rank_el`/`mhcflurry_mt_*` 等已经跑出来的真实预测值去算 `binding_evidence_score`/`presentation_evidence_score`（公式和第二步讲的一致，权重同样来自profile），事件层的安全性/克隆性/APPM等字段直接读 `raw_events.tsv` 里已有的值，然后调用neoag打分核心函数算出每条肽的 `efficacy_score` 并排序。免疫逃逸乘子（`escape_multiplier`）由于这条路径没有单独喂VEP/CNV/HLA-LOH文件，默认按1.0（不生效）处理——**如果用户的场景对免疫逃逸风险比较敏感（比如做过治疗、可能有HLA杂合性丢失），要提醒这一项没有被评估**。

---

## 第四步：解读排序表，评估可信度，给建议

打开 `ranked_peptides.tsv`，做以下几件事：

### 1. 看整体分布，判断这次打分是否"有区分度"
```python
import csv
from collections import Counter
rows = list(csv.DictReader(open("ranked_peptides.tsv"), delimiter="\t"))
print(Counter(r["final_priority"] for r in rows))
```
如果几乎所有候选都挤在同一档（比如99%都是D），说明这批候选整体证据都偏弱，排序更多是矬子里拔将军，不代表真的有"强候选"；如果A/B档有一批数量合理的候选，说明这批数据本身有信号。

### 2. 抽查关键字段是不是"占位符恒定值"——这是判断可信度最关键的一步
读几行 `ranked_events.tsv`，看 `driver_relevance`、`tumor_specificity` 这两列：
- 如果**全部事件的这两列都是同一个数字**（比如全是0.0，或全是0.7/0.3/0.5这种整齐的数），说明 `raw_events.tsv` 是用旧版本的适配器生成的，这两项实际上没有真实计算，只是占位符，**排序结果里事件层的证据强度部分参考价值有限**，需要提醒用户：想要真实的驱动基因/组织特异性证据，得用最新版适配器重新生成 `raw_events.tsv`（接入了OncoKB/GTEx真实数据源），而不是直接用手头这份文件。
- 如果这两列的数值有明显差异（不同基因给出不同分数），说明是真实计算出来的，可以正常参考。

同样检查肽段表里的 `self_similarity_score` 是否全部恒为0——如果是，说明"这条肽和自身正常蛋白有多像"这个安全检查实际上没有真实跑过，所有候选在这一项上都被默认判定为"安全"，**这不等于真的安全，只是没评估**，需要在解读时明确指出这个局限，别让用户误以为高排名=已验证安全。

### 3. 结合表达证据判断排序前几名是否"真的在表达"
对照 `event_expression`/`rna_alt_reads` 这两列，如果排名很靠前的候选这两项都是0或很低，要单独提醒用户："这条候选虽然结合/呈递预测分数很高，但没有RNA层面的直接表达证据，存在'预测出来的肽段其实根本没被转录'的风险"。

### 4. 给出下一步建议（按用户实际场景调整，不要机械套模板）
- 优先推荐排名靠前、**同时**表达证据扎实（`rna_alt_reads`>0）、安全性PASS、免疫逃逸未评估但对应HLA没有已知问题的候选，作为湿实验验证（ELISpot/四聚体）的第一批对象
- 如果A/B档候选数量太少，建议检查是不是筛选阈值（呈递门控里的EL rank/stabpan阈值）设得过严，可以回到第二步读到的TOML数值，跟用户讨论要不要针对这批数据放宽
- 如果发现`driver_relevance`/`tumor_specificity`/`self_similarity_score`是占位符恒定值，明确建议"用最新适配器重新生成一份 raw_events.tsv 再排序"，而不是直接拿现在这份排序结果做决策
- 如果用户的应用场景涉及免疫逃逸风险（既往治疗史、怀疑HLA丢失），建议补跑 `immune-escape` 相关步骤（需要VEP/CNV/HLA-LOH文件）后再排序，现在这份结果里 `escape_multiplier` 全部是1.0，等于没评估这个风险

---

## 常见追问处理

- **"为什么这条肽排在前面/后面"** → 打开 `ranked_peptides.tsv` 里这条肽对应的行，把 `l3_*_score` 各维度分、四个乘子（`safety_status`对应的`safety_multiplier`、`appm_multiplier`、`presentation_gate_multiplier`、`ccf_multiplier`）逐项列出来，对照第二步读到的真实权重逐项解释，不要只说"综合分较高/较低"。
- **"这个权重能不能改"** → 可以，改 `<neo仓库>/profiles/<profile名>.toml` 对应字段，然后重新跑第三步的脚本；提醒用户改的是"用户自定义profile"还是"default.toml"（后者会影响所有没显式覆盖该字段的profile，因为 `load_profile()` 会自动继承default.toml）。
- **"这个排序靠谱吗"** → 直接引用第四步第2、3点的检查结果回答，不要笼统地说"靠谱"或"不靠谱"。
