---
name: neoag-scoring
description: Run neoag_v03's peptide ranking/scoring pipeline on a variant_peptides_annotated.tsv (extract-variant-peptides wide annotation table, already carrying NetMHCpan/MHCflurry/PRIME/BigMHC/IEDB/NetMHCstabPan columns) plus a matching raw_events.tsv, and interpret the resulting ranked_peptides.tsv for the user. Use this whenever the user asks to rank, score, prioritize, or sort neoantigen/variant peptide candidates, mentions "variant_peptides_annotated.tsv" or "raw_events.tsv" by name, points at a neoag pipeline output directory, or asks "which of these peptides should I validate first" / "trust this ranking" / "why is this peptide ranked here". Always run the five-step workflow below in order (inventory the directory -> explain input fields -> explain weights from the live TOML -> run the scorer -> interpret the output) rather than jumping straight to running the scorer.
---

# neoag 排序与结果解读

## 总览

这个skill覆盖五步（第零到第四步），**必须按顺序执行，不要跳过前几步直接跑排序**：

0. 先清点用户给的目录/文件里实际有哪些辅助证据、这份neo仓库代码本身修复到了哪个版本，判断这次能不能完成"全套真实数据打分"，哪些会退化成默认值——**如果发现已经有现成的排序表，直接进入第四步解读，不要重新跑**
1. 解释 `variant_peptides_annotated.tsv` 里的字段属于哪个评价维度，让用户能自己判断要不要信任后续排序
2. 解释排序用的权重和阈值——**这些数字必须现场从用户实际使用的 profile TOML 文件里读出来再讲给用户**，禁止在对话里凭记忆背诵一份可能已经过时的数字
3. 运行 `scripts/score_from_annotated_table.py` 生成 `ranked_peptides.tsv` / `ranked_events.tsv`
4. 打开生成的排序表，评估这次排序结果本身有多可信，给出下一步建议

**极简版决策树**（想快速定位该看哪一节时用）：拿到路径 → 先查`provenance.v03.json`+清点目录(第零步) → 有现成排序表就直接第四步 → 没有就第一步开始走完整流程。用户问"为什么这条候选排这个名次" → 第二步(权重来源)+第四步。用户问"这两张表对不上" → 下面"事件排序 vs 肽段排序"这一节。

---

## 第零步：先清点目录——判断这次能补全到什么程度

### 0.1 先看有没有 `provenance.v03.json`，有就优先信它

真实的neoag产出目录根目录下通常有一个 `provenance.v03.json`，专门记录这次运行用的是哪个代码版本/commit、跑了哪些步骤、用了哪个profile。**如果这个文件存在，先读它**，这是判断"这份排序结果基于什么版本跑出来的"最直接、最权威的信息源，比后面靠猜数据模式准确得多：

```bash
cat <run_dir>/provenance.v03.json
```
读到里面记录的版本/commit信息后，直接告诉用户"这次运行是用XX版本的代码跑的"，不需要再靠下面0.5的间接推断。**如果这个文件不存在，才需要走0.5那套"猜"的方法。**

### 0.2 目录结构清点

真实的neoag产出目录远不止一张宽表，通常长这样（示例，实际子目录不一定全部存在）：

```
<run_dir>/
├── provenance.v03.json  ← 优先读这个，见0.1
├── scoring/            ← 如果这里已经有 ranked_peptides.v03.tsv，说明已经跑过完整排序了
│   ├── ranked_peptides.v03.tsv
│   └── ranked_events.v03.tsv
├── parsed/              ← 干净的 raw_events.tsv / raw_peptides.tsv（打分最基础的输入）
├── presentation/         ← presentation_evidence.tsv 等真实呈递证据（已经算好binding/presentation分数）
├── clonality/            ← ccf_2.tsv / ccf_lite.tsv 真实克隆性证据
├── immune_escape/        ← peptide_escape_flags.tsv 真实免疫逃逸证据
├── appm/                 ← appm_summary.tsv 等APPM通路完整性证据
├── safety/               ← event_safety.tsv / peptide_safety.tsv
└── upstream/tools/       ← variant_peptides.annotated.tsv 宽表、原始VEP/netmhcpan/mhcflurry输出等
```

**拿到用户的路径后，先用 `scripts/score_from_annotated_table.py --run-dir <路径> --outdir <任意临时目录>` 跑一次清点**（不加`--force`时，这一步只会打印清点结果，不会真的重新打分，即使发现已有排序表也只是提示不会覆盖）。脚本会打印类似：

```
=== 目录清点结果 ===
  [found]    已有排序表 scoring/ranked_peptides.v03.tsv
  [missing]  真实免疫逃逸证据 immune_escape/peptide_escape_flags.tsv
  ...
```

根据清点结果分两种情况：

### 0.3 情况A：`scoring/ranked_peptides.v03.tsv` 已经存在
说明用户之前已经跑过neoag的排序命令了。**直接跳到第四步，对这份现成的表做解读，不要重新跑排序**。可以顺带告诉用户："检测到你之前已经跑过排序了，我直接读现成的结果"；如果用户明确要求换一组参数重新算，再用 `--force` 重新跑（这时候正常走第一~三步流程）。

### 0.4 情况B：还没有排序表，需要新跑
把清点结果翻译成"能不能做全套真实打分"讲给用户，例如：

> 这批数据里，`parsed/raw_events.tsv` 和宽表都在，可以正常打分；`presentation/presentation_evidence.tsv` 和 `clonality/ccf_lite.tsv` 也在，说明结合呈递证据和克隆性都是真实计算出来的，不用退化成占位值；但没有找到 `immune_escape/peptide_escape_flags.tsv`，所以这次排序里**免疫逃逸这一项不会生效**（所有候选的escape_multiplier都按1.0算，不代表真的排除了逃逸风险）；也没有找到独立的 `appm/appm_summary.tsv`，不过 `raw_events.tsv` 本身如果已经带了 `appm_mhc_i_integrity`/`appm_mhc_ii_integrity` 这两列，脚本会直接用那份数据，这个不算缺失。

**必须在跑第三步之前，把"这次哪些维度是真实数据、哪些会退化成默认值"讲清楚，不要等用户看完排序结果自己发现。**

### 0.5 更可靠的补充检查：直接问代码，而不是猜数据模式

如果能拿到neo仓库本身的路径（不只是产出目录），**优先用这一步确认某个字段到底有没有真实计算逻辑，而不是靠后面第四步"看数值像不像占位符"这种间接推断**——数据模式有时会有歧义（比如"未收录基因的中性默认值"和"完全没修复的占位符"可能长得很像），但代码存在与否是确定性的：

```bash
# self_similarity_score 有没有真实计算(而不是恒为0的占位符)
test -f <neo仓库>/src/neoag_v03/self_similarity.py && echo "有真实自身相似性计算模块" || echo "没有，self_similarity_score仍是占位符"

# driver_relevance 有没有接入真实OncoKB查表
grep -q "_parse_oncokb_gene_list" <neo仓库>/src/neoag_v03/driver_gene_db.py && echo "支持OncoKB真实格式" || echo "可能还是旧版/占位符"

# tumor_specificity 有没有接入真实GTEx查表
grep -q "_parse_gtex_gct" <neo仓库>/src/neoag_v03/tumor_specificity.py && echo "支持GTEx官方GCT格式" || echo "可能还是旧版/占位符"

# L3权重里normal_tissue_safety/apm_integrity是否已经清零(避免重复计算)
grep -A1 "normal_tissue_safety" <neo仓库>/profiles/default.toml | head -2
```
把检测结果记在心里，后面第一步、第二步、第四步讲解字段含义和评估可信度时，直接依据这个确定性结论来说，不用每次都重新猜。

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

**优先用 `--run-dir` 模式**（如果第零步清点用的是一个完整目录）：

```bash
python3 scripts/score_from_annotated_table.py \
    --run-dir <neoag产出目录路径> \
    --profile <profile名，不指定则用default> \
    --neoag-root <neo仓库根目录路径> \
    --outdir <输出目录>
```
这个模式会自动发现 `parsed/raw_events.tsv`、`upstream/tools/variant_peptides*.tsv`、`presentation/presentation_evidence.tsv`、`clonality/ccf_2.tsv`(或`ccf_lite.tsv`)、`immune_escape/peptide_escape_flags.tsv`，**能找到的真实证据会直接使用，不会重算或退化成默认值**；找不到的会在第零步已经跟用户说明过了。

如果手头只有零散的两个文件（没有完整目录结构），退回基础模式：

```bash
python3 scripts/score_from_annotated_table.py \
    --annotated variant_peptides_annotated.tsv \
    --raw-events raw_events.tsv \
    [--presentation-evidence presentation_evidence.tsv]   # 有就加，没有就省略
    [--ccf ccf_lite.tsv]                                  # 同上
    [--escape-flags peptide_escape_flags.tsv]             # 同上
    --profile <profile名> \
    --neoag-root <neo仓库根目录路径> \
    --outdir <输出目录>
```

跑之前先检查两个文件的关联键对不对得上：

```bash
cut -f1 raw_events.tsv | tail -n +2 | sort -u > /tmp/eids.txt
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="variant_key") vk=i} NR>1{print $vk}' variant_peptides_annotated.tsv | sort -u > /tmp/vks.txt
comm -23 /tmp/vks.txt /tmp/eids.txt | wc -l   # 这个数字如果很大，说明两个文件对不上，见下面"数据不匹配"处理
```

脚本会打印类似：
```
=== 目录清点结果 ===
  ...
Using real clonality evidence from clonality/ccf_lite.tsv (overrides raw_events.tsv's own ccf columns).
Using real presentation evidence from presentation/presentation_evidence.tsv (skipping wide-table recompute where matched).
Using real immune-escape evidence from immune_escape/peptide_escape_flags.tsv.

Scored 1655 peptide x HLA rows across 45 events.
  1230/1655 peptides matched real presentation evidence (rest fell back to wide-table recompute).
  890/1655 peptides matched real immune-escape evidence.
WARNING: skipped N rows whose variant_key had no matching event_id in --raw-events.
Wrote: <outdir>/ranked_peptides.tsv
Wrote: <outdir>/ranked_events.tsv
```
**如果"matched real xxx evidence"的比例很低（比如远小于总数），要如实告诉用户"这份sidecar文件覆盖的候选没有想象中多，剩下的还是走了wide-table重算/默认值1.0"，不要因为脚本用了真实文件就默认全部候选都是真实数据。**

**数据不匹配处理**：如果 `skipped` 的数量占比很高（比如前置检查里`comm`命令数字接近variant_key总数），如实告诉用户"这两个文件对不上，raw_events.tsv 里的事件数明显少于宽表里的变异数，可能是 raw_events.tsv 经过了预筛选，或者两个文件不是同一批次生成的"，不要悄悄跳过不提。

这个脚本做的事（可以照实讲给用户，不用回避细节）：优先使用 `--run-dir`/`--presentation-evidence`/`--ccf`/`--escape-flags` 指向的真实sidecar证据；没有sidecar覆盖的候选，就对宽表每一行直接复用 `netmhcpan_mt_rank_ba`/`netmhcpan_mt_rank_el`/`mhcflurry_mt_*` 等已经跑出来的真实预测值去算 `binding_evidence_score`/`presentation_evidence_score`（公式和第二步讲的一致，权重同样来自profile），事件层的安全性/克隆性/APPM等字段直接读 `raw_events.tsv` 里已有的值，然后调用neoag打分核心函数算出每条肽的 `efficacy_score` 并排序。**免疫逃逸乘子（`escape_multiplier`）只有在提供了 `--escape-flags`（或`--run-dir`能自动发现`immune_escape/peptide_escape_flags.tsv`）时才会生效，否则默认按1.0（不生效）处理**——如果用户的场景对免疫逃逸风险比较敏感（比如做过治疗、可能有HLA杂合性丢失），要提醒这一项是否被真的评估过。

---

## 第四步：解读排序表，评估可信度，给建议

不管这份 `ranked_peptides.tsv`（或 `scoring/ranked_peptides.v03.tsv`）是第三步刚跑出来的，还是第零步发现用户已经有的现成结果，**下面的检查方法都一样适用**——列名是同一套打分函数产出的（`efficacy_score`/`final_priority`/`l3_*_score`/`safety_status`/`escape_multiplier`等），不需要区分来源。

打开表格，做以下几件事：

### 1. 看整体分布，判断这次打分是否"有区分度"
```python
import csv
from collections import Counter
rows = list(csv.DictReader(open("ranked_peptides.tsv"), delimiter="\t"))
print(Counter(r["final_priority"] for r in rows))
```
如果几乎所有候选都挤在同一档（比如99%都是D），说明这批候选整体证据都偏弱，排序更多是矬子里拔将军，不代表真的有"强候选"；如果A/B档有一批数量合理的候选，说明这批数据本身有信号。

### 2. 抽查关键字段是不是"占位符恒定值"——这是判断可信度最关键的一步

**优先用第零步0.5的代码检测结果**：如果已经确认过neo仓库有没有 `self_similarity.py`/`_parse_oncokb_gene_list`/`_parse_gtex_gct` 这些模块，直接依据那个确定性结论来判断，不需要再靠下面这种"猜数据模式"的方法；只有拿不到neo仓库路径、只能看产出数据本身时，才用下面的方法。

**driver_relevance / tumor_specificity**——不要只看"是否全部相同"，这个判断标准现在不够可靠了：即使是真实的OncoKB/GTEx查表，"参考表里没收录的基因"也会拿到统一的中性默认值（driver_relevance=0.3、tumor_specificity=0.5），如果这批变异恰好大多命中冷门基因，看起来会跟"完全没修复的占位符"很像。更可靠的做法是**挑几个明确的知名基因验证**：
- 表里如果有 KRAS、TP53、EGFR、PIK3CA 这类常见驱动基因，检查它们的 `driver_relevance` 是不是明显偏高（真实OncoKB查表算出来的典型值在0.9附近）——如果连这些耳熟能详的基因都只拿到0.3左右的中性值，基本可以确定没有真实查表在生效
- 如果找不到这类知名基因作参照，退而求其次看数值有没有明显差异（不同基因给出不同分数，至少说明不是写死的常量）

**self_similarity_score**——原来的判断方法（"全部恒为0就是没修复"）现在会误判，因为修复后的真实计算，对frameshift/融合/剪接这类没有对应野生型的肽段，本来就应该输出0（这是正确答案，不是bug）。正确的排查方法是**同时对照 `wildtype_peptide` 这一列**：
- 挑几行 `wildtype_peptide` **非空**（说明是missense类型，理应有真实野生型可比较）的候选，看它们的 `self_similarity_score` 是不是也恒为0——如果是，大概率没打上这个修复补丁，这条安全检查实质上没有真实跑过，需要提醒用户"高排名不等于已排除自身免疫风险，只是没评估"
- 如果 `wildtype_peptide` 为空（本来就没有野生型可比较，比如融合/移码肽），对应的 `self_similarity_score=0` 是**符合预期的正确结果**，不代表检查失效

**即使确认这个字段是真实计算的，也要提醒用户一个已知局限**：目前的计算**只覆盖"这条肽跟它自己对应的野生型像不像"**，不做全人类蛋白组范围的"这条肽跟任意一个不相关的自身蛋白像不像"这种更广的搜索（这部分功能还在方案阶段，尚未实现）。所以`self_similarity_score`低不代表这条肽在整个人体蛋白组里都找不到相似序列，只能说明它跟自己的野生型对照差异大。

### 3. 结合表达证据判断排序前几名是否"真的在表达"
对照 `event_expression`/`rna_alt_reads` 这两列，如果排名很靠前的候选这两项都是0或很低，要单独提醒用户："这条候选虽然结合/呈递预测分数很高，但没有RNA层面的直接表达证据，存在'预测出来的肽段其实根本没被转录'的风险"。

### 4. 给出下一步建议（按用户实际场景调整，不要机械套模板）
- 优先推荐排名靠前、**同时**表达证据扎实（`rna_alt_reads`>0）、安全性PASS、免疫逃逸未评估但对应HLA没有已知问题的候选，作为湿实验验证（ELISpot/四聚体）的第一批对象
- 如果A/B档候选数量太少，建议检查是不是筛选阈值（呈递门控里的EL rank/stabpan阈值）设得过严，可以回到第二步读到的TOML数值，跟用户讨论要不要针对这批数据放宽
- 如果发现`driver_relevance`/`tumor_specificity`是占位符恒定值，明确建议"用最新适配器重新生成一份 raw_events.tsv 再排序"，而不是直接拿现在这份排序结果做决策
- 如果发现`self_similarity_score`对有野生型对照的候选也恒为0，同样建议重新生成；如果这个字段已经是真实计算的，仍要提醒"目前只覆盖跟自身野生型的比较，还没有做全蛋白组搜索，自身免疫风险的排查不是100%完整"
- 如果用户的应用场景涉及免疫逃逸风险（既往治疗史、怀疑HLA丢失），建议补跑 `immune-escape` 相关步骤（需要VEP/CNV/HLA-LOH文件）后再排序，现在这份结果里 `escape_multiplier` 全部是1.0，等于没评估这个风险

---

## 事件排序 vs 肽段排序——两套独立并行的排名，不要混为一谈

`ranked_events.tsv` 和 `ranked_peptides.tsv` 衡量的是完全不同的问题，而且**互相不耦合**，解读时必须分开讲清楚，避免用户想当然地认为"这两张表应该一致"：

| | `ranked_events.tsv`（变异排序） | `ranked_peptides.tsv`（肽段排序） |
|---|---|---|
| 排的是什么 | **变异事件本身**（一个基因上的一个突变） | **具体的"肽段×HLA"候选组合** |
| 回答的问题 | 这个变异本身值不值得关注（检出可信度高不高、有没有表达、是不是打在驱动基因上、克隆性/持续性如何） | 这个具体候选能不能被呈递、会不会引发T细胞识别、打了安不安全 |
| 用的分数 | `event_score`（`[event_weights]`配置） | `efficacy_score`（`[l3_weights]` + 五个独立乘子） |
| 一对多关系 | 一个变异对应一行 | 一个变异可能衍生出**好几行**（不同肽段长度的滑窗、不同HLA分型分别是独立候选） |

**关键点：这两套分数是独立并行计算的，不是一个乘进另一个**。`event_score`高，不代表这个变异衍生出来的肽段候选在`ranked_peptides.tsv`里也会排前面；反过来，一条肽段候选`efficacy_score`很高，也不代表它对应的变异在`ranked_events.tsv`里排名靠前。这是有意的架构设计——避免"变异证据强"直接拉高"肽段呈递弱"的候选，也避免"某条肽段结合预测意外地强"就让本来证据一般的变异显得很重要。

**给用户解读时的正确姿势**：
- 如果用户问"这个变异重要吗" → 看 `ranked_events.tsv`
- 如果用户问"这条候选该不该优先做验证" → 看 `ranked_peptides.tsv`，**不要用变异层面的`event_score`去佐证肽段层面的结论**
- 如果用户发现"这个变异在事件表里排很前，但它衍生的肽段候选在肽段表里都排很后" → 这不是矛盾，如实解释："这个变异本身证据很扎实（比如高表达、打在KRAS这种驱动基因上），但它能滑窗出来的肽段，结合/呈递预测都不理想，或者安全性/呈递门控没通过——变异值得关注不等于它一定能产生好的治疗靶点"

---

---

## 权重/阈值怎么改——只讲方法，绝不代替用户动手改

**硬性规则：这个skill的范围内，不管用户怎么要求（包括"帮我改一下""你直接改就行""改完顺便跑一下"），都不要用任何工具（str_replace/bash/create_file等）去改 `profiles/*.toml` 里的任何字段。** 只负责：①告诉用户具体改哪个文件、哪一行、改成什么语法；②用户自己改完之后，可以帮忙重新跑第三步的脚本、或者帮忙确认改动是否生效。**改动动作本身必须由用户手动完成。**

为什么要这样限制：这些权重和阈值编码的是"呈递比结合重要""亚克隆突变要打几折"这类具体的科研/临床判断，属于会实际影响候选优先级、进而影响湿实验资源投入方向的决策参数。这类改动应该由用户自己明确做出、自己心里有数改了什么，而不是Claude在后台悄悄改掉——否则同一份数据前后两次跑出不同排序时，用户会搞不清楚是数据变了还是配置被谁动过了，可复现性和可追责性都会被破坏。

### 该怎么教用户改

先明确告诉用户两个选择，并说明各自影响范围：

1. **直接改 `profiles/default.toml`**：会影响**所有没有显式覆盖该字段的profile**（因为 `load_profile()` 会自动把default.toml当基础再合并其他profile，参见第二步）。适合"这个改动想对所有场景生效"的情况。
2. **改用户自己场景的profile**（比如 `profiles/leukemia.toml`）或新建一个自定义profile文件：只影响这一个profile，不动全局默认值。更安全，适合"只想给这批数据单独试一组参数"的情况。**新建自定义profile时，只需要在新文件里写你想覆盖的字段，其余会自动从default.toml继承，不用整份复制。**

然后按用户想改的内容，指给他们具体的TOML区块（**现场`cat`一下实际文件，把真实存在的区块和当前数值报给用户，不要凭记忆列，因为字段是否存在、当前值是多少，取决于这份仓库具体打过的补丁版本**）：

| 想改什么 | 去哪个区块 | 格式提示 |
|---|---|---|
| 十维综合分各维度的相对重要性 | `[l3_weights]` | `hla_presentation = 0.25` 这种 `key = 数字`，数字之间的比例才有意义（代码会自动除以总和归一化），不用凑到1.0 |
| 事件层排序（ranked_events.tsv）用的权重 | `[event_weights]` | 同上 |
| 结合/呈递证据里NetMHCpan、MHCflurry各占多少 | `[presentation_weights]` | 同上 |
| 呈递门控硬阈值（grade等级、EL rank上限、稳定性下限、不达标打几折） | `[gates]` | `max_el_rank = 2.0` / `failure_multiplier = 0.25` 这种 |
| 正常组织表达安全性判定的CAUTION/FAIL界限 | `[safety]` | 各种 `*_tpm` 结尾的阈值 |
| 免疫逃逸各机制对应打几折、封顶到哪一档 | `[immune_escape.multipliers]` / `[immune_escape.priority_caps]` | 乘子是 `0~1` 的小数，封顶档位是字符串比如 `"C"` |
| 驱动基因/肿瘤特异性查表用哪份参考文件、未收录基因给多少分 | `[driver_genes]` / `[tumor_specificity_gtex]` | `reference_path = "..."` 指向真实数据库文件路径 |

举例给用户的具体操作应该长这样（把方括号里的换成真实值，直接告诉用户改哪一行）：

> 打开 `profiles/default.toml`，找到 `[l3_weights]` 这一段，把 `hla_presentation = 0.25` 这一行的 `0.25` 改成你想要的数字，保存。不用改其他行，代码会自动重新按比例归一化。改完之后告诉我，我帮你重新跑一遍排序看看效果有没有变化。

### 用户改完之后

用户确认改完了，**这时候可以做的事**：`cat` 一下改动后的文件确认真的存的是新值（帮用户核对，不是帮用户改），然后重新跑第三步的打分脚本，对比改动前后 `ranked_peptides.tsv` 里同一批候选的排名/分数变化，讲清楚这次改动实际造成了什么影响。

---


---

## 常见追问处理

- **"为什么这条肽排在前面/后面"** → 打开 `ranked_peptides.tsv` 里这条肽对应的行，把 `l3_*_score` 各维度分、四个乘子（`safety_status`对应的`safety_multiplier`、`appm_multiplier`、`presentation_gate_multiplier`、`ccf_multiplier`）逐项列出来，对照第二步读到的真实权重逐项解释，不要只说"综合分较高/较低"。
- **"这个权重能不能改"** → 不要自己动手改，按下面"权重/阈值怎么改"这一节的方法教用户自己改。
- **"这个排序靠谱吗"** → 直接引用第四步第2、3点的检查结果回答，不要笼统地说"靠谱"或"不靠谱"。
- **"这个变异排名很靠前，为什么它的肽段候选反而排名靠后"（或反过来）** → 不是矛盾，参照上面"事件排序 vs 肽段排序"这一节解释：两套分数独立计算，变异本身证据强不代表它衍生出的具体肽段候选就一定呈递/结合得好。
