# NeoAg v0.5.1：双生成器、突变致剪接与正常背景三证据链

## 1. 版本目标

v0.5.1在v0.5.0正式Splice Provenance Layer之上新增三条互相独立、可以审计的证据链：

```text
RNA-driven
  精确RNA junction / event
  → ImmunoPepper
  → moPepGen
  → event + peptide 或 event + ORF 共识

DNA-causal
  精确DNA variant
  → splice2neo
  → 精确variant-junction关系
  → EasyQuant靶向重定量
  → pVACsplice变异致剪接肽-HLA结果

normal-background
  匹配正常或协议匹配正常的位点覆盖
  + 正常junction目录
  + k4neo大规模序列索引筛查
```

三条链只在最终共识层汇合。工具数量本身不是证据数量；来自同一BAM或同一上游模型的相关结果不会被重复计票。

## 2. 精确身份契约

### Junction

```text
SJ|BUILD|CHR|INTRON_START_1BASED|INTRON_END_1BASED|STRAND
```

### Variant

```text
VAR|BUILD|CHR|POS_1BASED|REF|ALT
```

### 外部查询

EasyQuant和k4neo只能使用第一遍建库生成的内容寻址query ID：

```text
SEQ|<stable digest>
```

禁止以下回填：

- 同基因最大reads；
- 最近基因组位置；
- 不带strand的模糊junction；
- 近似肽段；
- 只按gene或transcript名称合并；
- 使用另一轮临时event ID生成的query map。

## 3. RNA-driven链

### ImmunoPepper

作为主剪接图肽段生成器。短读长局部翻译仍标记为局部转录本/ORF假说，不能自动声称为完整全长ORF。

### moPepGen

作为第二个多组学肽段生成器。它的FASTA记录通常是非经典肽段产物，不一定是完整ORF，因此v0.5.1区分：

- `DUAL_GENERATOR_EXACT_PEPTIDE`：同一canonical event、完全相同肽段；
- `DUAL_GENERATOR_EXACT_ORF`：同一event、相同protein SHA-256、相同frame，且两者都提供可比较的有效ORF；
- `SINGLE_GENERATOR_RNA_SUPPORTED`：只有一个可解析生成器。

moPepGen只有通过明确的provenance map或带可解析junction的GVF，才可形成高等级共识。无法映射的FASTA记录保留在冲突表，不按gene猜测。

## 4. DNA-causal链

### splice2neo

负责统一DNA剪接效应预测与RNA junction证据。核心输出必须包括：

- `junc_id`，格式为`chr:pos1-pos2:strand`；
- 变异染色体、位置、REF、ALT；
- 对应转录本/CDS/肽段上下文；
- SpliceAI、Pangolin、MMSplice等存在时的分数。

### EasyQuant

第一遍建库输出：

```text
splice_easyquant_input.tsv
splice_easyquant_query_map.tsv
```

EasyQuant结果只能按`name == query_id`导回。未知或重复映射query ID进入冲突表。支持状态区分：

- `TARGETED_REQUANT_SUPPORTED`；
- `TARGETED_SPANNING_ONLY`；
- `TARGETED_REQUANT_NEGATIVE`。

### pVACsplice

pVACsplice结果只有同时解析到：

```text
exact variant + exact strand-aware junction + canonical event
```

才进入DNA-causal presentation。只给junction start/stop但没有strand，且无法通过显式map唯一解析时，结果被拒绝，不做最近位置或gene级匹配。

## 5. normal-background链

### Coverage-aware正常证据

只有正常样本/面板在位点上具有充分覆盖且未检出，才能形成：

```text
NOT_DETECTED_ADEQUATE_COVERAGE
```

低覆盖或未评估不能解释为安全。

### k4neo

k4neo结果仅按项目生成的`cts_id == query_id`导回。系统严格区分：

- `NOT_DETECTED_KMER_SCREEN`：大规模序列索引未检出；
- `NOT_DETECTED_ADEQUATE_COVERAGE`：具体位点有充分覆盖且未检出。

k4neo阴性不能单独证明某个正常组织位点具有充分覆盖，因此仅k4neo阴性会触发R3上限。k4neo正常关键组织阳性可形成强负面证据。

## 6. 输出表

v0.5.1新增：

```text
splice_variants.tsv
splice_causal_links.tsv
splice_sequence_queries.tsv
splice_targeted_quantification.tsv
splice_pvacsplice_predictions.tsv
splice_evidence_chains.tsv
splice_easyquant_input.tsv
splice_easyquant_query_map.tsv
splice_k4neo_input.tsv
splice_k4neo_query_map.tsv
```

继续保留v0.5.0的junction、event、transcript hypothesis、ORF、peptide origin、pVACbind、normal background、consensus、conflict和兼容投影表。

## 7. 推荐生产运行

```bash
bash scripts/run_splice_provenance_v051.sh \
  --sample-id PATIENT_001 \
  --outdir results/PATIENT_001/splice_v051 \
  --junctions results/regtools/junctions.tsv \
  --spladder-gff3 results/spladder/events.confirmed.gff3 \
  --irfinder results/irfinder/IRFinder-IR-nondir.txt \
  --immunopepper-meta results/immunopepper/peptides_meta.tsv \
  --mopepgen-fasta results/mopepgen/variant_peptides.fasta \
  --mopepgen-gvf results/mopepgen/alt_splice.gvf \
  --mopepgen-provenance-map results/mopepgen/neoag_map.tsv \
  --splice2neo results/splice2neo/export.tsv \
  --run-easyquant --easyquant-bam results/rna/tumor.bam \
  --normal-coverage results/normal/coverage.tsv \
  --run-k4neo --k4neo-database /refs/k4neo/db --k4neo-index /refs/k4neo/index.tsv \
  --k4neo-license-accepted \
  --run-pvacsplice --pvacsplice-junctions results/regtools/cis.tsv \
  --annotated-vcf results/vep/tumor.vep.vcf.gz \
  --ref-fasta /refs/GRCh38.fa --gtf /refs/gencode.gtf \
  --hla HLA-A*02:01,HLA-B*07:02 \
  --strict --overwrite
```

## 8. 版本边界

本版本完成的是：

- 外部工具输出适配；
- 内容寻址query map；
- 精确variant-junction-event-peptide-HLA回链；
- 三链状态、冲突和共识；
- 两遍/三遍式生产编排。

本版本不声称：

- 计算ORF一定被内源性翻译；
- HLA预测等于真实呈递；
- k4neo阴性等于所有正常组织绝对不存在；
- 单个工具的结果等于实验确认；
- 软件测试等于临床验证。
