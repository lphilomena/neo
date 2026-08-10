# NeoAg Event Pipeline v0.5.0 正式 Splice Provenance Layer

NeoAg Event Pipeline 是研究型肿瘤新抗原候选优先级分析流程。它将 SNV/InDel、fusion、splice、结构变异以及 peptide-only 候选统一转换为标准事件表和 peptide-HLA 表，并叠加 HLA 呈递、APPM、CCF、safety、immune escape、validation plan 和 report 证据层。

本包是轻量级 online release，包含源代码、CLI 入口、Nextflow workflow、测试、fixture、profile、安装脚本和文档。不包含大型参考库、授权工具、conda 环境、缓存 work 目录、真实患者数据或生产结果。

重要边界：本流程输出的是 computational triage 和实验验证规划结果，不构成临床诊断、临床耐药判定或已验证治疗建议。

## v0.5.0 正式 Splice Provenance Layer

v0.5.0 将 v0.4.4 的精确 junction 修复升级为具有完整参照完整性的正式剪接来源模型：

```text
junction → splice event → transcript hypothesis → ORF → peptide origin → peptide-HLA presentation
```

本版本新增 canonical 实体注册表、event-junction 与 peptide-origin 关系表、保守的正常背景状态机、按独立证据组计算的共识、pVACbind FASTA Index 严格回链，并继续投影生成兼容的 `raw_events.tsv`、`raw_peptides.tsv` 和 `rna_junction_evidence.tsv`。生产入口为 `scripts/run_splice_provenance_v050.sh`，Python CLI 为 `neoag-splice-layer`。详见 `docs/V050_SPLICE_PROVENANCE_LAYER.md` 与 `CHANGELOG_V050_SPLICE_PROVENANCE_LAYER.md`。

## v0.4.4 精确 Junction 与 Provenance 修复

v0.4.4 建立 canonical splice-junction ID，删除同基因和邻近位点 reads 回填，将上游声明值与精确验证值分离；生产合并不再丢弃后续来源记录，splice 跨域确认必须落在同一 canonical junction。详见 `CHANGELOG_V044_SPLICE_PROVENANCE.md` 与 `docs/V044_SPLICE_PROVENANCE.md`。

## 功能概览

本流程可以：

- 将 pVACtools-like SNV/fusion/splice 输出解析为 `raw_events.tsv` 和 `raw_peptides.tsv`。
- 从 VEP 注释 VCF 生成 sliding-window variant peptides；当 VCF 缺少 CSQ 注释时，可在配置完整时自动先运行 VEP 注释。
- 使用 NetMHCpan、MHCflurry 以及可选稳定性/免疫原性工具构建 MHC 呈递证据。
- 构建 APPM 2.0 证据，包括输入完整性、冲突、肽段修饰因子和 immune-context 注释。
- 基于 purity、CNV 和 VAF context 估计 CCF/clonality。
- 基于 normal expression、normal ligandome、normal junction、matched-normal 和 reference proteome context 构建 peptide safety 证据。
- 基于 HLA LOH、APPM、CCF、B2M/JAK/APM context 及相关证据表构建 immune escape 证据。
- 为 frameshift、splice、exon-junction、fusion 和 SV 候选生成 long peptide / minigene 验证设计。
- 生成患者沟通版和科研技术版 HTML 报告。
- 通过 CLI 或内置 Nextflow wrapper 运行 fixture workflow。

`.tsv` 后缀是 schema 兼容标签，不代表软件版本。当前 release 版本是 v0.5.0，但仍写出 schema-compatible 表格，便于旧下游脚本继续读取相同文件名。

## Agent Skills 和 Coordinator

本 release 增加了仓库内置的 agent skills pack，位于 `.agents/skills/`，并提供轻量 coordinator CLI：

```bash
neoag-agent --message "比较 recommendation 和 NetMHCpan42 排序差异" --result-dir results/sample --outdir work/agent_plan
```

默认模式只生成 dry-run plan。对支持的低风险 skill，可加 `--execute` 执行。skill 列表、输入、输出和解释边界见 `docs/AGENT_SKILLS_P0_P1.md`。

## 新机器安装和运行：三个宏 Skill

新机器部署和患者样本运行统一使用三个公开 Open-Neo 宏 Skill。README 只提供人类
导航，机器可读 manifest 仍是准绳：

- `.agents/skills/open-neo-install-check/SKILL.md`：Skill1，负责新机器安装、
  参考库/工具发现、经批准的安装或修复、Doctor、smoke test 和生产运行就绪检查。
- `.agents/skills/open-neo-run/SKILL.md`：Skill2，负责输入质控、路线选择、
  Gateway 受控执行/续跑、多工具证据生成、旧加权排序和证据共识排序。
- `.agents/skills/open-neo-review/SKILL.md`：Skill3，负责只读结果审阅、
  实验优先级表和患者版/技术版报告。

### Skill1：安装并验收新机器

在目标机器 clone 发布分支：

```bash
mkdir -p /home/na/project
git clone --branch na0707_upload_release \
  https://github.com/lphilomena/neo.git \
  /home/na/project/neo
cd /home/na/project/neo
```

如果新机器上 `open-neo` console script 还不在 `PATH`，可先使用项目环境里的
模块入口：

```bash
export PYTHONPATH="$PWD/src"
alias open-neo='/home/na/miniforge3/envs/neoag-tools/bin/python -m neoag.open_neo.cli'
```

运行安装检查宏。生产使用默认按 `all-open` 安装；它会在许可允许的前提下安装
开放生产工具集，同步全套生产参考资产，并在 `configs/local/` 写入本机配置后的
manifest。

```bash
open-neo install-check \
  --project-root "$PWD" \
  --deployment-tier full \
  --mode install \
  --installer-profile all-open \
  --asset-source-host na@10.200.50.134 \
  --asset-source-root /mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata4git \
  --tools-root /home/na/project/open-neo-deploy/env_tool \
  --reference-root /home/na/project/open-neo-deploy/refs \
  --conda-base /home/na/miniforge3 \
  --allow-download \
  --approved \
  --outdir work/install-check-full
```

`--installer-profile minimal` 适合 review/core 场景；
`--installer-profile standard` 是较轻的生产主路径安装。重复运行 Skill1 是安全
的：已经安装的工具、已经同步的资产和 PASS checkpoint 会在签名匹配时复用。
如果安装中断，用 `resume` 续跑：

```bash
open-neo install-check \
  --project-root "$PWD" \
  --deployment-tier full \
  --mode resume \
  --installer-profile all-open \
  --approved \
  --outdir work/install-check-full
```

完整安装默认没有墙钟超时。若操作者显式提供 `--install-timeout SECONDS`，
中断或超时时会先终止整个安装进程组，再写入可控 checkpoint。

### Skill2：通过 Gateway 运行患者样本

先在目标机器本机启动 NeoAg Gateway。除非经过审阅，不建议把 Gateway 暴露到
内网或公网。

```bash
cd /home/na/project/neo
source /home/na/miniforge3/bin/activate neoag-tools

mkdir -p work/neoag_gateway
nohup env PYTHONPATH="$PWD/src" \
  python -m neoag.controlled_execution.gateway \
  --host 127.0.0.1 \
  --port 8000 \
  --project-root "$PWD" \
  --outdir "$PWD/work/neoag_gateway" \
  --allowed-root "$PWD/work" \
  --allowed-root /mnt/zzbnew/Public/neoag_results \
  > work/neoag_gateway/gateway.log 2>&1 &

curl -s http://127.0.0.1:8000/health
```

有样本 manifest 时提交完整流程：

```bash
open-neo run \
  --sample-manifest configs/local/sample.yaml \
  --tools-manifest configs/local/tools_manifest.configured.yaml \
  --reference-manifest configs/local/reference_manifest.configured.yaml \
  --outdir /mnt/zzbnew/Public/neoag_results/CASE001 \
  --mode execute \
  --approved \
  --gateway-url http://127.0.0.1:8000 \
  --gateway-wait
```

没有 manifest 时，也可以直接给 BAM/VCF/RNA FASTQ：

```bash
open-neo run \
  --sample-id CASE001 \
  --tumor-dna-bam /path/to/tumor.bam \
  --normal-dna-bam /path/to/normal.bam \
  --somatic-vcf /path/to/somatic.pass.vcf.gz \
  --tumor-rna-fastq /path/to/R1.fq.gz /path/to/R2.fq.gz \
  --rna-threads 12 \
  --outdir /mnt/zzbnew/Public/neoag_results/CASE001 \
  --mode execute \
  --approved \
  --gateway-url http://127.0.0.1:8000 \
  --gateway-wait
```

Skill2 会在安全时并行运行同类工具。成对 GRCh38 DNA 输入下，纯度/CNV 证据
阶段可同时运行 FACETS、Sequenza 和 PURPLE，再生成跨工具 purity/ploidy 推荐。
HLA LOH 必须先等到非单工具纯度共识，再并行启动 LOHHLA 和 SpecHLA。两个
HLA LOH 工具都有可用输出时报告标记为 `dual_tool_consensus`；若其中一个失败
或无输出，流程继续使用另一个工具的单工具证据，并在
`hla_loh_tool_status.tsv`、`hla_loh_summary.json`、`hla_loh_review.md` 和
`recommended_hla_loh.tsv` 中明确标记 `single_tool_result`。

中断后的样本续跑：

```bash
open-neo run \
  --result-dir /mnt/zzbnew/Public/neoag_results/CASE001 \
  --mode resume \
  --approved \
  --gateway-url http://127.0.0.1:8000 \
  --gateway-wait
```

### Skill3：审阅结果并生成报告

Skill2 完成排序后，对结果目录只读运行 Skill3：

```bash
open-neo review \
  --result-dir /mnt/zzbnew/Public/neoag_results/CASE001 \
  --reports patient,technical,onepage \
  --outdir /mnt/zzbnew/Public/neoag_results/CASE001/review
```

Skill3 不重跑重型工具。它检查结果完整性，比较旧加权排序和证据共识排序，
生成事件级实验优先级表，并输出有边界的患者版/技术版报告。缺失证据或单工具
证据只能写成部分证据，不能解释成生物学阴性。

## 快速开始

从项目根目录运行：

```bash
python -m pip install -e .
neoag run-demo --outdir work/demo_v043 --sample-id DEMO001
```

重要 demo 输出包括：

- `work/demo_v043/scoring/ranked_peptides.tsv`
- `work/demo_v043/scoring/ranked_peptides.evidence_consensus.tsv`（并行证据共识排序，不替换原排序）
- `work/demo_v043/scoring/ranked_events.evidence_consensus.tsv`
- `work/demo_v043/scoring/evidence_states.tsv`
- `work/demo_v043/scoring/evidence_conflicts.tsv`
- `work/demo_v043/scoring/weighted_vs_consensus_comparison.tsv`
- `work/demo_v043/scoring/ranked_events.tsv`
- `work/demo_v043/scoring/validation_plan.tsv`
- `work/demo_v043/reports/evidence_report.html`
- `work/demo_v043/reports/evidence_report.patient.html`
- `work/demo_v043/reports/evidence_report.technical.html`
- `work/demo_v043/appm/appm_summary.tsv`
- `work/demo_v043/appm/appm_peptide_modifiers.tsv`
- `work/demo_v043/clonality/ccf_lite.tsv`
- `work/demo_v043/safety/peptide_safety.tsv`
- `work/demo_v043/immune_escape/peptide_escape_flags.tsv`

### 并行证据共识排序

每次流水线运行都会保留现有固定权重的 `ranked_peptides.tsv`。流程随后从
`comprehensive_peptide_evidence.tsv` 独立生成标准化证据状态、R1-R4 证据等级、
同事件来源赛道内的 Pareto 分层，以及肽段和事件共识排序。缺失证据会保持为
`MISSING`，不会解释为生物学阴性；硬失败和 priority cap 只影响共识分支。

对已有综合证据表运行：

```bash
neoag evidence-rank \
  --comprehensive-evidence results/sample/scoring/comprehensive_peptide_evidence.tsv \
  --weighted-baseline results/sample/scoring/ranked_peptides.tsv \
  --rules configs/ranking/sarcoma_evidence_consensus_v1.toml \
  --provenance results/sample/provenance.json \
  --outdir results/sample/scoring/evidence_consensus \
  --mode parallel --track all \
  --emit-event-ranking --compare-weighted --deterministic
```

编程 Agent 应通过轻量的对外 Skill2 包装器调用正式算法：

```bash
neoag-skill run open-neo-run \
  --outdir results/sample/scoring/evidence_consensus \
  --arg comprehensive_evidence=results/sample/scoring/comprehensive_peptide_evidence.tsv \
  --arg weighted_baseline=results/sample/scoring/ranked_peptides.tsv
```

`open-neo-run` 与兼容入口 `neoag-ranking` 都调用正式的
`neoag evidence-rank` CLI。Skill 层不再维护独立的 R1-R4、Pareto、hard fail、
priority cap 或 event 去重算法。

使用通用审计命令比较任意两个排序：

```bash
neoag-ranking-compare \
  --left results/sample/scoring/ranked_peptides.weighted_baseline.tsv \
  --left-name weighted_baseline \
  --right results/sample/scoring/ranked_peptides.evidence_consensus.tsv \
  --right-name evidence_consensus \
  --outdir results/sample/scoring/ranking_comparison
```

输出包括 Top10/20/50/100 overlap、Spearman 相关性、候选升降级、高位 hard-fail、
event type 与 HLA 构成、证据冲突比例、missing evidence 比例和人工复核候选。

第一阶段不提供 `--replace-primary-ranking`。逐行差异原因写在
`ranking_compare_weighted_vs_consensus.tsv`；旧文件名继续作为兼容别名。字段定义和算法边界见
`docs/EVIDENCE_CONSENSUS.md`。所有阈值均标记为
`PROVISIONAL_RESEARCH_ONLY`，第一阶段仅用于比较算法行为，不代表临床验证。

运行测试：

```bash
python -m pip install -e '.[test]'
pytest -q
```

默认测试命令会刻意跳过 integration、benchmark 和 external-tool 测试。

## 常用运行命令

### 准备环境

仅用于 fixture 开发：

```bash
python -m pip install -e '.[test]'
pytest -q
neoag run-demo --outdir work/demo_v043 --sample-id DEMO001
```

需要外部工具的运行：

```bash
bash scripts/setup_tools_env.sh
source conf/tools.env.sh
python -m pip install -e '.[test]'
neoag check-tools
```

较小的开发/测试环境：

```bash
NEOAG_TOOLS_LITE=1 bash scripts/setup_tools_env.sh
source conf/tools.env.sh
python -m pip install -e '.[test]'
pytest -q
```

### 完整生产流程

生产流程并行执行 HLA 分型、纯度/CNV/CCF、HLA LOH/APPM 和 RNA 表达分析，再合并 SNV/InDel、fusion、splice 三类候选肽，统一进行呈递预测和排序。Fusion 或 splice 事件只有在跨断点/异常 junction 肽经过 NetMHCpan、MHCflurry 以及共同的免疫原性和安全性证据层后，才算完成新抗原评估。

完整依赖关系和验收标准见 `docs/PRODUCTION_WORKFLOW.md`。完整阶段执行使用 `configs/workflows/production_workflow.example.toml` 和 `neoag-production-run`；已有分支结果时使用 `conf/run.production_multisource.example.toml` 合并。缺少 fusion/splice 来源时流程继续运行并标记 `LOW_CONFIDENCE`；缺少必需呈递预测结果时才会失败。

下面的脚本用于运行 manifest 驱动的完整生产流程。首次运行会创建 private manifest 并停止，必须先填入真实的样本路径、HLA 分型、BAM/FASTQ/VCF、参考数据和授权工具配置，才能再次执行。

```bash
#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-/root/neo/src/na0707_upload_release}"
SAMPLE_ID="${SAMPLE_ID:-SAMPLE001}"
MANIFEST="${MANIFEST:-${PROJECT_ROOT}/configs/workflows/production_workflow.private.toml}"
OUTDIR="${OUTDIR:-${PROJECT_ROOT}/results/${SAMPLE_ID}_production}"

cd "$PROJECT_ROOT"
source conf/tools.env.sh
python -m pip install --no-deps -e .

if [[ ! -f "$MANIFEST" ]]; then
  cp configs/workflows/production_workflow.example.toml "$MANIFEST"
  echo "已创建 $MANIFEST"
  echo "请填写所有 /path/... 路径和各阶段命令，然后重新运行本脚本。"
  exit 2
fi

# 环境/工具预检，然后生成只展示命令的生产执行计划。
neoag check-tools
neoag-production-run \
  --project-root "$PROJECT_ROOT" \
  --manifest "$MANIFEST" \
  --outdir "$OUTDIR"

# 执行缺失阶段；已经存在声明输出的阶段会自动复用。
neoag-production-run \
  --project-root "$PROJECT_ROOT" \
  --manifest "$MANIFEST" \
  --outdir "$OUTDIR" \
  --execute

# 完整生产结果的最低验收条件。
test -s "$OUTDIR/production_run_summary.json"
test -s "$OUTDIR/production_stage_status.tsv"
test -s "$OUTDIR/final/scoring/ranked_peptides.tsv"
test -s "$OUTDIR/final/scoring/comprehensive_peptide_evidence.tsv"
test -s "$OUTDIR/final/reports/evidence_report.html"

python -m json.tool "$OUTDIR/production_run_summary.json"
```

常规断点续跑不要添加 `--force`；只有明确需要重算已有输出的阶段时才使用。最终状态为 `LOW_CONFIDENCE` 表示排序已经完成，但 fusion、splice 等一个或多个可选来源缺失或失败；解释排序前应检查 `production_stage_status.tsv`。

### 从已有 pVAC-like 表运行

当你已有 pVACseq/pVACfuse/pVACsplice-like aggregated tables 时使用：

```bash
neoag run \
  --outdir results/sample \
  --sample-id SAMPLE001 \
  --profile default \
  --pvac data/fixtures/pvacseq_aggregated.tsv \
  --immunogenicity-stub
```

### 从预生成 raw intermediates 运行

当 `parsed/raw_events.tsv` 和 `parsed/raw_peptides.tsv` 已经存在时使用：

```bash
neoag run \
  --outdir results/sample \
  --sample-id SAMPLE001 \
  --profile default \
  --raw-events results/sample/parsed/raw_events.tsv \
  --raw-peptides results/sample/parsed/raw_peptides.tsv \
  --netmhcpan results/sample/presentation/netmhcpan.xls \
  --mhcflurry results/sample/presentation/mhcflurry.csv \
  --expression results/sample/parsed/expression.tsv \
  --hla-loh results/sample/tools/hla_loh.tsv \
  --purity results/sample/tools/purity.tsv \
  --cnv results/sample/tools/cnv_segments.tsv
```

### Sliding-window Variant Peptides 到排序结果

当你有 somatic SNV/InDel VCF，并希望通过 sliding window 生成 mutant peptides、预测 peptide-HLA 呈递并生成 event/peptide 排序结果时使用。

如果 VCF 已含 VEP `CSQ` 注释，pipeline 会直接使用。如果缺少 `CSQ`，在 VEP、cache、reference FASTA 和 plugins 配置完整时，`run-full` 会先自动运行 VEP 注释。

```bash
cat > conf/run.sliding.private.toml <<'TOML'
[sample]
id = "SAMPLE001"
profile = "default"

[tools]
stub = false
enabled = ["netmhcpan", "mhcflurry"]
immunogenicity_stub = false

[inputs]
entry_mode = "snv_indel"
variant_peptide_extraction = true
variants_vcf = "/path/to/sample.somatic.pass.vcf.gz"
tumor_sample_name = "TUMOR"
hla_alleles = ["HLA-A*02:01", "HLA-B*07:02", "HLA-C*07:02"]
extract_appm_from_vcf = false
normal_expression = "resources/normal_expression.example.tsv"
normal_hla_ligands = "resources/normal_hla_ligands.example.tsv"
TOML

neoag run-full \
  --config conf/run.sliding.private.toml \
  --outdir results/SAMPLE001_sliding
```

关键输出：

- `results/SAMPLE001_sliding/upstream/tools/variant_peptides.tsv`
- `results/SAMPLE001_sliding/upstream/tools/variant_peptides.annotated.tsv`
- `results/SAMPLE001_sliding/upstream/parsed/raw_events.tsv`
- `results/SAMPLE001_sliding/upstream/parsed/raw_peptides.tsv`
- `results/SAMPLE001_sliding/scoring/ranked_events.tsv`
- `results/SAMPLE001_sliding/scoring/ranked_peptides.tsv`
- `results/SAMPLE001_sliding/scoring/comprehensive_peptide_evidence.tsv`（合并全部注释、预测器、RNA、APPM、CCF、安全性、免疫逃逸、排序与验证证据）
- `results/SAMPLE001_sliding/scoring/validation_plan.tsv`
- `results/SAMPLE001_sliding/reports/evidence_report.html`
- `results/SAMPLE001_sliding/reports/evidence_report.patient.html`
- `results/SAMPLE001_sliding/reports/evidence_report.technical.html`

当 `inputs.purity_recommendation` 指向 `purity_recommendation.json` 时，流程优先使用多工具中位数；
CCF 保留所有工具值、范围和一致性状态，强冲突时自动降为低置信度。

手动调试 variant peptide extraction：

```bash
neoag extract-variant-peptides \
  --input-vcf /path/to/sample.vep.annotated.vcf.gz \
  --output results/SAMPLE001_sliding/upstream/tools/variant_peptides.tsv \
  --sample-id SAMPLE001 \
  --lengths 8,9,10,11 \
  --mini-len 27 \
  --hla-alleles HLA-A*02:01,HLA-B*07:02,HLA-C*07:02 \
  --tumor-sample-name TUMOR \
  --normal-proteome-fasta /path/to/Homo_sapiens.GRCh38.pep.all.fa \
  --filter-normal-proteome
```

没有授权 predictor 的 smoke test 可以在 TOML 中设置 `stub = true`，或在直接调用 `run` 时添加 `--immunogenicity-stub`。生产排序应使用真实 NetMHCpan/MHCflurry 输出以及真实 normal-expression/normal-ligand evidence，不应使用 fixture resources。

### 构建标准证据 sidecars

```bash
neoag build-evidence-layer \
  --outdir results/sample \
  --profile default \
  --sample-id SAMPLE001 \
  --raw-events results/sample/parsed/raw_events.tsv \
  --raw-peptides results/sample/parsed/raw_peptides.tsv \
  --expression results/sample/parsed/gene_expression.tsv \
  --rna-vaf results/sample/parsed/rna_vaf.tsv \
  --rna-junction results/sample/parsed/rna_junctions.tsv \
  --fusion-evidence results/sample/parsed/fusion_evidence.tsv \
  --normal-expression resources/normal_expression.example.tsv \
  --normal-hla-ligands resources/normal_hla_ligands.example.tsv
```

### HLA LOH 转换与交叉检查

```bash
neoag convert-lohhla \
  -i results/sample/tools/LOHHLA.HLAlossPrediction_CI.xls \
  -o results/sample/tools/lohhla.hla_loh.tsv

neoag convert-spechla \
  -i results/sample/tools/merge.hla.copy.txt \
  -o results/sample/tools/spechla.hla_loh.tsv

neoag crosscheck-hla-loh \
  --lohhla-hla-loh results/sample/tools/lohhla.hla_loh.tsv \
  --spechla-hla-loh results/sample/tools/spechla.hla_loh.tsv \
  --out results/sample/tools/hla_loh.crosscheck.tsv \
  --consensus-out results/sample/tools/hla_loh.consensus.tsv
```

### 生成报告

生成默认综合报告，同时输出患者沟通版和科研技术版报告：

```bash
neoag report \
  --profile default \
  --ranked-events results/sample/scoring/ranked_events.tsv \
  --ranked-peptides results/sample/scoring/ranked_peptides.tsv \
  --appm-summary results/sample/appm/appm_summary.tsv \
  --validation-plan results/sample/scoring/validation_plan.tsv \
  --outdir results/sample \
  --audience both \
  --out results/sample/reports/evidence_report.html
```

### Nextflow Fixture 运行

建议使用项目 wrapper，不要直接调用 `nextflow`。该 wrapper 会优先使用当前 checkout 的 `bin/neoag`，设置项目路径，并避免将 Nextflow metadata 写入 root-owned 位置。

```bash
export NXF_HOME=/path/to/writable/nextflow_cache
bin/neoag-nextflow -version
bin/neoag-nextflow run workflows/main.nf \
  -w /tmp/neoag_nf_work \
  --pvac_files data/fixtures/pvacseq_aggregated.tsv \
  --outdir results/demo_nf \
  --sample_id NF_DEMO
```

查看某个命令的完整参数：

```bash
neoag <command> --help
```

## Workflow 依赖矩阵

| Workflow / command | 最小输入 | 工具 | 参考库/数据 |
| --- | --- | --- | --- |
| Fixture demo: `neoag run-demo --outdir work/demo_v043 --sample-id DEMO001` | 内置 fixture | Python package 之外不需要 | 内置 fixtures/resources |
| Parsed pVAC results: `neoag run --outdir results/sample --sample-id SAMPLE001 --pvac data/fixtures/pvacseq_aggregated.tsv --immunogenicity-stub` | pVAC-like TSVs | 若输入已存在则不需外部工具 | 可选 normal expression/ligand tables |
| Raw intermediates: `neoag run --raw-events ... --raw-peptides ...` | `raw_events.tsv`, `raw_peptides.tsv` | 如提供则使用 NetMHCpan/MHCflurry 输出；可选 evidence tools | 可选 expression、LOH、purity、CNV、normal evidence |
| Full upstream run: `neoag run-full --config conf/run.sample.private.toml --outdir results/sample` | Run config | 取决于 enabled tools | 取决于 enabled tools |
| Binding prediction only: `peptide-predict` | Peptide/HLA table | 按选择使用 NetMHCpan、MHCflurry、PRIME/BigMHC/DeepImmuno | HLA alleles；predictor model data |
| VEP annotation: `vep-annotate` | VCF | VEP | VEP cache、reference FASTA、plugins |
| Variant peptide extraction: `extract-variant-peptides` | VEP-annotated VCF | Python；可选 VEP pre-step | Reference FASTA、可选 normal proteome |
| WES SNV calling: `snv-call-wes` | Tumor/normal BAM | GATK4 | GRCh38 FASTA、gnomAD AF VCF、PoN、intervals |
| WES SNV full: `snv-run-full-wes` | Somatic VCF 或 BAMs | BAM 模式需要 GATK；如启用则需 pVAC/binding tools | GRCh38 FASTA、HLA、可选 normal evidence |
| SV WGS raw build: `sv-build-raw` | SV VCF、FASTA、GTF、HLA | Python | Reference FASTA、GTF、HLA file |
| SV WES raw build: `sv-build-raw-wes` | SV VCF、FASTA、GTF、HLA、capture BED | Python | Reference FASTA、GTF、capture BED、HLA file |
| SV score: `sv-score` | Raw events/peptides | 除非 `--binding-stub`，否则需要 NetMHCpan/MHCflurry | HLA alleles、可选 evidence tables |
| Long-read SV wrapper | FASTQ/BAM 或 Sniffles2 VCF | 按选择使用 minimap2/samtools/Sniffles2 | Reference FASTA、GTF、HLA |
| Fusion discovery | FASTQ/BAM 或 caller outputs | 按选择使用 STAR-Fusion、FusionCatcher、Arriba、EasyFuse | CTAT/EasyFuse/fusion caller references |
| Immune escape evidence: `immune-escape` | Raw peptides、APPM/CCF/LOH evidence | 可选 LOHHLA/FACETS upstream | HLA LOH、CNV、VEP/APM/JAK/B2M evidence |
| Nextflow fixture | 内置 pVAC fixture | Java/Nextflow runtime | 内置 fixtures；可写 `NXF_HOME` |

## 测试

默认 pytest 只运行快速 unit tests：

```bash
pytest -q
```

显式运行更大的测试分组：

```bash
pytest -q --run-integration
pytest -q --run-benchmark
pytest -q --run-external
pytest -q --run-all
```

也支持 marker 形式：

```bash
pytest -q -m unit
pytest -q -m integration --run-integration
pytest -q -m benchmark --run-benchmark
pytest -q -m external --run-external
```

这样可以避免 lightweight release 用户在普通 `pytest` 下误运行耗时 Nextflow、benchmark 或 external-tool 测试。

## 常见错误与处理

| 现象 | 可能原因 | 处理 |
| --- | --- | --- |
| 找不到 `open-neo` | console script 还不在 `PATH` | 使用 Skill1 章节中的模块入口，然后重新运行 `open-neo install-check`。 |
| 安装、工具、参考库或 smoke 检查失败 | 机器还没有达到生产就绪 | 运行 `open-neo install-check --mode resume --approved`，查看 `deployment_report.md` 和 `production_run_readiness.tsv`。 |
| 样本执行中断或连接断开 | Gateway job 或生产阶段被中断 | 运行 `open-neo run --mode resume --approved --gateway-url http://127.0.0.1:8000 --gateway-wait`。 |
| 证据标记为 `UNASSESSED`、`PARTIAL` 或 `single_tool_result` | 某个输入、工具或参考库缺失，或只有一个工具产出可用证据 | 可以继续审阅，但报告和实验决策中必须保留部分证据标签。 |

## Release Boundary

不要提交或打包：

- `.git`, `.venv`, `.nextflow`, `.pytest_cache`
- `tools/`, `results/`, `work/`, `dist/`, `conda_packs/`
- `conf/tools.env.local.sh`
- `conf/site.config`
- `conf/private/*`
- `conf/*.private.toml`
- 真实患者数据或样本标识符
- 授权工具二进制文件
- `data/ref`、`data/vep` 等大型参考库

准备 online release 前，运行 `scripts/check_release_boundary.sh`。

## 更多文档

- `docs/V043_CCF21.md`：CCF 2.1 说明。
- `docs/V042_P1_APPM_EXPLAINABILITY.md`：APPM explainability 说明。
- `docs/V04_EVIDENCE_SAFETY_ESCAPE.md`：safety 和 immune-escape evidence 说明。
- `RELEASE.md`：release boundary 和测试总结。

## 解释边界

本流程用于研究分层和验证规划。排序候选需要结合 assay validation、疾病背景、HLA typing、tumor purity、expression/protein support、safety evidence、immune-escape context 以及适当的临床或科研治理流程进行复核。

### NetMHCpan 4.2c 容器运行环境

如果服务器因为缺少 `tcsh` 或 glibc 版本过低无法直接运行官方 NetMHCpan 4.2c，请使用 [docs/NETMHCPAN_CONTAINER.md](docs/NETMHCPAN_CONTAINER.md) 中的 Docker/Apptainer 运行方式。镜像只包含系统依赖，官方授权的 `tools/netMHCpan` 目录在运行时挂载。

### Priority tool containers

Docker/Apptainer runtimes for NetMHCpan, NetMHCstabpan, HLA-LA, SpecHLA, PURPLE/AMBER/COBALT, and EasyFuse are documented in [docs/PRIORITY_TOOL_CONTAINERS.md](docs/PRIORITY_TOOL_CONTAINERS.md). These images contain only runtime dependencies; licensed tools and large reference data are mounted from host paths.

## LLM-assisted Coordinator P1

This release adds an optional LLM-assisted Coordinator layer on top of the P0 Skills Pack. The default mode is dependency-free and rule-based; installing the optional `agent-llm` extra enables LiteLLM/LangGraph integration.

Plan only:

```bash
neoag-llm-agent --message "compare recommendation and NetMHCpan42 rankings" \
  --file ranked_peptides.recommendation.tsv \
  --file ranked_peptides.netmhcpan42.tsv \
  --outdir work/llm_plan --mode plan
```

Execute safe Skills:

```bash
neoag-llm-agent --message "compare recommendation and NetMHCpan42 rankings" \
  --file ranked_peptides.recommendation.tsv \
  --file ranked_peptides.netmhcpan42.tsv \
  --outdir work/llm_execute --mode execute-safe
```

Local Qwen/vLLM through LiteLLM/OpenAI-compatible API:

```bash
neoag-llm-agent --message "update patient report" \
  --file evidence_report.v04x_latest.html \
  --file ranked_peptides.recommendation.tsv \
  --file ranked_peptides.netmhcpan42.tsv \
  --outdir work/llm_report --mode execute-safe \
  --llm-provider litellm --model openai/qwen3-32b \
  --api-base http://localhost:8000/v1 --api-key-env LOCAL_VLLM_API_KEY
```

The Coordinator does not replace Project B CLI/Nextflow. It plans and calls registered Skills; high-impact operations such as HPC submission, installation, deletion, and overwrite require explicit approval.

See `docs/LLM_COORDINATOR_P1.md` and `docs/MODEL_API_AND_AGENT_FRAMEWORK_SELECTION.md`.

- [Tool inventory](docs/TOOL_INVENTORY.md): external tools, Docker images, environment variables, references, and licensing boundaries.

## Skills A/B/C/D 分层体系

本版本新增四类 Skills：

- **A 类入口适配型**：`neoag-vcf`、`neoag-fusion`、`neoag-splice`、`neoag-sv-wgs`、`neoag-sv-wes`、`neoag-peptide-csv`。
- **B 类公共证据分析型**：`neoag-hla-typing-loh`、`neoag-presentation`、`neoag-expression`、`neoag-rna-evidence`、`neoag-ccf`、`neoag-appm-escape`、`neoag-safety`、`neoag-ranking`。
- **C 类审阅/报告/实验设计型**：`neoag-ranking-compare`、`neoag-experiment-design`、`neoag-patient-report`、`neoag-technical-report`、`neoag-concept-explainer`。
- **D 类工程治理/执行控制型**：`neoag-input-qc`、`neoag-doctor`、`neoag-tool-reference-qc`、`neoag-run-demo-and-smoke`、`neoag-pipeline-full`、`neoag-release-qc`、`neoag-gateway-submit`、`neoag-hpc-runner`。

使用：

```bash
neoag-skill list
neoag-skill describe neoag-vcf
neoag-skill validate --root . --outdir work/skill_validate
neoag-skill run neoag-peptide-csv --outdir work/peptides --arg peptide_csv=peptides.tsv
```

Skills 是 SOP 封装，不承担临床决策，不包含患者 BAM/FASTQ/VCF 或大型参考库；HPC、安装、覆盖、删除等高风险路径默认 dry-run 或需要人工确认。

## Open-Neo 宏 Skills

新增三个稳定的对外入口，并继续复用现有 A/B/C/D Skills 与正式算法：

- `open-neo install-check`：迁移部署计划、经批准的安装/修复、Doctor、参考数据与 smoke 验收。
- `open-neo run`：manifest/CLI/目录输入质控、Gateway 受控生产执行或续跑、按证据域进行工具共识、综合证据表、旧加权基线与证据共识排序。
- `open-neo review`：事件级审阅、R1/R2 首批实验清单、R3 补证据队列及患者版/技术版报告。

当输入恰好是一对肿瘤 RNA FASTQ 时，`open-neo run` 现在会自动生成标准
`rna_fusion_splice_v1` production profile。通过 Doctor 检查并获得 Gateway
批准后，可依次运行 FASTQ QC、STAR、Salmon、EasyFuse、STAR-Fusion、Arriba、
RegTools、可选且经审核的 SNAF/SpliceMutr workflow、融合/剪接交叉验证、异常肽
生成及双排序。缺少可选工具或正常背景参考时只标记为 `UNASSESSED` 或
`SAFETY_PARTIAL`，不会解释为阴性证据。

`plan`、`verify` 和结果审阅不改动生产结果；安装、修复、执行和续跑必须显式批准。详见 [docs/OPEN_NEO_MACRO_SKILLS.md](docs/OPEN_NEO_MACRO_SKILLS.md)。
## 长读长 RNA splice 标准流程

长读长 RNA 的转录本重建、结构注释、ORF/蛋白翻译、SNAF/SpliceMutr 交叉验证、junction 肽生成和 MHC 预测统一使用机器可读配置 `configs/workflows/longrna_splice_profile.yaml`。大型参考数据和模型不提交到 Git，必须通过 reference/tools manifest 和资产同步流程解析。

初始化样本工作目录并检查流程计划：

```bash
bash scripts/run_longrna_splice.sh --dry-run \
  --sample-id SAMPLE \
  --input-dir /path/to/longrna_fastq \
  --workdir /path/to/work/SAMPLE/longrna_splice
```

正式阶段使用同一入口初始化标准目录，完成各阶段后运行验证：

```bash
bash scripts/run_longrna_splice.sh \
  --sample-id SAMPLE \
  --input-dir /path/to/longrna_fastq \
  --workdir /path/to/work/SAMPLE/longrna_splice
bash scripts/verify_longrna_splice.sh /path/to/work/SAMPLE/longrna_splice
```

最终证据必须区分“完整 ORF/蛋白已确认”和“仅有异常 junction 肽”。只有同时具备转录本 ID、SQANTI3/IsoQuant 结构、CDS 起止、frame、完整蛋白序列且蛋白序列跨越 junction，才允许标记为完整 ORF/蛋白确认。
