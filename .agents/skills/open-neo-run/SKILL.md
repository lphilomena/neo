---
name: open-neo-run
description: Public macro Skill2 that detects raw DNA/RNA and processed inputs, probes callable tools/references, generates an auditable production DAG, executes approved primary and cross-validation tools, and emits weighted plus evidence-consensus rankings.
---

# Open-Neo Run

The canonical stage order, critical/evidence-missing behavior, RNA FASTQ production profile, and resume contract are defined in `references/PIPELINE_STAGES.md`. Validate requests with `references/INPUT_SCHEMA.json` and preserve `references/OUTPUT_SCHEMA.json`.

`--mode resume` reads the prior `run_state.json`. A completed macro step is reusable only when its recorded input and output signatures still match; production-stage reuse remains delegated to the production runner's declared-output checks.

## Use when

- A code-capable agent receives a sample manifest or one or more supported input files.
- The user wants input checking, deterministic route selection, Pipeline execution, all-tool evidence integration, and two parallel rankings.
- Existing results need evidence-consensus ranking added without rerunning upstream tools.

## Inputs

Preferred: `sample_manifest.yaml`.

Supported direct entries include tumor/normal DNA BAM or FASTQ, tumor RNA BAM or FASTQ, somatic VCF, fusion caller table, splice junction table, WGS/WES SV VCF, peptide-HLA table, standard raw intermediates, a production manifest, an input directory, an existing result directory, or a case-root directory containing completed upstream tool results.

When raw DNA/RNA BAM or paired FASTQ inputs are supplied without an explicit
production manifest, Skill2 generates a capability-aware production profile.
For RNA it includes multi-batch paired FASTQ merging when more than one R1/R2
batch is supplied, FASTQ QC, STAR alignment, Salmon gene/transcript TPM plus
RSEM expression cross-check when a matching RSEM reference is available,
EasyFuse as the primary fusion meta-workflow while preserving EasyFuse's embedded caller outputs, configured standalone STAR-Fusion/FusionCatcher/Arriba/JAFFAL outputs and completed caller roots as a provenance-tagged union, RegTools,
SNAF and SpliceMutr, cross-tool splice normalization, fusion/splice peptide
generation, presentation, evidence integration, and dual ranking. For DNA it
can include BWA/samtools alignment, Mutect2, OptiType or command-template HLA
callers, FACETS/Sequenza plus configured PURPLE/ASCAT, LOHHLA, VEP peptide
generation and unified ranking.

For a tumor-normal BAM pair, run BAM-matcher sample-identity QC before paired
variant, purity/CNV and HLA-LOH analyses when both the isolated legacy tool and
a GRCh38-compatible identity SNP panel are declared. `MISMATCH` is a hard
paired-analysis failure; low comparable-site coverage is `INSUFFICIENT_DATA`
and requires review rather than being treated as a mismatch.

With the default `all-available` policy, paired GRCh38 DNA automatically uses
the robust FACETS omni2p5 profile (`CVAL_PRE=50`, `CVAL_PROC=300`,
`MIN_NHET=10`, `TARGET_ROWS=1000000`), resume-safe Sequenza
(`SEQUENZA_BIN_WINDOW`, default `500`, with support for existing chr/merged/binned
seqz inputs), and the repository-owned
PURPLE runner when their validated assets are present. The purity/CNV review
emits recommended purity/ploidy and normalized CNV segments, plus a per-tool
6p21 MHC-region CNV/LOH cross-check. ASCAT remains conditional on a validated
sample-level command template.

When LOHHLA and SpecHLA are available, both HLA LOH methods are run with the
same recommended purity/ploidy. LOHHLA calls require both BAF-informed copy
number `<0.5` and paired P value `<0.01`; the raw P values, copy numbers,
confidence intervals, coverage support, and SpecHLA call evidence are retained
in the two-field allele-level cross-check. Consensus states are
`CONSENSUS_LOST`, `CONSENSUS_RETAINED`, `DISCORDANT`, and `UNASSESSED`.


When a request supplies a completed `case_root` plus `somatic_vcf`, Skill2 uses
the repository-owned `scripts/run_production_case.sh` wrapper as the preferred
production-case entrypoint. This wrapper discovers standard completed outputs
under the case root and output tree, generates `manifest/production.results.toml`
with `scripts/generate_production_from_results_manifest.py`, then runs
`neoag.production_runner --execute` with the sarcoma RNA-supported v2 weighted
profile and v3 Evidence-consensus rules. It accepts explicit overrides for
Sequenza/PURPLE, gene TPM, transcript TPM/`quant.sf`, RNA FASTQ/BAM/VAF,
STAR/GTF/reference assets, fusion caller roots, normal read-through, normal expression/junction/ligandome/proteome assets through `--asset-root`, predictor dependency roots, NetMHCpan,
NetMHCstabpan, PRIME, BigMHC, DeepImmuno, and RNA threads. When overrides are absent it must reuse the
latest non-empty `gene_tpm.tsv`, `transcript_tpm.tsv`, Salmon `quant.sf`, RSEM
`*.genes.results`/`*.isoforms.results`, existing `rna_alt_vaf.tsv`, or a
non-empty STAR RNA BAM so patient reports can show transcript expression, RNA
site depth, RNA alt reads, and RNA VAF. The wrapper is for existing upstream
results and ranking/report production; it must not replace the raw-input Gateway
DAG for new heavy DNA/RNA tool execution.

## Modes

- `plan`: inspect and write a route/run plan.
- `dry-run`: add Doctor/preflight checks without heavy execution.
- `execute`: run the selected backend; explicit approval is required.
- `resume`: reuse completed stages and rerun missing/failed stages; explicit approval is required.
- `ranking-only`: reuse comprehensive evidence and the weighted baseline.

## Procedure

1. Detect inputs in deterministic order: manifest declarations, explicit CLI fill-ins, directory scanning, then extension/header inference. Validate non-empty files, HLA syntax, VCF samples/build, BAM indexes, capture BED and output writability.
2. Probe tool entrypoints, validated command templates, references, licenses and input compatibility. Write `capability_decisions.tsv`; PATH presence alone is not treated as a safe sample-level runner. For SNV/InDel candidate generation, require VEP plugin support for MT/WT extraction (`Wildtype.pm` and `Frameshift.pm` via `refs.vep_plugins` or `NEOAG_VEP_PLUGINS`) whenever variant candidates are expected.
3. For tumor-normal BAM input, verify genotype identity with BAM-matcher when configured. Never use the bundled hg19 panel with GRCh38 BAMs.
4. Route to the existing fine-grained internal Skills and generate `capability_aware.production.toml` for raw inputs.
5. Run Doctor/preflight.
6. Use `pipeline-full` for the dry-run DAG and submit approved execute/resume requests through NeoAg Gateway to the production runner. For existing completed case roots, prefer `scripts/run_production_case.sh` after explicit approval; it creates `manifest/production.results.toml` and invokes the production runner directly with fixed profiles and predictor environment pins. Raw heavy manifests must use the same final integration standards as the wrapper: `profiles/sarcoma_rna_supported_v2_provisional.toml`, `configs/ranking/sarcoma_evidence_consensus_v3_source_chain.toml`, required NetMHCpan/MHCflurry/NetMHCstabpan/NetChop presentation predictors, optional PRIME/MixMHCpred/BigMHC/DeepImmuno immunogenicity support when configured, and `patient,technical` reports.
For normal-background safety, prefer explicit inputs and otherwise pass through the fixed local asset root: indexed normal junctions, normal expression including GTEx/HSPC where available, normal HLA ligandome, and normal/reference proteome. Build or reuse the normal-junction sqlite index before splice review so background checks are streaming/pre-indexed rather than repeatedly scanning large tables. For normal-proteome safety, prefer explicit `NEOAG_NORMAL_PROTEOME_FASTA`/`NEOAG_NORMAL_PROTEOME`, then canonical asset names, then a non-empty FASTA discovered under `data/normal/proteome`; never treat a missing preferred filename as absence of the whole background.
7. Reuse existing gene/transcript TPM and RNA alt/VAF tables, or plan/run Salmon/RSEM gene plus transcript quantification from tumor RNA FASTQ and RNA ref/alt counting from tumor RNA BAM plus somatic VCF. When multiple paired RNA FASTQ batches are supplied, merge R1 files and R2 files first and pass the merged pair to all downstream RNA tools. Retain fusion/splice junction read evidence. Enforce the wrapper-compatible RNA allele evidence rule for raw runs too: use only one of RNA FASTQ, RNA BAM, or an existing RNA VAF table. Before final ranking, verify that `ranked_peptides.evidence_consensus.tsv` carries concrete `gene_expression_tpm`, `transcript_expression_tpm`, `rna_depth`, `rna_alt_reads`, and `rna_vaf` values when the corresponding upstream evidence exists; do not let reports silently fall back to generic “未提供/未计算” text when a usable Salmon/RSEM table or RNA BAM/VAF table is present. For SNV/InDel rows, also verify `wildtype_peptide`, WT binding predictions and `mutant_specificity_status`/`mutant_specificity_state`; if these are absent, rerun VEP/peptide extraction with Wildtype/Frameshift plugins and do not treat the final report/Excel as complete.
   For the automatic RNA profile, require HLA, FASTA/GTF, STAR index,
   EasyFuse reference, Salmon index plus tx2gene or RSEM reference before execute. EasyFuse remains the primary meta-workflow, but final fusion evidence must be built from the union of EasyFuse pass/unfiltered tables, EasyFuse/standalone STAR-Fusion, Arriba, FusionCatcher, JAFFAL and any completed caller roots that are present, with caller provenance and single-caller/targeted-rescue status retained.
   When Salmon and RSEM are both configured in auto mode, Salmon remains the
   primary expression handoff and RSEM is scheduled as an expression cross-check.
   SNAF and SpliceMutr are default splice stages. Their database/workflow
   assets must be present before execution; missing assets block the splice
   branch and are reported explicitly.
8. Cross-check HLA typing, LOHHLA/SpecHLA HLA LOH, fusion, splice, presentation and FACETS/Sequenza/PURPLE/ASCAT purity/CNV/CCF evidence by domain; missing evidence remains `UNASSESSED`. Fusion consensus is a union with provenance, not a mutually exclusive caller choice. When using the production-case wrapper, preserve its single-RNA-input-mode rule: choose one of RNA FASTQ, RNA BAM, or existing RNA VAF.
9. Build `all_tool_results.tsv`, long-form tool evidence and explicit consensus/conflict outputs.
10. Preserve the weighted baseline, generate independent Evidence consensus rankings, compare both rankings, and write run/audit manifests. Treat `ranked_peptides.evidence_consensus.tsv` and `ranked_events.evidence_consensus.tsv` as the final report/Excel ranking inputs; use `ranked_peptides.tsv` only as the legacy weighted baseline or compatibility source.
11. Generate the technical Pipeline report by default. Do not generate a patient-facing report unless explicitly requested; the final patient report belongs to `open-neo-review`.

## Required outputs

- `input_status.json`, `route_plan.json`
- `capability_plan.json`, `capability_decisions.tsv`, `capability_aware.production.toml`
- `manifest/production.results.toml` when using `scripts/run_production_case.sh`
- `rna_preprocessing_status.tsv`, `rna_preprocessing_summary.json`
- `gene_tpm.tsv`, `transcript_tpm.tsv`, `rna_alt_vaf.tsv` when generated or discovered; final ranking tables must retain the joined expression and RNA VAF fields used by reports
- `wildtype_peptide`, WT binding columns and mutant-specificity status for SNV/InDel candidates in `ranked_peptides.evidence_consensus.tsv`; missing MT/WT evidence is a required gate failure, not a silent `UNASSESSED` final state
- `all_tool_results.tsv`
- `tool_run_status.tsv`, `tool_consensus_summary.tsv`, `tool_evidence.long.tsv`
- `hla_typing_consensus.tsv`, `hla_loh_consensus.tsv`, `fusion_consensus.tsv`, `splice_consensus.tsv`
- `presentation_consensus.tsv`, `purity_cnv_consensus.tsv`, `ccf_consensus.tsv`
- `sample_identity_consensus.tsv` when tumor-normal genotype identity is assessed
- `evidence_conflicts.tsv`, `evidence_source_conflicts.tsv`
- `ranked_peptides.evidence_consensus.tsv` as the final report/Excel ranking table with concrete `evidence_grade`/R grade fields
- `ranked_events.evidence_consensus.tsv` as the final event ranking table
- `ranked_peptides.weighted_baseline.tsv` as the preserved legacy weighted baseline for audit/comparison only
- `ranking_compare_weighted_vs_consensus.md`
- `run_manifest.json`, `audit_log.jsonl`
- `reports/evidence_report.technical.html` as the default Pipeline report; an explicitly requested patient report is a non-final Pipeline snapshot
- `manifests/rna_fusion_splice.production.toml` and
  `manifests/rna_fusion_splice.requirements.tsv` for automatic RNA FASTQ runs

## Safety boundary

`execute` and `resume` require approval. Raw heavy execution requires Gateway dispatch. The existing-results production-case path may use the repository-owned `scripts/run_production_case.sh` wrapper directly after approval because it consumes completed upstream outputs, validates required files, and calls the production runner with pinned profiles. The Skill must not silently convert missing evidence to a negative result or overwrite the weighted baseline. Generated commands are restricted to repository-owned runners or administrator-reviewed `command_template` entries in the tools manifest.

The automatic RNA FASTQ profile is itself the reviewed repository-owned
production manifest generator. Execution still requires Gateway approval and
passes only when all required sample/reference assets exist. User-supplied
SNAF/SpliceMutr workflows are never invented or downloaded at run time.
