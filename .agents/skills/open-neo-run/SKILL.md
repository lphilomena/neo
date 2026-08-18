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

Supported direct entries include tumor/normal DNA BAM or FASTQ, tumor RNA BAM or FASTQ, somatic VCF, fusion caller table, splice junction table, WGS/WES SV VCF, peptide-HLA table, standard raw intermediates, a production manifest, an input directory, or an existing result directory.

When raw DNA/RNA BAM or paired FASTQ inputs are supplied without an explicit
production manifest, Skill2 generates a capability-aware production profile.
For RNA it includes multi-batch paired FASTQ merging when more than one R1/R2
batch is supplied, FASTQ QC, STAR alignment, Salmon gene/transcript TPM plus
RSEM expression cross-check when a matching RSEM reference is available,
EasyFuse, STAR-Fusion, Arriba, RegTools,
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
`MIN_NHET=10`, `TARGET_ROWS=1000000`), Sequenza, and the repository-owned
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

## Modes

- `plan`: inspect and write a route/run plan.
- `dry-run`: add Doctor/preflight checks without heavy execution.
- `execute`: run the selected backend; explicit approval is required.
- `resume`: reuse completed stages and rerun missing/failed stages; explicit approval is required.
- `ranking-only`: reuse comprehensive evidence and the weighted baseline.

## Procedure

1. Detect inputs in deterministic order: manifest declarations, explicit CLI fill-ins, directory scanning, then extension/header inference. Validate non-empty files, HLA syntax, VCF samples/build, BAM indexes, capture BED and output writability.
2. Probe tool entrypoints, validated command templates, references, licenses and input compatibility. Write `capability_decisions.tsv`; PATH presence alone is not treated as a safe sample-level runner.
3. For tumor-normal BAM input, verify genotype identity with BAM-matcher when configured. Never use the bundled hg19 panel with GRCh38 BAMs.
4. Route to the existing fine-grained internal Skills and generate `capability_aware.production.toml` for raw inputs.
5. Run Doctor/preflight.
6. Use `pipeline-full` for the dry-run DAG and submit approved execute/resume requests through NeoAg Gateway to the production runner.
7. Reuse existing gene/transcript TPM and RNA alt/VAF tables, or plan/run Salmon/RSEM gene plus transcript quantification from tumor RNA FASTQ and RNA ref/alt counting from tumor RNA BAM plus somatic VCF. When multiple paired RNA FASTQ batches are supplied, merge R1 files and R2 files first and pass the merged pair to all downstream RNA tools. Retain fusion/splice junction read evidence.
   For the automatic RNA profile, require HLA, FASTA/GTF, STAR index,
   EasyFuse reference, CTAT library, Salmon index plus tx2gene or RSEM reference before execute.
   When Salmon and RSEM are both configured in auto mode, Salmon remains the
   primary expression handoff and RSEM is scheduled as an expression cross-check.
   SNAF and SpliceMutr are default splice stages. Their database/workflow
   assets must be present before execution; missing assets block the splice
   branch and are reported explicitly.
8. Cross-check HLA typing, LOHHLA/SpecHLA HLA LOH, fusion, splice, presentation and FACETS/Sequenza/PURPLE/ASCAT purity/CNV/CCF evidence by domain; missing evidence remains `UNASSESSED`.
9. Build `all_tool_results.tsv`, long-form tool evidence and explicit consensus/conflict outputs.
10. Preserve the weighted baseline, generate independent Evidence consensus rankings, compare both rankings, and write run/audit manifests.
11. Generate the technical Pipeline report by default. Do not generate a patient-facing report unless explicitly requested; the final patient report belongs to `open-neo-review`.

## Required outputs

- `input_status.json`, `route_plan.json`
- `capability_plan.json`, `capability_decisions.tsv`, `capability_aware.production.toml`
- `rna_preprocessing_status.tsv`, `rna_preprocessing_summary.json`
- `gene_tpm.tsv`, `transcript_tpm.tsv`, `rna_alt_vaf.tsv` when generated
- `all_tool_results.tsv`
- `tool_run_status.tsv`, `tool_consensus_summary.tsv`, `tool_evidence.long.tsv`
- `hla_typing_consensus.tsv`, `hla_loh_consensus.tsv`, `fusion_consensus.tsv`, `splice_consensus.tsv`
- `presentation_consensus.tsv`, `purity_cnv_consensus.tsv`, `ccf_consensus.tsv`
- `sample_identity_consensus.tsv` when tumor-normal genotype identity is assessed
- `evidence_conflicts.tsv`, `evidence_source_conflicts.tsv`
- `ranked_peptides.weighted_baseline.tsv`
- `ranked_peptides.evidence_consensus.tsv`
- `ranked_events.evidence_consensus.tsv`
- `ranking_compare_weighted_vs_consensus.md`
- `run_manifest.json`, `audit_log.jsonl`
- `reports/evidence_report.technical.html` as the default Pipeline report; an explicitly requested patient report is a non-final Pipeline snapshot
- `manifests/rna_fusion_splice.production.toml` and
  `manifests/rna_fusion_splice.requirements.tsv` for automatic RNA FASTQ runs

## Safety boundary

`execute` and `resume` require approval and Gateway dispatch. The Skill must not silently convert missing evidence to a negative result or overwrite the weighted baseline. Generated commands are restricted to repository-owned runners or administrator-reviewed `command_template` entries in the tools manifest.

The automatic RNA FASTQ profile is itself the reviewed repository-owned
production manifest generator. Execution still requires Gateway approval and
passes only when all required sample/reference assets exist. User-supplied
SNAF/SpliceMutr workflows are never invented or downloaded at run time.
