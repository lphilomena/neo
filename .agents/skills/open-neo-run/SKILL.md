---
name: open-neo-run
description: Public macro Skill2 that detects Open-Neo inputs, routes VCF/fusion/splice/SV/peptide/result modes, plans or executes the pipeline, builds the all-tool evidence matrix, and emits weighted plus evidence-consensus rankings.
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

When exactly one paired tumor RNA FASTQ input is supplied without DNA inputs
or an explicit production manifest, Skill2 automatically generates the
`rna_fusion_splice_v1` production profile. The profile includes FASTQ QC, STAR
alignment, Salmon gene/transcript TPM, EasyFuse, STAR-Fusion, Arriba, RegTools,
optional SNAF/SpliceMutr, cross-tool splice normalization, fusion/splice peptide
generation, presentation, evidence integration, and dual ranking.

## Modes

- `plan`: inspect and write a route/run plan.
- `dry-run`: add Doctor/preflight checks without heavy execution.
- `execute`: run the selected backend; explicit approval is required.
- `resume`: reuse completed stages and rerun missing/failed stages; explicit approval is required.
- `ranking-only`: reuse comprehensive evidence and the weighted baseline.

## Procedure

1. Detect inputs in deterministic order: manifest declarations, explicit CLI fill-ins, directory scanning, then extension/header inference. Validate non-empty files, HLA syntax, VCF samples/build, BAM indexes, capture BED and output writability.
2. Route to the existing fine-grained internal Skills.
3. Run Doctor/preflight.
4. Use `pipeline-full` for the dry-run DAG and submit approved execute/resume requests through NeoAg Gateway to the production runner.
5. Reuse existing gene/transcript TPM and RNA alt/VAF tables, or plan/run Salmon/RSEM gene plus transcript quantification from tumor RNA FASTQ and RNA ref/alt counting from tumor RNA BAM plus somatic VCF. Retain fusion/splice junction read evidence.
   For the automatic RNA profile, require HLA, FASTA/GTF, STAR index,
   EasyFuse reference, CTAT library, Salmon index, and tx2gene before execute.
   SNAF and SpliceMutr remain optional cohort/workflow stages; when their
   reviewed workflow files are absent they are `UNASSESSED`, not negative.
6. Cross-check HLA typing, HLA LOH, fusion, splice, presentation and purity/CNV/CCF evidence by domain; missing evidence remains `UNASSESSED`.
7. Build `all_tool_results.tsv`, long-form tool evidence and explicit consensus/conflict outputs.
8. Preserve the weighted baseline.
9. Generate independent peptide- and event-level Evidence consensus rankings.
10. Compare both rankings and write run/audit manifests.

## Required outputs

- `input_status.json`, `route_plan.json`
- `rna_preprocessing_status.tsv`, `rna_preprocessing_summary.json`
- `gene_tpm.tsv`, `transcript_tpm.tsv`, `rna_alt_vaf.tsv` when generated
- `all_tool_results.tsv`
- `tool_run_status.tsv`, `tool_consensus_summary.tsv`, `tool_evidence.long.tsv`
- `hla_typing_consensus.tsv`, `hla_loh_consensus.tsv`, `fusion_consensus.tsv`, `splice_consensus.tsv`
- `presentation_consensus.tsv`, `purity_cnv_consensus.tsv`, `ccf_consensus.tsv`
- `evidence_conflicts.tsv`, `evidence_source_conflicts.tsv`
- `ranked_peptides.weighted_baseline.tsv`
- `ranked_peptides.evidence_consensus.tsv`
- `ranked_events.evidence_consensus.tsv`
- `ranking_compare_weighted_vs_consensus.md`
- `run_manifest.json`, `audit_log.jsonl`
- `manifests/rna_fusion_splice.production.toml` and
  `manifests/rna_fusion_splice.requirements.tsv` for automatic RNA FASTQ runs

## Safety boundary

`execute` and `resume` require approval and Gateway dispatch. The Skill must not silently convert missing evidence to a negative result or overwrite the weighted baseline. Raw candidate generation from BAM/FASTQ requires either an explicit production manifest with reviewed stage commands or the repository-owned automatic RNA FASTQ profile; the declared RNA quantification and RNA allele-count stages are separately controlled by Gateway.

The automatic RNA FASTQ profile is itself the reviewed repository-owned
production manifest generator. Execution still requires Gateway approval and
passes only when all required sample/reference assets exist. User-supplied
SNAF/SpliceMutr workflows are never invented or downloaded at run time.
