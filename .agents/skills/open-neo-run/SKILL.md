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
EasyFuse as the only scheduled fusion meta-workflow while preserving its embedded caller evidence, RegTools,
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
NetMHCstabpan, PRIME, BigMHC, DeepImmuno, RNA threads, and the global CPU,
memory and concurrent-stage budgets. When overrides are absent it must reuse the
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

## Controlled parallel scheduling

Raw production DAGs use dependency-aware, resource-budgeted parallel execution.
Set `total_cpus`, `total_memory_gb`, and `max_parallel_stages` in the request or
sample manifest. Skill2 writes those limits under `[run]` and writes `cpus` plus
`memory_gb` for every generated stage. Defaults are the requested RNA/thread
budget, 80% of currently available memory, and three concurrent stages.

Only dependency-ready stages whose combined reservations fit both global
budgets may start. Generated commands reduce explicit `--threads`/`--cores`
values to the stage CPU reservation, and common numerical-library thread
variables receive the same limit. Independent RNA expression, alignment and
fusion stages may therefore overlap without oversubscribing the machine.
Legacy manifests without scheduling fields remain serial. On SIGINT/SIGTERM,
the runner terminates every active stage process group; resume then reuses only
declared complete outputs.

## Fixed assets

Skill2 defaults to `open-neo/open-neo-public-assets` for redistributable fixed
assets. In `plan`/`dry-run` it only records whether the local Dataset marker is
ready. In an approved `execute`/`resume` (inside the Gateway for raw heavy runs)
it uses the official `hf download` CLI with a persistent local cache, resumes
missing chunks across interrupted runs, verifies archive checksums, extracts only when the
published archive identity changed, and then reuses the local tree. Select the
deployment location with `--public-asset-root` or
`OPEN_NEO_PUBLIC_ASSET_ROOT`; `--asset-root` and
`OPEN_NEO_REFERENCE_ROOT` are accepted as existing deployment locations. Use
`--public-asset-cache` for the resumable archive cache and
`--no-sync-public-assets` to require a pre-provisioned tree.
Large Dataset assets must not fall back to `curl`; a missing `hf` command is a
prerequisite error that should be repaired by Skill1 before resuming Skill2.

The Dataset excludes licensed HLA/LOH/presentation/predictor assets. Skill2
continues to resolve those from the machine-local licensed/tool configuration
created by Skill1 and reports them missing instead of attempting public
download.

## Procedure

1. Detect inputs in deterministic order: manifest declarations, explicit CLI fill-ins, directory scanning, then extension/header inference. Validate non-empty files, HLA syntax, VCF samples/build, BAM indexes, capture BED and output writability.
2. Probe tool entrypoints, validated command templates, references, licenses and input compatibility. Write `capability_decisions.tsv`; PATH presence alone is not treated as a safe sample-level runner. For SNV/InDel candidate generation, require VEP plugin support for MT/WT extraction (`Wildtype.pm` and `Frameshift.pm` via `refs.vep_plugins` or `NEOAG_VEP_PLUGINS`) whenever variant candidates are expected.
   Resolve and pass the deployed VEP cache root explicitly to every candidate-generation command. Accept a cache only when it contains `homo_sapiens/<version>_GRCh38`; do not allow an empty `~/.vep` to shadow it. For tumor-normal identity, use the repository `run_bam_matcher_pair.sh`, pin its Java runtime, verify BAM/FASTA/loci-VCF chromosome naming before launch, and reuse a completed non-empty identity TSV on resume.
3. For tumor-normal BAM input, verify genotype identity with BAM-matcher when configured. Never use the bundled hg19 panel with GRCh38 BAMs.
   Treat an explicit BAM-matcher mismatch or execution failure as a blocking sample-integrity error. If BAM-matcher, its GRCh38 loci panel, or a compatible Java runtime is unavailable, continue only with `sample_identity_status=UNASSESSED` and surface that limitation in evidence, ranking and both reports; never describe unavailable identity evidence as a match.
   For WGS/WES/PANEL DNA input, run the production DNA-SV branch before fusion integration. Require explicit tumor and normal VCF sample names for multisample VCFs and paired BAM calling; never infer sample order. Accept only strict `PASS` SV records by default. WES/PANEL runs must provide the assay-matched capture BED or stop before calling. Preserve caller-native Manta/SvABA/GRIDSS outputs in separate directories and normalize accepted records to an `adjacency_key` containing genome build, both breakpoints and orientations.
   Keep `DNA_SV_event` separate from `expressed_rearrangement_product`. DNA-only events may remain `ORF_HYPOTHESIS_UNRESOLVED`, but they must not directly generate neoantigen peptides. Link DNA and RNA evidence only by an exact orientation-compatible adjacency or by a unique, orientation-compatible transcript projection. Gene-pair-only matching is forbidden; ambiguous same-gene-pair adjacencies remain unlinked. Persist the DNA-RNA link table, link method, reason, caller support and split/discordant evidence into final evidence, ranking and Skill3 reports.
4. Route to the existing fine-grained internal Skills and generate `capability_aware.production.toml` for raw inputs.
5. Run Doctor/preflight.
6. Use `pipeline-full` for the dry-run DAG and submit approved execute/resume requests through NeoAg Gateway to the production runner. For existing completed case roots, prefer `scripts/run_production_case.sh` after explicit approval; it creates `manifest/production.results.toml` and invokes the production runner directly with fixed profiles and predictor environment pins. Raw heavy manifests must use the same final integration standards as the wrapper: `profiles/sarcoma_rna_supported_v2_provisional.toml`, `configs/ranking/sarcoma_evidence_consensus_v3_source_chain.toml`, required NetMHCpan/MHCflurry/NetMHCstabpan/NetChop presentation predictors, optional PRIME/MixMHCpred/BigMHC/DeepImmuno immunogenicity support when configured, and `patient,technical` reports.
For normal-background safety, prefer explicit inputs and otherwise pass through the fixed local asset root: indexed normal junctions, normal expression including GTEx/HSPC where available, normal HLA ligandome, and normal/reference proteome. Build or reuse the normal-junction sqlite index before splice review so background checks are streaming/pre-indexed rather than repeatedly scanning large tables. For normal-proteome safety, prefer explicit `NEOAG_NORMAL_PROTEOME_FASTA`/`NEOAG_NORMAL_PROTEOME`, then canonical asset names, then a non-empty FASTA discovered under `data/normal/proteome`; never treat a missing preferred filename as absence of the whole background.
For Fusion and Splice candidates, keep direct peptide/junction safety evidence separate from source-gene expression context. Direct safety evidence consists of the full-peptide normal-proteome exact match, normal transcript/junction match, normal immunopeptidome match, and similar-self-peptide cross-reactivity assessment. Normal-tissue, critical-organ, and HSPC expression of a fusion partner or splice-source gene are auxiliary context only: they must not independently reject a junction peptide or prove it safe. `peptide_safety.tsv` must carry the separate status fields plus `final_safety_conclusion` and direct `safety_evidence_completeness`. On resume, if an existing safety table lacks these fields, rerun the safety, Evidence consensus ranking, and report stages while reusing valid upstream caller outputs.
7. Reuse existing gene/transcript TPM and RNA alt/VAF tables, or plan/run Salmon/RSEM gene plus transcript quantification from tumor RNA FASTQ and RNA ref/alt counting from tumor RNA BAM plus somatic VCF. When multiple paired RNA FASTQ batches are supplied, merge R1 files and R2 files first and pass the merged pair to all downstream RNA tools. Retain fusion/splice junction read evidence. Enforce the wrapper-compatible RNA allele evidence rule for raw runs too: use only one of RNA FASTQ, RNA BAM, or an existing RNA VAF table. Before final ranking, verify that `ranked_peptides.evidence_consensus.tsv` carries concrete `gene_expression_tpm`, `transcript_expression_tpm`, `rna_depth`, `rna_alt_reads`, and `rna_vaf` values when the corresponding upstream evidence exists; do not let reports silently fall back to generic “未提供/未计算” text when a usable Salmon/RSEM table or RNA BAM/VAF table is present. For SNV/InDel ALT=0 rows, preserve coverage power as three distinct states: depth absent/0=`RNA_NO_COVERAGE_UNKNOWN`, 0<depth<10=`RNA_NO_ALT_LOW_COVERAGE`, and depth>=10=`RNA_NO_ALT_ADEQUATE_COVERAGE`. Only the last is significant negative RNA evidence; it must cap progression and must not be described as merely unassessed. For SNV/InDel rows, also verify `wildtype_peptide`, WT binding predictions and `mutant_specificity_status`/`mutant_specificity_state`; if these are absent, rerun VEP/peptide extraction with Wildtype/Frameshift plugins and do not treat the final report/Excel as complete.
   Before running presentation predictors on splice candidates, preserve the complete raw junction pool and apply the generic splice prefilter. Audit coordinate/alignment QC, explicit unique junction reads, total junction coverage, PSI, matched-normal support, GTEx/normal-cohort junction support, annotated normal isoforms, ORF/frame, NMD, true junction-spanning peptide generation and normal-proteome exact matches. Do not substitute gene TPM or caller aggregate reads for unique reads, total coverage or PSI. Reject explicit failures; retain missing-evidence candidates in a bounded `REVIEW` lane rather than silently treating them as PASS. Write `parsed/splice_prefilter_funnel.tsv`, `parsed/splice_prefilter_decisions.tsv`, and the unfiltered archive before NetMHCpan/MHCflurry/auxiliary prediction.
   Keep structural novelty, normal-catalog membership and tumor specificity as separate fields. An exact altered junction-spanning sequence may establish `ALTERED_JUNCTION_SPANNING_SEQUENCE`; absence from the recount3/GTEx presence catalog is only `NOT_LISTED_IN_NORMAL_CATALOG` and is never an adequate-coverage normal negative. If no assay-compatible normal RNA cohort or explicit local normal coverage table exists, set `UNASSESSED_NO_COMPATIBLE_NORMAL_RNA_COHORT`, limit normal safety to N1 and cap the final candidate at R3. Preserve the normal-catalog dataset/release provenance and observed read/sample denominators.
   For SpliceMutr with one tumor RNA sample and no matched/local normal RNA cohort, use the bundled GTEx cohort builder with at least 12 distinct-donor, assay-compatible RNA-seq samples selected from the installed GTEx junction matrix. Label every such normal as `normal_source=GTEx_V11` and `normal_match=PUBLIC_PROXY`, retain tissue, donor and RIN provenance, and run LeafCutter in `analysis_mode=outlier`. Do not run the replicated-group differential model with a one-sample tumor group and do not manufacture replicates from an aggregate catalog. Use `analysis_mode=differential` only when both biological groups have adequate independent replicates. Public-normal non-detection remains coverage-qualified proxy evidence, never a matched-normal negative.
   Treat the splice funnel as an ordered formal gate: exact junction -> explicit unique reads -> total junction coverage -> PSI -> coverage-qualified matched normal -> coverage-aware GTEx/normal cohort -> annotated normal isoform exclusion -> ORF/frame -> NMD -> residue-level junction-spanning peptide -> normal proteome and transcriptome exclusion -> HLA presentation. Only rows with `splice_formal_gate_pass=yes` enter `FORMAL_SPLICE_CANDIDATE`. Missing prerequisites may be predicted only in `EXPLORATION_EVIDENCE_INCOMPLETE`, remain capped at R3 and must not be counted as formal HLA-presentation passes.
   For the automatic RNA profile, require HLA, FASTA/GTF, STAR index,
   EasyFuse reference, Salmon index plus tx2gene or RSEM reference before execute. Fusion discovery schedules EasyFuse only; do not run duplicate standalone STAR-Fusion, Arriba or FusionCatcher jobs in addition to EasyFuse. For every patient, use the versioned fixed short-read panel `OPEN_NEO_SHORT_READ_FUSION_V1`: Arriba, STAR-Fusion and FusionCatcher as embedded caller identities, with EasyFuse PASS reported separately as the aggregator decision. A multi-caller conclusion requires the same normalized, orientation-aware adjacency to be supported by at least two fixed callers; gene-pair-only matching is forbidden. Emit `fusion_caller_availability.tsv` even when a caller produces zero calls, so a completed zero-call run is distinguishable from a missing caller. JAFFAL/other long-read calls are optional orthogonal support and never change the fixed-panel denominator. Build the formal event set as EasyFuse `fusions.pass.csv` union high-confidence rows from its complete `fusions.csv`. Internal or diagnostic rescue requires a valid frame and neo-peptide, excludes `cis_near`, remains explicitly `TARGETED_RESCUE`, and never silently becomes EasyFuse PASS or ordinary cross-validated consensus. Preserve `candidate_union_source`, fixed caller list/count, aggregator status, long-read support and rescue reason. `INTERNAL_CALLER_HIGH_CONFIDENCE` candidates trigger `CAP_CANDIDATE_UNION_INTERNAL_RESCUE` and cannot rank above R3 without independent or orthogonal confirmation.
   Pass the run-level recommended HLA file into fusion intermediate generation so junction-crossing 8-11mer candidates expand into peptide-HLA rows. The stage must fail when events exist but no peptide-HLA data rows are produced. A header-only peptide table is not successful and must not be reused during resume.
   A fusion peptide is eligible for the fusion-neoantigen lane only when residue-level mapping proves that it contains at least one amino acid from each fusion partner. Persist `junction_position_in_peptide_1based`, left/right gene and peptide segments, `fusion_junction_display`, transcript/ORF identifiers and `fusion_peptide_classification`. A caller label alone must not set `crosses_junction=yes`. One-sided peptides are excluded from the fusion-specific lane; missing boundary mapping is a blocking provenance failure. When one event yields near-identical peptide sequences, retain them as alternative transcript/ORF hypotheses and require exact patient transcript/frame adjudication rather than counting them as independent confirmed fusion neoantigens.
   Apply the Fusion ORF gate before formal candidate ranking. `ORF_SUPPORTED` requires an exact transcript and breakpoint pair, feasible translation start, resolved frame, translated ORF, exact peptide-in-ORF placement and NMD assessment. Preserve an otherwise plausible cross-junction candidate with an incomplete chain as `EXPLORATION_ORF_REQUIRED`, cap it at R3 and keep it out of the ORF-supported candidate layer. Invalid or contradictory ORF placement is `REJECTED_ORF_INVALID` and cannot be rescued by presentation scores.
   When Salmon and RSEM are both configured in auto mode, Salmon remains the
   primary expression handoff and RSEM is scheduled as an expression cross-check.
   SNAF and SpliceMutr are default splice stages. Their database/workflow
   assets must be present before execution; missing assets block the splice
   branch and are reported explicitly. Existing-results manifest generation
   must always write `manifest/splice_source_status.tsv` with separate rows for
   primary junction evidence, sNAF and SpliceMutr. A missing source is never
   treated as a negative result or silently omitted from the final status.
8. Cross-check HLA typing, LOHHLA/SpecHLA HLA LOH, fusion, splice, presentation and FACETS/Sequenza/PURPLE/ASCAT purity/CNV/CCF evidence by domain; missing evidence remains `UNASSESSED`. For every DNA-supported SNV/InDel, recompute CCF from tumor VAF or ALT/depth, consensus purity, the matched local total and allele-specific copy numbers, and biologically feasible mutation multiplicities. Propagate a Wilson 95% VAF interval through the CCF formula, retain multiplicity ambiguity and matched-normal ALT warnings, and label diploid fallback as low confidence. A high VAF or point CCF may be described only as compatible with clonality, never as proof. Fusion consensus is one authoritative adjacency-level union with provenance, not a mutually exclusive caller choice or a second gene-pair-only summary. Link Manta/SvABA/GRIDSS evidence after RNA consensus and write the DNA-SV status back into the published `fusion_consensus.tsv`; keep DNA-SV, long-read and fixed short-read support as separate columns. When using the production-case wrapper, preserve its single-RNA-input-mode rule: choose one of RNA FASTQ, RNA BAM, or existing RNA VAF.
9. Build `all_tool_results.tsv`, long-form tool evidence and explicit consensus/conflict outputs.
   Before accepting reused core presentation evidence, validate its schema and
   exact peptide-HLA key coverage. Stop when NetMHCpan/MHCflurry coverage
   collapses or supplied evidence would make most rows fall back to sentinel
   values (`rank=99`, `score=0`). Never substitute a final normalized table for
   a raw caller input unless the receiving field explicitly accepts normalized
   evidence.
   Compute auxiliary predictor coverage on unique applicable combinations.
   For NetMHCstabpan this means supported HLA-I 8-11mers; retain unsupported
   combinations as NOT_APPLICABLE, emit the missing applicable combinations
   for incremental resume, and merge newly completed rows with existing rows.
10. Preserve the weighted baseline, generate independent Evidence consensus rankings, compare both rankings, and write run/audit manifests. Treat `ranked_peptides.evidence_consensus.tsv` and `ranked_events.evidence_consensus.tsv` as the final report/Excel ranking inputs; use `ranked_peptides.tsv` only as the legacy weighted baseline or compatibility source.
    Preserve quantitative presentation fields in the final peptide table: NetMHCpan MT/WT EL and BA percentile ranks plus IC50 when emitted, MHCflurry MT/WT affinity percentile and presentation score, NetMHCstabpan MT/WT stability output, PRIME/BigMHC/DeepImmuno values, peptide length, mutation/junction position and computed MT/WT rank deltas or ratios. Record predictor allele support and extrapolation metadata when the tool provides it; otherwise keep these fields `UNASSESSED` rather than inferring support from a returned score.
    For every paired SNV/InDel peptide, also derive `mutation_position_role`, `mutation_position_interpretation`, `wt_self_reactivity_risk_status`, and `wt_self_reactivity_risk_reason`. Treat conventional MHC-I P2/P-omega changes as primary-anchor changes, internal non-primary-anchor changes only as *putative* TCR-facing positions, and keep MHC-II roles unresolved unless a binding register is available. A strong-binding WT peptide is an explicit self-reactivity/tolerance review signal and caps evidence-based progression; MT/WT rank difference, fold change or DAI must never independently establish immunogenicity.
11. Generate the technical Pipeline report by default. Do not generate a patient-facing report unless explicitly requested; the final patient report belongs to `open-neo-review`.
    When a patient snapshot is explicitly requested, default to 20 event-level representatives per applicable track and 100 cross-track peptide candidates. Preserve `event_top_n` and `candidate_top_n` as overridable run parameters and record them in the production manifest.
12. Run final regression gates: core predictor coverage must be non-zero, mass
    sentinel fallback is forbidden, and an unexpected all-R4 distribution is
    review-blocking until its evidence coverage is explained. Materialize only
    tracks present in canonical events. Preserve the Splice funnel and label
    PASS, REVIEW, FAIL and UNASSESSED separately; top-N REVIEW retention is not
    strict funnel validation.

## Required outputs

- `input_status.json`, `route_plan.json`
- `capability_plan.json`, `capability_decisions.tsv`, `capability_aware.production.toml`
- `manifest/production.results.toml` when using `scripts/run_production_case.sh`
- `rna_preprocessing_status.tsv`, `rna_preprocessing_summary.json`
- `manifest/splice_source_status.tsv`, including explicit AVAILABLE/MISSING status for primary junctions, sNAF and SpliceMutr
- `gene_tpm.tsv`, `transcript_tpm.tsv`, `rna_alt_vaf.tsv` when generated or discovered; final ranking tables must retain the joined expression and RNA VAF fields used by reports
- `wildtype_peptide`, WT binding columns and mutant-specificity status for SNV/InDel candidates in `ranked_peptides.evidence_consensus.tsv`; missing MT/WT evidence is a required gate failure, not a silent `UNASSESSED` final state
- `peptide_safety.tsv` with separate normal-proteome, normal transcript/junction, normal immunopeptidome, similar-peptide cross-reactivity, source-gene expression, critical-organ expression, hematopoietic expression, final safety conclusion and evidence-completeness fields
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
- `run_issue_log.json`, persisted beside the run state and appended on every execute/resume. Record the stage, symptom, root cause, failed command/log, temporary recovery, generic code/config fix, affected Skill1/2/3, validation result, and whether upstream outputs were safely reused.
- `production_resource_summary.json` and `production_resource_schedule.tsv`, recording global budgets, per-stage reservations and actual start/finish times.
- `reports/evidence_report.technical.html` as the default Pipeline report; an explicitly requested patient report is a non-final Pipeline snapshot
- `manifests/rna_fusion_splice.production.toml` and
  `manifests/rna_fusion_splice.requirements.tsv` for automatic RNA FASTQ runs

## Safety boundary

`execute` and `resume` require approval. Raw heavy execution requires Gateway dispatch. The existing-results production-case path may use the repository-owned `scripts/run_production_case.sh` wrapper directly after approval because it consumes completed upstream outputs, validates required files, and calls the production runner with pinned profiles. The Skill must not silently convert missing evidence to a negative result or overwrite the weighted baseline. Generated commands are restricted to repository-owned runners or administrator-reviewed `command_template` entries in the tools manifest.

Resume must reuse only outputs that pass declared existence and semantic checks. A prior `FAILED`, `BLOCKED`, or `LOW_CONFIDENCE` stage is not reusable merely because a partial file exists; retry it when the capability becomes available, write large outputs through temporary files plus atomic rename, and preserve the prior failure in `run_issue_log.json`.

The automatic RNA FASTQ profile is itself the reviewed repository-owned
production manifest generator. Execution still requires Gateway approval and
passes only when all required sample/reference assets exist. User-supplied
SNAF/SpliceMutr workflows are never invented or downloaded at run time.
