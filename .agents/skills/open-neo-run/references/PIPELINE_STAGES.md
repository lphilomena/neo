# Open-Neo run stages

Skill2 is an orchestration layer. Production algorithms remain in their
existing modules and fine-grained Skills.

1. Input inventory, schema validation and deterministic routing.
   Plan/dry-run records size and nanosecond mtime instead of hashing files at
   least 50 MiB, so WGS BAMs and large references are never fully read merely
   to construct a plan.
2. Approval/Gateway boundary for execute or resume.
3. Optional automatic `rna_fusion_splice_v1` production manifest generation.
4. Doctor and tool/reference preflight.
5. RNA preprocessing or reuse: gene TPM, transcript TPM and RNA alt/VAF.
6. Production DAG planning and execution with per-stage checkpoints. Raw DNA
   resume regenerates the same capability-aware manifest and reuses only
   signature-matching completed stages.
7. Candidate normalization for VCF, fusion, splice, SV or peptide entries.
8. Presentation, expression, RNA, HLA/APPM, CCF and safety evidence assembly.
9. Canonical `all_tool_results.tsv` and explicit conflicts/consensus tables.
10. Immutable weighted baseline plus independent Evidence-consensus ranking.
11. Output manifest, run manifest, `run_state.json` and audit log.

## Resource scheduler

Generated production manifests reserve `cpus` and `memory_gb` per stage and
define global `total_cpus`, `total_memory_gb`, and `max_parallel_stages`
budgets. A dependency-ready stage starts only when all three limits permit it.
The schedule is recorded in `production_resource_schedule.tsv`; active process
groups are terminated together when the parent run is interrupted, so resume
can safely reuse completed outputs and restart incomplete stages.

## Production Case Wrapper

For a completed case root plus somatic VCF, Skill2 should use
`scripts/run_production_case.sh` to convert standard upstream results into a
production manifest and run the final production runner. The wrapper requires
`--sample-id`, `--case-root`, `--outdir`, `--somatic-vcf`, and an approved local
NetMHCstabpan tree, then applies the sarcoma RNA-supported v2 weighted profile
and v3 Evidence-consensus rules by default. Optional overrides include
Sequenza/PURPLE paths, RNA FASTQ/BAM/VAF, STAR index, GTF, reference FASTA,
fusion caller roots, normal read-through/background assets, predictor dependency roots, NetMHCpan, NetMHCstabpan, PRIME/BigMHC/DeepImmuno evidence, and `--rna-threads`.
Exactly one RNA allele-evidence mode may be used: FASTQ pair, RNA BAM, or
existing RNA VAF. This path is for reusing completed upstream results and
producing final rankings/reports; raw heavy tool execution remains in the
Gateway-backed production DAG.

## RNA FASTQ profile

The built-in profile includes multi-batch paired FASTQ merging, FASTQ QC, STAR, Salmon gene/transcript TPM,
RSEM expression cross-check when a matching RSEM reference is available,
EasyFuse as the primary fusion meta-workflow, plus a provenance-tagged union of EasyFuse pass/unfiltered evidence, EasyFuse/standalone STAR-Fusion, Arriba, FusionCatcher, JAFFAL and completed caller roots when present, RegTools, optional reviewed SNAF/SpliceMutr workflows, fusion/splice
cross-validation, candidate peptide generation and downstream ranking.
Missing optional evidence remains `UNASSESSED` or `SAFETY_PARTIAL`.

## Purity/CNV

FACETS, Sequenza, PURPLE and ASCAT can run in parallel after BAM identity.
Sequenza is resume-safe: it can restart from an existing binned seqz, an
existing merged seqz, or completed per-chromosome seqz files before rerunning
only the remaining merge/binning/R-fit steps. `SEQUENZA_BIN_WINDOW` or
`BIN_WINDOW` controls the binning scale; the production default is 500.

## Resume

`run_state.json` records step status plus input and output signatures. A step is
reusable only when its previous status is reusable and all declared file
signatures still match. Large files use size and nanosecond mtime instead of rehashing.
Production stages additionally verify their declared outputs before reuse.

## Raw Heavy Finalization

Raw DNA/RNA Gateway DAGs must converge on the same final production standards as
`scripts/run_production_case.sh`: sarcoma RNA-supported v2 weighted profile,
v3 Evidence-consensus rules, required NetMHCpan/MHCflurry/NetMHCstabpan/NetChop
presentation predictors, optional PRIME/BigMHC/DeepImmuno evidence when configured, normal expression, normal ligandome, normal/reference proteome and pre-indexed normal junction backgrounds, and `patient,technical` reports. The raw DAG still
owns heavy upstream execution, but the generated production manifest should be
usable as the final production-results manifest shape once upstream stages have
completed. RNA allele evidence remains mutually exclusive across FASTQ-derived
STAR pileup, supplied RNA BAM pileup, and an existing RNA VAF table. Fusion evidence remains cumulative: preserve raw caller outputs and report whether each event is multi-caller, single-caller, or targeted-rescue supported.
