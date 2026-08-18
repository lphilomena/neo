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

## RNA FASTQ profile

The built-in profile includes multi-batch paired FASTQ merging, FASTQ QC, STAR, Salmon gene/transcript TPM,
RSEM expression cross-check when a matching RSEM reference is available,
EasyFuse, STAR-Fusion, FusionCatcher, Arriba, RegTools, optional reviewed SNAF/SpliceMutr workflows, fusion/splice
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
