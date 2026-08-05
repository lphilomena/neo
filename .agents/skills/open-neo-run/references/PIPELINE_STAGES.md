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

The built-in profile includes FASTQ QC, STAR, Salmon, EasyFuse, STAR-Fusion,
Arriba, RegTools, optional reviewed SNAF/SpliceMutr workflows, fusion/splice
cross-validation, candidate peptide generation and downstream ranking.
Missing optional evidence remains `UNASSESSED` or `SAFETY_PARTIAL`.

## Resume

`run_state.json` records step status plus input and output signatures. A step is
reusable only when its previous status is reusable and all declared file
signatures still match. Large files use size and nanosecond mtime instead of rehashing.
Production stages additionally verify their declared outputs before reuse.
