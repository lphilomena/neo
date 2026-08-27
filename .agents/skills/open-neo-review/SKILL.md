---
name: open-neo-review
description: Public macro Skill3 for read-only event-level integrity review, evidence-driven experiment priority, validation design, ranking comparison, and bounded patient/technical reporting.
---

# Open-Neo Review

The event-first selection rules, experiment recommendations, report boundaries, and `NEEDS_RANKING` behavior are defined in `references/REVIEW_RULES.md`. Public inputs and outputs must conform to the JSON schemas in this Skill directory.

## Use when

- `open-neo-run` has completed Evidence consensus ranking.
- The user asks which events deserve validation, why rankings differ, or which assay should be run first.
- Patient-facing and technical review artifacts are required.

## Required input

The only CLI entry is `result_dir`, but the directory must contain:

- `run_manifest.json`
- `ranked_events.evidence_consensus.tsv`
- `ranked_peptides.evidence_consensus.tsv`
- `ranked_peptides.weighted_baseline.tsv`
- `all_tool_results.tsv`
- `validation_plan.tsv`

Optional evidence includes APPM, HLA LOH, CCF, purity, peptide safety, conflicts, ranking comparison, fusion/splice consensus with caller provenance, RNA expression/RNA VAF provenance fields, normal-background safety status, `clinical_context.yaml`, and `disease_profile.yaml`.

Use `--reports patient,technical,onepage` to select report artifacts. Use
`--reports none` for event review and experiment tables without document/PPT
generation. Report generation is MEDIUM risk; table-only review is LOW risk.
Use `--event-top-n` and `--candidate-top-n` independently from `--top-n`.
The defaults are 20 event representatives per applicable track and 100
cross-track peptide candidates; `--top-n` continues to control the first-batch
experiment set only.

## Procedure

1. Verify required files, run identity, available input/reference hashes, Evidence consensus completion, event-peptide mapping, hard-fail propagation, and missing-evidence semantics.
2. Return `NEEDS_RANKING` when event-level consensus is absent. Never substitute weighted Top20.
3. Preserve `pipeline_r_grade` and `pipeline_event_rank`; write independent `review_status`, `review_reason`, and `experiment_priority` fields.
4. Review event-level representatives, with at most two peptide-HLA pairs per event and explicit phase/redundancy handling.
5. Invoke ranking comparison, experiment design, HLA-LOH/APPM review, CCF/clonality review, the shared formal patient-report renderer, technical report, and bounded concept explanations. For patient reports, require event rows to carry concrete RNA measurements from the canonical evidence tables: `gene_expression_tpm`, `transcript_expression_tpm`, `rna_depth`, `rna_alt_reads`, and `rna_vaf`. If any value is absent, report the upstream evidence gap or ID-mapping gap from `expression_source`, `transcript_expression_source`, `expression_evidence_status`, `rna_vaf_source`, and `rna_support_status` instead of a generic omission. Fusion and splice sections must consume Skill2/production-wrapper consensus outputs and caller-provenance columns; Skill3 must not independently rescan EasyFuse, STAR-Fusion, Arriba, FusionCatcher, JAFFAL, SNAF or SpliceMutr raw outputs. Render safety as separate full-peptide normal-proteome, normal transcript/junction, normal immunopeptidome, similar-peptide cross-reactivity, source-gene normal-expression, critical-organ expression and hematopoietic-expression dimensions, followed by a final conclusion and evidence completeness. For Fusion/Splice, label source-gene expression as auxiliary context and never present fusion-partner TPM/HSPC expression as evidence that the exact junction peptide occurs in normal cells.
   Add an event-level splice filtering funnel from Skill2 artifacts. Show raw/unified junction count, alignment/coordinate QC, explicit unique reads, total coverage, PSI, matched-normal and normal-cohort comparisons, annotated-normal-isoform exclusion, ORF, NMD, junction-spanning peptide, normal-proteome exclusion and HLA presentation. For each stage show entered, assessed, passed, failed and unassessed counts plus the possible survivor range; never infer missing stages from Top candidates. Also show CCF coverage by event track, separating reliable, low-confidence, unresolved and RNA-only-not-applicable events, and state how missing DNA-event CCF limits clonality and tumor-coverage interpretation. For SNV/InDel, display the CCF point estimate, 95% interval, local total/major/minor CN, mutation multiplicity assumption, purity source, matched-normal contamination status and confidence. Use `与克隆性相容` only when the interval and evidence quality support it; never convert high VAF alone into a clonal call.
   For every Fusion peptide, show `5′ partner residues | fusion junction | 3′ partner residues`, the junction position within the peptide, exact transcript and ORF identifiers, and whether both sides are present. Also show normal-proteome exact-match status. If boundary mapping is absent, state that fusion specificity is unproven; if the peptide is one-sided, do not call it a fusion neoantigen. Explain near-identical peptide alternatives as transcript/ORF/frame hypotheses requiring patient-specific adjudication.
   Add a quantitative presentation table for patient and technical reports. Show NetMHCpan MT/WT EL and BA percentile ranks and IC50, MHCflurry MT/WT affinity percentile and presentation score, NetMHCstabpan stability outputs, MT/WT deltas or ratios, PRIME/BigMHC/DeepImmuno values, peptide length and mutation/junction position. Show predictor allele support and extrapolation status only when explicitly recorded. State that lower percentile rank is stronger, and that a ≤1% rank is a research prioritization reference rather than proof of immunogenicity or in-vivo presentation.
   For SNV/InDel candidates, display the mutation's structural-role interpretation separately from MT/WT numeric differences: conventional MHC-I primary anchor, putative TCR-facing, mixed, or unresolved. Explicitly flag a WT peptide that still predicts as a strong binder to the same HLA as a self-reactivity and immune-tolerance review requirement. Do not describe marginal MT advantage, agretopicity or DAI as positive immunogenicity evidence by itself; state that position labels are sequence-based hypotheses unless supported by binding-register, structural or functional evidence.
6. Build a deterministic first-batch research set considering grade, RNA, safety, HLA diversity, clonality, event type, phase and redundancy. It is not a vaccine optimizer.
7. Generate short-peptide, long-peptide, minigene, targeted-RNA, and manual-review lanes.

## Outputs

- `review_integrity.json`, `review_integrity_checks.tsv`, `review_blocking_issues.tsv`
- `run_issue_log.json`, copied forward from Skill2 and appended with report/evidence defects, their source layer, generic renderer or ranking fix, whether Skill2 must be resumed, and post-regeneration validation status
- `candidate_review.tsv`
- `first_batch_experiment_set.tsv`
- `evidence_completion_queue.tsv`
- `manual_review_candidates.tsv`
- `experiment_candidates.tsv`
- `short_peptide_pool.tsv`
- `long_peptide_design.tsv`
- `minigene_design.tsv`
- `targeted_rna_validation_plan.tsv`
- APPM/HLA-LOH and CCF review files
- weighted-vs-consensus comparison files
- one formal patient report (`reports/patient_report.html`), rendered by `neoag.reports_dual.make_patient_report`, including per-candidate transcript expression, RNA depth/alt/VAF and event-type-aware layered safety evidence when available
- technical report (`md/html/docx` when available)
- `onepage_summary.pptx` when `python-pptx` is available

## Clinical boundary

Allowed wording: computational candidate, experiment priority, missing evidence, suggested research validation route.

Forbidden wording: confirmed neoantigen, guaranteed benefit, clinical resistance, ineffective immunotherapy, drug recommendation, or established vaccine/treatment plan.

Skill3 is read-only with respect to Skill2 outputs. It uses `ranked_peptides.evidence_consensus.tsv`, `ranked_events.evidence_consensus.tsv`, `fusion_consensus.tsv`, `splice_consensus.tsv`, and `all_tool_results.tsv` as the authoritative source of report content.
It must not create a second `production_patient/` report tree or maintain an independent simplified patient-report template.
When a report defect originates in missing or invalid upstream evidence, record it and return it to Skill2 instead of fabricating a display value. When the defect is renderer-only, regenerate Skill3 without altering Skill2 outputs and record the before/after integrity check.

Before selecting candidates or rendering reports, Skill3 must also:

- inspect R1-R4 counts and block review when all-R4 coincides with absent core
  NetMHCpan/MHCflurry coverage;
- calculate tool coverage over unique applicable peptide-HLA combinations and
  distinguish missing, unsupported and failed results;
- show only event tracks present in the canonical event ranking;
- read the Splice prefilter funnel when Splice is present and explicitly
  distinguish strict PASS from bounded REVIEW and UNASSESSED stages.
