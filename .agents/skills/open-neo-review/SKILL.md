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

## Procedure

1. Verify required files, run identity, available input/reference hashes, Evidence consensus completion, event-peptide mapping, hard-fail propagation, and missing-evidence semantics.
2. Return `NEEDS_RANKING` when event-level consensus is absent. Never substitute weighted Top20.
3. Preserve `pipeline_r_grade` and `pipeline_event_rank`; write independent `review_status`, `review_reason`, and `experiment_priority` fields.
4. Review event-level representatives, with at most two peptide-HLA pairs per event and explicit phase/redundancy handling.
5. Invoke ranking comparison, experiment design, HLA-LOH/APPM review, CCF/clonality review, the shared formal patient-report renderer, technical report, and bounded concept explanations. For patient reports, require event rows to carry concrete RNA measurements from the canonical evidence tables: `gene_expression_tpm`, `transcript_expression_tpm`, `rna_depth`, `rna_alt_reads`, and `rna_vaf`. If any value is absent, report the upstream evidence gap or ID-mapping gap from `expression_source`, `transcript_expression_source`, `expression_evidence_status`, `rna_vaf_source`, and `rna_support_status` instead of a generic omission. Fusion and splice sections must consume Skill2/production-wrapper consensus outputs and caller-provenance columns; Skill3 must not independently rescan EasyFuse, STAR-Fusion, Arriba, FusionCatcher, JAFFAL, SNAF or SpliceMutr raw outputs.
6. Build a deterministic first-batch research set considering grade, RNA, safety, HLA diversity, clonality, event type, phase and redundancy. It is not a vaccine optimizer.
7. Generate short-peptide, long-peptide, minigene, targeted-RNA, and manual-review lanes.

## Outputs

- `review_integrity.json`, `review_integrity_checks.tsv`, `review_blocking_issues.tsv`
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
- one formal patient report (`reports/patient_report.html`), rendered by `neoag.reports_dual.make_patient_report`, including per-candidate transcript expression and RNA depth/alt/VAF when available
- technical report (`md/html/docx` when available)
- `onepage_summary.pptx` when `python-pptx` is available

## Clinical boundary

Allowed wording: computational candidate, experiment priority, missing evidence, suggested research validation route.

Forbidden wording: confirmed neoantigen, guaranteed benefit, clinical resistance, ineffective immunotherapy, drug recommendation, or established vaccine/treatment plan.

Skill3 is read-only with respect to Skill2 outputs. It uses `ranked_peptides.evidence_consensus.tsv`, `ranked_events.evidence_consensus.tsv`, `fusion_consensus.tsv`, `splice_consensus.tsv`, and `all_tool_results.tsv` as the authoritative source of report content.
It must not create a second `production_patient/` report tree or maintain an independent simplified patient-report template.
