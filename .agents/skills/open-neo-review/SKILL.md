---
name: open-neo-review
description: Public macro Skill3 that reads event-level Evidence consensus results, compares rankings, creates an event-deduplicated experimental priority set, and generates patient and technical reports without changing Pipeline ranks.
---

# Open-Neo Review

## Use when

- `open-neo-run` has produced peptide- and event-level Evidence consensus results.
- The user wants to know which events should be validated first and why.
- A patient-facing or technical review report is required.

## Required input

- `result_dir` containing `ranked_events.evidence_consensus.tsv` and `ranked_peptides.evidence_consensus.tsv`. The weighted baseline is optional and enables comparison output.

## Procedure

1. Verify that all files refer to a coherent result set.
2. Use the event-level consensus table as the primary review input.
3. Keep no more than one or two representative peptide-HLA pairs per event/phase group.
4. Map R1-R4 and missing evidence to an independent experiment-priority state.
5. Invoke the internal experiment-design and ranking-compare Skills.
6. Create a deterministic R1/R2 first-batch set with event, HLA, and assay diversity; keep R3 in a separate evidence-completion queue.
7. Generate patient and technical reports while preserving research boundaries.

## Outputs

- `candidate_review.tsv`
- `first_batch_experiment_set.tsv`
- `experiment_design/experiment_candidates.tsv`
- short-peptide, long-peptide, minigene, and targeted-RNA plans
- patient report (`md/html/docx` when python-docx is available)
- technical report (`md/html`)

## Boundary

The first-batch set is a transparent heuristic for wet-lab planning, not a validated vaccine-set optimizer. This Skill must never rewrite the Pipeline's R grade, weighted rank, or Evidence consensus rank.
