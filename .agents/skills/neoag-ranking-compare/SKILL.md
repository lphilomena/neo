---
name: neoag-ranking-compare
description: Compare any two peptide rankings and audit rank shifts and evidence quality.
category: C - Review/report/design
risk_level: LOW
approval_required: false
---

# neoag-ranking-compare

## Required inputs

- `left`
- `right`

## Optional inputs

- `left_name`
- `right_name`

## Example

```bash
neoag-ranking-compare \
  --left ranked_peptides.weighted_baseline.tsv \
  --left-name weighted_baseline \
  --right ranked_peptides.evidence_consensus.tsv \
  --right-name evidence_consensus \
  --outdir comparison
```

Legacy `--netmhcpan42` and `--recommendation` arguments remain supported.

## Outputs

- `ranking_compare_report.md`
- `ranking_comparison_summary.json`
- `topn_overlap.tsv`
- `candidate_rank_changes.tsv`
- `high_rank_hard_fail.tsv`
- `top_composition.tsv`
- `evidence_qc_summary.tsv`
- `manual_review_candidates.tsv`

The comparison is an audit of ranking behavior, not a treatment decision.
