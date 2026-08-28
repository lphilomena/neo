---
name: neoag-ranking
description: Thin Skill2 wrapper around the production neoag evidence-rank CLI.
category: B - Public evidence analysis
risk_level: LOW
approval_required: false
---

# neoag-ranking

This compatibility Skill delegates to `neoag evidence-rank`. It does not
implement scoring, R1-R4 assignment, Pareto ranking, or event deduplication.

## Required inputs

- `comprehensive_evidence`
- `weighted_baseline`

## Optional inputs

- `rules`
- `provenance`
- `track`

## Outputs

- `all_tool_results.tsv`
- `ranked_peptides.weighted_baseline.tsv`
- `ranked_peptides.evidence_consensus.tsv`
- `ranked_events.evidence_consensus.tsv`
- `ranking_compare_weighted_vs_consensus.md`

## Run

```bash
neoag-skill run neoag-ranking \
  --outdir work/neoag-ranking \
  --arg comprehensive_evidence=results/scoring/comprehensive_peptide_evidence.tsv \
  --arg weighted_baseline=results/scoring/ranked_peptides.tsv
```

The weighted baseline remains the primary compatibility ranking. The evidence
consensus is a parallel, research-only candidate-review output.
