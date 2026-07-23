---
name: open-neo-run
description: Public macro Skill2 wrapper for the production evidence-consensus ranking CLI.
category: B - Public evidence analysis
risk_level: LOW
approval_required: false
---

# open-neo-run

## Purpose

Run the protected parallel weighted-baseline and evidence-consensus ranking.
This Skill is an SOP/adapter only. It calls `neoag evidence-rank`; all R1-R4,
hard-fail, priority-cap, Pareto, tie-break, and event-deduplication logic remains
in `src/neoag/evidence_consensus.py`.

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
neoag-skill run open-neo-run \
  --outdir work/open-neo-run \
  --arg comprehensive_evidence=results/scoring/comprehensive_peptide_evidence.tsv \
  --arg weighted_baseline=results/scoring/ranked_peptides.tsv
```

The consensus ranking is research-only and does not replace the primary
weighted ranking.
