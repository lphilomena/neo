---
name: open-neo-run
description: Public macro Skill2 that detects Open-Neo inputs, routes VCF/fusion/splice/SV/peptide/result modes, plans or executes the pipeline, builds the all-tool evidence matrix, and emits weighted plus evidence-consensus rankings.
---

# Open-Neo Run

## Use when

- A code-capable agent receives a sample manifest or one or more supported input files.
- The user wants input checking, deterministic route selection, Pipeline execution, all-tool evidence integration, and two parallel rankings.
- Existing results need evidence-consensus ranking added without rerunning upstream tools.

## Inputs

Preferred: `sample_manifest.yaml`.

Supported direct entries include somatic VCF, fusion caller table, splice junction table, WGS/WES SV VCF, peptide-HLA table, standard raw intermediates, a production manifest, or an existing result directory.

## Modes

- `plan`: inspect and write a route/run plan.
- `dry-run`: add Doctor/preflight checks without heavy execution.
- `execute`: run the selected backend; explicit approval is required.
- `resume`: reuse completed stages and rerun missing/failed stages; explicit approval is required.
- `ranking-only`: reuse comprehensive evidence and the weighted baseline.

## Procedure

1. Detect and validate inputs; never rely on filename alone when table headers or manifests are available.
2. Route to the existing fine-grained internal Skills.
3. Run Doctor/preflight.
4. Use production runner/Nextflow/NeoAg CLI as the actual execution layer.
5. Build `all_tool_results.tsv` and tool-consensus/conflict outputs.
6. Preserve the weighted baseline.
7. Generate independent peptide- and event-level Evidence consensus rankings.
8. Compare both rankings and write run/audit manifests.

## Required outputs

- `input_status.json`, `route_plan.json`
- `all_tool_results.tsv`
- `ranked_peptides.weighted_baseline.tsv`
- `ranked_peptides.evidence_consensus.tsv`
- `ranked_events.evidence_consensus.tsv`
- `ranking_compare_weighted_vs_consensus.md`
- `run_manifest.json`, `audit_log.jsonl`

## Safety boundary

`execute` and `resume` require approval. The Skill must not silently convert missing evidence to a negative result, overwrite the weighted baseline, or combine raw DNA SV with other branches without a production manifest or prebuilt standard SV intermediates.
