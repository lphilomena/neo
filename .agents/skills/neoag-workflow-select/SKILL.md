---
name: neoag-workflow-select
description: Select a safe NeoAg workflow from sample, tool, and reference manifests. Use it before BAM/VCF/intermediate/result workflows or QC consensus runs. It writes a dry-run plan only and never invokes external tools.
---

# NeoAg Workflow Select

## Purpose

Turn available inputs into a deterministic, auditable workflow plan before any
heavy tool runs. This skill replaces ad hoc shell routing with manifest-driven
selection and explicit readiness, risk, and missing-evidence states.

## Inputs

- required `sample_manifest.yaml`;
- recommended `tools_manifest.yaml` and `reference_manifest.yaml`;
- output directory for the local plan.

## Outputs

- `workflow_selection.json`;
- `workflow_selection.md`.

## Procedure

1. Run the selector without execution:

   ```bash
   neoag-workflow-select \
     --sample-manifest configs/local/sample_manifest.yaml \
     --tools-manifest configs/local/tools_manifest.yaml \
     --reference-manifest configs/local/reference_manifest.yaml \
     --outdir work/workflow_selection
   ```

2. Review the selected entry mode and every stage status.
3. Run Doctor before any stage marked `MEDIUM` or `HIGH`.
4. Use `neoag-production-run` or the named high-level skill only after reviewing
   the generated plan. Execution remains a separate, approval-controlled step.

## Routing Rules

- existing ranked peptides: result review;
- canonical raw events and raw peptides: intermediates-to-ranking;
- somatic VCF plus HLA: VCF-to-ranking;
- somatic VCF without HLA but with a suitable normal/RNA source: HLA typing then VCF-to-ranking;
- paired tumor/normal BAM: paired-BAM production plan;
- insufficient required inputs: `BLOCKED`.

Purity/CNV and HLA LOH are planned as consensus evidence stages. Missing or
failed tools must remain `UNASSESSED`, `PARTIAL`, or tool-specific failure; they
must not be replaced with synthetic negative or successful outputs.

## Boundaries

- This skill never executes shell commands or external tools.
- Paths are redacted to basenames by default; use `--show-paths` only for a
  private local plan.
- BAM/FASTQ execution, full pipeline execution, downloads, installations, HPC
  submission, and overwrites require explicit approval.
- NetMHCpan and other licensed assets are never bundled into plans or releases.
