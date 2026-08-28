# Workflow selection and QC orchestration

`neoag-workflow-select` is the read-only entry point for choosing between
result review, canonical-intermediate ranking, somatic-VCF processing, and
paired-BAM production. It reads manifests and writes a plan; it does not launch
Nextflow, external tools, containers, or HPC jobs.

## Why the deploy-dev QC modules were not copied directly

The experimental FACETS, PURPLE, LOHHLA, and SpecHLA modules suppressed tool
failures and emitted placeholder TSV files. A successful workflow status could
therefore mean that no biological result was produced. The release branch keeps
the existing production scripts and high-level consensus skills instead:

- `neoag-purity-cnv-run-and-review` for FACETS, PURPLE, Sequenza, and ASCAT;
- `neoag-hla-typing-run-and-compare` for OptiType, SpecHLA, and HLA-LA;
- `neoag-hla-loh-appm-review` for LOHHLA/SpecHLA and APPM interpretation;
- `neoag-production-run` for dependency ordering and declared-output checks.

Every external-tool stage must preserve its own exit status and declared output
status. Missing evidence is `UNASSESSED`; a tool failure is a failure or partial
run, not a negative biological result.

## Example

```bash
neoag-workflow-select \
  --sample-manifest configs/local/sample_manifest.yaml \
  --tools-manifest configs/local/tools_manifest.yaml \
  --reference-manifest configs/local/reference_manifest.yaml \
  --outdir work/workflow_selection
```

Review `workflow_selection.md`, run Doctor, then use the recommended production
runner or skill through the Gateway when execution approval is required.
