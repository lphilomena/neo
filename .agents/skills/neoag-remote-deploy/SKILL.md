---
name: neoag-remote-deploy
description: Use this skill when deploying or migrating Project B / NeoAg Event Pipeline to a new machine with a programming agent. It provides a fixed SOP for release checksum verification, unpacking, preflight checks, core Python installation, runtime smoke tests, local manifest generation, read-only Doctor, deployment reporting, and safe tiered external-tool readiness checks. Do not use it for interpreting existing neoantigen results only.
---

# NeoAg Remote Deploy

## Purpose

Deploy Project B on a new machine with a programming agent using a fixed,
reproducible, safety-first SOP. Make the target machine ready for Tier 0/Tier 1
result review first, then use Doctor to decide whether Tier 2/Tier 3
external-tool deployment is needed.

## When To Use

Use this skill when the user asks to migrate, deploy, install, validate a release
tarball, configure local manifests, run deployment smoke tests, run Doctor, or
prepare a deployment report for a new server, workstation, or HPC node.

Do not use this skill when the task is only to review existing ranked peptides,
compare reports, explain HLA results, or analyze patient outputs.

## Inputs

- project checkout path or release tarball;
- optional `.sha256` file for the release tarball;
- optional deployment tier: `tier0`, `tier1`, `tier2`, or `tier3`;
- optional local `tools_manifest.yaml`, `reference_manifest.yaml`, and
  `sample_manifest.yaml`;
- optional container image manifests.

## Outputs

Write deployment outputs under a work directory, usually `work/remote_deploy/`:

- `preflight_report.md` and `preflight_status.tsv`;
- local manifests under `configs/local/` by default;
- `smoke_test_report.md`;
- `doctor/doctor_summary.md` and Doctor TSV/JSON outputs;
- `deployment_report.md`;
- `audit_log.jsonl`.

## Core Procedure

Follow this sequence. Do not skip ahead to full pipeline execution.

1. Checksum and unpack if given a release tarball:
   `scripts/01_unpack_release.sh`.
2. Preflight target machine:
   `scripts/00_preflight.sh --project-root <root> --outdir <outdir>`.
3. Install core Python entry points only:
   `scripts/02_install_core.sh --project-root <root> --python <python>`.
4. Check runtime entry points:
   `scripts/03_check_runtime.sh --project-root <root> --outdir <outdir>`.
5. Generate local manifests:
   `scripts/04_configure_manifests.py --project-root <root> --outdir configs/local`.
6. Run smoke tests for the requested tier:
   `scripts/05_run_smoke_tests.sh --project-root <root> --outdir <outdir> --tier tier0`.
7. Run read-only Doctor:
   `scripts/06_run_doctor.sh --project-root <root> --manifest-dir configs/local --outdir <outdir>/doctor`.
8. Write deployment report:
   `scripts/07_write_deploy_report.py --project-root <root> --workdir <outdir>`.

Fast path from an existing checkout:

```bash
bash scripts/bootstrap_agent_deploy.sh --outdir work/agent_deploy
```

The bundled deployment scripts are more explicit and should be used when another
agent needs a step-by-step audit trail.

## Deployment Tiers

Read `references/DEPLOYMENT_TIERS.md` before choosing a tier.

- Tier 0: result review and report generation. No heavy external tools.
- Tier 1: core scoring from existing intermediate tables.
- Tier 2: containerized prediction/tool subset with Doctor mini smoke.
- Tier 3: full production/HPC from BAM/FASTQ/VCF after approval.

Default to Tier 0 or Tier 1 unless the user explicitly needs full external-tool
execution.

## Safety Rules

- Do not copy patient BAM/FASTQ/VCF into the repository.
- Do not write `/home`, `/mnt`, `/root`, HPC private paths, or patient paths into
  tracked source files.
- Do not package licensed tools, license files, or controlled references.
- Do not download large references, install system packages, submit HPC jobs,
  delete files, or overwrite results without explicit approval.
- `check-tools` or `which tool` is not enough for production readiness. Use
  Doctor mini smoke or mark the tool `PARTIAL`.
- Missing expression, HLA LOH, CCF, reference, or tool evidence must be reported
  as `MISSING`, `PARTIAL`, or `UNASSESSED`, not as biological negative evidence.

## References

Read these only as needed:

- `references/DEPLOYMENT_TIERS.md`: tier selection and acceptance criteria.
- `references/TOOLS_KNOWN_ISSUES.md`: common migration failures and fixes.
- `references/REFERENCE_MANIFEST_SCHEMA.md`: local reference manifest guidance.
- `references/TOOLS_MANIFEST_SCHEMA.md`: local tools manifest guidance.

## Failure Codes

Use stable failure codes in reports:

`CHECKSUM_FAILED`, `PYTHON_ENV_FAILED`, `CORE_INSTALL_FAILED`,
`RUN_DEMO_FAILED`, `PYTEST_FAILED`, `NEXTFLOW_PERMISSION_DENIED`,
`NEXTFLOW_FAILED`, `TOOLS_MANIFEST_MISSING`, `REFERENCE_MANIFEST_MISSING`,
`LICENSED_TOOL_MISSING`, `REFERENCE_PATH_MISSING`,
`TOOL_PATH_OK_BUT_SMOKE_FAILED`, `PRIVATE_PATH_DETECTED`,
`PATIENT_DATA_IN_RELEASE`.
