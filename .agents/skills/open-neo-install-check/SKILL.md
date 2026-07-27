---
name: open-neo-install-check
description: Public macro Skill1 for Open-Neo installation, new-machine migration, environment/reference/tool validation, Doctor and smoke-test readiness. Use it before running analysis on a new or changed machine.
---

# Open-Neo Install Check

## Use when

- Deploying or migrating Open-Neo to another machine.
- Determining whether the machine is ready for `review`, `core`, `prediction`, or `full` use.
- Checking tools, models, caches, references, licensed assets, release boundaries, and minimal smoke tests.

## Do not use when

- The user only wants to interpret existing results; use `open-neo-review`.
- The environment is already verified and the user wants to run a case; use `open-neo-run`.

## Modes

- `plan`: collect inventory and prepare local manifest templates without executing smoke commands.
- `verify`: run read-only Doctor and selected smoke tests.
- `repair` / `install`: require explicit human approval; the macro Skill itself never bypasses license terms.

## Required input

- `project_root`, or a verified release tarball supplied by the user.

## Procedure

1. Verify checksum when supplied.
2. Record Python, Java, Docker/Apptainer, Nextflow, disk, and platform information.
3. Generate machine-local `tools_manifest`, `reference_manifest`, and `paths.env` templates.
4. For approved `repair`/`install`, delegate to the portable `neoag-remote-deploy` new-machine installer; downloads remain opt-in and licensed assets remain external.
5. Run NeoAg Doctor and optional mini smoke tests.
6. Evaluate readiness against the requested deployment tier.
7. Write `deployment_report.md`, machine-readable status files, and an audit log.

## Outputs

- `environment_inventory.tsv`
- `doctor/doctor_status.json`
- `deployment_status.tsv`
- `deployment_report.md`
- `manifests/tools_manifest.local.yaml`
- `manifests/reference_manifest.local.yaml`
- `skill_result.json`

## Safety boundary

Do not install or redistribute licensed tools, download large references, overwrite production settings, delete files, or submit HPC jobs unless the user explicitly approves the operation.

## Contracts and failure handling

- Validate public inputs against `references/INPUT_SCHEMA.json` before execution.
- Emit the stable result contract described by `references/OUTPUT_SCHEMA.json`.
- Use the canonical failure codes and remediation in `references/FAILURE_CODES.md`.
- Every invocation writes `skill_result.json` and a sibling `run_state.json`; install and repair actions remain approval gated.
