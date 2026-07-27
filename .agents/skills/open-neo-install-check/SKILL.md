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
- `resume`: rerun an interrupted idempotent deployment, or reuse a matching PASS checkpoint; approval remains required.

## Required input

- `project_root`, or a verified release tarball supplied by the user.

## Procedure

1. Require and verify the checksum for a release archive, reject traversal/link/device members, safely stage it under the output directory, and identify one project root.
2. Record Python, Java, Docker/Apptainer, Nextflow, disk, and platform information.
3. Generate comprehensive machine-local `tools_manifest`, `reference_manifest`, `paths.env`, and `production_assets.local.tsv`. Generic `/srv` targets are rewritten to the selected machine roots.
4. For approved `repair`/`install`, delegate to the portable `neoag-remote-deploy` new-machine installer; downloads remain opt-in and licensed assets remain external.
5. For mutating modes, run a read-only Doctor before installation and a full Doctor after installation.
6. Evaluate required tools, alternative capability groups, critical references and sidecars against the requested tier. Missing required references can never produce READY.
7. Write timeout-protected deployment checkpoints and an installation status delta for safe resume/review.
8. Write `deployment_report.md`, machine-readable status files, and an audit log.

For asset synchronization, set `--asset-source-root /mounted/assets` for a
local/mounted asset tree or `--asset-source-host user@host` for remote source
paths. Approved execution stops before installation when required sources are
not reachable.

## Outputs

- `environment_inventory.tsv`
- `doctor/doctor_status.json`
- `deployment_status.tsv`
- `tier_requirements.tsv`, `deployment_delta.tsv`
- `deployment_checkpoint.json` for mutating modes
- `deployment_report.md`
- `manifests/tools_manifest.local.yaml`
- `manifests/reference_manifest.local.yaml`
- `skill_result.json`
- `run_state.json`

## Safety boundary

Do not install or redistribute licensed tools, download large references, overwrite production settings, delete files, or submit HPC jobs unless the user explicitly approves the operation.

## Contracts and failure handling

- Validate public inputs against `references/INPUT_SCHEMA.json` before execution.
- Emit the stable result contract described by `references/OUTPUT_SCHEMA.json`.
- Use the canonical failure codes and remediation in `references/FAILURE_CODES.md`.
- Every invocation writes `skill_result.json` and a sibling `run_state.json`; install and repair actions remain approval gated.
