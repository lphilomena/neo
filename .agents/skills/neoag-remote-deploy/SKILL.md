---
name: neoag-remote-deploy
description: Use this skill when deploying or migrating Project B / NeoAg Event Pipeline to a new machine with a programming agent. It provides a fixed SOP for release checksum verification, unpacking, preflight checks, core Python installation, optional env_tool/reference/licensed-tool asset migration, portable activation rewriting, runtime smoke tests, local manifest generation, read-only Doctor, deployment reporting, and safe tiered external-tool readiness checks. Do not use it for interpreting existing neoantigen results only.
---

# NeoAg Remote Deploy

## Purpose

Deploy Project B on a new machine with a programming agent using a fixed,
reproducible, safety-first SOP. Make the target machine ready for Tier 0/Tier 1
result review first, then use Doctor to decide whether Tier 2/Tier 3
external-tool deployment is needed. When full production execution is requested,
use the production-asset extension in this skill to migrate or rebuild
`env_tool`, references, and licensed-tool wrappers before running real data.

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
- optional container image manifests;
- optional source asset paths or old-machine SSH host for migration:
  `env_tool`, reference data root, and licensed-tool root;
- explicit human approval for any `--execute`, `--write`, rsync, large download,
  licensed-tool staging, or real pipeline run.

## Outputs

Write deployment outputs under a work directory, usually `work/remote_deploy/`:

- `preflight_report.md` and `preflight_status.tsv`;
- local manifests under `configs/local/` by default;
- `smoke_test_report.md`;
- `doctor/doctor_summary.md` and Doctor TSV/JSON outputs;
- `deployment_report.md`;
- `audit_log.jsonl`;
- optional `asset_migration_plan.tsv`, `asset_migration_report.md`,
  `production_runtime_status.tsv`, and rewritten local activation files when
  production assets are approved.

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

For Tier 3 or real VCF/BAM/FASTQ execution, continue with the production asset
extension before running real data:

9. Plan migration/rebuild of machine assets without copying anything:
   `scripts/08_plan_asset_migration.sh --project-root <root> --outdir <outdir> --tools-root <target-env_tool> --reference-root <target-reference-root> --licensed-root <target-licensed-root>`.
10. After approval, sync approved assets from the old machine or staged local
    paths. Default is dry-run; use `--execute` only with approval:
    `scripts/09_sync_production_assets.sh --old-host <user@host> --old-env-tool <path> --old-reference-root <path> --old-licensed-root <path> --tools-root <target-env_tool> --reference-root <target-reference-root> --licensed-root <target-licensed-root> --outdir <outdir> --execute`.
11. When README-listed external tools should be installed or rebuilt on the
    target machine, use the consolidated installer. It defaults to Miniforge3 at
    `<target-env_tool>/miniforge3`, supports open conda/git tools by group, and
    delegates licensed/local archives to `12_install_local_licensed_tools.sh`:
    `scripts/13_install_readme_tools.sh --project-root <root> --tools-root <target-env_tool> --licensed-root <target-licensed-root> --reference-root <target-reference-root> --core-env --vep --gatk --immunogenicity --allow-download --execute`.
12. When licensed tools are available as files or directories already visible
    on the target machine, install them into the target licensed-tool root without
    creating `/mnt`, `/home`, or old-machine symlinks:
    `scripts/12_install_local_licensed_tools.sh --licensed-root <target-licensed-root> --netmhcpan-tar <netMHCpan.tar.gz> --mixmhcpred-dir <MixMHCpred_install> --execute`.
13. Rewrite local activation and wrappers so the new machine uses portable
    paths, not old `/home`, `/mnt`, or stale conda prefixes:
    `scripts/10_rewrite_production_activation.sh --project-root <root> --tools-root <target-env_tool> --reference-root <target-reference-root> --licensed-root <target-licensed-root> --write`.
14. Validate the production runtime before real data:
    `scripts/11_validate_production_runtime.sh --project-root <root> --tools-root <target-env_tool> --outdir <outdir>/production_runtime --mini-prime`.

Do not run `run-full`, `pipeline-full --execute`, or any patient workflow until
step 14 shows that VEP, reference FASTA, NetMHCpan, PRIME/MixMHCpred, and the
enabled Python dependencies are available or explicitly waived.

Fast path from an existing checkout:

```bash
bash scripts/bootstrap_agent_deploy.sh --outdir work/agent_deploy
```

The bundled deployment scripts are more explicit and should be used when another
agent needs a step-by-step audit trail.

When local licensed-tool installers are missing, first retrieve an official or
user-approved download URL. Then use `12_install_local_licensed_tools.sh` with
`--allow-download`; do not download from third-party mirrors or bypass license,
registration, login, or click-through controls.

For new-machine migrations, prefer the consolidated entrypoint when the user
wants one command for asset sync, tool installation, activation rewrite,
runtime validation, and optional real VCF smoke:

```bash
bash .agents/skills/neoag-remote-deploy/scripts/16_install_new_machine.sh \
  --asset-source-host na@10.200.50.134 \
  --allow-download \
  --execute
```

For README-listed open/conda tools, prefer `13_install_readme_tools.sh` over
running many installer scripts manually. It defaults to Miniforge3 under
`/root/neo/env_tool/miniforge3` or `<tools-root>/miniforge3`; use
`--no-install-miniforge` only when a site-managed conda must be used instead.
Use `--allow-download` for Miniforge, conda packages, git clones, VEP cache, or
approved tool URLs.

`--all-open` and `--all` also enable the production-adjacent open/containerized
tooling used by copy-number and HLA workflows: SpecHLA, HLA-LA, Sequenza, and
HMF PURPLE/AMBER/COBALT. The installer expects their large databases and
container image tarballs to come from `configs/assets/production_assets.tsv`
when `--sync-assets` is used. Sequenza is installed as a conda env from
`conda/env.neoag-sequenza.yml`; splice helpers are installed from `conda/env.neoag-splice.yml`; SpecHLA, HLA-LA, and HMF PURPLE are registered
by loading staged container images when Docker is available and by writing
portable wrappers/environment variables into the production activation script.
Because `--all-open` includes BigMHC, it also installs/repairs torch by default. SHERPA-Presentation is not publicly auto-downloadable; install it only with an authorized `--sherpa-source`, `--sherpa-archive`, or `--sherpa-container-image`.
Use `--skip-torch-install` only for a deliberately lightweight install, and then
run real VCF smoke with `--skip-real-vcf-bigmhc` until torch is installed.

VEP is pinned to release 105 by default (`--vep-version 105`) and must match the bundled `homo_sapiens/105_GRCh38` cache.

The core environment installer must keep MHCflurry compatible with modern
TensorFlow/Keras by installing the matching `tf-keras` shim and exporting
`TF_USE_LEGACY_KERAS=1` in generated activation files.

NetMHCpan repair must rewrite any copied frontend that still defaults to an old
conda prefix such as `/home/na/miniforge3`; after repair, `netMHCpan -h` must
work with `CONDA_BASE=<target-env_tool>/miniforge3`.

`13_install_readme_tools.sh --run-real-vcf-smoke` runs the default
M1ML150017383 VCF smoke test after installation. The smoke test runs
MHCflurry by default, skips NetMHCstabpan by default because it is slow,
and accepts `--real-vcf-smoke-top-n <N>` for a smaller or larger test.
Use `--skip-real-vcf-mhcflurry` only as a temporary fallback on hosts with
unresolved TensorFlow/Keras compatibility issues.
Use `--skip-real-vcf-bigmhc` only as a temporary fallback on hosts where torch
or BigMHC models are not ready; this is not a full production predictor smoke.

External assets are not bundled in the skill, but the installer accepts staged
asset locations so a new machine can prepare itself reproducibly:

- Real VCF smoke inputs: `--real-vcf`, `--real-annotated-vcf`,
  `--real-vcf-hla-alleles`, and `--real-vcf-hla-file`.
- BigMHC models: `--bigmhc-models-dir <dir>` copies a local model directory, or
  combine it with `--bigmhc-models-host <user@host>` to copy from a source
  server using rsync.
- Large model/reference sets can also be synchronized from a TSV manifest with
  `--asset-manifest configs/assets/production_assets.tsv --sync-assets
  --asset-source-host <user@host>`. The manifest stores paths and markers only,
  not the asset payloads.
- SpecHLA, HLA-LA, Sequenza, and HMF PURPLE/AMBER/COBALT assets are part of the
  production asset manifest. SHERPA-Presentation is registered only from authorized local source/archive/container assets. Asset sync dereferences source symlinks (`rsync
  -L`), so stable `/mnt` links can point at real source directories while the
  target receives concrete files/directories. When a `/mnt` symlink only
  resolves on the source host, use `--asset-source-host` so rsync resolves it
  there.
- Licensed tools remain explicit inputs via `--netmhcpan-tar`,
  `--netmhcpan-dir`, `--netmhcpan-url`, `--mixmhcpred-dir`,
  `--mixmhcpred-archive`, and `--mixmhcpred-url`; do not bundle or download
  them unless the user has rights and approves the source.

Production asset fast path after explicit approval:

```bash
bash .agents/skills/neoag-remote-deploy/scripts/08_plan_asset_migration.sh \
  --project-root . \
  --tools-root /root/neo/env_tool \
  --reference-root /root/neo/neodata4git \
  --licensed-root /root/neo/licensed_tools \
  --outdir work/agent_deploy

# Review work/agent_deploy/asset_migration_report.md first, then run approved sync/rewrite.
```

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
- Do not copy or package licensed tools unless the user confirms the target use
  is licensed and the transfer is allowed. Prefer wrappers that point to a
  locally staged official installation.
- Do not migrate raw old-machine activation files verbatim. Always rewrite
  `activate_neoag_production_refs.sh`, `common.sh` defaults, VEP, NetMHCpan,
  PRIME, and MixMHCpred wrappers for the target machine.
- Do not trust long-running PRIME jobs until PRIME/MixMHCpred smoke tests have
  produced non-empty output. A running `PRIME.x` process with 100% CPU and a
  one-byte output file is not sufficient readiness evidence.
- `check-tools` or `which tool` is not enough for production readiness. Use
  Doctor mini smoke or mark the tool `PARTIAL`.
- Missing expression, HLA LOH, CCF, reference, or tool evidence must be reported
  as `MISSING`, `PARTIAL`, or `UNASSESSED`, not as biological negative evidence.

## References

Read these only as needed:

- `references/DEPLOYMENT_TIERS.md`: tier selection and acceptance criteria.
- `references/TOOLS_KNOWN_ISSUES.md`: common migration failures and fixes.
- `references/PRODUCTION_ASSET_MIGRATION.md`: env_tool, references, licensed
  tools, activation rewrites, and smoke-test acceptance criteria.
- `references/REFERENCE_MANIFEST_SCHEMA.md`: local reference manifest guidance.
- `references/TOOLS_MANIFEST_SCHEMA.md`: local tools manifest guidance.

## Failure Codes

Use stable failure codes in reports:

`CHECKSUM_FAILED`, `PYTHON_ENV_FAILED`, `CORE_INSTALL_FAILED`,
`RUN_DEMO_FAILED`, `PYTEST_FAILED`, `NEXTFLOW_PERMISSION_DENIED`,
`NEXTFLOW_FAILED`, `TOOLS_MANIFEST_MISSING`, `REFERENCE_MANIFEST_MISSING`,
`LICENSED_TOOL_MISSING`, `REFERENCE_PATH_MISSING`,
`TOOL_PATH_OK_BUT_SMOKE_FAILED`, `PRIVATE_PATH_DETECTED`,
`PATIENT_DATA_IN_RELEASE`, `ASSET_SYNC_NOT_APPROVED`,
`ASSET_SOURCE_MISSING`, `ACTIVATION_REWRITE_REQUIRED`,
`PRIME_SMOKE_FAILED`, `MIXMHCPRED_SMOKE_FAILED`, `VEP_CACHE_MISSING`,
`NETMHCPAN_LICENSE_TOOL_MISSING`.
