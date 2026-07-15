# Tool migration known issues

## Nextflow permission denied

Fix:

```bash
find bin -maxdepth 1 -type f -exec chmod +x {} \;
```

## VEP path/cache incomplete

Declare `vep.executable` and `vep_cache` in local manifests. Do not copy old
server paths into tracked files.

## MHCflurry path OK but model load fails

Mark `PARTIAL`; install/fetch models in the intended environment, then rerun
Doctor mini smoke.

## NetMHCpan / NetMHCstabpan license boundary

Do not redistribute licensed binaries, data directories, or license files. The
user must stage the official install locally and configure `tools_manifest.yaml`.

## LOHHLA / Polysolver / Novoalign

`which LOHHLA` is not sufficient. Require reference/config smoke or report
`PARTIAL`.

## FACETS / ASCAT / PURPLE

If wrapper exists but reference paths are absent, report
`REFERENCE_PATH_MISSING`.

## PRIME / BigMHC / MixMHCpred

If entrypoint exists but smoke fails, report
`TOOL_PATH_OK_BUT_SMOKE_FAILED`.

## Private paths in release

If release audit finds `/home`, `/mnt`, `/root`, patient IDs, site mount points,
or license files, mark `UNSAFE` and run release cleanup before publishing.

## Real Migration Failures Fixed In 2026-07 Deployment Test

### `neoag-doctor: command not found`

Cause: the project package source was present but the console script entry point
was not installed or wrapped in the active environment.

Fix:

- install editable entry points with `python -m pip install -e .`; or
- add a project-local `bin/neoag-doctor` wrapper that runs
  `python -m neoag_v03.controlled_execution.doctor` with `PYTHONPATH=src`.

### `FileNotFoundError: vep`

Cause: VEP existed in a separate env but no `vep` command was visible to the
runtime process, or `NEOAG_VEP_BIN` was overwritten by a stale old-machine path.

Fix:

- write a target-machine `env_tool/bin/vep` wrapper;
- export `NEOAG_VEP_BIN` to that wrapper;
- validate `vep --help` after sourcing production activation.

### Missing `inputs.reference_fasta` / `NEOAG_REFERENCE_FASTA`

Cause: reference manifests pointed to old paths or activation did not export the
FASTA path required for automatic VEP annotation.

Fix:

- stage GRCh38 FASTA and `.fai` on the target machine;
- export `NEOAG_REFERENCE_FASTA` and `NEOAG_GENCODE_GTF`; or
- declare them in local manifests/run configs.

### NetMHCpan Exists But Is Not Runnable

Cause: copied wrappers may hard-code an old conda sysroot or old host paths.

Fix:

- create a target-machine wrapper under `env_tool/bin/netMHCpan`;
- set `NEOAG_NETMHCPAN_HOME` to the licensed local install;
- set temp dir to a writable location, usually `/tmp`;
- validate with `netMHCpan -h` and a small prediction smoke test.

### PRIME Appears To Run Forever And Writes Only A One-Byte Output

Cause candidates observed during migration:

- PRIME entry point recursively called itself instead of the official script;
- PRIME was running from a root-owned or non-writable licensed-tool directory;
- MixMHCpred wrapper was missing or used the wrong Python environment;
- `lib/PRIME.x` was copied from an incompatible build or used stale temp paths;
- MixMHCpred was fixed after PRIME processes had already started.

Fix:

- stop the bad run before trusting outputs;
- stage PRIME under the migration `env_tool/tools/prime` or another writable
  target-machine path;
- ensure `PRIME` is the official entry point and `lib/PRIME.x` is compiled on
  the target machine;
- use a target-machine MixMHCpred wrapper and set `MIXMHCPRED_REAL_BIN`;
- run `11_validate_production_runtime.sh --mini-prime` and require non-empty
  output before real VCF/BAM/FASTQ execution.

### BigMHC Python Dependency Drift

Cause: BigMHC and immunogenicity code run in the runtime Python env, not always
in the tool install env.

Fix: validate these imports in the env used by `neoag-v03`:

```bash
python -c 'import torch,numpy,pandas,scipy,sklearn,psutil'
```

