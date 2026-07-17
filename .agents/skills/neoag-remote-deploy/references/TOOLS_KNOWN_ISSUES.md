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
Observed error:

```text
netMHCpan: conda sysroot loader missing at /home/na/miniforge3/envs/neoag-tools/x86_64-conda-linux-gnu/sysroot/lib/ld-linux-x86-64.so.2
```

Fix:

- create a target-machine wrapper under `env_tool/bin/netMHCpan`;
- set `NEOAG_NETMHCPAN_HOME` to the licensed local install;
- rewrite copied licensed frontend defaults so `CONDA_BASE` points to the target
  machine, usually `/root/neo/env_tool/miniforge3`;
- set temp dir to a writable location, usually `/tmp`;
- validate with `netMHCpan -h` and a small prediction smoke test.

`scripts/install_netmhcpan.sh --repair` must rewrite tcsh launchers and any
frontend containing stale `/home/na/miniforge3` or other old-machine conda
defaults.

### MHCflurry Fails With Keras Cannot Be Imported

Observed error:

```text
ImportError: Keras cannot be imported. Check that it is installed.
```

Cause: MHCflurry 2.x runs against modern TensorFlow/Keras and needs the matching
legacy `tf-keras` shim with `TF_USE_LEGACY_KERAS=1`.

Fix:

```bash
source /root/neo/env_tool/miniforge3/etc/profile.d/conda.sh
conda activate neoag-tools
TF_KERAS_SPEC="$(python - <<'PY'
import tensorflow as tf
major, minor, *_ = tf.__version__.split(".")
print(f"tf-keras>={major}.{minor},<{major}.{int(minor) + 1}")
PY
)"
pip install "$TF_KERAS_SPEC"
```

The consolidated installer now performs this repair automatically after the core
environment exists.

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

Observed error during real VCF smoke:

```text
ModuleNotFoundError: No module named 'torch'
```

Fix:

- full `--all-open` / `--all` installs must install torch by default because
  BigMHC is included by default;
- prefer CPU torch on new machines:
  `pip install --index-url https://download.pytorch.org/whl/cpu torch`;
- if an approved local wheel cache is provided, pass `--torch-wheel-dir <dir>`;
- if a CUDA torch wheel is used, install all matching `nvidia-*-cu12`,
  `nvidia-nvjitlink-cu12`, `triton`, `filelock`, and `sympy==1.13.1`;
- use `--skip-real-vcf-bigmhc` only as a temporary smoke-test fallback, not as
  production readiness.

### Asset Symlinks Resolve On Source But Not Target

Observed with VEP plugins and SpecHLA DB: `/mnt/.../neodata4git` entries were
absolute symlinks to `/home/na/...`. They resolved on the source host but not on
the new machine.

Fix:

- `15_sync_asset_manifest.sh` uses `rsync -aL` so source symlinks are
  dereferenced and the target receives real files/directories;
- use `--asset-source-host` when the symlink target only exists on the source
  host;
- use real marker files, such as `ref/hla.ref.extend.fa` for SpecHLA DB, rather
  than a bare directory marker.

### VEP Plugin Directory Present But Plugin Files Missing

Observed error:

```text
VEP Wildtype plugin missing
VEP Frameshift plugin missing
```

Cause: the manifest path existed but resolved to an empty or broken symlink on
the target.

Fix: copy real `Wildtype.pm` and `Frameshift.pm` into
`<reference-root>/work/vep_plugins`, and keep `env_tool/work/vep_plugins` linked
there.
