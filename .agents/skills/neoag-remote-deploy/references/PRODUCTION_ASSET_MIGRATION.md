# Production Asset Migration

Use this reference only for Tier 3 or real VCF/BAM/FASTQ execution on a new
machine. The source checkout is not enough for production execution: the target
machine must also have a working tool environment, references, and licensed or
restricted companion tools.

## Assets

Recommended target layout:

```text
/root/neo/
  neoag_migration_bundle/ or src/na0707_upload_release/
  env_tool/
  neodata4git/ or refs/
  licensed_tools/
```

`env_tool` should contain the conda/miniforge base, activation scripts, wrappers,
PRIME, MixMHCpred wrappers, NetMHCpan wrappers, and other local tool shims.

Reference roots should contain at minimum:

```text
GRCh38 FASTA and .fai
GENCODE GTF
VEP cache root with homo_sapiens/105_GRCh38 or the configured version
normal proteome / ligandome when enabled
HLA, FACETS, ASCAT, LOHHLA, fusion, and other workflow-specific references
```

Licensed/restricted tools may include:

```text
NetMHCpan
NetMHCstabpan
MixMHCpred
PRIME
Novoalign / LOHHLA companion resources, when applicable
```

## Required Sequence

1. Run `08_plan_asset_migration.sh`; do not copy anything yet.
2. Review the generated plan with the user, especially licensed-tool and large
   reference entries.
3. Run `09_sync_production_assets.sh --execute` only after approval when
   assets must be copied from an old machine.
4. If README-listed open/conda tools must be rebuilt locally, run
   `13_install_readme_tools.sh --execute` with the needed tool-group flags. The
   default conda base is `/root/neo/env_tool/miniforge3` or
   `<tools-root>/miniforge3`; add `--allow-download` when Miniforge, conda
   packages, git clones, VEP cache, or approved URLs are needed.
5. If licensed tools are already available on the target machine as archives or
   install directories, run `12_install_local_licensed_tools.sh --execute` to
   install them into `/root/neo/licensed_tools` without target symlinks.
6. Run `10_rewrite_production_activation.sh --write` on the target machine.
7. Run `11_validate_production_runtime.sh --mini-prime`.
8. Run Doctor with mini smoke.
9. Only then run real data.

## README Tool Installer

Use the consolidated README tool installer when a new machine should build tools
instead of linking/copying an old `env_tool` tree:

```bash
bash .agents/skills/neoag-remote-deploy/scripts/13_install_readme_tools.sh \
  --project-root /root/neo/src/na0707_upload_release \
  --tools-root /root/neo/env_tool \
  --licensed-root /root/neo/licensed_tools \
  --reference-root /root/neo/neodata4git \
  --core-env --vep --gatk --immunogenicity --optitype \
  --allow-download \
  --execute
```

The installer defaults to Miniforge3 at `/root/neo/env_tool/miniforge3` and can
create it automatically when `--allow-download --execute` is used. Add
`--no-install-miniforge` only when a site-managed conda installation must be
used. Heavy or workflow-specific groups such as `--vep-cache`, `--fusion`,
`--facets`, and `--ascat-pyclone` should be enabled only when needed.

## Local Download Fallback

If the target machine does not already have the required installer archive or
install directory, the agent may search the web for the official or user-approved
download location for the requested version. For licensed tools such as
NetMHCpan, NetMHCstabpan, and MixMHCpred, use only vendor/lab/project URLs or a
URL explicitly provided by the user, and respect all license, registration, and
institutional access requirements.

After approval, download and install with:

```bash
bash .agents/skills/neoag-remote-deploy/scripts/12_install_local_licensed_tools.sh \
  --licensed-root /root/neo/licensed_tools \
  --netmhcpan-url <official_or_user_approved_archive_url> \
  --mixmhcpred-url <official_or_user_approved_archive_url> \
  --allow-download \
  --execute
```

The script stores downloads under `<outdir>/downloads`, extracts/copies into
`/root/neo/licensed_tools`, and intentionally avoids target symlinks to `/mnt`,
`/home`, or old-machine paths.

## Path Rewrite Rules

Do not use raw activation files copied from another server. Rewrite these for
the target machine:

- `activate_neoag_production_refs.sh`
- `scripts/common.sh` default tool root, when using a migration bundle
- VEP wrapper and `NEOAG_VEP_BIN`
- NetMHCpan wrappers and temp directory
- PRIME path and temp directory
- MixMHCpred wrapper and Python environment
- `MIXMHCPRED_REAL_BIN`, `PRIME_HOME`, `MIXMHCPRED_BIN`, `BIGMHC_DIR`

Private old paths such as `/home/na`, stale `/mnt/...` mounts, and old conda
prefixes are allowed only in local untracked configuration when they are valid on
the target machine. They must not be committed to source.

## Acceptance Checks

Production readiness requires more than `which`:

```bash
source /root/neo/env_tool/activate_neoag_production_refs.sh
neoag-doctor --help
vep --help
netMHCpan -h
/root/neo/env_tool/wrappers/mixMHCpred_install/MixMHCpred -h
/root/neo/env_tool/tools/prime/PRIME -h
python -c 'import torch,numpy,pandas,scipy,sklearn,psutil'
```

A PRIME smoke test must produce a non-empty output file. A long-running
`PRIME.x` process with nearly 100% CPU and a one-byte output file is not a pass;
it usually means the wrapper/core/temp-path state needs repair or the job should
be stopped and rerun after repair.

## Known Fixes From Migration Testing

- `neoag-doctor` may be missing as an entry point. Add a wrapper or install the
  package entry points.
- VEP may exist but not be on `PATH`; create a `vep` wrapper and set
  `NEOAG_VEP_BIN`.
- `reference_fasta` must be exported as `NEOAG_REFERENCE_FASTA` or declared in
  the run config/reference manifest.
- NetMHCpan wrappers copied from another server may hard-code a conda sysroot.
  Prefer a target-machine wrapper and a writable temp directory.
- PRIME wrappers must call the official PRIME entry point, not recursively call
  themselves.
- PRIME must use a writable local temp directory and a C++ core compiled for the
  target machine (`lib/PRIME.x`).
- MixMHCpred should use a target-machine wrapper and the intended migration
  Python environment.
- BigMHC and immunogenicity paths require Python packages in the runtime env:
  `torch`, `numpy`, `pandas`, `scipy`, `sklearn`, and `psutil`.


## Asset Manifest Sync

Large assets such as BigMHC models, EasyFuse references, CTAT libraries, and
other indexed references should stay outside Git. Use
`configs/assets/production_assets.tsv` to declare source and target paths, then
sync them during installation:

```bash
bash .agents/skills/neoag-remote-deploy/scripts/13_install_readme_tools.sh \
  --asset-manifest configs/assets/production_assets.tsv \
  --sync-assets \
  --asset-source-host na@10.200.50.134 \
  --execute
```

The installer calls `15_sync_asset_manifest.sh`. Default mode is dry-run; with
`--execute`, directories and files are copied with rsync. Required assets fail
the install if sync or marker/checksum verification fails. Optional assets are
reported but do not stop the install.
