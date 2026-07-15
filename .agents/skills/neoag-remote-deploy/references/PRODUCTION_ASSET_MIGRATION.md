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
3. Run `09_sync_production_assets.sh --execute` only after approval.
4. Run `10_rewrite_production_activation.sh --write` on the target machine.
5. Run `11_validate_production_runtime.sh --mini-prime`.
6. Run Doctor with mini smoke.
7. Only then run real data.

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
