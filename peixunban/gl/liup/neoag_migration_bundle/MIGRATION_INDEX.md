# NeoAg migration package index

This directory is organized as a three-part migration kit for a fresh machine.

## 1. Core migration package

Archive:

- /mnt/zjl-bgi-zzb/peixunban/gl/liup/neoag_migration_bundle_20260714.tar.gz
- /mnt/zjl-bgi-zzb/peixunban/gl/liup/neoag_migration_bundle_20260714.tar.gz.sha256

Purpose:

- install scripts
- doctor scripts
- release tarball
- test manifests
- production-style manifests
- README and migration reports

Start here on a new machine.

## 2. Companion references and authorized tools

Directory:

- /mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata4git

Purpose:

- real GRCh38/CTAT/VEP/FACETS/ASCAT/HLA references
- restricted or license-sensitive tools and predictors
- container image tar files
- environment-variable file: neodata4git.env.sh

Do not publish this companion directory to unauthorized users. Keep it on an internal file share or provide it only to users who are allowed to use the licensed tools.

## 3. test_env_tools reference patch

Archive:

- /mnt/zjl-bgi-zzb/peixunban/gl/liup/test_env_tools_refs_neodata4git_patch_20260714.tar.gz
- /mnt/zjl-bgi-zzb/peixunban/gl/liup/test_env_tools_refs_neodata4git_patch_20260714.tar.gz.sha256

Purpose:

- production activation helper
- one-command doctor helper
- normal ligandome TSV
- Arriba reference entrypoint
- Arriba official v2.5.1 GRCh38 database

## New-machine validation

After installing the fresh conda/tool environments under test_env_tools and making neodata4git available at the expected path, run:

```bash
source /mnt/zjl-bgi-zzb/peixunban/gl/liup/test_env_tools/activate_neoag_production_refs.sh
/mnt/zjl-bgi-zzb/peixunban/gl/liup/test_env_tools/run_doctor_neodata4git.sh
```

Expected status: PARTIAL with no blocking issues. PARTIAL remains because some optional/container tools and sample input placeholders are intentionally not required for a dry-run migration doctor.

## Root bootstrap entry

For a new machine targeting `/root/neo/env_tool`, run:

```bash
sudo bash bootstrap/install_to_root_neo_env_tool.sh
```

This clones the release branches, extracts the migration/patch bundles, installs Miniforge under `/root/neo/env_tool/miniforge3`, then runs tier1/tier2 installers.
