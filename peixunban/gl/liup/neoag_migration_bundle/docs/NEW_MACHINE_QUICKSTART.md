# New machine quickstart

1. Unpack the core migration package.

```bash
tar -xzf neoag_migration_bundle_20260714.tar.gz
cd neoag_migration_bundle
```

2. Install non-licensed tools into the fresh test_env_tools root.

```bash
bash scripts/install_tier1_core.sh
bash scripts/install_tier2_tools.sh
bash scripts/install_gatk_vep_facets_ascat.sh
```

3. Make the companion directory available.

Expected path:

```text
/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata4git
```

If a different mount path is used, update:

- test_env_tools/activate_neoag_production_refs.sh
- refs/reference_manifest.neodata4git.yaml
- configs/local/tools_manifest.neodata4git.yaml

4. Apply the reference patch if it was distributed separately.

```bash
tar -xzf test_env_tools_refs_neodata4git_patch_20260714.tar.gz -C /mnt/zjl-bgi-zzb/peixunban/gl/liup
```

5. Run production-style doctor.

```bash
source /mnt/zjl-bgi-zzb/peixunban/gl/liup/test_env_tools/activate_neoag_production_refs.sh
/mnt/zjl-bgi-zzb/peixunban/gl/liup/test_env_tools/run_doctor_neodata4git.sh
```

Pass criterion for migration validation:

- blocking_issues.tsv contains only the header line.
