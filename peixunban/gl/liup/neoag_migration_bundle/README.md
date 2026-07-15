# NeoAg migration bundle

This bundle is prepared for migrating `neo-na0707_upload_release-3` to a fresh Linux machine.

## Contents

- `release/`: source release tarball and SHA256 checksum.
- `env/`: conda environment YAML and activation example.
- `configs/local/`: local tool/sample manifest templates.
- `refs/`: reference manifest template, tiny test reference set, and reference staging notes.
- `licensed/`: licensed tool staging notes.
- `scripts/`: tiered installation wrappers around the project install skill/scripts.
- `reports/`: smoke test and Doctor outputs from validation.
- `logs/`: selected installation logs for troubleshooting.

## Minimal fresh migration

```bash
tar -xzf neoag_migration_bundle_20260714.tar.gz
cd neoag_migration_bundle
(cd release && sha256sum -c neo-na0707_upload_release-3.sha256)
export NEOAG_TOOLS_ROOT=/path/to/test_env_tools
bash scripts/install_tier1_core.sh
bash scripts/install_tier2_tools.sh
bash scripts/run_doctor_bundle_test.sh
```

Expected result: Doctor exits successfully with no blocking issues. The status may be `PARTIAL` because optional external tools and licensed tools are intentionally not bundled.

## Optional installation layers

```bash
bash scripts/install_gatk_vep_facets_ascat.sh
bash scripts/install_fusion_tools.sh
bash scripts/install_optional_immuno_hla_tools.sh
```

Licensed tools:

```bash
NETMHCPAN_TARBALL=/path/to/netMHCpan.tar.gz \
NETMHCSTABPAN_TARBALL=/path/to/netMHCstabpan.tar.gz \
bash scripts/install_licensed_tools.sh
```

Reference staging:

```bash
NEOAG_REF_ROOT=/path/to/neoag_refs bash scripts/stage_references_template.sh
REFERENCE_MANIFEST=refs/reference_manifest.local_template.yaml bash scripts/run_full_doctor.sh
```

This bundle does not include patient data, large production reference data, or licensed NetMHCpan/NetMHCstabpan packages.

## Root bootstrap

Use `bootstrap/install_to_root_neo_env_tool.sh` as the first-step installer on a new machine when `/root/neo/env_tool` is the target environment/tool root.
