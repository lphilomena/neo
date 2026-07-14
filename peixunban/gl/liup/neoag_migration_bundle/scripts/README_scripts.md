# Migration install scripts

Run scripts from inside the unpacked migration bundle. Set `NEOAG_TOOLS_ROOT` to the target install root.

Recommended order:

```bash
export NEOAG_TOOLS_ROOT=/path/to/test_env_tools
bash scripts/install_tier1_core.sh
bash scripts/install_tier2_tools.sh
bash scripts/run_doctor_bundle_test.sh
```

Optional heavier layers:

```bash
bash scripts/install_gatk_vep_facets_ascat.sh
bash scripts/install_fusion_tools.sh
bash scripts/install_optional_immuno_hla_tools.sh
```

Licensed tools require user-provided packages:

```bash
NETMHCPAN_TARBALL=/path/to/netMHCpan.tar.gz \
NETMHCSTABPAN_TARBALL=/path/to/netMHCstabpan.tar.gz \
bash scripts/install_licensed_tools.sh
```

Reference staging template:

```bash
NEOAG_REF_ROOT=/path/to/neoag_refs bash scripts/stage_references_template.sh
REFERENCE_MANIFEST=refs/reference_manifest.local_template.yaml bash scripts/run_full_doctor.sh
```

Notes:

- `install_all_nonlicensed.sh` can download/install many large packages. Use it only on a machine intended for full deployment.
- VEP cache is not downloaded unless `INSTALL_VEP_CACHE=1` is set.
- NetMHCpan, NetMHCstabpan, Polysolver/Novoalign and similar tools must follow local license rules.
