# Tools manifest guidance

Recommended local path: `configs/local/tools_manifest.yaml`.

Example:

```yaml
tools:
  vep:
    executable: /opt/conda/envs/neoag-vep/bin/vep
    mode: conda_or_container
    required: false
  netmhcpan:
    executable: /opt/netMHCpan/bin/netMHCpan
    mode: local_license
    license_required: true
    required: false
  mhcflurry:
    executable: mhcflurry-predict
    mode: conda
    required: false
  lohhla:
    executable: LOHHLA
    mode: local_or_container
    required: false
  facets:
    executable: runFACETS.R
    mode: conda_or_container
    required: false
```

Rules:

- Prefer containers for fragile external tools.
- Mark licensed tools with `license_required: true`.
- Keep local tool paths in untracked local manifests, not source code.
- If a tool exists but mini smoke fails, report `PARTIAL`, not `READY`.
