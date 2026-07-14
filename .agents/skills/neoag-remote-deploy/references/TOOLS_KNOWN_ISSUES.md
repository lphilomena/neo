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
