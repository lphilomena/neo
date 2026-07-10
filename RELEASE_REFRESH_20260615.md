# NeoAg Event Pipeline v0.3-rc — Server Refresh 2026-06-15

This file records the server refresh because `RELEASE.md` is owned by root in this checkout and cannot be edited by user `na`.

## Verification

- Command: `.venv/bin/python -m pytest -q`
- Result: `115 passed`
- Pytest collection is limited to project tests via `pyproject.toml`; bundled third-party tool tests under `tools/` are excluded.
- `tests/conftest.py` forces the current checkout `src/` path ahead of a stale root-owned editable `.venv` path.

## Fixes Applied

- Fixed NetMHCpan classic allele formatting for `-a` style input: `HLA-A*02:01` -> `HLA-A0201`.
- Kept NetMHCpan PEPTIDEMHC allele formatting separate: `HLA-A02:01`.
- Disabled pVACseq enrich in generated SNV WES stub configs, because fixture VCFs may intentionally lack CSQ annotations.
- Added README guidance to load `conf/tools.env.sh` before Nextflow runs so Java is visible.

## Environment Check

After `source conf/tools.env.sh`:

- Java: `${NEOAG_CONDA_BASE}/envs/neoag-fusion/bin/java`
- Nextflow: `bin/nextflow -version` reports `26.04.3 build 12259`

## Remaining Note

The directory `.venv/` is owned by root, so user `na` cannot rewrite the editable install pointer directly. The project-level test bootstrap and `bin/neoag-v03` wrapper both force the current checkout source path. For a fully clean environment, recreate `.venv/` as user `na` or change ownership before running `pip install -e ".[test]"`.

## Follow-up Optimization Applied

- Added `bin/neoag-nextflow`, which sources `conf/tools.env.sh` before launching Nextflow.
- Added VCF preflight checks for pVACseq enrich readiness: file exists, VEP CSQ header, and sample GT format.
- Added `scripts/write_release_manifest.py` to emit `release_manifest.json` with pytest, Java, Nextflow, tool-check, git status, and release file hashes.
- Replaced `datetime.utcnow()` calls with timezone-aware UTC timestamps.

## V04 Evidence/Safety/Escape Patch Merge - 2026-06-16

- Merged `/tmp/neoag_v04_evidence_safety_escape_patch_only_20260616.tar` into the server checkout.
- Added v04 evidence, peptide safety, immune escape, SV/WES phase 1.5, docs, fixtures, and tests from the patch bundle.
- Repaired schema compatibility after the patch overwrite so existing provenance, APPM, presentation, CCF, HLA LOH, and immune escape outputs retain their expected columns.
- Restored the VEP compatibility package required by the CLI import path.
- Extended ranked peptide output to carry `escape_flag`, `resistance_risk`, and related immune-escape fields through scoring.
- Verification: `.venv/bin/python -m pytest -q` -> `80 passed`.
- Release manifest refreshed: `release_manifest.json`.

