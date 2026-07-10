# Installation README update — 2026-07-01

Updated according to migration tests on the 169 machine.

## Main changes

- Rewrote `README.md` / `README_zh.md` with validated installation sequence, tested demo commands, Nextflow permission fix, and tool-specific troubleshooting.
- Rewrote `docs/TOOLS_SETUP.md` with per-tool installation commands and known failure fixes.
- Fixed `scripts/setup_tools_env.sh` to prefer conda C++ runtime for MHCflurry model download.
- Added `libstdcxx-ng`, `sysroot_linux-64`, and `patchelf` to `conda/env.neoag-tools.yml`; added `libstdcxx-ng` to lite env.
- Fixed VEP install script to export `NEOAG_VEP_BIN` and add the VEP env to PATH.
- Fixed NetMHCpan installer to use `NEOAG_CONDA_BASE` / `conda info --base` instead of hardcoded `${NEOAG_CONDA_BASE}`, and to fall back to direct binary execution if conda loader is unavailable.
- Rewrote ASCAT/PyClone installer to use conda by default and create `bin/ascat.R` / `bin/pyclone` wrappers.
- Rewrote PRIME/MixMHCpred/BigMHC installer to set required variables, compile `PRIME.x`, install Python dependencies, and create `bin/bigmhc_predict`.
- Added `scripts/install_facets.sh` and `scripts/install_lohhla.sh`.
- Updated `conf/tools.env.sh` to resolve conda base dynamically and detect VEP env.

## Validation

- `bash -n` passed for updated shell scripts.
- `python -m compileall -q src` passed.
- `PYTHONPATH=src pytest -q` passed: `177 passed, 95 skipped`.
- `PYTHONPATH=src python -m neoag_v03.cli run-demo --outdir ... --sample-id READMEDEMO` completed successfully.
