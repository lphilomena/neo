# Tool Setup for NeoAg v0.4.3

This document records the external-tool installation procedure validated during the 169-machine migration test and the fixes made after that test.

The lightweight package does **not** bundle large third-party tools, references, VEP cache, NetMHCpan license package, patient data, or conda environments. Install only the tools required by your workflow.

## 1. Base system packages

Ubuntu/Debian:

```bash
sudo apt-get update
sudo apt-get install -y \
  bash coreutils curl wget git tar gzip unzip bzip2 xz-utils \
  ca-certificates build-essential openjdk-17-jre-headless rsync file
```

## 2. Basic package check

```bash
find bin -maxdepth 1 -type f -exec chmod +x {} \;
find scripts -maxdepth 1 -type f -name '*.sh' -exec chmod +x {} \;
python -m pip install -e '.[test]'
pytest -q
neoag-v03 run-demo --outdir work/demo_v043 --sample-id DEMO001
```

Migration-test reference result: `175 passed, 95 skipped` for the light test suite. Skips are expected for external tools and benchmark tests.

## 3. Environment entry point

Create or update the primary tool env:

```bash
bash scripts/setup_tools_env.sh
source conf/tools.env.sh
neoag-v03 check-tools
```

If `mhcflurry-downloads fetch` fails with a `CXXABI` or `libstdc++` error:

```bash
conda activate neoag-tools
conda install -n neoag-tools -c conda-forge -y 'libstdcxx-ng>=13'
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
mhcflurry-downloads fetch
```

`check-tools` confirms executables, not necessarily model/reference completeness. For MHCflurry, also smoke-test the model command.

## 4. Tool-specific installation

### 4.1 GATK

```bash
bash scripts/install_gatk.sh
source conf/tools.env.sh
neoag-v03 check-tools
```

### 4.2 VEP

```bash
bash scripts/install_vep.sh
source conf/tools.env.sh
neoag-v03 check-tools
```

VEP cache may be slow to download. Either run:

```bash
bash scripts/install_vep_cache.sh
```

or set an existing cache:

```bash
export NEOAG_VEP_CACHE=/path/to/vep_cache
export NEOAG_VEP_CACHE_VERSION=105
source conf/tools.env.sh
```

The migration test found that VEP could be installed but not detected until the VEP env was added to `PATH`. The updated script writes `NEOAG_VEP_BIN` and the env `bin` path into `conf/tools.env.sh`.

### 4.3 NetMHCpan

NetMHCpan requires a DTU academic license tarball.

```bash
mkdir -p vendor
cp /path/to/netMHCpan-4.2c.Linux.tar.gz vendor/
export NEOAG_CONDA_BASE="$(conda info --base)"
bash scripts/install_netmhcpan.sh vendor/netMHCpan-4.2c.Linux.tar.gz
source conf/tools.env.sh
neoag-v03 check-tools
netMHCpan -h | head
```

If an existing installation needs only wrapper repair:

```bash
bash scripts/install_netmhcpan.sh --repair
```

The previous hardcoded `/home/na/miniforge3` path was removed. Use `NEOAG_CONDA_BASE` to override conda location.

### 4.4 NetMHCstabpan

```bash
bash scripts/install_netmhcstabpan.sh --iedb
source conf/tools.env.sh
neoag-v03 check-tools
```

### 4.5 DeepImmuno

```bash
bash scripts/install_deepimmuno.sh
source conf/tools.env.sh
neoag-v03 check-tools
```

### 4.6 PRIME / MixMHCpred / BigMHC

```bash
bash scripts/install_immunogenicity_tools.sh
source conf/tools.env.sh
neoag-v03 check-tools
```

Fixes incorporated after migration testing:

- sets `NEOAG_PRIME_BIN`, `MIXMHCPRED_BIN`, and `BIGMHC_DIR`;
- compiles `tools/prime/lib/PRIME.x` rather than `PRIME.x.bin`;
- installs Python packages required by MixMHCpred/BigMHC: `numpy`, `pandas`, `psutil`, and `torch`;
- creates `bin/bigmhc_predict` wrapper.

If BigMHC cloning fails due to network interruption, rerun the script or pre-stage `tools/bigmhc` manually.

### 4.7 FACETS

```bash
bash scripts/install_facets.sh
source conf/tools.env.sh
neoag-v03 check-tools
```

This installs the R package/wrapper. Real FACETS analysis still needs pileup input, common-SNP VCF, and sample-specific fit/export settings.

### 4.8 ASCAT / PyClone-VI

```bash
bash scripts/install_ascat_pyclone.sh
source conf/tools.env.sh
neoag-v03 check-tools
```

Do **not** run `pip install mamba` to satisfy this script. That installs a Python test framework named `mamba`, not the conda-forge solver. The updated script uses `conda` by default.

### 4.9 LOHHLA

```bash
bash scripts/install_lohhla.sh
source conf/tools.env.sh
neoag-v03 check-tools
```

For real LOHHLA runs, also configure in `conf/tools.env.local.sh`:

```bash
export POLYSOLVER_HOME=/path/to/polysolver
export NOVOALIGN_LICENSE_FILE=/path/to/novoalign.lic
```

LOHHLA also requires patient HLA calls, HLA FASTA/resources, tumor/normal BAM, and purity/ploidy information.

### 4.10 Fusion tools

The release includes wrappers for EasyFuse/STAR-Fusion/FusionCatcher when their environments and references are already prepared. `scripts/seed_easyfuse_conda_envs.sh` is **not** a full fresh installer; it seeds from an existing Nextflow conda cache. If `work/.nextflow_conda` does not exist, first prepare the fusion environments/references or use site-provided tool bundles.

Arriba is optional and may be installed separately via conda or upstream binaries. Put `arriba` on `PATH` before checking.

## 5. Expected `check-tools` interpretation

A minimal scoring run may not require every tool. Example statuses:

- Required for common presentation: `netmhcpan`, `mhcflurry`.
- Recommended for SNV/InDel annotation: `vep`, optionally `gatk` for upstream calling.
- Recommended for CCF/APPM/escape: `facets`, `pyclone`, `lohhla` depending on available data.
- Optional immunogenicity: `prime`, `bigmhc_im`, `deepimmuno`.
- Optional fusion: `easyfuse`, `star_fusion`, `fusioncatcher`, `arriba`.

If a tool is `MISSING`, the corresponding evidence layer should be marked missing/unassessed, not negative.

## 6. Production strict mode

For production-like runs:

```bash
export NEOAG_STRICT_MODE=1
```

Strict mode forbids stub tool outputs. Use demo/stub mode only for smoke testing and software validation.

## 7. Common failure table

| Failure | Fix |
|---|---|
| `Permission denied` on `bin/neoag-nextflow` | `find bin -maxdepth 1 -type f -exec chmod +x {} \;` |
| MHCflurry `CXXABI_1.3.15` error | install `libstdcxx-ng>=13`, export conda lib in `LD_LIBRARY_PATH` |
| VEP installed but missing | run updated `scripts/install_vep.sh`, then `source conf/tools.env.sh` |
| NetMHCpan uses wrong conda path | set `NEOAG_CONDA_BASE=$(conda info --base)`, rerun install or `--repair` |
| `mamba: unrecognized arguments` | do not use pip mamba; updated script uses conda |
| PRIME syntax error in `PRIME.x` | updated script compiles `lib/PRIME.x` |
| BigMHC lacks torch/pandas/psutil | updated script installs required Python packages |
| LOHHLA missing | run `scripts/install_lohhla.sh` and configure Polysolver/Novoalign |
| FACETS missing | run `scripts/install_facets.sh` |
