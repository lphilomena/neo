#!/usr/bin/env bash
# Install splice-derived neoantigen tools in isolated environments.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONDA_BASE="${NEOAG_CONDA_BASE:-$(conda info --base 2>/dev/null)}"
ENV_ROOT="${NEOAG_ENV_ROOT:-${CONDA_BASE}/envs}"
TOOLS_ROOT="${NEOAG_SPLICE_TOOLS_ROOT:-${NEOAG_TOOLS_ROOT:-${ROOT}}/tools/splice}"
ARCHIVE_ROOT="${NEOAG_CONTAINER_ARCHIVE_ROOT:-${ROOT}/work/container_images}"
SNAF_VERSION="0.7.0"
DO_CORE=0 DO_SNAF=0 DO_ASNEO=0 DO_NEOSPLICE=0 DO_SPLICE2NEO=0 DO_SPLICEMUTR=0

usage() {
  cat <<'EOF'
Usage: bash scripts/install_splice_neoantigen_tools.sh [targets]

Targets:
  --core         RegTools + samtools + bcftools; verify pVACsplice
  --snaf         SNAF in a pinned Python 3.7 environment
  --asneo        ASNEO source and Python dependencies (GRCh37 only)
  --neosplice    Stage NeoSplice 0.0.3 source (matched tumor/normal RNA required)
  --splice2neo   Pull splice2neo v0.6.14 and export its image archive
  --splicemutr   Build the supplied SpliceMutr and LeafCutter environments
  --all          Run every target

Environment:
  NEOAG_ENV_ROOT=/path/to/envs
  NEOAG_SPLICE_TOOLS_ROOT=/path/to/tool_sources
  NEOAG_CONTAINER_ARCHIVE_ROOT=/path/to/container_archives
  NEOAG_SNAF_SOURCE=/path/to/pinned/SNAF/source
  NEOAG_SPLICEMUTR_SOURCE=/path/to/SpliceMutr/source
  NEOAG_PIP_INDEX_URL=https://your.site/pypi/simple/  Optional site mirror
EOF
}

if [[ $# -eq 0 ]]; then DO_CORE=1; fi
while [[ $# -gt 0 ]]; do
  case "$1" in
    --core) DO_CORE=1 ;;
    --snaf) DO_SNAF=1 ;;
    --asneo) DO_ASNEO=1 ;;
    --neosplice) DO_NEOSPLICE=1 ;;
    --splice2neo) DO_SPLICE2NEO=1 ;;
    --splicemutr) DO_SPLICEMUTR=1 ;;
    --all) DO_CORE=1; DO_SNAF=1; DO_ASNEO=1; DO_NEOSPLICE=1; DO_SPLICE2NEO=1; DO_SPLICEMUTR=1 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown target: $1" >&2; usage >&2; exit 2 ;;
  esac
  shift
done

mkdir -p "${ENV_ROOT}" "${TOOLS_ROOT}" "${ARCHIVE_ROOT}"
MAMBA="${CONDA_BASE}/bin/mamba"
[[ -x "${MAMBA}" ]] || MAMBA="${CONDA_BASE}/bin/conda"

download_valid_tar() {
  local url="$1" output="$2" attempt tmp
  if tar -tzf "${output}" >/dev/null 2>&1; then return 0; fi
  tmp="${output}.partial"
  for attempt in 1 2 3 4 5 6 7 8; do
    echo "==> Download ${output} (attempt ${attempt}/8)"
    if command -v wget >/dev/null 2>&1; then
      wget --continue --tries=8 --timeout=30 -O "${tmp}" "${url}" || true
    else
      curl -fL --retry 8 --retry-delay 5 --connect-timeout 20 --max-time 7200 \
        -C - -o "${tmp}" "${url}" || true
    fi
    if tar -tzf "${tmp}" >/dev/null 2>&1; then
      mv "${tmp}" "${output}"
      return 0
    fi
  done
  echo "ERROR: unable to download a valid archive: ${url}" >&2
  return 1
}

if [[ "${DO_CORE}" == 1 ]]; then
  env="${ENV_ROOT}/neoag-splice"
  if [[ ! -x "${env}/bin/regtools" ]]; then
    "${MAMBA}" create -y -p "${env}" -c conda-forge -c bioconda regtools=1.0.0 samtools bcftools
  fi
  "${env}/bin/regtools" --version
  pvac=""
  for candidate in \
    "${NEOAG_PVACSPLICE_BIN:-}" \
    "$(command -v pvacsplice 2>/dev/null || true)" \
    "${CONDA_BASE}/envs/neoag-tools/bin/pvacsplice" \
    "${HOME}/miniforge3/envs/neoag-tools/bin/pvacsplice" \
    "${HOME}/miniconda3/envs/neoag-tools/bin/pvacsplice"
  do
    if [[ -n "${candidate}" && -x "${candidate}" ]]; then
      pvac="${candidate}"
      break
    fi
  done
  [[ -n "${pvac}" ]] || { echo "ERROR: pVACsplice missing; set NEOAG_PVACSPLICE_BIN" >&2; exit 1; }
  "${pvac}" run -h >/dev/null
fi

if [[ "${DO_SNAF}" == 1 ]]; then
  env="${ENV_ROOT}/neoag-snaf"
  pip_args=()
  [[ -n "${NEOAG_PIP_INDEX_URL:-}" ]] && pip_args+=(--index-url "${NEOAG_PIP_INDEX_URL}")
  [[ -x "${env}/bin/python" ]] || "${MAMBA}" create -y -p "${env}" -c conda-forge python=3.7 pip
  snaf_source="${NEOAG_SNAF_SOURCE:-${TOOLS_ROOT}/SNAF-e23ce39512a1a7f58c74e59b4b7cedc89248b908}"
  if [[ -f "${snaf_source}/setup.py" ]]; then
    "${env}/bin/pip" install "${pip_args[@]}" "${snaf_source}"
  else
    "${env}/bin/pip" install "${pip_args[@]}" "snaf==${SNAF_VERSION}"
  fi
  "${env}/bin/python" -c 'import snaf; print("SNAF", getattr(snaf, "__version__", "installed"))'
fi

if [[ "${DO_SPLICEMUTR}" == 1 ]]; then
  src="${NEOAG_SPLICEMUTR_SOURCE:-${TOOLS_ROOT}/SpliceMutr-main}"
  main_yml="${src}/envs/splicemutr_packages.yml"
  leaf_yml="${src}/envs/leafcutter_package.yml"
  [[ -f "${main_yml}" && -f "${leaf_yml}" ]] || {
    echo "ERROR: SpliceMutr source or environment files missing: ${src}" >&2
    exit 1
  }
  [[ -x "${ENV_ROOT}/neoag-splicemutr/bin/Rscript" ]] || \
    "${MAMBA}" env create -y -p "${ENV_ROOT}/neoag-splicemutr" -f "${main_yml}"
  [[ -x "${ENV_ROOT}/neoag-splicemutr-leafcutter/bin/Rscript" ]] || \
    "${MAMBA}" env create -y -p "${ENV_ROOT}/neoag-splicemutr-leafcutter" -f "${leaf_yml}"
  "${ENV_ROOT}/neoag-splicemutr/bin/python" -c 'import numpy, pandas; print("SpliceMutr Python dependencies ready")'
  "${ENV_ROOT}/neoag-splicemutr/bin/Rscript" -e 'suppressPackageStartupMessages(library(GenomicFeatures)); cat("SpliceMutr R dependencies ready\n")'
fi

if [[ "${DO_ASNEO}" == 1 ]]; then
  env="${ENV_ROOT}/neoag-asneo"
  archive="${TOOLS_ROOT}/ASNEO-master.tar.gz"
  src="${TOOLS_ROOT}/ASNEO"
  [[ -x "${env}/bin/python" ]] || "${MAMBA}" create -y -p "${env}" -c conda-forge -c bioconda python=3.9 pip bedtools samtools
  "${env}/bin/pip" install pandas numpy scipy pysam biopython scikit-learn 'xgboost==1.7.6'
  "${env}/bin/pip" install --no-build-isolation sj2psi
  download_valid_tar "https://codeload.github.com/bm2-lab/ASNEO/tar.gz/refs/heads/master" "${archive}"
  if [[ ! -f "${src}/ASNEO.py" ]]; then mkdir -p "${src}"; tar -xzf "${archive}" --strip-components=1 -C "${src}"; fi
  echo "ASNEO staged at ${src}; upstream supports GRCh37/hg19, not GRCh38."
fi

if [[ "${DO_NEOSPLICE}" == 1 ]]; then
  archive="${TOOLS_ROOT}/NeoSplice-0.0.3.tar.gz"
  src="${TOOLS_ROOT}/NeoSplice-0.0.3"
  download_valid_tar "https://codeload.github.com/pirl-unc/NeoSplice/tar.gz/refs/tags/0.0.3" "${archive}"
  if [[ ! -f "${src}/README.md" ]]; then mkdir -p "${src}"; tar -xzf "${archive}" --strip-components=1 -C "${src}"; fi
  echo "NeoSplice source staged at ${src}; production runs require matched tumor/normal RNA inputs."
fi

if [[ "${DO_SPLICE2NEO}" == 1 ]]; then
  command -v docker >/dev/null 2>&1 || { echo "ERROR: docker is required for splice2neo" >&2; exit 1; }
  image="ghcr.io/tron-bioinformatics/splice2neo:v0.6.14"
  docker pull "${image}"
  docker save -o "${ARCHIVE_ROOT}/splice2neo_v0.6.14.tar" "${image}"
  sha256sum "${ARCHIVE_ROOT}/splice2neo_v0.6.14.tar" > "${ARCHIVE_ROOT}/splice2neo_v0.6.14.tar.sha256"
fi

echo "==> Splice neoantigen tool installation completed"
