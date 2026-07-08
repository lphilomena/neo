#!/usr/bin/env bash
# Install Ensembl VEP via conda and optionally download cache for offline use.
#
# Usage:
#   bash scripts/install_vep.sh              # install vep into neoag-tools
#   bash scripts/install_vep.sh --cache      # also install homo_sapiens cache (~10GB+)
#
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
ENV_NAME="${NEOAG_VEP_ENV:-neoag-vep}"
INSTALL_CACHE=false
for arg in "$@"; do
  case "$arg" in
    --cache) INSTALL_CACHE=true ;;
  esac
done

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda required" >&2
  exit 1
fi

# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"

if ! conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  echo "==> Creating ${ENV_NAME} from conda/env.neoag-vep.yml ..."
  conda env create -n "${ENV_NAME}" -f "${ROOT}/conda/env.neoag-vep.yml"
fi

conda activate "${ENV_NAME}"

if ! command -v vep >/dev/null 2>&1; then
  echo "==> Installing ensembl-vep into ${ENV_NAME} ..."
  conda install -y -c bioconda -c conda-forge ensembl-vep
fi

echo "==> VEP version:"
vep --help 2>&1 | head -3 || vep -h 2>&1 | head -3

if [[ "${INSTALL_CACHE}" == "true" ]]; then
  echo "==> Installing VEP cache (homo_sapiens, can take long and use >10GB) ..."
  vep_install -a cf -s homo_sapiens -y GRCh38 -n
  echo "Cache installed. Pipeline can use: vep --cache --offline ..."
else
  echo ""
  echo "NOTE: neoag upstream uses --cache --offline by default."
  echo "Run cache install when ready:"
  echo "  conda activate ${ENV_NAME}"
  echo "  vep_install -a cf -s homo_sapiens -y"
  echo ""
  echo "Or use online VEP (set NEOAG_VEP_ONLINE=1 in run config / see docs/TOOLS_SETUP.md)."
fi

PREFIX="$(conda env list | awk -v n="${ENV_NAME}" '$1==n {print $NF}')"
TOOLS_ENV="${ROOT}/conf/tools.env.sh"
mkdir -p "${ROOT}/conf"
if [[ ! -f "${TOOLS_ENV}" ]]; then
  cat > "${TOOLS_ENV}" <<EOF
export NEOAG_PROJECT_ROOT="${ROOT}"
export NEOAG_TOOLS_ROOT="${ROOT}"
export NEOAG_CONDA_BASE="$(conda info --base)"
export NEOAG_CONDA_ENV="neoag-tools"
EOF
fi
if ! grep -q 'VEP — installed via scripts/install_vep.sh' "${TOOLS_ENV}"; then
  cat >> "${TOOLS_ENV}" <<EOF

# VEP — installed via scripts/install_vep.sh
export NEOAG_VEP_ENV="${ENV_NAME}"
export NEOAG_VEP_BIN="${PREFIX}/bin/vep"
if [[ -d "${PREFIX}/bin" ]]; then
  export PATH="${PREFIX}/bin:\${PATH}"
fi
EOF
else
  echo "==> conf/tools.env.sh already contains a VEP install block; check NEOAG_VEP_BIN if needed."
fi

echo "==> Done. Test: conda activate ${ENV_NAME} && vep --help | head"
