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
VEP_VERSION="${NEOAG_VEP_VERSION:-105}"
INSTALL_CACHE=false
for arg in "$@"; do
  case "$arg" in
    --cache) INSTALL_CACHE=true ;;
  esac
done

if [[ -n "${NEOAG_CONDA_BASE:-}" ]]; then
  export PATH="$NEOAG_CONDA_BASE/bin:$PATH"
fi

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda required" >&2
  exit 1
fi

CONDA_BASE="${NEOAG_CONDA_BASE:-$(conda info --base)}"
TOOLS_ROOT="${NEOAG_TOOLS_ROOT:-$(dirname "$CONDA_BASE")}"
CONDA_PKGS_DIR="${NEOAG_CONDA_PKGS_DIR:-$TOOLS_ROOT/conda_pkgs}"
mkdir -p "$CONDA_PKGS_DIR"
conda config --remove-key pkgs_dirs >/dev/null 2>&1 || true
conda config --add pkgs_dirs "$CONDA_PKGS_DIR" >/dev/null 2>&1 || true

# shellcheck disable=SC1091
source "$CONDA_BASE/etc/profile.d/conda.sh"

if [[ ! -x "$CONDA_BASE/envs/${ENV_NAME}/bin/vep" ]]; then
  echo "==> Creating ${ENV_NAME} from conda/env.neoag-vep.yml ..."
  conda create -n "${ENV_NAME}" --override-channels -c conda-forge -c bioconda -y "ensembl-vep=${VEP_VERSION}.*"
fi

conda activate "${ENV_NAME}"

CURRENT_VEP_VERSION="$(conda list -n "${ENV_NAME}" ensembl-vep 2>/dev/null | awk '$1=="ensembl-vep" {print $2; exit}')"
if ! command -v vep >/dev/null 2>&1 || [[ "${CURRENT_VEP_VERSION}" != ${VEP_VERSION}* ]]; then
  echo "==> Installing ensembl-vep ${VEP_VERSION}.* into ${ENV_NAME} ..."
  conda install --override-channels -c conda-forge -c bioconda -y "ensembl-vep=${VEP_VERSION}.*"
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
  echo "  export NEOAG_VEP_VERSION=${VEP_VERSION}"
  echo "  conda activate ${ENV_NAME}"
  echo "  vep_install -a cf -s homo_sapiens -y"
  echo ""
  echo "Or use online VEP (set NEOAG_VEP_ONLINE=1 in run config / see docs/TOOLS_SETUP.md)."
fi

PREFIX="${CONDA_BASE}/envs/${ENV_NAME}"
WRAPPER_DIR="${TOOLS_ROOT}/tools/bin"
mkdir -p "${WRAPPER_DIR}"
cat > "${WRAPPER_DIR}/vep" <<EOF
#!/usr/bin/env bash
exec "${CONDA_BASE}/bin/conda" run -n "${ENV_NAME}" vep "\$@"
EOF
chmod +x "${WRAPPER_DIR}/vep"
VEP_BIN="${WRAPPER_DIR}/vep"
TOOLS_ENV="${ROOT}/conf/tools.env.sh"
mkdir -p "${ROOT}/conf"
if [[ ! -f "${TOOLS_ENV}" ]]; then
  cat > "${TOOLS_ENV}" <<EOF
export NEOAG_PROJECT_ROOT="${ROOT}"
export NEOAG_TOOLS_ROOT="${TOOLS_ROOT}"
export NEOAG_CONDA_BASE="${NEOAG_CONDA_BASE:-$(conda info --base)}"
export NEOAG_CONDA_ENV="neoag-tools"
EOF
fi
if ! grep -q 'VEP — installed via scripts/install_vep.sh' "${TOOLS_ENV}"; then
  cat >> "${TOOLS_ENV}" <<EOF

# VEP — installed via scripts/install_vep.sh
export PATH="${WRAPPER_DIR}:\${PATH}"
export NEOAG_VEP_ENV="${ENV_NAME}"
export NEOAG_VEP_VERSION="${VEP_VERSION}"
export NEOAG_VEP_BIN="${VEP_BIN}"
EOF
else
  echo "==> conf/tools.env.sh already contains a VEP install block; check NEOAG_VEP_BIN if needed."
fi

echo "==> Done. Test: conda activate ${ENV_NAME} && vep --help | head"
