#!/usr/bin/env bash
# Install Ensembl VEP via conda and optionally download cache for offline use.
#
# Usage:
#   bash scripts/install_vep.sh              # install VEP into neoag-vep
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

CONDA_BASE="${NEOAG_CONDA_BASE:-$(command conda info --base)}"
TOOLS_ROOT="${NEOAG_TOOLS_ROOT:-$(dirname "$CONDA_BASE")}"
CONDA_PKGS_DIR="${NEOAG_CONDA_PKGS_DIR:-$TOOLS_ROOT/conda_pkgs}"
mkdir -p "$CONDA_PKGS_DIR"
conda config --remove-key pkgs_dirs >/dev/null 2>&1 || true
conda config --add pkgs_dirs "$CONDA_PKGS_DIR" >/dev/null 2>&1 || true

# shellcheck disable=SC1091
source "$CONDA_BASE/etc/profile.d/conda.sh"

conda_safe() {
  set +u
  "$CONDA_BASE/bin/conda" "$@"
  local rc=$?
  set -u
  return "$rc"
}

resolve_env_prefix() {
  conda_safe env list | awk -v n="${ENV_NAME}" '$1==n {print $NF; exit}'
}

ENV_PREFIX="$(resolve_env_prefix)"
if [[ -z "${ENV_PREFIX}" || ! -x "${ENV_PREFIX}/bin/vep" ]]; then
  echo "==> Creating ${ENV_NAME} from conda/env.neoag-vep.yml ..."
  conda_safe create -n "${ENV_NAME}" --override-channels -c conda-forge -c bioconda -y "ensembl-vep=${VEP_VERSION}.*"
fi

ENV_PREFIX="$(resolve_env_prefix)"
if [[ -z "${ENV_PREFIX}" || ! -x "${ENV_PREFIX}/bin/perl" ]]; then
  echo "ERROR: could not resolve Perl for conda env ${ENV_NAME}" >&2
  exit 1
fi

run_in_vep_env() (
  unset PERL5LIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT
  export PATH="${ENV_PREFIX}/bin:${PATH}"
  "$@"
)

CURRENT_VEP_VERSION="$(conda_safe list -p "${ENV_PREFIX}" ensembl-vep 2>/dev/null | awk '$1=="ensembl-vep" {print $2; exit}')"
if [[ ! -x "${ENV_PREFIX}/bin/vep" || "${CURRENT_VEP_VERSION}" != ${VEP_VERSION}* ]]; then
  echo "==> Installing ensembl-vep ${VEP_VERSION}.* into ${ENV_NAME} ..."
  conda_safe install -p "${ENV_PREFIX}" --override-channels -c conda-forge -c bioconda -y "ensembl-vep=${VEP_VERSION}.*"
fi

if ! run_in_vep_env "${ENV_PREFIX}/bin/perl" -MDBI -e 'exit 0'; then
  echo "==> Installing missing Perl DBI into ${ENV_NAME} ..."
  conda_safe install -p "${ENV_PREFIX}" --override-channels -c conda-forge -c bioconda -y perl-dbi
fi
run_in_vep_env "${ENV_PREFIX}/bin/perl" -MDBI -e 'print "DBI $DBI::VERSION\n"'

echo "==> VEP version:"
VEP_HELP="$(mktemp)"
if ! run_in_vep_env "${ENV_PREFIX}/bin/vep" --help >"${VEP_HELP}" 2>&1; then
  cat "${VEP_HELP}" >&2
  rm -f "${VEP_HELP}"
  echo "ERROR: VEP smoke test failed with ${ENV_PREFIX}/bin/perl" >&2
  exit 1
fi
head -3 "${VEP_HELP}"
rm -f "${VEP_HELP}"

if [[ "${INSTALL_CACHE}" == "true" ]]; then
  echo "==> Installing VEP cache (homo_sapiens, can take long and use >10GB) ..."
  run_in_vep_env "${ENV_PREFIX}/bin/vep_install" -a cf -s homo_sapiens -y GRCh38 -n
  echo "Cache installed. Pipeline can use: vep --cache --offline ..."
else
  echo ""
  echo "NOTE: neoag upstream uses --cache --offline by default."
  echo "Run cache install when ready:"
  echo "  NEOAG_VEP_ENV=${ENV_NAME} NEOAG_VEP_VERSION=${VEP_VERSION} bash scripts/install_vep.sh --cache"
  echo ""
  echo "Or use online VEP (set NEOAG_VEP_ONLINE=1 in run config / see docs/TOOLS_SETUP.md)."
fi

PREFIX="${ENV_PREFIX}"
WRAPPER_DIR="${TOOLS_ROOT}/tools/bin"
mkdir -p "${WRAPPER_DIR}"
cat > "${WRAPPER_DIR}/vep" <<EOF
#!/usr/bin/env bash
unset PERL5LIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT
export PATH="${PREFIX}/bin:\${PATH}"
exec "${PREFIX}/bin/vep" "\$@"
EOF
chmod +x "${WRAPPER_DIR}/vep"
VEP_BIN="${WRAPPER_DIR}/vep"
TOOLS_ENV="${NEOAG_TOOLS_ENV:-${ROOT}/conf/tools.env.local.sh}"
mkdir -p "${ROOT}/conf"
if [[ ! -f "${TOOLS_ENV}" ]]; then
  cat > "${TOOLS_ENV}" <<EOF
export NEOAG_PROJECT_ROOT="${ROOT}"
export NEOAG_TOOLS_ROOT="${TOOLS_ROOT}"
export NEOAG_CONDA_BASE="${CONDA_BASE}"
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

echo "==> Done. Test: ${VEP_BIN} --help"
