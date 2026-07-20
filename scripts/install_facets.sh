#!/usr/bin/env bash
# Install FACETS into a conda/R environment and expose a runFACETS.R wrapper for check-tools.
#
# Usage:
#   bash scripts/install_facets.sh
#   source conf/tools.env.sh
#   neoag check-tools | grep facets
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CONDA_BASE="${NEOAG_CONDA_BASE:-$(command conda info --base)}"
export PATH="${CONDA_BASE}/bin:${PATH}"
ENV_NAME="${NEOAG_FACETS_ENV:-neoag-facets}"
TOOLS_ENV="${ROOT}/conf/tools.env.sh"
BIN_DIR="${ROOT}/bin"
mkdir -p "${BIN_DIR}"

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda not found" >&2
  exit 1
fi
conda_safe() {
  set +u
  conda "$@"
  local rc=$?
  set -u
  return "$rc"
}

# shellcheck disable=SC1091
set +u
source "${CONDA_BASE}/etc/profile.d/conda.sh"
set -u

env_exists() { conda_safe env list | awk '{print $1}' | grep -qx "$1"; }
env_has_facets() {
  conda_safe run -n "$1" Rscript -e 'quit(status=ifelse(requireNamespace("facets", quietly=TRUE),0,1))' >/dev/null 2>&1
}

if [[ -x "${BIN_DIR}/runFACETS.R" ]] && "${BIN_DIR}/runFACETS.R" --version >/dev/null 2>&1; then
  echo "==> FACETS wrapper already present: ${BIN_DIR}/runFACETS.R"
  "${BIN_DIR}/runFACETS.R" --version || true
  echo "==> Done. Run: source conf/tools.env.sh && neoag check-tools | grep facets"
  exit 0
fi

FACETS_ENV="${ENV_NAME}"
for candidate in "${ENV_NAME}" neoag-fusion-r36 neoag-fusion; do
  if env_exists "${candidate}" && env_has_facets "${candidate}"; then
    FACETS_ENV="${candidate}"
    echo "==> Reusing existing facets R package in conda env: ${FACETS_ENV}"
    break
  fi
done

if ! env_has_facets "${FACETS_ENV}"; then
  if env_exists "${ENV_NAME}"; then
    conda_safe install -n "${ENV_NAME}" -c conda-forge -c bioconda -y r-base r-optparse r-devtools r-remotes htslib || true
  else
    conda_safe create -n "${ENV_NAME}" -c conda-forge -c bioconda -y r-base r-optparse r-devtools r-remotes htslib
  fi
  FACETS_ENV="${ENV_NAME}"
  if ! env_has_facets "${FACETS_ENV}"; then
    if ! conda_safe install -n "${FACETS_ENV}" -c conda-forge -c bioconda -y r-facets; then
      SRC_ROOT="${NEOAG_TOOLS_SRC_ROOT:-$(dirname "${CONDA_BASE}")/tools/src}"
      mkdir -p "${SRC_ROOT}"
      if [[ ! -d "${SRC_ROOT}/pctGCdata/.git" ]]; then
        git clone --depth 1 https://github.com/mskcc/pctGCdata.git "${SRC_ROOT}/pctGCdata"
      else
        git -C "${SRC_ROOT}/pctGCdata" pull --ff-only
      fi
      if [[ ! -d "${SRC_ROOT}/facets/.git" ]]; then
        git clone --depth 1 https://github.com/mskcc/facets.git "${SRC_ROOT}/facets"
      else
        git -C "${SRC_ROOT}/facets" pull --ff-only
      fi
      "${CONDA_BASE}/bin/conda" run -n "${FACETS_ENV}" R CMD INSTALL "${SRC_ROOT}/pctGCdata"
      "${CONDA_BASE}/bin/conda" run -n "${FACETS_ENV}" R CMD INSTALL "${SRC_ROOT}/facets"
    fi
  fi
fi
ENV_NAME="${FACETS_ENV}"

cat > "${BIN_DIR}/runFACETS.R" <<EOF
#!/usr/bin/env bash
set -euo pipefail
if [[ "\${1:-}" == "--version" || "\${1:-}" == "-v" ]]; then
  ${CONDA_BASE}/bin/conda run -n "${ENV_NAME}" Rscript -e 'cat(as.character(utils::packageVersion("facets")), "\\n")'
  exit 0
fi
if [[ "\$#" -eq 0 ]]; then
  echo "FACETS wrapper. Use scripts/facets_fit_from_pileup.R or run: conda activate ${ENV_NAME}; Rscript your_facets_script.R" >&2
  exit 0
fi
${CONDA_BASE}/bin/conda run -n "${ENV_NAME}" Rscript "\$@"
EOF
chmod +x "${BIN_DIR}/runFACETS.R"

if [[ ! -f "${TOOLS_ENV}" ]]; then
  cat > "${TOOLS_ENV}" <<EOF
export NEOAG_PROJECT_ROOT="${ROOT}"
export NEOAG_TOOLS_ROOT="$(dirname "${CONDA_BASE}")"
export NEOAG_CONDA_BASE="${CONDA_BASE}"
export NEOAG_CONDA_ENV="neoag-tools"
EOF
fi
if ! grep -q 'FACETS — installed via scripts/install_facets.sh' "${TOOLS_ENV}"; then
  cat >> "${TOOLS_ENV}" <<EOF

# FACETS — installed via scripts/install_facets.sh
export NEOAG_FACETS_ENV="${ENV_NAME}"
export FACETS_HOME="${BIN_DIR}"
export PATH="${BIN_DIR}:\${PATH}"
EOF
fi

"${BIN_DIR}/runFACETS.R" --version || true
echo "==> Done. Run: source conf/tools.env.sh && neoag check-tools | grep facets"
