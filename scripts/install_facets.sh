#!/usr/bin/env bash
# Install FACETS into a conda/R environment and expose a runFACETS.R wrapper for check-tools.
#
# Usage:
#   bash scripts/install_facets.sh
#   source conf/tools.env.sh
#   neoag-v03 check-tools | grep facets
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CONDA_BASE="${NEOAG_CONDA_BASE:-$(conda info --base)}"
ENV_NAME="${NEOAG_FACETS_ENV:-neoag-facets}"
TOOLS_ENV="${ROOT}/conf/tools.env.sh"
BIN_DIR="${ROOT}/bin"
mkdir -p "${BIN_DIR}"

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda not found" >&2
  exit 1
fi
# shellcheck disable=SC1091
source "${CONDA_BASE}/etc/profile.d/conda.sh"

env_exists() { conda env list | awk '{print $1}' | grep -qx "$1"; }
env_has_facets() {
  conda run -n "$1" Rscript -e 'quit(status=ifelse(requireNamespace("facets", quietly=TRUE),0,1))' >/dev/null 2>&1
}

if [[ -x "${BIN_DIR}/runFACETS.R" ]] && "${BIN_DIR}/runFACETS.R" --version >/dev/null 2>&1; then
  echo "==> FACETS wrapper already present: ${BIN_DIR}/runFACETS.R"
  "${BIN_DIR}/runFACETS.R" --version || true
  echo "==> Done. Run: source conf/tools.env.sh && neoag-v03 check-tools | grep facets"
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
    conda install -n "${ENV_NAME}" -c conda-forge -c bioconda -y r-base r-optparse r-devtools r-remotes htslib || true
  else
    conda create -n "${ENV_NAME}" -c conda-forge -c bioconda -y r-base r-optparse r-devtools r-remotes htslib
  fi
  FACETS_ENV="${ENV_NAME}"
  if ! env_has_facets "${FACETS_ENV}"; then
    conda install -n "${FACETS_ENV}" -c conda-forge -c bioconda -y r-facets || \
    conda run -n "${FACETS_ENV}" Rscript -e 'if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes", repos="https://cloud.r-project.org"); remotes::install_github("mskcc/facets")'
  fi
fi
ENV_NAME="${FACETS_ENV}"

cat > "${BIN_DIR}/runFACETS.R" <<EOF
#!/usr/bin/env bash
set -euo pipefail
if [[ "\${1:-}" == "--version" || "\${1:-}" == "-v" ]]; then
  conda run -n "${ENV_NAME}" Rscript -e 'cat(as.character(utils::packageVersion("facets")), "\\n")'
  exit 0
fi
if [[ "\$#" -eq 0 ]]; then
  echo "FACETS wrapper. Use scripts/facets_fit_from_pileup.R or run: conda activate ${ENV_NAME}; Rscript your_facets_script.R" >&2
  exit 0
fi
conda run -n "${ENV_NAME}" Rscript "\$@"
EOF
chmod +x "${BIN_DIR}/runFACETS.R"

if [[ ! -f "${TOOLS_ENV}" ]]; then
  cat > "${TOOLS_ENV}" <<EOF
export NEOAG_PROJECT_ROOT="${ROOT}"
export NEOAG_TOOLS_ROOT="${ROOT}"
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
echo "==> Done. Run: source conf/tools.env.sh && neoag-v03 check-tools | grep facets"
