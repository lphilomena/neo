#!/usr/bin/env bash
# Install LOHHLA source tree and expose a LOHHLA launcher for check-tools.
# LOHHLA itself requires external dependencies for production runs: patient HLA calls,
# HLA FASTA, tumor purity/ploidy, and typically Novoalign/Polysolver resources.
#
# Usage:
#   bash scripts/install_lohhla.sh
#   export POLYSOLVER_HOME=/path/to/polysolver
#   export NOVOALIGN_LICENSE_FILE=/path/to/novoalign.lic
#   source conf/tools.env.sh
#   neoag-v03 check-tools | grep lohhla
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TOOLS_ROOT="${NEOAG_TOOLS_ROOT:-${ROOT}}"
TARGET="${LOHHLA_HOME:-${TOOLS_ROOT}/tools/lohhla}"
TOOLS_ENV="${ROOT}/conf/tools.env.sh"
BIN_DIR="${ROOT}/bin"
REPO="${LOHHLA_GIT_URL:-https://bitbucket.org/mcferrine/lohhla.git}"
mkdir -p "$(dirname "${TARGET}")" "${BIN_DIR}"

if [[ ! -f "${TARGET}/LOHHLAscript.R" ]]; then
  git clone "${REPO}" "${TARGET}"
fi

cat > "${BIN_DIR}/LOHHLA" <<EOF
#!/usr/bin/env bash
set -euo pipefail
if [[ "\${1:-}" == "--version" || "\${1:-}" == "-v" || "\${1:-}" == "-h" || "\${1:-}" == "--help" ]]; then
  echo "LOHHLA wrapper for ${TARGET}/LOHHLAscript.R"
  exit 0
fi
exec Rscript "${TARGET}/LOHHLAscript.R" "\$@"
EOF
chmod +x "${BIN_DIR}/LOHHLA"

if [[ ! -f "${TOOLS_ENV}" ]]; then
  cat > "${TOOLS_ENV}" <<EOF
export NEOAG_PROJECT_ROOT="${ROOT}"
export NEOAG_TOOLS_ROOT="${TOOLS_ROOT}"
export NEOAG_CONDA_ENV="neoag-tools"
EOF
fi
if ! grep -q 'LOHHLA — installed via scripts/install_lohhla.sh' "${TOOLS_ENV}"; then
  cat >> "${TOOLS_ENV}" <<EOF

# LOHHLA — installed via scripts/install_lohhla.sh
export LOHHLA_HOME="${TARGET}"
export PATH="${BIN_DIR}:${TARGET}:\${PATH}"
# Set these in conf/tools.env.local.sh for real runs:
export POLYSOLVER_HOME="\${POLYSOLVER_HOME:-}"
export NOVOALIGN_LICENSE_FILE="\${NOVOALIGN_LICENSE_FILE:-}"
EOF
fi

echo "==> LOHHLA source installed at ${TARGET}. Production runs still need Polysolver/Novoalign/HLA references."
echo "==> Run: source conf/tools.env.sh && neoag-v03 check-tools | grep lohhla"
