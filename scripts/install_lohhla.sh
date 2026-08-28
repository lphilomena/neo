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
#   neoag check-tools | grep lohhla
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TOOLS_ROOT="${NEOAG_TOOLS_ROOT:-${ROOT}}"
TARGET="${LOHHLA_HOME:-${TOOLS_ROOT}/tools/lohhla}"
TOOLS_ENV="${NEOAG_TOOLS_ENV:-${ROOT}/conf/tools.env.local.sh}"
BIN_DIR="${ROOT}/bin"
CONDA_BASE="${NEOAG_CONDA_BASE:-${TOOLS_ROOT}/miniforge3}"
RSCRIPT_BIN="${NEOAG_LOHHLA_RSCRIPT:-}"
if [[ -z "$RSCRIPT_BIN" && -x "${CONDA_BASE}/envs/neoag-facets/bin/Rscript" ]]; then
  RSCRIPT_BIN="${CONDA_BASE}/envs/neoag-facets/bin/Rscript"
fi
if [[ -z "$RSCRIPT_BIN" ]]; then
  RSCRIPT_BIN="$(command -v Rscript || true)"
fi
# The historical Bitbucket clone URL no longer supports anonymous git access.
# Use the public mirror by default while retaining an override for institutions.
REPO="${LOHHLA_GIT_URL:-https://github.com/slagtermaarten/LOHHLA.git}"
REF="${NEOAG_LOHHLA_REF:-b38c4770995b24628a4e038fccb1a9cd57c4f305}"
GITHUB_PROXY_PREFIX="${NEOAG_GITHUB_PROXY_PREFIX:-https://ghproxy.net/}"
mkdir -p "$(dirname "${TARGET}")" "${BIN_DIR}"

curl_supports_retry_all_errors() {
  command -v curl >/dev/null 2>&1 || return 1
  { curl --help all 2>/dev/null || curl --help 2>/dev/null; } | grep -q -- '--retry-all-errors'
}

download_file() {
  local url="$1" destination="$2"
  local -a curl_args=(-fL --retry 5 --connect-timeout 30)
  if curl_supports_retry_all_errors; then
    curl_args+=(--retry-all-errors)
  else
    echo "INFO: curl lacks --retry-all-errors; using portable retry options." >&2
  fi
  curl "${curl_args[@]}" -o "$destination" "$url"
}

if [[ ! -f "${TARGET}/LOHHLAscript.R" ]]; then
  tmp="$(mktemp -d)"
  archive="$tmp/source.tar.gz"
  direct="https://github.com/slagtermaarten/LOHHLA/archive/${REF}.tar.gz"
  if [[ "$REPO" == "https://github.com/slagtermaarten/LOHHLA.git" ]]; then
    for url in "${GITHUB_PROXY_PREFIX}${direct}" "$direct"; do
      download_file "$url" "$archive" && break
      rm -f "$archive"
    done
  fi
  if [[ -s "$archive" ]]; then
    mkdir -p "$TARGET"
    tar -xzf "$archive" -C "$tmp"
    source_dir="$(find "$tmp" -mindepth 1 -maxdepth 1 -type d | head -1)"
    cp -a "$source_dir/." "$TARGET/"
  else
    rm -rf "$TARGET" 2>/dev/null || true
    git clone --depth 1 "$REPO" "$TARGET"
  fi
  rm -rf "$tmp"
fi

cat > "${BIN_DIR}/LOHHLA" <<EOF
#!/usr/bin/env bash
set -euo pipefail
if [[ "\${1:-}" == "--version" || "\${1:-}" == "-v" || "\${1:-}" == "-h" || "\${1:-}" == "--help" ]]; then
  echo "LOHHLA wrapper for ${TARGET}/LOHHLAscript.R"
  exit 0
fi
exec "${RSCRIPT_BIN}" "${TARGET}/LOHHLAscript.R" "\$@"
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
export NEOAG_LOHHLA_RSCRIPT="${RSCRIPT_BIN}"
export PATH="${BIN_DIR}:${TARGET}:\${PATH}"
# Set these in conf/tools.env.local.sh for real runs:
export POLYSOLVER_HOME="\${POLYSOLVER_HOME:-}"
export NOVOALIGN_LICENSE_FILE="\${NOVOALIGN_LICENSE_FILE:-}"
EOF
fi

echo "==> LOHHLA source installed at ${TARGET}. Production runs still need Polysolver/Novoalign/HLA references."
echo "==> Run: source conf/tools.env.sh && neoag check-tools | grep lohhla"
