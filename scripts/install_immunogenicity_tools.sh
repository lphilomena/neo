#!/usr/bin/env bash
# Install PRIME + MixMHCpred + BigMHC_IM for neoag immunogenicity evidence.
#
# Usage:
#   bash scripts/install_immunogenicity_tools.sh
#   source conf/tools.env.sh
#   neoag check-tools | grep -E 'prime|bigmhc'
#
# Notes from deployment tests:
# - PRIME must compile lib/PRIME.x, not a separate PRIME.x.bin, because the PRIME wrapper calls PRIME.x.
# - MixMHCpred and BigMHC require Python packages (numpy, pandas, psutil, torch).
# - Network cloning of BigMHC can fail on slow links; rerun the script or pre-stage tools/bigmhc.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TOOLS="${NEOAG_TOOLS_ROOT:-${ROOT}}/tools"
BIN_DIR="${ROOT}/bin"
PRIME_DIR="${TOOLS}/prime"
MIX_DIR="${TOOLS}/mixMHCpred_install"
BIGMHC_DIR="${TOOLS}/bigmhc"
TOOLS_ENV="${ROOT}/conf/tools.env.sh"
PYTHON_BIN="${NEOAG_IMMUNO_PYTHON:-python3}"
GITHUB_PROXY_PREFIX="${NEOAG_GITHUB_PROXY_PREFIX:-https://ghproxy.net/}"
PRIME_REF="${NEOAG_PRIME_REF:-7b18d4e11042141e7102f7c69be2b0e03d138dab}"
MIX_REF="${NEOAG_MIXMHCPRED_REF:-c29e4db17abe6266bfee72750efb713459540d18}"
BIGMHC_REF="${NEOAG_BIGMHC_REF:-c7e37a249317704bf96a1e3881a7ece3c3c977a6}"

mkdir -p "${TOOLS}" "${BIN_DIR}"

install_github_snapshot() {
  local repo="$1" ref="$2" target="$3"
  local archive_url="https://github.com/${repo}/archive/${ref}.tar.gz"
  local tmp archive extracted
  tmp="$(mktemp -d)"
  archive="$tmp/source.tar.gz"
  for url in "${GITHUB_PROXY_PREFIX}${archive_url}" "$archive_url"; do
    echo "Downloading pinned ${repo}@${ref} from ${url}"
    curl -fL --retry 5 --retry-all-errors --connect-timeout 30 -o "$archive" "$url" && break
    rm -f "$archive"
  done
  [[ -s "$archive" ]] || { rm -rf "$tmp"; echo "ERROR: failed to download ${repo}@${ref}" >&2; return 1; }
  tar -xzf "$archive" -C "$tmp"
  extracted="$(find "$tmp" -mindepth 1 -maxdepth 1 -type d | head -1)"
  [[ -n "$extracted" ]] || { rm -rf "$tmp"; echo "ERROR: invalid snapshot for ${repo}" >&2; return 1; }
  mkdir -p "$target"
  cp -a "$extracted/." "$target/"
  rm -rf "$tmp"
}

echo "[0/4] Python dependencies for MixMHCpred/BigMHC"
"${PYTHON_BIN}" -m pip install numpy pandas psutil || true
if [[ "${NEOAG_SKIP_TORCH_INSTALL:-0}" != "1" ]]; then
  if "${PYTHON_BIN}" -c 'import torch' >/dev/null 2>&1; then
    echo "torch already available in ${PYTHON_BIN}; skip reinstall"
  else
    "${PYTHON_BIN}" -m pip install torch || echo "WARN: torch install failed; BigMHC smoke test may fail. Install torch manually or set NEOAG_SKIP_TORCH_INSTALL=1 if already available." >&2
  fi
fi

echo "[1/4] PRIME"
if [[ ! -f "${PRIME_DIR}/lib/run_PRIME.pl" ]]; then
  install_github_snapshot GfellerLab/PRIME "$PRIME_REF" "$PRIME_DIR"
fi
if [[ ! -x "${PRIME_DIR}/PRIME" ]]; then
  curl -fsSL https://raw.githubusercontent.com/GfellerLab/PRIME/master/PRIME -o "${PRIME_DIR}/PRIME"
fi
chmod +x "${PRIME_DIR}/PRIME"
PRIME_TEMP_DIR="${PRIME_DIR}/temp"
mkdir -p "${PRIME_TEMP_DIR}"
if [[ ! -d "${PRIME_TEMP_DIR}" || ! -w "${PRIME_TEMP_DIR}" ]]; then
  echo "ERROR: PRIME runtime temp directory is not writable: ${PRIME_TEMP_DIR}" >&2
  echo "Run: NEOAG_PRIME_RUNTIME_USER=<runtime-user> bash scripts/fix_prime_temp.sh" >&2
  exit 1
fi
echo "PRIME runtime temp: ${PRIME_TEMP_DIR} (writable)"
if [[ -f "${PRIME_DIR}/lib/PRIME.cc" ]]; then
  echo "[1b/4] Compile PRIME lib/PRIME.x"
  (cd "${PRIME_DIR}/lib" && g++ -O3 PRIME.cc -o PRIME.x)
  chmod +x "${PRIME_DIR}/lib/PRIME.x"
fi

echo "[2/4] MixMHCpred"
if [[ ! -x "${MIX_DIR}/MixMHCpred" ]]; then
  install_github_snapshot GfellerLab/MixMHCpred "$MIX_REF" "$MIX_DIR"
fi
chmod +x "${MIX_DIR}/MixMHCpred" 2>/dev/null || true

# MixMHCpred invokes `python3` internally. Keep it on the same Python runtime
# where its dependencies were installed, including when PRIME calls it.
PYTHON_DIR="$(dirname "$(readlink -f "${PYTHON_BIN}")")"
cat > "${BIN_DIR}/MixMHCpred" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export PATH="${PYTHON_DIR}:\${PATH}"
exec "${MIX_DIR}/MixMHCpred" "\$@"
EOF
chmod +x "${BIN_DIR}/MixMHCpred"

echo "[3/4] BigMHC"
if [[ ! -f "${BIGMHC_DIR}/src/predict.py" ]]; then
  # Merge source into an existing asset directory so pre-staged model weights survive.
  install_github_snapshot KarchinLab/bigmhc "$BIGMHC_REF" "$BIGMHC_DIR"
fi

cat > "${BIN_DIR}/bigmhc_predict" <<EOF
#!/usr/bin/env bash
set -euo pipefail
cd "${BIGMHC_DIR}/src"
exec "${PYTHON_BIN}" predict.py "\$@"
EOF
chmod +x "${BIN_DIR}/bigmhc_predict"

mkdir -p "${ROOT}/conf"
if [[ ! -f "${TOOLS_ENV}" ]]; then
  cat > "${TOOLS_ENV}" <<EOF
export NEOAG_PROJECT_ROOT="${ROOT}"
export NEOAG_TOOLS_ROOT="${ROOT}"
export NEOAG_CONDA_ENV="neoag-tools"
EOF
fi
if ! grep -q 'PRIME / MixMHCpred / BigMHC — installed via scripts/install_immunogenicity_tools.sh' "${TOOLS_ENV}"; then
  cat >> "${TOOLS_ENV}" <<EOF

# PRIME / MixMHCpred / BigMHC — installed via scripts/install_immunogenicity_tools.sh
export PRIME_HOME="${PRIME_DIR}"
export MIXMHCPRED_HOME="${MIX_DIR}"
export BIGMHC_DIR="${BIGMHC_DIR}"
export NEOAG_PRIME_BIN="${PRIME_DIR}/PRIME"
export MIXMHCPRED_BIN="${BIN_DIR}/MixMHCpred"
export PATH="${PRIME_DIR}:${MIX_DIR}:${BIN_DIR}:\${PATH}"
EOF
fi

# Export for this shell too.
export PRIME_HOME="${PRIME_DIR}"
export MIXMHCPRED_HOME="${MIX_DIR}"
export BIGMHC_DIR="${BIGMHC_DIR}"
export NEOAG_PRIME_BIN="${PRIME_DIR}/PRIME"
export MIXMHCPRED_BIN="${BIN_DIR}/MixMHCpred"
export PATH="${PRIME_DIR}:${MIX_DIR}:${BIN_DIR}:${PATH}"

echo "[4/4] Smoke tests"
if [[ -f "${PRIME_DIR}/test/test.txt" && -x "${NEOAG_PRIME_BIN}" && -x "${MIXMHCPRED_BIN}" ]]; then
  "${NEOAG_PRIME_BIN}" -i "${PRIME_DIR}/test/test.txt" \
    -o /tmp/prime_smoke.tsv \
    -a A0101,A2501,B0801,B1801 \
    -mix "${MIXMHCPRED_BIN}" >/tmp/prime_smoke.log 2>&1 || {
      cat /tmp/prime_smoke.log >&2
      echo "WARN: PRIME smoke failed; inspect PRIME/MixMHCpred dependencies." >&2
    }
  head -3 /tmp/prime_smoke.tsv 2>/dev/null || true
fi
if [[ -f "${BIGMHC_DIR}/data/example1.csv" ]]; then
  (cd "${BIGMHC_DIR}/src" && "${PYTHON_BIN}" predict.py -i=../data/example1.csv -m=im -d=cpu -a=0 -p=1 -c=1 >/tmp/bigmhc_smoke.log 2>&1) || {
    cat /tmp/bigmhc_smoke.log >&2
    echo "WARN: BigMHC smoke failed; ensure torch/pandas/psutil are installed in ${PYTHON_BIN}." >&2
  }
fi

echo "OK if check-tools reports prime and bigmhc_im as OK. Run: source conf/tools.env.sh && neoag check-tools | grep -E 'prime|bigmhc'"
