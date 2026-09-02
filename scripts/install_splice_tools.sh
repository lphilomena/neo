#!/usr/bin/env bash
# Install/register splice-junction neoantigen helper tools.
#
# Installs:
#   - RegTools in a dedicated neoag-splice conda env
#   - pVACsplice wrapper from the existing neoag-tools/pVACtools env
#   - SNAF and SpliceMutr from pinned approved source snapshots (default on)
# Optional:
#   - ASNEO / NeoSplice / splice2neo source directories as registered wrappers
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CONDA_BASE="${NEOAG_CONDA_BASE:-$(command conda info --base)}"
TOOLS_ROOT="${NEOAG_TOOLS_ROOT:-${ROOT}}"
PORTABLE_ENVS="${NEOAG_CONDA_ENVS_PATH:-${CONDA_ENVS_PATH:-}}"
[[ -z "$PORTABLE_ENVS" ]] || export CONDA_ENVS_PATH="$PORTABLE_ENVS"
ENV_NAME="${NEOAG_SPLICE_ENV:-neoag-splice}"
SNAF_ENV_NAME="${NEOAG_SNAF_ENV:-neoag-snaf}"
TOOLS_ENV="${NEOAG_TOOLS_ENV:-${ROOT}/conf/tools.env.local.sh}"
BIN_DIR="${ROOT}/bin"
YML="${ROOT}/conda/env.neoag-splice.yml"
export CONDA_CHANNEL_ALIAS="${NEOAG_CONDA_CHANNEL_ALIAS:-${CONDA_CHANNEL_ALIAS:-https://conda.anaconda.org}}"

INSTALL_SNAF="${NEOAG_INSTALL_SNAF:-1}"
INSTALL_SPLICEMUTR="${NEOAG_INSTALL_SPLICEMUTR:-1}"
SNAF_SOURCE="${SNAF_SOURCE:-${NEOAG_SNAF_SOURCE:-}}"
SNAF_GIT_URL="${SNAF_GIT_URL:-${NEOAG_SNAF_GIT_URL:-https://github.com/frankligy/SNAF.git}}"
SNAF_GIT_REF="${SNAF_GIT_REF:-${NEOAG_SNAF_GIT_REF:-e23ce39512a1a7f58c74e59b4b7cedc89248b908}}"
SNAF_PACKAGE_URL="${SNAF_PACKAGE_URL:-${NEOAG_SNAF_PACKAGE_URL:-https://gh-proxy.com/https://github.com/frankligy/SNAF/archive/${SNAF_GIT_REF}.tar.gz}}"
SNAF_PIP_INDEX_URL="${SNAF_PIP_INDEX_URL:-${NEOAG_SNAF_PIP_INDEX_URL:-https://mirrors.aliyun.com/pypi/simple}}"
SNAF_ARCHIVE_CACHE="${SNAF_ARCHIVE_CACHE:-${NEOAG_SNAF_ARCHIVE_CACHE:-${TOOLS_ROOT}/sources/SNAF-${SNAF_GIT_REF}.tar.gz}}"
ASNEO_SOURCE="${ASNEO_SOURCE:-${NEOAG_ASNEO_SOURCE:-}}"
NEOSPLICE_SOURCE="${NEOSPLICE_SOURCE:-${NEOAG_NEOSPLICE_SOURCE:-}}"
SPLICE2NEO_SOURCE="${SPLICE2NEO_SOURCE:-${NEOAG_SPLICE2NEO_SOURCE:-}}"

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

set +u
# shellcheck disable=SC1091
source "${CONDA_BASE}/etc/profile.d/conda.sh"
set -u

if [[ "${NEOAG_USE_MAMBA:-0}" == "1" ]] && command -v mamba >/dev/null 2>&1; then
  CONDA_RUNNER=mamba
else
  CONDA_RUNNER=conda
fi

env_exists() {
  local portable_root="${PORTABLE_ENVS%%:*}"
  if [[ -n "${portable_root}" && -f "${portable_root}/$1/conda-meta/history" ]]; then
    return 0
  fi
  conda_safe env list | awk '{print $1}' | grep -qx "$1"
}
env_has_regtools() { env_exists "$1" && conda_safe run -n "$1" regtools junctions extract -h >/dev/null 2>&1; }

# SNAF 0.7.0 (frankligy/SNAF@e23ce395) pins tensorflow==2.3.0, which only has
# wheels for Python 3.5-3.8. Installing SNAF as one pip transaction lets the
# resolver backtrack that stack onto TensorFlow 2.12+ wheels, or conda run can
# inherit an outer .venv/base interpreter. Skill1 therefore uses the env prefix
# python, pre-pins TF 2.3 / protobuf 3.20.3, then installs SNAF with --no-deps.
SNAF_ENV_PREFIX=""
SNAF_PYTHON_BIN=""

resolve_named_env_prefix() {
  local env_name="$1"
  local portable_root="${PORTABLE_ENVS%%:*}"
  if [[ -n "${portable_root}" && ( -x "${portable_root}/${env_name}/bin/python" || -f "${portable_root}/${env_name}/conda-meta/history" ) ]]; then
    printf '%s\n' "${portable_root}/${env_name}"
    return
  fi
  printf '%s\n' "${CONDA_BASE}/envs/${env_name}"
}

python_is_38() {
  [[ -x "$1" ]] || return 1
  "$1" -c 'import sys; raise SystemExit(0 if sys.version_info[:2] == (3, 8) else 1)'
}

ensure_snaf_env() {
  local prefix python_bin found_version
  prefix="$(resolve_named_env_prefix "${SNAF_ENV_NAME}")"
  python_bin="${prefix}/bin/python"
  if env_exists "${SNAF_ENV_NAME}" && python_is_38 "${python_bin}"; then
    echo "==> Reusing ${SNAF_ENV_NAME} at ${prefix} (Python 3.8)"
  else
    if env_exists "${SNAF_ENV_NAME}"; then
      found_version="$("${python_bin}" -c 'import sys; print(sys.version.split()[0])' 2>/dev/null || echo missing)"
      echo "==> Recreating ${SNAF_ENV_NAME}; SNAF TensorFlow 2.3 requires Python 3.8, found ${found_version}"
      if [[ "$(basename "${prefix}")" == "${SNAF_ENV_NAME}" && -d "${prefix}/conda-meta" ]]; then
        rm -rf "${prefix}"
      fi
      "${CONDA_RUNNER}" env remove -n "${SNAF_ENV_NAME}" -y >/dev/null 2>&1 || true
    fi
    echo "==> Creating ${SNAF_ENV_NAME} with Python 3.8 for SNAF's pinned TensorFlow 2.3 dependency"
    "${CONDA_RUNNER}" create -n "${SNAF_ENV_NAME}" --override-channels \
      -c conda-forge -c bioconda python=3.8 pip -y
    prefix="$(resolve_named_env_prefix "${SNAF_ENV_NAME}")"
    python_bin="${prefix}/bin/python"
  fi
  [[ -x "${python_bin}" ]] || { echo "ERROR: SNAF environment Python is missing: ${python_bin}" >&2; exit 45; }
  python_is_38 "${python_bin}" || { echo "ERROR: ${python_bin} is not Python 3.8" >&2; exit 45; }
  SNAF_ENV_PREFIX="${prefix}"
  SNAF_PYTHON_BIN="${SNAF_ENV_PREFIX}/bin/python"
}

snaf_pip() {
  "${SNAF_PYTHON_BIN}" -m pip install --index-url "${SNAF_PIP_INDEX_URL}" "$@"
}

install_snaf_package() {
  local package="$1"
  local constraints
  constraints="$(mktemp "${TMPDIR:-/tmp}/snaf-tf.XXXXXX")"
  cat > "${constraints}" <<'EOF'
tensorflow==2.3.0
protobuf==3.20.3
numpy==1.18.5
h5py==2.10.0
scipy==1.4.1
EOF
  echo "==> Pre-installing tensorflow==2.3.0 and protobuf==3.20.3 before SNAF so pip cannot backtrack onto TensorFlow 2.12+"
  snaf_pip -c "${constraints}" --prefer-binary "tensorflow==2.3.0" "protobuf==3.20.3"
  echo "==> Installing remaining SNAF 0.7.0 runtime pins with TensorFlow held at 2.3.0"
  snaf_pip -c "${constraints}" \
    "pandas==1.3.4" \
    "numpy==1.18.5" \
    "numba==0.53.0" \
    "mhcflurry==2.0.5" \
    "h5py==2.10.0" \
    "anndata==0.7.6" \
    "seaborn==0.11.2" \
    "biopython==1.79" \
    "requests==2.26.0" \
    "xmltodict==0.12.0" \
    "xmltramp2==3.1.1" \
    "tqdm==4.62.3" \
    "scipy==1.4.1" \
    "statsmodels==0.13.1" \
    "lifelines==0.26.4" \
    "umap-learn==0.5.2" \
    "plotly==5.4.0" \
    "Werkzeug==2.0.2" \
    "flask==2.0.2" \
    "dash==2.0.0" \
    "dash-dangerously-set-inner-html==0.0.2" \
    "mygene==3.2.2" \
    "adjustText==0.8"
  echo "==> Installing SNAF with --no-deps to keep the TensorFlow 2.3 stack"
  snaf_pip --no-deps "${package}"
  snaf_pip -c "${constraints}" "protobuf==3.20.3"
  rm -f "${constraints}"
}

write_snaf_wrapper() {
  cat > "${BIN_DIR}/snaf-neoag" <<EOF
#!/usr/bin/env bash
set -euo pipefail
exec "${SNAF_PYTHON_BIN}" -m snaf "\$@"
EOF
  chmod +x "${BIN_DIR}/snaf-neoag"
}

if [[ "${NEOAG_FORCE_ENV_UPDATE:-0}" == "1" ]]; then
  if env_exists "${ENV_NAME}"; then
    "${CONDA_RUNNER}" env update -n "${ENV_NAME}" -f "${YML}" --prune
  else
    "${CONDA_RUNNER}" env create -n "${ENV_NAME}" -f "${YML}" -y
  fi
elif env_has_regtools "${ENV_NAME}"; then
  echo "==> ${ENV_NAME} already has RegTools; skipping env update"
elif env_exists "${ENV_NAME}"; then
  "${CONDA_RUNNER}" env update -n "${ENV_NAME}" -f "${YML}" --prune
else
  "${CONDA_RUNNER}" env create -n "${ENV_NAME}" -f "${YML}" -y
fi

mkdir -p "${BIN_DIR}" "${ROOT}/conf" "${TOOLS_ROOT}/tools"

cat > "${BIN_DIR}/regtools-neoag" <<EOF
#!/usr/bin/env bash
set -euo pipefail
exec "${CONDA_BASE}/bin/conda" run -n "${ENV_NAME}" regtools "\$@"
EOF
[[ -z "$PORTABLE_ENVS" ]] || sed -i "2a export CONDA_ENVS_PATH=\"${PORTABLE_ENVS}\"" "${BIN_DIR}/regtools-neoag"
chmod +x "${BIN_DIR}/regtools-neoag"

PVACSPLICE_BIN="${CONDA_BASE}/envs/neoag-tools/bin/pvacsplice"
if [[ -x "${PVACSPLICE_BIN}" ]]; then
  cat > "${BIN_DIR}/pvacsplice-neoag" <<EOF
#!/usr/bin/env bash
set -euo pipefail
exec "${PVACSPLICE_BIN}" "\$@"
EOF
  chmod +x "${BIN_DIR}/pvacsplice-neoag"
else
  echo "WARN: pvacsplice not found in neoag-tools; install core env first." >&2
fi

install_python_source() {
  local label="$1" source_dir="$2" target_dir="$3" wrapper="$4" module="$5"
  local install_env="${6:-${ENV_NAME}}"
  [[ -n "${source_dir}" ]] || return 0
  [[ -e "${source_dir}" ]] || { echo "ERROR: ${label} source missing: ${source_dir}" >&2; exit 42; }
  mkdir -p "$(dirname "${target_dir}")"
  rsync -a --delete "${source_dir}/" "${target_dir}/"
  if [[ -f "${target_dir}/pyproject.toml" || -f "${target_dir}/setup.py" ]]; then
    "${CONDA_BASE}/bin/conda" run -n "${install_env}" python -m pip install "${target_dir}"
  fi
  cat > "${BIN_DIR}/${wrapper}" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export PYTHONPATH="${target_dir}:\${PYTHONPATH:-}"
exec "${CONDA_BASE}/bin/conda" run -n "${install_env}" python -m ${module} "\$@"
EOF
  chmod +x "${BIN_DIR}/${wrapper}"
}


download_verified_tarball() {
  local label="$1" url="$2" dest="$3"
  local attempts="${NEOAG_DOWNLOAD_RETRIES:-5}"
  local attempt part
  mkdir -p "$(dirname "${dest}")"
  for attempt in $(seq 1 "${attempts}"); do
    part="${dest}.part.${attempt}"
    rm -f "${part}"
    echo "==> Downloading ${label} snapshot (attempt ${attempt}/${attempts})"
    if curl -fL --retry 3 --connect-timeout 20 --max-time "${NEOAG_DOWNLOAD_MAX_TIME:-1800}" -o "${part}" "${url}" \
      && [[ -s "${part}" ]] \
      && tar -tzf "${part}" >/dev/null 2>&1; then
      mv "${part}" "${dest}"
      return 0
    fi
    rm -f "${part}"
    sleep "${NEOAG_DOWNLOAD_RETRY_SLEEP:-10}"
  done
  echo "ERROR: ${label} source snapshot download failed or produced an incomplete tarball: ${url}" >&2
  exit 44
}

if [[ "${INSTALL_SNAF}" == "1" ]]; then
  SNAF_HOME="${TOOLS_ROOT}/tools/SNAF"
  ensure_snaf_env
  if [[ -n "${SNAF_SOURCE}" ]]; then
    [[ -e "${SNAF_SOURCE}" ]] || { echo "ERROR: SNAF source missing: ${SNAF_SOURCE}" >&2; exit 42; }
    mkdir -p "$(dirname "${SNAF_HOME}")"
    rsync -a --delete "${SNAF_SOURCE}/" "${SNAF_HOME}/"
    install_snaf_package "${SNAF_HOME}"
    write_snaf_wrapper
  else
    echo "==> Installing SNAF from pinned source snapshot ${SNAF_PACKAGE_URL}"
    if [[ ! -s "${SNAF_ARCHIVE_CACHE}" ]]; then
      command -v curl >/dev/null 2>&1 || { echo "ERROR: curl is required to download SNAF" >&2; exit 43; }
      mkdir -p "$(dirname "${SNAF_ARCHIVE_CACHE}")"
      download_verified_tarball "SNAF" "${SNAF_PACKAGE_URL}" "${SNAF_ARCHIVE_CACHE}"
    else
      echo "==> Reusing cached SNAF snapshot ${SNAF_ARCHIVE_CACHE}"
    fi
    install_snaf_package "${SNAF_ARCHIVE_CACHE}"
    write_snaf_wrapper
  fi
else
  echo "==> Skipping SNAF because NEOAG_INSTALL_SNAF=${INSTALL_SNAF}"
fi

if [[ "${INSTALL_SPLICEMUTR}" == "1" ]]; then
  NEOAG_CONDA_BASE="${CONDA_BASE}" NEOAG_TOOLS_ROOT="${TOOLS_ROOT}" \
    bash "${ROOT}/scripts/install_splicemutr.sh"
else
  echo "==> Skipping SpliceMutr because NEOAG_INSTALL_SPLICEMUTR=${INSTALL_SPLICEMUTR}"
fi

install_python_source "ASNEO" "${ASNEO_SOURCE}" "${TOOLS_ROOT}/tools/ASNEO" "asneo-neoag" "ASNEO"
install_python_source "NeoSplice" "${NEOSPLICE_SOURCE}" "${TOOLS_ROOT}/tools/NeoSplice" "neosplice-neoag" "NeoSplice"
install_python_source "splice2neo" "${SPLICE2NEO_SOURCE}" "${TOOLS_ROOT}/tools/splice2neo" "splice2neo-neoag" "splice2neo"

if [[ -f "${TOOLS_ENV}" ]]; then
  if grep -q 'NEOAG_SPLICE_ENV' "${TOOLS_ENV}"; then
    sed -i "s|^export NEOAG_SPLICE_ENV=.*|export NEOAG_SPLICE_ENV="${ENV_NAME}"|" "${TOOLS_ENV}"
    sed -i "s|^export NEOAG_REGTOOLS_BIN=.*|export NEOAG_REGTOOLS_BIN="${BIN_DIR}/regtools-neoag"|" "${TOOLS_ENV}"
    sed -i "s|^export NEOAG_PVACSPLICE_BIN=.*|export NEOAG_PVACSPLICE_BIN="${BIN_DIR}/pvacsplice-neoag"|" "${TOOLS_ENV}"
  else
    cat >> "${TOOLS_ENV}" <<EOF

# Splice neoantigen tools
export NEOAG_SPLICE_ENV="${ENV_NAME}"
export NEOAG_REGTOOLS_BIN="${BIN_DIR}/regtools-neoag"
export NEOAG_PVACSPLICE_BIN="${BIN_DIR}/pvacsplice-neoag"
export PATH="${BIN_DIR}:${CONDA_BASE}/envs/${ENV_NAME}/bin:\$PATH"
EOF
  fi
else
  cat > "${TOOLS_ENV}" <<EOF
export NEOAG_SPLICE_ENV="${ENV_NAME}"
export NEOAG_REGTOOLS_BIN="${BIN_DIR}/regtools-neoag"
export NEOAG_PVACSPLICE_BIN="${BIN_DIR}/pvacsplice-neoag"
export PATH="${BIN_DIR}:${CONDA_BASE}/envs/${ENV_NAME}/bin:\$PATH"
EOF
fi

if [[ -f "${TOOLS_ENV}" ]]; then
  if grep -q '^export NEOAG_SNAF_ENV=' "${TOOLS_ENV}"; then
    sed -i "s|^export NEOAG_SNAF_ENV=.*|export NEOAG_SNAF_ENV=\"${SNAF_ENV_NAME}\"|" "${TOOLS_ENV}"
  else
    echo "export NEOAG_SNAF_ENV=\"${SNAF_ENV_NAME}\"" >> "${TOOLS_ENV}"
  fi
  if grep -q '^export NEOAG_SNAF_BIN=' "${TOOLS_ENV}"; then
    sed -i "s|^export NEOAG_SNAF_BIN=.*|export NEOAG_SNAF_BIN=\"${BIN_DIR}/snaf-neoag\"|" "${TOOLS_ENV}"
  else
    echo "export NEOAG_SNAF_BIN=\"${BIN_DIR}/snaf-neoag\"" >> "${TOOLS_ENV}"
  fi
  snaf_python="${SNAF_PYTHON_BIN:-}"
  if [[ -z "${snaf_python}" ]]; then
    snaf_python="$(resolve_named_env_prefix "${SNAF_ENV_NAME}")/bin/python"
  fi
  [[ -x "${snaf_python}" ]] || {
    echo "ERROR: SNAF environment Python is missing: ${snaf_python}" >&2
    exit 45
  }
  if grep -q '^export SNAF_PYTHON=' "${TOOLS_ENV}"; then
    sed -i "s|^export SNAF_PYTHON=.*|export SNAF_PYTHON=\"${snaf_python}\"|" "${TOOLS_ENV}"
  else
    echo "export SNAF_PYTHON=\"${snaf_python}\"" >> "${TOOLS_ENV}"
  fi
fi

echo "==> Splice tools smoke"
"${BIN_DIR}/regtools-neoag" junctions extract -h | head -8 || true
if [[ -x "${BIN_DIR}/pvacsplice-neoag" ]]; then
  "${BIN_DIR}/pvacsplice-neoag" --help | head -8 || true
fi
if [[ -x "${BIN_DIR}/snaf-neoag" ]]; then
  snaf_smoke_dir="$(mktemp -d)"
  mkdir -p "${snaf_smoke_dir}/assets"
  (
    cd "${snaf_smoke_dir}"
    "${SNAF_PYTHON_BIN}" -c \
      'import snaf, tensorflow as tf; assert tf.__version__.startswith("2.3"), tf.__version__; print("SNAF import OK; TensorFlow " + tf.__version__)'
  )
  rm -rf "${snaf_smoke_dir}"
fi
if [[ -x "${BIN_DIR}/splicemutr-neoag" ]]; then
  "${BIN_DIR}/splicemutr-neoag" doctor
fi

echo "==> Done. Run: bash scripts/run_splice_tool_smoke.sh"
