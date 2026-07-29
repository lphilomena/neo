#!/usr/bin/env bash
# Install the current official SplAdder, IRFinder-S, ImmunoPepper and pVACbind releases.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONDA_BASE="${NEOAG_CONDA_BASE:-${CONDA_PREFIX:-${HOME}/miniforge3}}"
ENV_ROOT="${NEOAG_ENV_ROOT:-${CONDA_BASE}/envs}"
SPLADDER_ENV="${NEOAG_SPLADDER_ENV:-${ENV_ROOT}/neoag-spladder}"
IMMUNOPEPPER_ENV="${NEOAG_IMMUNOPEPPER_ENV:-${ENV_ROOT}/neoag-immunopepper}"
TOOLS_ROOT="${NEOAG_SPLICE_TOOLS_ROOT:-${NEOAG_TOOLS_ROOT:-${ROOT}}/tools/splice}"
BIN_DIR="${NEOAG_BIN_DIR:-${ROOT}/bin}"

SPLADDER_VERSION="${NEOAG_SPLADDER_VERSION:-3.1.1}"
IMMUNOPEPPER_VERSION="${NEOAG_IMMUNOPEPPER_VERSION:-2.0.0}"
IMMUNOPEPPER_REF="${NEOAG_IMMUNOPEPPER_REF:-aedd6ed6af1a0baf7d30f21111e628e84e167cdd}"
PVACTOOLS_VERSION="${NEOAG_PVACTOOLS_VERSION:-7.1.1}"
PVACTOOLS_BASE_ENV="${NEOAG_PVACTOOLS_BASE_ENV:-${CONDA_BASE}/envs/neoag-tools}"
IRFINDER_VERSION="${NEOAG_IRFINDER_VERSION:-2.0.1}"
IRFINDER_IMAGE="${NEOAG_IRFINDER_IMAGE:-cloxd/irfinder:${IRFINDER_VERSION}}"

DO_SPLADDER=0 DO_IRFINDER=0 DO_IMMUNOPEPPER=0 DO_PVACBIND=0

usage() {
  cat <<'EOF'
Usage: bash scripts/install_extended_splice_tools.sh [targets]

Targets:
  --spladder       Install SplAdder 3.1.1 in an isolated Python environment
  --irfinder-s     Pull the official IRFinder-S 2.0.1 Docker image
  --immunopepper   Install ImmunoPepper 2.0.0 from a pinned official commit
  --pvacbind       Install pVACtools 7.1.1, including pVACbind, side by side
  --all            Install all four tools (default)

Portable paths:
  NEOAG_CONDA_BASE=/path/to/miniforge
  NEOAG_ENV_ROOT=/large/path/to/envs
  NEOAG_SPLADDER_ENV=/path/to/spladder/env
  NEOAG_IMMUNOPEPPER_ENV=/path/to/immunopepper/env
  NEOAG_SPLICE_TOOLS_ROOT=/large/path/to/tool_sources
  NEOAG_BIN_DIR=/path/to/project/bin

IRFinder-S wrapper boundary:
  Set NEOAG_IRFINDER_WORKDIR to a directory containing inputs/references and
  pass container paths below /work. The wrapper does not mount the host root.
EOF
}

if [[ $# -eq 0 ]]; then set -- --all; fi
while [[ $# -gt 0 ]]; do
  case "$1" in
    --spladder) DO_SPLADDER=1 ;;
    --irfinder-s) DO_IRFINDER=1 ;;
    --immunopepper) DO_IMMUNOPEPPER=1 ;;
    --pvacbind) DO_PVACBIND=1 ;;
    --all) DO_SPLADDER=1; DO_IRFINDER=1; DO_IMMUNOPEPPER=1; DO_PVACBIND=1 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown target: $1" >&2; usage >&2; exit 2 ;;
  esac
  shift
done

MAMBA="${CONDA_BASE}/bin/mamba"
[[ -x "${MAMBA}" ]] || MAMBA="${CONDA_BASE}/bin/conda"
[[ -x "${MAMBA}" ]] || { echo "ERROR: conda/mamba missing under ${CONDA_BASE}" >&2; exit 3; }
mkdir -p "${ENV_ROOT}" "${TOOLS_ROOT}" "${BIN_DIR}"

create_env() {
  local env="$1" python_version="$2"
  if [[ ! -x "${env}/bin/python" ]]; then
    "${MAMBA}" create -y -p "${env}" -c conda-forge "python=${python_version}" pip
  fi
}

write_python_wrapper() {
  local output="$1" executable="$2"
  cat > "${output}" <<EOF
#!/usr/bin/env bash
set -euo pipefail
exec "${executable}" "\$@"
EOF
  chmod +x "${output}"
}

write_pvactools_overlay_wrapper() {
  local output="$1" executable="$2" overlay="$3"
  cat > "${output}" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export PYTHONPATH="${overlay}${PYTHONPATH:+:${PYTHONPATH}}"
exec "${executable}" "\$@"
EOF
  chmod +x "${output}"
}

if [[ "${DO_SPLADDER}" == 1 ]]; then
  env="${SPLADDER_ENV}"
  create_env "${env}" 3.11
  "${env}/bin/python" -m pip install --upgrade "spladder==${SPLADDER_VERSION}"
  write_python_wrapper "${BIN_DIR}/spladder-neoag" "${env}/bin/spladder"
  "${env}/bin/python" -c \
    "import spladder; assert spladder.__version__ == '${SPLADDER_VERSION}'"
fi

if [[ "${DO_IMMUNOPEPPER}" == 1 ]]; then
  env="${IMMUNOPEPPER_ENV}"
  source_dir="${TOOLS_ROOT}/ImmunoPepper-${IMMUNOPEPPER_VERSION}"
  create_env "${env}" 3.9
  if [[ -e "${source_dir}" ]]; then
    if [[ ! -d "${source_dir}/.git" ]] || \
       [[ "$(git -C "${source_dir}" rev-parse HEAD 2>/dev/null || true)" != "${IMMUNOPEPPER_REF}" ]]; then
      echo "ERROR: existing ImmunoPepper source is not the pinned sparse checkout: ${source_dir}" >&2
      echo "Move it aside or set NEOAG_SPLICE_TOOLS_ROOT to another directory." >&2
      exit 5
    fi
  else
    mkdir -p "${source_dir}"
    git -C "${source_dir}" init
    git -C "${source_dir}" remote add origin https://github.com/ratschlab/immunopepper.git
    git -C "${source_dir}" config remote.origin.promisor true
    git -C "${source_dir}" config remote.origin.partialclonefilter blob:none
    git -C "${source_dir}" sparse-checkout init --cone
    git -C "${source_dir}" sparse-checkout set immunopepper
    git -C "${source_dir}" fetch --depth=1 --filter=blob:none origin "${IMMUNOPEPPER_REF}"
    git -C "${source_dir}" checkout --detach FETCH_HEAD
  fi
  "${env}/bin/python" -m pip install --upgrade pip setuptools wheel "cython==3.0.11"
  "${env}/bin/python" -m pip install \
    "numpy==1.26.4" "pandas==2.2.3" "scipy>=1.11,<1.14" h5py \
    "psutil==6.1.0" "pyarrow==14.0.2" "pyspark==3.5.1" \
    "tblib==3.0.0" "pyyaml==6.0.2" biopython mhctools pyensembl \
    "spladder==3.1.0"
  "${env}/bin/python" -m pip install --no-deps "${source_dir}"
  write_python_wrapper "${BIN_DIR}/immunopepper-neoag" "${env}/bin/immunopepper"
  "${BIN_DIR}/immunopepper-neoag" --help >/dev/null
fi

if [[ "${DO_PVACBIND}" == 1 ]]; then
  overlay="${TOOLS_ROOT}/pvactools-${PVACTOOLS_VERSION}-overlay"
  [[ -x "${PVACTOOLS_BASE_ENV}/bin/python" ]] || {
    echo "ERROR: tested pVACtools base environment is missing: ${PVACTOOLS_BASE_ENV}" >&2
    echo "Install the standard neoag-tools environment or set NEOAG_PVACTOOLS_BASE_ENV." >&2
    exit 6
  }
  "${PVACTOOLS_BASE_ENV}/bin/python" -m pip install --no-deps --upgrade \
    --target "${overlay}" "pvactools==${PVACTOOLS_VERSION}"
  for tool in pvacbind pvacseq pvacfuse pvacsplice; do
    write_pvactools_overlay_wrapper \
      "${BIN_DIR}/${tool}-neoag7" "${PVACTOOLS_BASE_ENV}/bin/${tool}" "${overlay}"
    "${BIN_DIR}/${tool}-neoag7" --help >/dev/null
  done
  PYTHONPATH="${overlay}${PYTHONPATH:+:${PYTHONPATH}}" \
    "${PVACTOOLS_BASE_ENV}/bin/python" -c \
    "import importlib.metadata as m; assert m.version('pvactools') == '${PVACTOOLS_VERSION}'"
fi

if [[ "${DO_IRFINDER}" == 1 ]]; then
  command -v docker >/dev/null 2>&1 || { echo "ERROR: Docker is required for IRFinder-S" >&2; exit 4; }
  docker pull "${IRFINDER_IMAGE}"
  cat > "${BIN_DIR}/irfinder-s-neoag" <<EOF
#!/usr/bin/env bash
set -euo pipefail
workdir="\${NEOAG_IRFINDER_WORKDIR:-\$PWD}"
[[ -d "\${workdir}" ]] || { echo "ERROR: NEOAG_IRFINDER_WORKDIR is not a directory: \${workdir}" >&2; exit 2; }
exec docker run --rm --user "\$(id -u):\$(id -g)" \
  -v "\${workdir}:/work" -w /work "${IRFINDER_IMAGE}" "\$@"
EOF
  chmod +x "${BIN_DIR}/irfinder-s-neoag"
  irfinder_help="$("${BIN_DIR}/irfinder-s-neoag" --help 2>&1 || true)"
  grep -q 'IRFinder version: 2.0.1' <<<"${irfinder_help}"
fi

versions="${TOOLS_ROOT}/extended_splice_versions.tsv"
printf 'tool\tversion\tinstall_mode\tlocation\n' > "${versions}"
[[ -x "${SPLADDER_ENV}/bin/spladder" ]] && printf 'SplAdder\t%s\tconda_pip\t%s\n' "${SPLADDER_VERSION}" "${SPLADDER_ENV}" >> "${versions}"
[[ -x "${IMMUNOPEPPER_ENV}/bin/immunopepper" ]] && printf 'ImmunoPepper\t%s@%s\tofficial_source\t%s\n' "${IMMUNOPEPPER_VERSION}" "${IMMUNOPEPPER_REF}" "${IMMUNOPEPPER_ENV}" >> "${versions}"
[[ -d "${TOOLS_ROOT}/pvactools-${PVACTOOLS_VERSION}-overlay/pvactools" ]] && printf 'pVACtools/pVACbind\t%s\tisolated_overlay\t%s\n' "${PVACTOOLS_VERSION}" "${TOOLS_ROOT}/pvactools-${PVACTOOLS_VERSION}-overlay" >> "${versions}"
docker image inspect "${IRFINDER_IMAGE}" >/dev/null 2>&1 && printf 'IRFinder-S\t%s\tdocker\t%s\n' "${IRFINDER_VERSION}" "${IRFINDER_IMAGE}" >> "${versions}"

echo "Installed versions:"
cat "${versions}"
echo "Run: bash scripts/verify_extended_splice_tools.sh"
