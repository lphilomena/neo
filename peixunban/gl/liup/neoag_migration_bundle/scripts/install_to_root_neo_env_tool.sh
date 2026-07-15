#!/usr/bin/env bash
set -euo pipefail

if [[ "$(id -u)" -ne 0 ]]; then
  echo "ERROR: run this script as root; target install root defaults to /root/neo/env_tool" >&2
  exit 1
fi

REPO_URL="${REPO_URL:-https://github.com/lphilomena/neo.git}"
TARGET_ROOT="${TARGET_ROOT:-/root/neo}"
TOOLS_ROOT="${NEOAG_TOOLS_ROOT:-${TARGET_ROOT}/env_tool}"
SRC_ROOT="${SRC_ROOT:-${TARGET_ROOT}/src}"
PRIMARY_BRANCH="${PRIMARY_BRANCH:-neo-na0707_upload_release}"
SECONDARY_BRANCH="${SECONDARY_BRANCH:-na0707_upload_release}"
MINIFORGE_URL="${MINIFORGE_URL:-https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh}"

log() { printf "[neo-root-install] %s\n" "$*"; }

mkdir -p "${TARGET_ROOT}" "${TOOLS_ROOT}" "${SRC_ROOT}" "${TOOLS_ROOT}/logs"

clone_or_update() {
  local branch="$1"
  local dest="${SRC_ROOT}/${branch}"
  if [[ -d "${dest}/.git" ]]; then
    log "updating ${branch} in ${dest}"
    git -C "${dest}" fetch origin "${branch}"
    git -C "${dest}" checkout "${branch}"
    git -C "${dest}" reset --hard "origin/${branch}"
  else
    log "cloning ${branch} to ${dest}"
    git clone -b "${branch}" --single-branch "${REPO_URL}" "${dest}"
  fi
}

clone_or_update "${PRIMARY_BRANCH}"
clone_or_update "${SECONDARY_BRANCH}" || log "WARN: secondary branch ${SECONDARY_BRANCH} is optional and was not cloned"

PRIMARY_REPO="${SRC_ROOT}/${PRIMARY_BRANCH}"
BUNDLE_TGZ="${PRIMARY_REPO}/peixunban/gl/liup/neoag_migration_bundle_20260714.tar.gz"
PATCH_TGZ="${PRIMARY_REPO}/peixunban/gl/liup/test_env_tools_refs_neodata4git_patch_20260714.tar.gz"
[[ -f "${BUNDLE_TGZ}" ]] || { echo "ERROR: missing ${BUNDLE_TGZ}" >&2; exit 2; }
[[ -f "${PATCH_TGZ}" ]] || { echo "ERROR: missing ${PATCH_TGZ}" >&2; exit 2; }

log "extracting migration bundle"
tar -xzf "${BUNDLE_TGZ}" -C "${TARGET_ROOT}"

log "extracting env/tool patch into ${TOOLS_ROOT}"
TMP_PATCH="${TARGET_ROOT}/.patch_extract.$$"
rm -rf "${TMP_PATCH}"
mkdir -p "${TMP_PATCH}"
tar -xzf "${PATCH_TGZ}" -C "${TMP_PATCH}"
cp -a "${TMP_PATCH}/test_env_tools/." "${TOOLS_ROOT}/"
rm -rf "${TMP_PATCH}"

if [[ -n "${NEODATA_TARBALL:-}" ]]; then
  log "extracting neodata companion from ${NEODATA_TARBALL}"
  tar -xzf "${NEODATA_TARBALL}" -C "${TARGET_ROOT}"
fi

if [[ ! -x "${TOOLS_ROOT}/miniforge3/bin/conda" ]]; then
  log "installing Miniforge under ${TOOLS_ROOT}/miniforge3"
  if [[ -n "${MINIFORGE_INSTALLER:-}" && -f "${MINIFORGE_INSTALLER}" ]]; then
    bash "${MINIFORGE_INSTALLER}" -b -p "${TOOLS_ROOT}/miniforge3"
  else
    curl -L "${MINIFORGE_URL}" -o "${TOOLS_ROOT}/logs/Miniforge3-Linux-x86_64.sh"
    bash "${TOOLS_ROOT}/logs/Miniforge3-Linux-x86_64.sh" -b -p "${TOOLS_ROOT}/miniforge3"
  fi
fi

export NEOAG_TOOLS_ROOT="${TOOLS_ROOT}"
export NEOAG_CONDA_BASE="${TOOLS_ROOT}/miniforge3"
export NEOAG_BUNDLE_ROOT="${TARGET_ROOT}/neoag_migration_bundle"
export PROJECT_ROOT="${NEOAG_BUNDLE_ROOT}/project/neo-na0707_upload_release-3"
export NEODATA_ROOT="${NEODATA_ROOT:-${TARGET_ROOT}/neodata4git}"
export PATH="${NEOAG_CONDA_BASE}/bin:${PATH}"

log "installing tier1 core"
bash "${NEOAG_BUNDLE_ROOT}/scripts/install_tier1_core.sh"
log "installing tier2 tools"
bash "${NEOAG_BUNDLE_ROOT}/scripts/install_tier2_tools.sh"

log "activation file: ${TOOLS_ROOT}/activate_neoag_production_refs.sh"
log "run: source ${TOOLS_ROOT}/activate_neoag_production_refs.sh && neoag-v03 --help"
if [[ -d "${NEODATA_ROOT}" ]]; then
  log "running doctor with NEODATA_ROOT=${NEODATA_ROOT}"
  bash "${TOOLS_ROOT}/run_doctor_neodata4git.sh" || true
else
  log "NEODATA_ROOT not found at ${NEODATA_ROOT}; copy/mount companion neodata4git before real production runs"
fi
