#!/usr/bin/env bash
# Install the EasyFuse-native fusion stack by default.
#
# EasyFuse owns STAR, Arriba, STAR-Fusion, and FusionCatcher module definitions.
# A standalone caller stack remains available only as an explicit compatibility
# fallback: NEOAG_FUSION_INSTALL_MODE=standalone.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CONDA_BASE="${NEOAG_CONDA_BASE:-$(command conda info --base)}"
TOOLS_ROOT="${NEOAG_TOOLS_ROOT:-${ROOT}}"
INSTALL_MODE="${NEOAG_FUSION_INSTALL_MODE:-easyfuse}"
PROJECT_BIN="${ROOT}/bin"
TOOLS_ENV="${NEOAG_TOOLS_ENV:-${ROOT}/conf/tools.env.local.sh}"
BIOC_CACHE_HELPER="${ROOT}/.agents/skills/neoag-remote-deploy/scripts/with_bioc_data_cache.sh"
CACHE_ROOT="${NEOAG_INSTALL_CACHE_ROOT:-${TOOLS_ROOT}/install_cache}"

[[ -x "${CONDA_BASE}/bin/conda" ]] || {
  echo "ERROR: conda not found at ${CONDA_BASE}/bin/conda" >&2
  exit 1
}

conda_safe() {
  set +u
  "${CONDA_BASE}/bin/conda" "$@"
  local rc=$?
  set -u
  return "$rc"
}

install_github_snapshot() {
  local label="$1" url="$2" target="$3"
  local work archive
  work="$(mktemp -d)"
  archive="${work}/source.tar.gz"
  echo "==> Downloading ${label} source snapshot"
  if command -v aria2c >/dev/null 2>&1; then
    aria2c -x 8 -s 8 -k 1M --allow-overwrite=true --file-allocation=none \
      -d "${work}" -o source.tar.gz "${url}"
  else
    curl -fL --retry 5 --retry-all-errors --continue-at - -o "${archive}" "${url}"
  fi
  mkdir -p "${work}/source" "${target}"
  tar -xzf "${archive}" -C "${work}/source" --strip-components=1
  rsync -a --delete "${work}/source/" "${target}/"
  rm -rf "${work}"
}

write_easyfuse_environment() {
  local easyfuse_home="$1" env_prefix="$2"
  mkdir -p "${PROJECT_BIN}" "${ROOT}/conf"
  cat > "${PROJECT_BIN}/nextflow" <<EOF
#!/usr/bin/env bash
set -euo pipefail
exec "${env_prefix}/bin/nextflow" "\$@"
EOF
  chmod +x "${PROJECT_BIN}/nextflow"

  touch "${TOOLS_ENV}"
  sed -i '/^# BEGIN OPEN-NEO FUSION$/,/^# END OPEN-NEO FUSION$/d' "${TOOLS_ENV}"
  cat >> "${TOOLS_ENV}" <<EOF

# BEGIN OPEN-NEO FUSION
# EasyFuse-native fusion stack
export NEOAG_FUSION_INSTALL_MODE="easyfuse"
export NEOAG_EASYFUSE_HOME="${easyfuse_home}"
export NEOAG_EASYFUSE_ENV_PREFIX="${env_prefix}"
export NEOAG_FUSION_ENV="$(basename "${env_prefix}")"
export NEOAG_NEXTFLOW="${env_prefix}/bin/nextflow"
export PATH="${PROJECT_BIN}:${env_prefix}/bin:\$PATH"
# END OPEN-NEO FUSION
EOF
}

install_easyfuse_native() {
  local version="${NEOAG_EASYFUSE_VERSION:-2.2.1}"
  local easyfuse_home="${NEOAG_EASYFUSE_HOME:-${TOOLS_ROOT}/tools/EasyFuse}"
  if [[ -z "${NEOAG_EASYFUSE_HOME:-}" \
    && -f "${easyfuse_home}/main.nf" \
    && ! -f "${easyfuse_home}/modules/arriba/environment.yml" ]]; then
    easyfuse_home="${TOOLS_ROOT}/tools/EasyFuse-v${version}"
    echo "==> Preserving legacy EasyFuse source; installing module-native release at ${easyfuse_home}"
  fi
  local env_prefix="${NEOAG_EASYFUSE_ENV_PREFIX:-${CONDA_BASE}/envs/neoag-easyfuse}"
  local env_yml="${easyfuse_home}/environment.yml"

  if [[ ! -f "${easyfuse_home}/main.nf" ]]; then
    mkdir -p "$(dirname "${easyfuse_home}")"
    install_github_snapshot "EasyFuse v${version}" \
      "${NEOAG_EASYFUSE_SNAPSHOT_URL:-https://github.com/TRON-Bioinformatics/EasyFuse/archive/refs/tags/v${version}.tar.gz}" \
      "${easyfuse_home}"
  else
    echo "==> EasyFuse source already present: ${easyfuse_home}"
  fi

  [[ -f "${env_yml}" ]] || {
    echo "ERROR: EasyFuse driver environment is missing: ${env_yml}" >&2
    exit 2
  }

  if [[ "${NEOAG_FORCE_ENV_UPDATE:-0}" == "1" && -x "${env_prefix}/bin/nextflow" ]]; then
    echo "==> Updating EasyFuse driver environment: ${env_prefix}"
    conda_safe env update -p "${env_prefix}" -f "${env_yml}" --prune
  elif [[ ! -x "${env_prefix}/bin/nextflow" ]]; then
    echo "==> Creating EasyFuse driver environment: ${env_prefix}"
    conda_safe env create -p "${env_prefix}" -f "${env_yml}" -y
  else
    echo "==> EasyFuse driver environment already ready: ${env_prefix}"
  fi

  local marker missing=0
  for marker in \
    modules/arriba/environment.yml \
    modules/starfusion/starfusion/environment.yml \
    modules/fusioncatcher/environment.yml; do
    if [[ ! -f "${easyfuse_home}/${marker}" ]]; then
      echo "ERROR: EasyFuse internal caller definition missing: ${marker}" >&2
      missing=1
    fi
  done
  [[ "${missing}" == "0" ]] || exit 3

  "${env_prefix}/bin/nextflow" -version | head -3
  write_easyfuse_environment "${easyfuse_home}" "${env_prefix}"
  echo "==> EasyFuse owns STAR, Arriba, STAR-Fusion, and FusionCatcher; no standalone caller copies were installed."
}

install_standalone_fallback() {
  local env_name="${NEOAG_FUSION_ENV:-neoag-fusion}"
  local env_prefix="${CONDA_BASE}/envs/${env_name}"
  local fusion_yml="${ROOT}/conda/env.neoag-fusion.yml"
  local star_home="${NEOAG_STAR_FUSION_HOME:-${TOOLS_ROOT}/tools/STAR-Fusion}"
  local fc_home="${NEOAG_FUSIONCATCHER_HOME:-${TOOLS_ROOT}/tools/fusioncatcher}"
  local package_key="${NEOAG_FUSION_BIOC_PACKAGE_KEY:-genomeinfodbdata-1.2.11}"
  local runner="${CONDA_BASE}/bin/conda"

  [[ -x "${BIOC_CACHE_HELPER}" ]] || {
    echo "ERROR: Bioconductor cache helper missing: ${BIOC_CACHE_HELPER}" >&2
    exit 4
  }
  if [[ "${NEOAG_USE_MAMBA:-0}" == "1" && -x "${CONDA_BASE}/bin/mamba" ]]; then
    runner="${CONDA_BASE}/bin/mamba"
  fi

  local transaction=(
    "${BIOC_CACHE_HELPER}"
    --conda-base "${CONDA_BASE}"
    --cache-root "${CACHE_ROOT}"
    --package-key "${package_key}"
    -- "${runner}"
  )
  if [[ -x "${env_prefix}/bin/arriba" ]]; then
    echo "==> Standalone fusion environment already ready: ${env_prefix}"
  elif [[ -d "${env_prefix}/conda-meta" ]]; then
    "${transaction[@]}" env update -p "${env_prefix}" -f "${fusion_yml}" --prune
  else
    "${transaction[@]}" env create -p "${env_prefix}" -f "${fusion_yml}" -y
  fi

  if [[ ! -x "${star_home}/STAR-Fusion" ]]; then
    install_github_snapshot "STAR-Fusion standalone fallback" \
      "${NEOAG_STAR_FUSION_SNAPSHOT_URL:-https://github.com/STAR-Fusion/STAR-Fusion/archive/refs/heads/master.tar.gz}" \
      "${star_home}"
    chmod +x "${star_home}/STAR-Fusion" 2>/dev/null || true
  fi
  if [[ ! -d "${fc_home}/bin" ]]; then
    install_github_snapshot "FusionCatcher standalone fallback" \
      "${NEOAG_FUSIONCATCHER_SNAPSHOT_URL:-https://github.com/ndaniel/fusioncatcher/archive/refs/heads/master.tar.gz}" \
      "${fc_home}"
  fi

  conda_safe run -p "${env_prefix}" arriba -h | head -6
  touch "${TOOLS_ENV}"
  sed -i '/^# BEGIN OPEN-NEO FUSION$/,/^# END OPEN-NEO FUSION$/d' "${TOOLS_ENV}"
  cat >> "${TOOLS_ENV}" <<EOF

# BEGIN OPEN-NEO FUSION
# Standalone fusion compatibility stack
export NEOAG_FUSION_INSTALL_MODE="standalone"
export NEOAG_FUSION_ENV="${env_name}"
export NEOAG_ARRIBA_BIN="${env_prefix}/bin/arriba"
export NEOAG_STAR_FUSION_HOME="${star_home}"
export NEOAG_FUSIONCATCHER_HOME="${fc_home}"
export PATH="${PROJECT_BIN}:${env_prefix}/bin:\$PATH"
# END OPEN-NEO FUSION
EOF
}

case "${INSTALL_MODE}" in
  easyfuse)
    install_easyfuse_native
    ;;
  standalone)
    echo "WARN: installing explicit standalone fusion compatibility stack" >&2
    install_standalone_fallback
    ;;
  *)
    echo "ERROR: unsupported NEOAG_FUSION_INSTALL_MODE=${INSTALL_MODE}; use easyfuse or standalone" >&2
    exit 64
    ;;
esac

echo "==> Fusion installation complete (${INSTALL_MODE})"
