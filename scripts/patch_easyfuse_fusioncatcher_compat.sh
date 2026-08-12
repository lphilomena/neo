#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CONDA_CACHE="${EASYFUSE_NXF_CONDA_CACHEDIR:-${NXF_CONDA_CACHEDIR:-${ROOT}/work/.nextflow_conda}}"
STAR252="${NEOAG_FUSIONCATCHER_STAR252:-${ROOT}/../open-neo-deploy/env_tool/conda_pkgs/star-2.5.2b-0/bin/STAR}"

echo "==> patch_easyfuse_fusioncatcher_compat $(date -Is)"

if [[ ! -x "${STAR252}" ]]; then
  echo "WARN: FusionCatcher STAR 2.5.2b not found: ${STAR252}" >&2
else
  star_version="$("${STAR252}" --version 2>/dev/null | head -1 || true)"
  if [[ "${star_version}" != "STAR_2.5.2b" ]]; then
    echo "WARN: FusionCatcher STAR candidate has unexpected version: ${STAR252} -> ${star_version}" >&2
  fi
fi

patch_configuration_cfg() {
  local cfg="$1"
  [[ -f "${cfg}" ]] || return 0
  if [[ -x "${STAR252}" ]]; then
    sed -i "s|^star *=.*|star =  $(dirname "${STAR252}")/|" "${cfg}"
  fi
}

patch_fusioncatcher_py() {
  local py="$1"
  [[ -f "${py}" ]] || return 0
  if grep -q "build_version\\[1\\].*pipeline_version\\[1\\]" "${py}" \
    && ! grep -q "startswith('1.33')" "${py}"; then
    perl -0pi -e "s/build_version\\[1\\]\\.strip\\(\\) == pipeline_version\\[1\\]\\.strip\\(\\) or old_build_version\\[1\\]\\.strip\\(\\) == build_version\\[1\\]\\.strip\\(\\)/build_version[1].strip() == pipeline_version[1].strip() or old_build_version[1].strip() == build_version[1].strip() or build_version[1].strip().startswith('1.33')/" "${py}"
  fi
  sed -i "s/startswith(1\\.33)/startswith('1.33')/" "${py}"
}

patch_configuration_cfg "${ROOT}/../open-neo-deploy/env_tool/conda_pkgs/fusioncatcher-1.00-py27h8c6ebc1_1/etc/configuration.cfg"
patch_configuration_cfg "${ROOT}/../open-neo-deploy/env_tool/tools/fusioncatcher/etc/configuration.cfg"

if [[ -d "${CONDA_CACHE}" ]]; then
  while IFS= read -r cfg; do
    patch_configuration_cfg "${cfg}"
  done < <(find "${CONDA_CACHE}" -path "*/etc/configuration.cfg" -type f 2>/dev/null)

  while IFS= read -r py; do
    patch_fusioncatcher_py "${py}"
  done < <(find "${CONDA_CACHE}" -path "*/bin/fusioncatcher.py" -type f 2>/dev/null)
fi

echo "==> patch_easyfuse_fusioncatcher_compat done"
