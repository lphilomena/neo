#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CONDA_CACHE="${EASYFUSE_NXF_CONDA_CACHEDIR:-${NXF_CONDA_CACHEDIR:-${ROOT}/work/.nextflow_conda}}"
STAR252="${NEOAG_FUSIONCATCHER_STAR252:-${ROOT}/../open-neo-deploy/env_tool/conda_pkgs/star-2.5.2b-0/bin/STAR}"
FUSIONCATCHER_REF="${NEOAG_FUSIONCATCHER_REF:-${ROOT}/../open-neo-deploy/refs/data/easyfuse/easyfuse_ref_v4/fusioncatcher_index}"

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

find_executable() {
  local name="$1"
  local found=""
  if command -v "${name}" >/dev/null 2>&1; then
    command -v "${name}"
    return 0
  fi
  if [[ -d "${CONDA_CACHE}" ]]; then
    found="$(find "${CONDA_CACHE}" -path "*/bin/${name}" -type f -perm -111 2>/dev/null | head -1 || true)"
    if [[ -n "${found}" ]]; then
      echo "${found}"
      return 0
    fi
  fi
  return 1
}

ensure_bowtie1_genome_index() {
  local ref_dir="${FUSIONCATCHER_REF}"
  local genome_index="${ref_dir}/genome_index"
  local genome_index2="${ref_dir}/genome_index2/index"
  local bowtie_build=""
  local bowtie2_inspect=""
  local tmp_dir=""
  local fasta=""

  [[ -d "${ref_dir}" ]] || return 0
  if compgen -G "${genome_index}/.1.ebwt*" >/dev/null; then
    return 0
  fi
  if [[ ! -f "${genome_index2}.1.bt2" && ! -f "${genome_index2}.1.bt2l" ]]; then
    echo "WARN: FusionCatcher Bowtie2 genome index not found: ${genome_index2}" >&2
    return 0
  fi

  bowtie_build="$(find_executable bowtie-build || true)"
  bowtie2_inspect="$(find_executable bowtie2-inspect || true)"
  if [[ -z "${bowtie_build}" || -z "${bowtie2_inspect}" ]]; then
    echo "WARN: cannot build FusionCatcher Bowtie1 genome_index; bowtie-build=${bowtie_build:-missing}, bowtie2-inspect=${bowtie2_inspect:-missing}" >&2
    return 0
  fi

  echo "==> building missing FusionCatcher Bowtie1 genome_index from genome_index2"
  tmp_dir="$(mktemp -d "${ref_dir}/.genome_index_build.XXXXXX")"
  fasta="${tmp_dir}/genome.fa"
  "${bowtie2_inspect}" "${genome_index2}" > "${fasta}"
  mkdir -p "${tmp_dir}/genome_index" "${genome_index}"
  "${bowtie_build}" "${fasta}" "${tmp_dir}/genome_index/"
  find "${genome_index}" -maxdepth 1 -type f -name "*.ebwt*" -delete
  find "${tmp_dir}/genome_index" -maxdepth 1 -type f -name "*.ebwt*" -exec mv {} "${genome_index}/" \;
  rm -rf "${tmp_dir}"
}

ensure_reference_aliases() {
  local ref_dir="${FUSIONCATCHER_REF}"
  [[ -d "${ref_dir}" ]] || return 0
  if [[ ! -e "${ref_dir}/lincrnas.txt" && -f "${ref_dir}/lncrnas.txt" ]]; then
    ln -s "lncrnas.txt" "${ref_dir}/lincrnas.txt" 2>/dev/null || cp -p "${ref_dir}/lncrnas.txt" "${ref_dir}/lincrnas.txt"
  fi
}

patch_configuration_cfg "${ROOT}/../open-neo-deploy/env_tool/conda_pkgs/fusioncatcher-1.00-py27h8c6ebc1_1/etc/configuration.cfg"
patch_configuration_cfg "${ROOT}/../open-neo-deploy/env_tool/tools/fusioncatcher/etc/configuration.cfg"
ensure_bowtie1_genome_index
ensure_reference_aliases

if [[ -d "${CONDA_CACHE}" ]]; then
  while IFS= read -r cfg; do
    patch_configuration_cfg "${cfg}"
  done < <(find "${CONDA_CACHE}" -path "*/etc/configuration.cfg" -type f 2>/dev/null)

  while IFS= read -r py; do
    patch_fusioncatcher_py "${py}"
  done < <(find "${CONDA_CACHE}" -path "*/bin/fusioncatcher.py" -type f 2>/dev/null)
fi

echo "==> patch_easyfuse_fusioncatcher_compat done"
