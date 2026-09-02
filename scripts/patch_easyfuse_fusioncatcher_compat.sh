#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CONDA_CACHE="${EASYFUSE_NXF_CONDA_CACHEDIR:-${NXF_CONDA_CACHEDIR:-${ROOT}/work/.nextflow_conda}}"
STAR_OVERRIDE="${NEOAG_FUSIONCATCHER_STAR:-${NEOAG_FUSIONCATCHER_STAR252:-}}"
FUSIONCATCHER_REF="${NEOAG_FUSIONCATCHER_REF:-${ROOT}/../open-neo-deploy/refs/data/easyfuse/easyfuse_ref_v4/fusioncatcher_index}"

echo "==> patch_easyfuse_fusioncatcher_compat $(date -Is)"

star_version() {
  local star="$1"
  "${star}" --version 2>/dev/null | head -1 | sed 's/^STAR_//'
}

fusioncatcher_required_star() {
  local cfg="$1" prefix py
  prefix="$(cd "$(dirname "${cfg}")/.." && pwd -P)"
  py="${prefix}/bin/fusioncatcher.py"
  [[ -f "${py}" ]] || return 1
  sed -nE "s/^[[:space:]]*correct_version[[:space:]]*=[[:space:]]*['\"]([^'\"]+)['\"].*/\1/p" "${py}" | head -1
}

find_star_for_version() {
  local required="$1" candidate version
  local -a candidates=()
  [[ -n "${STAR_OVERRIDE}" ]] && candidates+=("${STAR_OVERRIDE}")
  for candidate in \
    "${NEOAG_CONDA_BASE:-${HOME}/miniforge3}"/pkgs/star-"${required}"-*/bin/STAR \
    "${ROOT}/../open-neo-deploy/env_tool/conda_pkgs"/star-"${required}"-*/bin/STAR \
    "${CONDA_CACHE}"/env-*/bin/STAR; do
    [[ -x "${candidate}" ]] && candidates+=("${candidate}")
  done
  if command -v STAR >/dev/null 2>&1; then
    candidates+=("$(command -v STAR)")
  fi
  for candidate in "${candidates[@]}"; do
    version="$(star_version "${candidate}" || true)"
    if [[ "${version}" == "${required}" ]]; then
      printf '%s\n' "${candidate}"
      return 0
    fi
  done
  return 1
}

patch_configuration_cfg() {
  local cfg="$1" required star prefix launcher star_dir env_star
  [[ -f "${cfg}" ]] || return 0
  required="$(fusioncatcher_required_star "${cfg}" || true)"
  if [[ -z "${required}" ]]; then
    echo "WARN: cannot determine FusionCatcher STAR requirement for ${cfg}" >&2
    return 0
  fi
  star="$(find_star_for_version "${required}" || true)"
  if [[ -z "${star}" ]]; then
    echo "WARN: FusionCatcher requires STAR ${required}, but no matching binary was found for ${cfg}" >&2
    return 0
  fi
  sed -i "s|^star *=.*|star =  $(dirname "${star}")/|" "${cfg}"
  prefix="$(cd "$(dirname "${cfg}")/.." && pwd -P)"
  launcher="${prefix}/bin/fusioncatcher"
  star_dir="$(dirname "${star}")"
  if [[ -f "${launcher}" ]] && grep -q '^export PATH=\$fbin:\$PATH' "${launcher}"; then
    sed -i "s|^export PATH=.*|export PATH=${star_dir}:\$fbin:\$PATH|" "${launcher}"
  fi
  env_star="${prefix}/bin/STAR"
  if [[ "$(star_version "${env_star}" || true)" != "${required}" ]]; then
    if [[ -e "${env_star}" && ! -e "${env_star}.openneo-original" ]]; then
      mv "${env_star}" "${env_star}.openneo-original"
    else
      rm -f "${env_star}"
    fi
    ln -s "${star}" "${env_star}"
  fi
  echo "==> FusionCatcher STAR ${required}: ${star}"
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
  local required_index_files=(.1.ebwt .2.ebwt .3.ebwt .4.ebwt .rev.1.ebwt .rev.2.ebwt)
  local bowtie_build=""
  local bowtie2_inspect=""
  local tmp_dir=""
  local fasta=""
  local complete="yes"

  [[ -d "${ref_dir}" ]] || return 0
  for idx_name in "${required_index_files[@]}"; do
    if [[ ! -s "${genome_index}/${idx_name}" ]]; then
      complete="no"
      break
    fi
  done
  if [[ "${complete}" == "yes" ]]; then
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
  local optional_filters=(
    antisenses.txt
    rp11.txt
    cta.txt
    ctb.txt
    ctd.txt
    ctc.txt
    rp.txt
    celllines.txt
  )
  [[ -d "${ref_dir}" ]] || return 0
  if [[ ! -e "${ref_dir}/lincrnas.txt" && -f "${ref_dir}/lncrnas.txt" ]]; then
    ln -s "lncrnas.txt" "${ref_dir}/lincrnas.txt" 2>/dev/null || cp -p "${ref_dir}/lncrnas.txt" "${ref_dir}/lincrnas.txt"
  fi
  if [[ ! -e "${ref_dir}/prostates.txt" && -f "${ref_dir}/prostate_cancer.txt" ]]; then
    ln -s "prostate_cancer.txt" "${ref_dir}/prostates.txt" 2>/dev/null || cp -p "${ref_dir}/prostate_cancer.txt" "${ref_dir}/prostates.txt"
  fi
  if [[ ! -e "${ref_dir}/chimerdb3kb.txt" && -f "${ref_dir}/chimerdb4kb.txt" ]]; then
    ln -s "chimerdb4kb.txt" "${ref_dir}/chimerdb3kb.txt" 2>/dev/null || cp -p "${ref_dir}/chimerdb4kb.txt" "${ref_dir}/chimerdb3kb.txt"
  fi
  if [[ ! -e "${ref_dir}/chimerdb3pub.txt" && -f "${ref_dir}/chimerdb4pub.txt" ]]; then
    ln -s "chimerdb4pub.txt" "${ref_dir}/chimerdb3pub.txt" 2>/dev/null || cp -p "${ref_dir}/chimerdb4pub.txt" "${ref_dir}/chimerdb3pub.txt"
  fi
  if [[ ! -e "${ref_dir}/chimerdb3seq.txt" && -f "${ref_dir}/chimerdb4seq.txt" ]]; then
    ln -s "chimerdb4seq.txt" "${ref_dir}/chimerdb3seq.txt" 2>/dev/null || cp -p "${ref_dir}/chimerdb4seq.txt" "${ref_dir}/chimerdb3seq.txt"
  fi
  for filter_name in "${optional_filters[@]}"; do
    if [[ ! -e "${ref_dir}/${filter_name}" ]]; then
      : > "${ref_dir}/${filter_name}"
    fi
  done
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
