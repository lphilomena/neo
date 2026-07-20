#!/usr/bin/env bash
# Install a pinned SpliceMutr source snapshot and its portable core runtime.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CONDA_BASE="${NEOAG_CONDA_BASE:-$(conda info --base)}"
TOOLS_ROOT="${NEOAG_TOOLS_ROOT:-${ROOT}}"
ENV_NAME="${NEOAG_SPLICEMUTR_ENV:-neoag-splicemutr}"
REF="${NEOAG_SPLICEMUTR_REF:-ac0d17005cb37810bc1e6c9a50d7707f8bd3ae66}"
PACKAGE_URL="${NEOAG_SPLICEMUTR_PACKAGE_URL:-https://gh-proxy.com/https://github.com/FertigLab/splicemutr/archive/${REF}.tar.gz}"
ARCHIVE="${NEOAG_SPLICEMUTR_ARCHIVE:-${TOOLS_ROOT}/sources/SpliceMutr-${REF}.tar.gz}"
GENOMEINFO_ARCHIVE="${NEOAG_GENOMEINFODBDATA_ARCHIVE:-${TOOLS_ROOT}/sources/GenomeInfoDbData_1.2.11.tar.gz}"
GENOMEINFO_URL="${NEOAG_GENOMEINFODBDATA_URL:-https://mghp.osn.xsede.org/bir190004-bucket01/archive.bioconductor.org/packages/3.18/data/annotation/src/contrib/GenomeInfoDbData_1.2.11.tar.gz}"
GENOMEINFO_MD5="2a4cbfc2031992fed3c9445f450890a2"
HOME_DIR="${NEOAG_SPLICEMUTR_HOME:-${TOOLS_ROOT}/tools/SpliceMutr}"
YML="${ROOT}/conda/env.neoag-splicemutr.yml"
BIN_DIR="${ROOT}/bin"
TOOLS_ENV="${ROOT}/conf/tools.env.sh"
export CONDA_CHANNEL_ALIAS="${NEOAG_CONDA_CHANNEL_ALIAS:-${CONDA_CHANNEL_ALIAS:-https://conda.anaconda.org}}"

# Bioconda annotation packages download their payloads in post-link scripts.
# Resume interrupted transfers instead of restarting large files from zero.
CURL_RETRY_HOME="$(mktemp -d)"
tmp_dir=""
cleanup() {
  rm -rf "${CURL_RETRY_HOME}"
  [[ -z "${tmp_dir}" ]] || rm -rf "${tmp_dir}"
}
trap cleanup EXIT
cat > "${CURL_RETRY_HOME}/.curlrc" <<'EOF'
retry = 10
retry-all-errors
retry-delay = 2
connect-timeout = 30
continue-at = -
EOF
export CURL_HOME="${CURL_RETRY_HOME}"

# Bioconda's post-link helper streams this payload through stdout and cannot
# resume safely. Cache and verify it first, then serve it locally to that helper.
mkdir -p "$(dirname "${GENOMEINFO_ARCHIVE}")"
if [[ ! -s "${GENOMEINFO_ARCHIVE}" ]] || ! echo "${GENOMEINFO_MD5}  ${GENOMEINFO_ARCHIVE}" | md5sum -c - >/dev/null 2>&1; then
  curl -fL --retry 10 --retry-all-errors -C - -o "${GENOMEINFO_ARCHIVE}" "${GENOMEINFO_URL}"
fi
echo "${GENOMEINFO_MD5}  ${GENOMEINFO_ARCHIVE}" | md5sum -c -
export NEOAG_GENOMEINFODBDATA_ARCHIVE="${GENOMEINFO_ARCHIVE}"
cat > "${CURL_RETRY_HOME}/bash_env" <<'EOF'
curl() {
  local arg
  for arg in "$@"; do
    if [[ "${arg}" == *GenomeInfoDbData_1.2.11.tar.gz ]]; then
      cat "${NEOAG_GENOMEINFODBDATA_ARCHIVE}"
      return
    fi
  done
  command curl "$@"
}
EOF
export BASH_ENV="${CURL_RETRY_HOME}/bash_env"

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

env_exists() { conda_safe env list | awk '{print $1}' | grep -qx "$1"; }
env_is_ready() {
  env_exists "$1" && conda_safe run -n "$1" python -c 'import Bio, numpy, pandas, sklearn' >/dev/null 2>&1 \
    && conda_safe run -n "$1" Rscript -e 'library(BSgenome); library(GenomicFeatures); library(optparse)' >/dev/null 2>&1 \
    && conda_safe run -n "$1" snakemake --version >/dev/null 2>&1
}

if [[ "${NEOAG_FORCE_ENV_UPDATE:-0}" == "1" ]]; then
  if env_exists "${ENV_NAME}"; then
    "${CONDA_RUNNER}" env update -n "${ENV_NAME}" -f "${YML}" --prune
  else
    "${CONDA_RUNNER}" env create -n "${ENV_NAME}" -f "${YML}" -y
  fi
elif env_is_ready "${ENV_NAME}"; then
  echo "==> ${ENV_NAME} core runtime is already ready"
elif env_exists "${ENV_NAME}"; then
  "${CONDA_RUNNER}" env update -n "${ENV_NAME}" -f "${YML}"
else
  "${CONDA_RUNNER}" env create -n "${ENV_NAME}" -f "${YML}" -y
fi

mkdir -p "$(dirname "${ARCHIVE}")" "$(dirname "${HOME_DIR}")" "${BIN_DIR}" "${ROOT}/conf"
if [[ ! -s "${ARCHIVE}" ]]; then
  command -v curl >/dev/null 2>&1 || { echo "ERROR: curl is required to download SpliceMutr" >&2; exit 43; }
  curl -fL --retry 3 --connect-timeout 20 -o "${ARCHIVE}.part" "${PACKAGE_URL}"
  mv "${ARCHIVE}.part" "${ARCHIVE}"
else
  echo "==> Reusing cached SpliceMutr snapshot ${ARCHIVE}"
fi

tmp_dir="$(mktemp -d)"
tar -xzf "${ARCHIVE}" -C "${tmp_dir}"
source_dir="$(find "${tmp_dir}" -mindepth 1 -maxdepth 1 -type d | head -1)"
[[ -n "${source_dir}" ]] || { echo "ERROR: SpliceMutr archive has no source directory" >&2; exit 44; }
rm -rf "${HOME_DIR}.new"
mv "${source_dir}" "${HOME_DIR}.new"
printf '%s\n' "${REF}" > "${HOME_DIR}.new/NEOAG_PINNED_REF"
rm -rf "${HOME_DIR}"
mv "${HOME_DIR}.new" "${HOME_DIR}"

cat > "${BIN_DIR}/splicemutr-neoag" <<EOF
#!/usr/bin/env bash
set -euo pipefail
CONDA="${CONDA_BASE}/bin/conda"
ENV_NAME="${ENV_NAME}"
HOME_DIR="${HOME_DIR}"
case "\${1:-doctor}" in
  -h|--help)
    echo "Usage: splicemutr-neoag doctor | workflow WORKFLOW [snakemake args] | r SCRIPT [args] | python SCRIPT [args]"
    echo "WORKFLOW may be a .smk path or one of: prep, leafcutter, run, genotype"
    ;;
  doctor)
    test -f "\${HOME_DIR}/NEOAG_PINNED_REF"
    "\${CONDA}" run -n "\${ENV_NAME}" python -c 'import Bio, numpy, pandas, sklearn; print("SpliceMutr Python runtime OK")'
    "\${CONDA}" run -n "\${ENV_NAME}" Rscript -e 'suppressPackageStartupMessages({library(BSgenome); library(GenomicFeatures); library(optparse)}); cat("SpliceMutr R runtime OK\\n")'
    "\${CONDA}" run -n "\${ENV_NAME}" snakemake --version
    ;;
  workflow)
    shift
    workflow="\${1:?workflow name or .smk path is required}"; shift
    case "\${workflow}" in
      prep) workflow="\${HOME_DIR}/simulation/prep_references/prep_ref.smk" ;;
      leafcutter) workflow="\${HOME_DIR}/simulation/running_leafcutter/running_leafcutter.smk" ;;
      run) workflow="\${HOME_DIR}/simulation/running_splicemutr/run_splicemutr.smk" ;;
      genotype) workflow="\${HOME_DIR}/simulation/genotyping_samples/genotype_samples.smk" ;;
    esac
    exec "\${CONDA}" run -n "\${ENV_NAME}" snakemake -s "\${workflow}" "\$@"
    ;;
  r)
    shift; script="\${1:?R script is required}"; shift
    [[ "\${script}" = /* ]] || script="\${HOME_DIR}/Rscripts/\${script}"
    exec "\${CONDA}" run -n "\${ENV_NAME}" Rscript "\${script}" "\$@"
    ;;
  python)
    shift; script="\${1:?Python script is required}"; shift
    [[ "\${script}" = /* ]] || script="\${HOME_DIR}/python_scripts/\${script}"
    exec "\${CONDA}" run -n "\${ENV_NAME}" python "\${script}" "\$@"
    ;;
  *) echo "ERROR: unknown command: \$1" >&2; exit 2 ;;
esac
EOF
chmod +x "${BIN_DIR}/splicemutr-neoag"

for assignment in \
  "NEOAG_SPLICEMUTR_ENV=${ENV_NAME}" \
  "NEOAG_SPLICEMUTR_HOME=${HOME_DIR}" \
  "NEOAG_SPLICEMUTR_BIN=${BIN_DIR}/splicemutr-neoag"; do
  key="${assignment%%=*}"; value="${assignment#*=}"
  if grep -q "^export ${key}=" "${TOOLS_ENV}" 2>/dev/null; then
    sed -i "s|^export ${key}=.*|export ${key}=\"${value}\"|" "${TOOLS_ENV}"
  else
    printf 'export %s="%s"\n' "${key}" "${value}" >> "${TOOLS_ENV}"
  fi
done

echo "==> SpliceMutr smoke"
"${BIN_DIR}/splicemutr-neoag" doctor
echo "==> SpliceMutr installed at ${HOME_DIR} (${REF})"
echo "NOTE: a genome-specific BSgenome package/2bit, GTF/TxDb, STAR index, and cohort junction inputs are runtime assets."
