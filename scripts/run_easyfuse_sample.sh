#!/usr/bin/env bash
# Generic EasyFuse runner for one paired-end RNA-seq sample.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
# shellcheck source=/dev/null
source "${ROOT}/conf/tools.env.sh"
: "${NEOAG_CONDA_BASE:?ERROR: set NEOAG_CONDA_BASE to your conda/mamba installation root}"
if [[ ! -f "${NEOAG_EASYFUSE_HOME:-}/main.nf" ]]; then
  for easyfuse_candidate in \
    "${ROOT}/tools/EasyFuse" \
    "${ROOT}/../open-neo-deploy/env_tool/tools/EasyFuse" \
    "${ROOT}/../neoantigen/neoag_event_pipeline_v03_rc/tools/EasyFuse"; do
    if [[ -f "${easyfuse_candidate}/main.nf" ]]; then
      export NEOAG_EASYFUSE_HOME="${easyfuse_candidate}"
      break
    fi
  done
  unset easyfuse_candidate
fi
[[ -f "${NEOAG_EASYFUSE_HOME:-}/main.nf" ]] || {
  echo "ERROR: EasyFuse main.nf not found; set NEOAG_EASYFUSE_HOME to the deployed EasyFuse directory" >&2
  exit 1
}
if [[ -f "${NEOAG_EASYFUSE_HOME}/modules/arriba/environment.yml" \
  && -f "${NEOAG_EASYFUSE_HOME}/modules/starfusion/starfusion/environment.yml" \
  && -f "${NEOAG_EASYFUSE_HOME}/modules/fusioncatcher/environment.yml" ]]; then
  EASYFUSE_LAYOUT="module-native"
else
  EASYFUSE_LAYOUT="legacy"
fi
export PATH="${NEOAG_CONDA_BASE}/bin:${PATH}"
# EasyFuse Nextflow activates starfusion.yml; keep repo STAR-Fusion off PATH.
export PATH="$(echo "${PATH}" | tr ':' '\n' | grep -vE "${NEOAG_STAR_FUSION_HOME}\$" | paste -sd: -)"

SAMPLE_ID="${EASYFUSE_SAMPLE_ID:-${SAMPLE_ID:-sample}}"
FQ1="${EASYFUSE_FQ1:?ERROR: set EASYFUSE_FQ1=/path/sample_R1.fq.gz}"
FQ2="${EASYFUSE_FQ2:?ERROR: set EASYFUSE_FQ2=/path/sample_R2.fq.gz}"
REF="${NEOAG_EASYFUSE_REF:?ERROR: set NEOAG_EASYFUSE_REF=/path/to/easyfuse_ref_v4}"
STAR_INDEX="${EASYFUSE_STAR_INDEX:-${REF}/starfusion_index/ref_genome.fa.star.idx}"
STARFUSION_INDEX="${EASYFUSE_STARFUSION_INDEX:-${REF}/starfusion_index}"
OUT="${OUTDIR:-${ROOT}/results/easyfuse}"
LOG="${LOG:-${ROOT}/work/run_easyfuse_${SAMPLE_ID}.log}"
RUNTIME_DIR="${EASYFUSE_RUNTIME_DIR:-${ROOT}/work}"
INPUT="${EASYFUSE_INPUT_TSV:-${RUNTIME_DIR}/easyfuse_${SAMPLE_ID}_input.tsv}"
PREBUILD_LOG="${RUNTIME_DIR}/easyfuse_conda_prebuild.log"
PREBUILD_PID_FILE="${RUNTIME_DIR}/easyfuse_conda_prebuild.pid"
NXF_RUN_NAME="${EASYFUSE_RUN_NAME:-easyfuse_${SAMPLE_ID}}"
NXF_STEM="${NXF_RUN_NAME//[^A-Za-z0-9_.-]/_}"

ensure_input_tsv() {
  if [[ "${EASYFUSE_LAYOUT}" == "module-native" ]]; then
    printf 'sample\tfastq_1\tfastq_2\n%s\t%s\t%s\n' "${SAMPLE_ID}" "${FQ1}" "${FQ2}" > "${INPUT}"
  else
    printf '%s\t%s\t%s\n' "${SAMPLE_ID}" "${FQ1}" "${FQ2}" > "${INPUT}"
  fi
  if ! awk -F'\t' 'NF==3 {found=1} END{exit !found}' "${INPUT}"; then
    echo "ERROR: input TSV must have 3 tab-separated columns: ${INPUT}" >&2
    exit 1
  fi
}

STAR_TMP="${EASYFUSE_TMPDIR:-${RUNTIME_DIR}/star_tmp_${NXF_STEM}}"
export NXF_HOME="${EASYFUSE_NXF_HOME:-${RUNTIME_DIR}/.nextflow_home_${NXF_STEM}}"
export NXF_WORK="${EASYFUSE_NXF_WORK:-${RUNTIME_DIR}/.nextflow_work_${NXF_STEM}}"
mkdir -p "${OUT}" "$(dirname "${LOG}")" "${NXF_HOME}" "${NXF_WORK}" "${STAR_TMP}"
export TMPDIR="${STAR_TMP}"
ensure_input_tsv

# Each run receives its own Nextflow and STAR directories below. Do not kill
# processes by a global EasyFuse command pattern: those may belong to another
# approved case running on the same machine.

export NXF_DISABLE_CHECK_TTY=true
export CONDA_ALWAYS_YES=true
export MAMBA_ALWAYS_YES=true
export NEOAG_REAL_MAMBA="${NEOAG_CONDA_BASE}/bin/mamba"
export NEOAG_BIOC_CACHE_HELPER="${ROOT}/.agents/skills/neoag-remote-deploy/scripts/with_bioc_data_cache.sh"
export NEOAG_INSTALL_CACHE_ROOT="${NEOAG_INSTALL_CACHE_ROOT:-${RUNTIME_DIR}/install_cache}"
export NEOAG_EASYFUSE_BIOC_PACKAGE_KEY="${NEOAG_EASYFUSE_BIOC_PACKAGE_KEY:-genomeinfodbdata-1.2.11}"
mkdir -p "${ROOT}/work/easyfuse_bin"
cat > "${ROOT}/work/easyfuse_bin/mamba" <<'MAMBA_WRAPPER'
#!/usr/bin/env bash
set -euo pipefail
real="${NEOAG_REAL_MAMBA:-}"
if [[ -z "$real" || ! -x "$real" ]]; then
  real="$(command -v mamba.real || true)"
fi
if [[ -z "$real" || ! -x "$real" ]]; then
  for candidate in "${NEOAG_CONDA_BASE:-${HOME}/miniforge3}/bin/mamba" "${HOME}/miniforge3/bin/mamba"; do
    [[ -x "$candidate" && "$candidate" != "$0" ]] && real="$candidate" && break
  done
fi
[[ -n "$real" && -x "$real" ]] || { echo "ERROR: real mamba not found" >&2; exit 127; }
args=("$@")
case " $* " in
  *" -y "*|*" --yes "*) ;;
  *) args=(--yes "${args[@]}") ;;
esac
env_file=""
for ((i=0; i<${#args[@]}; i++)); do
  if [[ "${args[$i]}" == "--file" || "${args[$i]}" == "-f" ]]; then
    env_file="${args[$((i + 1))]:-}"
    break
  fi
done
if [[ -n "$env_file" && -f "$env_file" ]] \
  && grep -Eiq '(^|[=[:space:]-])(arriba|bioconductor-genomeinfodbdata)([=[:space:]]|$)' "$env_file" \
  && [[ -x "${NEOAG_BIOC_CACHE_HELPER:-}" ]]; then
  exec "${NEOAG_BIOC_CACHE_HELPER}" \
    --conda-base "${NEOAG_CONDA_BASE}" \
    --cache-root "${NEOAG_INSTALL_CACHE_ROOT}" \
    --package-key "${NEOAG_EASYFUSE_BIOC_PACKAGE_KEY}" \
    -- "$real" "${args[@]}"
fi
exec "$real" "${args[@]}"
MAMBA_WRAPPER
chmod +x "${ROOT}/work/easyfuse_bin/mamba"
resolve_easyfuse_java_home() {
  local candidate java_bin
  for candidate in \
    "${NEOAG_EASYFUSE_JAVA_HOME:-}" \
    "${JAVA_HOME:-}" \
    "${NEOAG_EASYFUSE_ENV_PREFIX:-}" \
    "${NEOAG_CONDA_BASE}/envs/${NEOAG_FUSION_ENV}" \
    "${NEOAG_CONDA_BASE}/envs/${NEOAG_GATK_ENV:-neoag-gatk}" \
    "${NEOAG_CONDA_BASE}/envs/neoag-runtime"; do
    if [[ -n "${candidate}" && -x "${candidate}/bin/java" ]]; then
      printf '%s\n' "${candidate}"
      return 0
    fi
  done
  java_bin="$(command -v java 2>/dev/null || true)"
  [[ -n "${java_bin}" ]] || return 1
  cd "$(dirname "${java_bin}")/.." && pwd -P
}
export JAVA_HOME="$(resolve_easyfuse_java_home)"
export PATH="${ROOT}/work/easyfuse_bin:${NEOAG_CONDA_BASE}/bin:${JAVA_HOME}/bin:${PATH}"
# EasyFuse 2.x uses legacy Groovy config variables. Nextflow 26 defaults to
# the strict v2 parser, so retain v1 unless a deployment explicitly overrides it.
export NXF_SYNTAX_PARSER="${NXF_SYNTAX_PARSER:-v1}"
resolve_nextflow() {
  local candidate
  for candidate in \
    "${NEOAG_NEXTFLOW:-}" \
    "${ROOT}/bin/nextflow" \
    "${JAVA_HOME}/bin/nextflow" \
    "${NEOAG_CONDA_BASE}/bin/nextflow"; do
    if [[ -n "${candidate}" && -x "${candidate}" ]]; then
      printf '%s\n' "${candidate}"
      return 0
    fi
  done
  if command -v nextflow >/dev/null 2>&1; then
    command -v nextflow
    return 0
  fi
  echo "ERROR: nextflow executable was not found" >&2
  return 1
}
export NEOAG_NEXTFLOW="$(resolve_nextflow)"

exec > >(tee -a "${LOG}") 2>&1
echo "==> run_easyfuse_sample $(date -Is)"
echo "    sample=${SAMPLE_ID}"
echo "    fq1=${FQ1}"
echo "    fq2=${FQ2}"
echo "    input=${INPUT}"
echo "    reference=${REF}"
echo "    star_index=${STAR_INDEX}"
echo "    starfusion_index=${STARFUSION_INDEX}"
echo "    output=${OUT}"
echo "    easyfuse_layout=${EASYFUSE_LAYOUT}"
echo "    nxf_home=${NXF_HOME}"
echo "    nxf_work=${NXF_WORK}"

[[ -f "${FQ1}" && -f "${FQ2}" ]] || {
  echo "ERROR: FASTQ not found:" >&2
  echo "  ${FQ1}" >&2
  echo "  ${FQ2}" >&2
  exit 1
}

[[ -f "${REF}/BEFORE_EXECUTING_EASYFUSE" ]] || {
  echo "ERROR: EasyFuse reference missing: ${REF}" >&2
  exit 1
}

[[ -f "${STAR_INDEX}/genomeParameters.txt" && -f "${STAR_INDEX}/Genome" ]] || {
  echo "ERROR: EasyFuse STAR index is incomplete: ${STAR_INDEX}" >&2
  echo "       expected genomeParameters.txt and Genome" >&2
  exit 1
}

[[ -d "${STARFUSION_INDEX}" ]] || {
  echo "ERROR: EasyFuse STAR-Fusion index missing: ${STARFUSION_INDEX}" >&2
  exit 1
}

CONDA_CACHE="${EASYFUSE_NXF_CONDA_CACHEDIR:-${RUNTIME_DIR}/.nextflow_conda}"
mkdir -p "${CONDA_CACHE}"
export NXF_CONDA_CACHEDIR="${EASYFUSE_NXF_CONDA_CACHEDIR:-${CONDA_CACHE}}"
echo "    nxf_conda_cachedir=${NXF_CONDA_CACHEDIR}"

ensure_easyfuse_env_compat_files() {
  local env_dir="${NEOAG_EASYFUSE_HOME}/environments"
  [[ -d "${env_dir}" ]] || { echo "ERROR: EasyFuse environments dir missing: ${env_dir}" >&2; exit 3; }

  if [[ ! -f "${env_dir}/easyfuse_src.yml" ]]; then
    echo "==> Creating EasyFuse v1 compatibility env: ${env_dir}/easyfuse_src.yml"
    cat > "${env_dir}/easyfuse_src.yml" <<'YAML'
name: easyfuse_src
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - bioconda::pyeasyfuse=2.0.3
  - bioconda::skewer
YAML
  fi

  if [[ ! -f "${env_dir}/requantification_wo_easyfuse.yml" ]]; then
    echo "==> Creating EasyFuse v1 compatibility env: ${env_dir}/requantification_wo_easyfuse.yml"
    if [[ -f "${env_dir}/alignment.yml" ]]; then
      cp "${env_dir}/alignment.yml" "${env_dir}/requantification_wo_easyfuse.yml"
    else
      cat > "${env_dir}/requantification_wo_easyfuse.yml" <<'YAML'
name: requantification_wo_easyfuse
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - bioconda::star=2.6.1d
  - bioconda::samtools=1.9.0
YAML
    fi
  fi
}


ensure_easyfuse_entrypoints() {
  local source=""
  local candidate
  for candidate in "${CONDA_CACHE}"/env-*/bin/easy-fuse; do
    [[ -x "${candidate}" ]] || continue
    if "${candidate}" --help >/dev/null 2>&1; then
      source="${candidate}"
      break
    fi
  done
  [[ -n "${source}" ]] || {
    echo "ERROR: no working easy-fuse entrypoint found in ${CONDA_CACHE}" >&2
    exit 1
  }

  local prefix target
  for prefix in "${CONDA_CACHE}"/env-*; do
    [[ -d "${prefix}/bin" ]] || continue
    target="${prefix}/bin/easy-fuse"
    [[ -e "${target}" ]] && continue
    [[ -x "${prefix}/bin/STAR" || -x "${prefix}/bin/bowtie-build" ]] || continue
    cat > "${target}" <<EOF
#!/usr/bin/env bash
exec ${source} "\$@"
EOF
    chmod +x "${target}"
    echo "    installed easy-fuse shim ${target} -> ${source}"
  done
}

wait_for_mamba_free() {
  while pgrep -f 'mamba env create' >/dev/null 2>&1; do
    sleep 15
  done
  rm -f "${NEOAG_CONDA_BASE}/pkgs/pkgs.lock" 2>/dev/null || true
}

select_easyfuse_env_yml() {
  local preferred="$1"
  shift
  if [[ -f "${NEOAG_EASYFUSE_HOME}/environments/${preferred}" ]]; then
    printf '%s' "$preferred"
    return 0
  fi
  local alt
  for alt in "$@"; do
    if [[ -f "${NEOAG_EASYFUSE_HOME}/environments/${alt}" ]]; then
      printf '%s' "$alt"
      return 0
    fi
  done
  echo "ERROR: none of the EasyFuse env files exist under ${NEOAG_EASYFUSE_HOME}/environments: ${preferred} $*" >&2
  exit 3
}

if [[ "${EASYFUSE_LAYOUT}" == "legacy" ]]; then
ensure_easyfuse_env_compat_files

prebuild_conda_env() {
  local env_id="$1"
  local yml="$2"
  local check_bin="$3"
  local prefix="${CONDA_CACHE}/env-${env_id}"

  if [[ -x "${prefix}/bin/${check_bin}" ]]; then
    bash "${ROOT}/scripts/fix_easyfuse_pyeasyfuse_env.sh" >/dev/null 2>&1 || true
    return 0
  fi

  wait_for_mamba_free
  echo "==> Pre-building EasyFuse ${yml} ..."
  rm -rf "${prefix}"
  mamba env create -y \
    --prefix "${prefix}" \
    --file "${NEOAG_EASYFUSE_HOME}/environments/${yml}"
  bash "${ROOT}/scripts/fix_easyfuse_pyeasyfuse_env.sh"
}

QC_ENV="${CONDA_CACHE}/env-574d468f667e5ead-1f348f31c1e78ea89e97e435a63f0c7d"
SRC_ENV="${CONDA_CACHE}/env-adab1ef12c1f56bf-14649bb80e8151aa81731d54781c13cc"

if [[ ! -x "${QC_ENV}/bin/fastp" || ! -x "${SRC_ENV}/bin/skewer" ]]; then
  wait_for_mamba_free
fi

prebuild_conda_env \
  "574d468f667e5ead-1f348f31c1e78ea89e97e435a63f0c7d" \
  "qc.yml" \
  "fastp"

if [[ -f "${NEOAG_EASYFUSE_HOME}/environments/easyfuse_src.yml" ]]; then
  prebuild_conda_env \
    "adab1ef12c1f56bf-14649bb80e8151aa81731d54781c13cc" \
    "easyfuse_src.yml" \
    "skewer"
else
  echo "==> EasyFuse easyfuse_src.yml not present; using v2 environment layout"
fi

prebuild_worker_running() {
  [[ -s "${PREBUILD_PID_FILE}" ]] || return 1
  local pid
  pid="$(cat "${PREBUILD_PID_FILE}")"
  [[ "${pid}" =~ ^[0-9]+$ ]] && kill -0 "${pid}" 2>/dev/null
}

start_prebuild_worker() {
  if prebuild_worker_running; then
    echo "==> background conda prebuild worker already running (PID=$(cat "${PREBUILD_PID_FILE}"))"
    return 0
  fi
  rm -f "${PREBUILD_PID_FILE}"
  (
    set +e
    EASYFUSE_RUNTIME_DIR="${RUNTIME_DIR}" \
    EASYFUSE_NXF_CONDA_CACHEDIR="${CONDA_CACHE}" \
    EASYFUSE_PREBUILD_LOG="${PREBUILD_LOG}" \
      bash "${ROOT}/scripts/easyfuse_prebuild_remaining_envs.sh"
    rc=$?
    rm -f "${PREBUILD_PID_FILE}"
    exit "${rc}"
  ) >>"${PREBUILD_LOG}" 2>&1 &
  local pid=$!
  printf '%s\n' "${pid}" > "${PREBUILD_PID_FILE}"
  echo "==> background conda prebuild worker PID=${pid} (log: ${PREBUILD_LOG})"
}

start_prebuild_worker

ALIGN_ENV="${CONDA_CACHE}/env-6f2b394c864eeaa5-8f88fe4572f59d9bb818f7644ca8f1fa"
echo "==> waiting for alignment env (STAR) before Nextflow ..."
while [[ ! -x "${ALIGN_ENV}/bin/STAR" ]]; do
  if ! prebuild_worker_running; then
    echo "ERROR: alignment env build failed; see ${PREBUILD_LOG}" >&2
    exit 1
  fi
  sleep 20
done
echo "==> alignment env ready"

bash "${ROOT}/scripts/fix_easyfuse_pyeasyfuse_env.sh"
bash "${ROOT}/scripts/seed_easyfuse_conda_envs.sh"
ensure_easyfuse_entrypoints

REQ_WO_ENV="${CONDA_CACHE}/env-requantification_wo_easyfuse"
if [[ ! -x "${REQ_WO_ENV}/bin/STAR" ]]; then
  wait_for_mamba_free
  echo "==> Pre-building requantification_wo_easyfuse.yml ..."
  rm -rf "${REQ_WO_ENV}"
  mamba env create -y \
    --prefix "${REQ_WO_ENV}" \
    --file "${NEOAG_EASYFUSE_HOME}/environments/$(select_easyfuse_env_yml requantification_wo_easyfuse.yml requantification.yml)"
fi

bash "${ROOT}/scripts/patch_easyfuse_star_avx2.sh"
bash "${ROOT}/scripts/patch_easyfuse_requant_star_index_cleanup.sh"
ensure_easyfuse_entrypoints
bash "${ROOT}/scripts/patch_easyfuse_fusioncatcher_compat.sh"
bash "${ROOT}/scripts/fix_easyfuse_pyeasyfuse_env.sh"
else
  echo "==> EasyFuse module-native layout: internal caller environments are managed by Nextflow"
fi

# FusionCatcher validates an exact historical STAR version even in the
# module-native EasyFuse layout. Patch every discovered Nextflow environment
# after it exists, regardless of which EasyFuse directory layout is in use.
bash "${ROOT}/scripts/patch_easyfuse_fusioncatcher_compat.sh"

export PATH="$(echo "${PATH}" | tr ':' '\n' | grep -vE '/envs/neoag-tools/bin$|/tools/fusioncatcher/bin$' | paste -sd: -)"

cd "${ROOT}/work"

# Separate run name avoids resume/lock collision with easyfuse_cfrna_test (session 9d03387c...).
NXF_HISTORY="${NXF_HOME:-${ROOT}/work/.nextflow_home}/history"
if [[ ! -f "${NXF_HISTORY}" ]]; then
  NXF_HISTORY="${ROOT}/work/.nextflow/history"
fi
if grep -qF "${NXF_RUN_NAME}" "${NXF_HISTORY}" 2>/dev/null; then
  NXF_RESUME_ARGS=(-resume "${NXF_RUN_NAME}")
else
  NXF_RESUME_ARGS=(-name "${NXF_RUN_NAME}")
fi

run_nextflow() {
  if [[ "${EASYFUSE_LAYOUT}" == "module-native" ]]; then
    "${NEOAG_NEXTFLOW}" run "${NEOAG_EASYFUSE_HOME}/main.nf" \
      "${NXF_RESUME_ARGS[@]}" \
      -c "${ROOT}/conf/easyfuse.nextflow.config" \
      -profile conda \
      -w "${NXF_WORK}" \
      --output "${OUT}" \
      --input_files "${INPUT}" \
      --reference "${REF}" </dev/null
  else
    "${NEOAG_NEXTFLOW}" run "${NEOAG_EASYFUSE_HOME}/main.nf" \
      "${NXF_RESUME_ARGS[@]}" \
      -c "${ROOT}/conf/easyfuse.nextflow.config" \
      -profile conda \
      -w "${NXF_WORK}" \
      --output "${OUT}" \
      --input_files "${INPUT}" \
      --reference "${REF}" \
      --star_index "${STAR_INDEX}" \
      --starfusion_index "${STARFUSION_INDEX}" \
      --annotation_db "${REF}/Homo_sapiens.GRCh38.110.gff3.db" \
      --reference_tsl "${REF}/Homo_sapiens.GRCh38.110.gtf.tsl" </dev/null
  fi
}

run_nextflow || {
  echo "==> Nextflow failed; retrying once with -resume ..."
  if [[ "${EASYFUSE_LAYOUT}" == "legacy" ]]; then
    bash "${ROOT}/scripts/patch_easyfuse_star_avx2.sh"
  fi
  bash "${ROOT}/scripts/patch_easyfuse_fusioncatcher_compat.sh"
  NXF_RESUME_ARGS=(-resume "${NXF_RUN_NAME}")
  run_nextflow
}

[[ "${EASYFUSE_LAYOUT}" != "legacy" ]] || bash "${ROOT}/scripts/patch_easyfuse_star_avx2.sh"

PASS_CSV="${OUT}/${SAMPLE_ID}/fusions.pass.csv"
echo ""
echo "==> Done. Check:"
echo "    ${PASS_CSV}"
echo ""
echo "Next (pipeline adapter):"
echo "  neoag build-intermediates --entry-mode fusion \\"
echo "    --easyfuse-tsv ${PASS_CSV} \\"
echo "    --sample-id ${SAMPLE_ID} \\"
echo "    --outdir ${OUT}/intermediates"
