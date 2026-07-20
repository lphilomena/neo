#!/usr/bin/env bash
# Generic EasyFuse runner for one paired-end RNA-seq sample.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
# shellcheck source=/dev/null
source "${ROOT}/conf/tools.env.sh"
: "${NEOAG_CONDA_BASE:?ERROR: set NEOAG_CONDA_BASE to your conda/mamba installation root}"
export PATH="${NEOAG_CONDA_BASE}/bin:${PATH}"
# EasyFuse Nextflow activates starfusion.yml; keep repo STAR-Fusion off PATH.
export PATH="$(echo "${PATH}" | tr ':' '\n' | grep -vE "${NEOAG_STAR_FUSION_HOME}\$" | paste -sd: -)"

SAMPLE_ID="${EASYFUSE_SAMPLE_ID:-${SAMPLE_ID:-sample}}"
FQ1="${EASYFUSE_FQ1:?ERROR: set EASYFUSE_FQ1=/path/sample_R1.fq.gz}"
FQ2="${EASYFUSE_FQ2:?ERROR: set EASYFUSE_FQ2=/path/sample_R2.fq.gz}"
REF="${NEOAG_EASYFUSE_REF:?ERROR: set NEOAG_EASYFUSE_REF=/path/to/easyfuse_ref_v4}"
INPUT="${EASYFUSE_INPUT_TSV:-${ROOT}/work/easyfuse_${SAMPLE_ID}_input.tsv}"
OUT="${OUTDIR:-${ROOT}/results/easyfuse}"
LOG="${LOG:-${ROOT}/work/run_easyfuse_${SAMPLE_ID}.log}"
PREBUILD_LOG="${ROOT}/work/easyfuse_conda_prebuild.log"

ensure_input_tsv() {
  printf '%s\t%s\t%s\n' "${SAMPLE_ID}" "${FQ1}" "${FQ2}" > "${INPUT}"
  if ! awk -F'\t' 'NF==3 {found=1} END{exit !found}' "${INPUT}"; then
    echo "ERROR: input TSV must have 3 tab-separated columns: ${INPUT}" >&2
    exit 1
  fi
}

STAR_TMP="${ROOT}/work/star_tmp"
mkdir -p "${OUT}" "$(dirname "${LOG}")" "${ROOT}/work/.nextflow_home" "${ROOT}/work/.nextflow_work" "${STAR_TMP}"
export TMPDIR="${STAR_TMP}"
ensure_input_tsv

# Avoid Nextflow session lock / STAR temp collisions with other EasyFuse runs.
for stale_pid in $(pgrep -f 'run_easyfuse_cfrna_test\.sh' 2>/dev/null || true); do
  echo "==> stopping stale easyfuse_cfrna_test PID=${stale_pid}"
  kill "${stale_pid}" 2>/dev/null || true
done
while pgrep -f 'nextflow.*EasyFuse/main.nf' >/dev/null 2>&1; do
  echo "==> waiting for other EasyFuse Nextflow runs to finish ..."
  sleep 30
done

export NXF_HOME="${ROOT}/work/.nextflow_home"
export NXF_WORK="${ROOT}/work/.nextflow_work"
export NXF_DISABLE_CHECK_TTY=true
export CONDA_ALWAYS_YES=true
export MAMBA_YES=1
export NEOAG_REAL_MAMBA="${NEOAG_CONDA_BASE}/bin/mamba"
export JAVA_HOME="${NEOAG_CONDA_BASE}/envs/${NEOAG_FUSION_ENV}"
export PATH="${ROOT}/work/easyfuse_bin:${NEOAG_CONDA_BASE}/bin:${JAVA_HOME}/bin:${PATH}"
export NEOAG_NEXTFLOW="${JAVA_HOME}/bin/nextflow"

exec > >(tee -a "${LOG}") 2>&1
echo "==> run_easyfuse_sample $(date -Is)"
echo "    sample=${SAMPLE_ID}"
echo "    fq1=${FQ1}"
echo "    fq2=${FQ2}"
echo "    input=${INPUT}"
echo "    reference=${REF}"
echo "    output=${OUT}"

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

CONDA_CACHE="${ROOT}/work/.nextflow_conda"
mkdir -p "${CONDA_CACHE}"

wait_for_mamba_free() {
  while pgrep -f 'mamba env create' >/dev/null 2>&1; do
    sleep 15
  done
  rm -f "${NEOAG_CONDA_BASE}/pkgs/pkgs.lock" 2>/dev/null || true
}

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

QC_ENV="${CONDA_CACHE}/env-749ebc089f673418-1f348f31c1e78ea89e97e435a63f0c7d"
SRC_ENV="${CONDA_CACHE}/env-e6082ee0f0a13e81-c203347e504f3b4d10ed96fdd01318ce"

if [[ ! -x "${QC_ENV}/bin/fastqc" || ! -x "${SRC_ENV}/bin/skewer" ]]; then
  wait_for_mamba_free
fi

prebuild_conda_env \
  "749ebc089f673418-1f348f31c1e78ea89e97e435a63f0c7d" \
  "qc.yml" \
  "fastqc"

prebuild_conda_env \
  "e6082ee0f0a13e81-c203347e504f3b4d10ed96fdd01318ce" \
  "easyfuse_src.yml" \
  "skewer"

if ! pgrep -f 'easyfuse_prebuild_remaining_envs\.sh' >/dev/null 2>&1; then
  nohup bash "${ROOT}/scripts/easyfuse_prebuild_remaining_envs.sh" >/dev/null 2>&1 &
  echo "==> background conda prebuild worker PID=$! (log: ${PREBUILD_LOG})"
else
  echo "==> background conda prebuild worker already running"
fi

ALIGN_ENV="${CONDA_CACHE}/env-6f2b394c864eeaa5-8f88fe4572f59d9bb818f7644ca8f1fa"
echo "==> waiting for alignment env (STAR) before Nextflow ..."
while [[ ! -x "${ALIGN_ENV}/bin/STAR" ]]; do
  if ! pgrep -f 'easyfuse_prebuild_remaining_envs\.sh' >/dev/null 2>&1 \
    && ! pgrep -f 'mamba env create.*6f2b394c864eeaa5' >/dev/null 2>&1; then
    echo "ERROR: alignment env build failed; see ${PREBUILD_LOG}" >&2
    exit 1
  fi
  sleep 20
done
echo "==> alignment env ready"

bash "${ROOT}/scripts/fix_easyfuse_pyeasyfuse_env.sh"
bash "${ROOT}/scripts/seed_easyfuse_conda_envs.sh"

REQ_WO_ENV="${CONDA_CACHE}/env-requantification_wo_easyfuse"
if [[ ! -x "${REQ_WO_ENV}/bin/STAR" ]]; then
  wait_for_mamba_free
  echo "==> Pre-building requantification_wo_easyfuse.yml ..."
  rm -rf "${REQ_WO_ENV}"
  mamba env create -y \
    --prefix "${REQ_WO_ENV}" \
    --file "${NEOAG_EASYFUSE_HOME}/environments/requantification_wo_easyfuse.yml"
fi

bash "${ROOT}/scripts/patch_easyfuse_star_avx2.sh"
bash "${ROOT}/scripts/fix_easyfuse_pyeasyfuse_env.sh"

export PATH="$(echo "${PATH}" | tr ':' '\n' | grep -vE '/envs/neoag-tools/bin$|/tools/fusioncatcher/bin$' | paste -sd: -)"

cd "${ROOT}/work"

# Separate run name avoids resume/lock collision with easyfuse_cfrna_test (session 9d03387c...).
NXF_RUN_NAME="${EASYFUSE_RUN_NAME:-easyfuse_${SAMPLE_ID}}"
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
  "${NEOAG_NEXTFLOW}" run "${NEOAG_EASYFUSE_HOME}/main.nf" \
    "${NXF_RESUME_ARGS[@]}" \
    -c "${ROOT}/conf/easyfuse.nextflow.config" \
    -profile conda \
    -w "${NXF_WORK}" \
    --output "${OUT}" \
    --input_files "${INPUT}" \
    --reference "${REF}" \
    --annotation_db "${REF}/Homo_sapiens.GRCh38.110.gff3.db" \
    --reference_tsl "${REF}/Homo_sapiens.GRCh38.110.gtf.tsl" </dev/null
}

run_nextflow || {
  echo "==> Nextflow failed; patching STAR and retrying once with -resume ..."
  bash "${ROOT}/scripts/patch_easyfuse_star_avx2.sh"
  NXF_RESUME_ARGS=(-resume "${NXF_RUN_NAME}")
  run_nextflow
}

bash "${ROOT}/scripts/patch_easyfuse_star_avx2.sh"

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
