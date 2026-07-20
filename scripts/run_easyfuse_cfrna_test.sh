#!/usr/bin/env bash
# EasyFuse smoke test on NAS cfrna_paper FASTQ (PRJNA1199029).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
# shellcheck source=/dev/null
source "${ROOT}/conf/tools.env.sh"
# EasyFuse Nextflow activates starfusion.yml; keep repo STAR-Fusion off PATH.
export PATH="$(echo "${PATH}" | tr ':' '\n' | grep -vE "${NEOAG_STAR_FUSION_HOME}\$" | paste -sd: -)"

REF="${NEOAG_EASYFUSE_REF:?ERROR: set NEOAG_EASYFUSE_REF=/path/to/easyfuse_ref_v4}"
INPUT="${ROOT}/work/easyfuse_cfrna_input.tsv"
OUT="${ROOT}/results/easyfuse_cfrna_test"
LOG="${ROOT}/work/run_easyfuse_cfrna_test.log"
PREBUILD_LOG="${ROOT}/work/easyfuse_conda_prebuild.log"

mkdir -p "${OUT}" "$(dirname "${LOG}")" "${ROOT}/work/.nextflow_home" "${ROOT}/work/.nextflow_work"

export NXF_HOME="${ROOT}/work/.nextflow_home"
export NXF_WORK="${ROOT}/work/.nextflow_work"
export NXF_DISABLE_CHECK_TTY=true
export CONDA_ALWAYS_YES=true
export MAMBA_YES=1
export NEOAG_REAL_MAMBA="${NEOAG_CONDA_BASE}/bin/mamba"
export JAVA_HOME="${NEOAG_CONDA_BASE}/envs/${NEOAG_FUSION_ENV}"
export PATH="${ROOT}/work/easyfuse_bin:${NEOAG_CONDA_BASE}/bin:${JAVA_HOME}/bin:${PATH}"
# Do not prepend repo STAR-Fusion — EasyFuse Nextflow uses its own starfusion conda env.
export NEOAG_NEXTFLOW="${JAVA_HOME}/bin/nextflow"

exec > >(tee -a "${LOG}") 2>&1
echo "==> run_easyfuse_cfrna_test $(date -Is)"
echo "    input=${INPUT}"
echo "    reference=${REF}"
echo "    output=${OUT}"

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

# Build remaining envs with mamba -y (Nextflow subprocess omits -y and hangs).
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

# Prebuild STAR-only env used by STAR_CUSTOM (may be created on-the-fly by Nextflow).
REQ_WO_ENV="${CONDA_CACHE}/env-requantification_wo_easyfuse"
if [[ ! -x "${REQ_WO_ENV}/bin/STAR" ]]; then
  wait_for_mamba_free
  echo "==> Pre-building requantification_wo_easyfuse.yml ..."
  rm -rf "${REQ_WO_ENV}"
  mamba env create -y \
    --prefix "${REQ_WO_ENV}" \
    --file "${NEOAG_EASYFUSE_HOME}/environments/requantification_wo_easyfuse.yml"
fi

# neoag-tools python shadows conda env `#!/usr/bin/env python` entrypoints during Nextflow tasks.
export PATH="$(echo "${PATH}" | tr ':' '\n' | grep -vE '/envs/neoag-tools/bin$' | paste -sd: -)"

cd "${ROOT}/work"

"${NEOAG_NEXTFLOW}" run "${NEOAG_EASYFUSE_HOME}/main.nf" \
  -c "${ROOT}/conf/easyfuse.nextflow.config" \
  -profile conda \
  -w "${NXF_WORK}" \
  --output "${OUT}" \
  --input_files "${INPUT}" \
  --reference "${REF}" \
  --annotation_db "${REF}/Homo_sapiens.GRCh38.110.gff3.db" \
  --reference_tsl "${REF}/Homo_sapiens.GRCh38.110.gtf.tsl" \
  -resume </dev/null

echo ""
echo "==> Done. Check:"
echo "    ${OUT}/SRR31741871/fusions.pass.csv"
echo ""
echo "Next (pipeline adapter):"
echo "  neoag build-intermediates --entry-mode fusion \\"
echo "    --easyfuse-tsv ${OUT}/SRR31741871/fusions.pass.csv \\"
echo "    --sample-id SRR31741871 --outdir results/easyfuse_cfrna_test/intermediates"
