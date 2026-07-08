#!/usr/bin/env bash
# LOHHLA validation run on bundled example BAMs.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${ROOT}/conf/tools.env.sh"

export PATH="${NEOAG_CONDA_BASE}/envs/neoag-fusion/bin:${NEOAG_CONDA_BASE}/envs/neoag-tools/bin:${PATH}"

PATIENT="example"
WORKDIR="${ROOT}/work/lohhla_example/out"
BAMDIR="${ROOT}/tools/lohhla/example-file/bam"
LOG="${ROOT}/work/run_lohhla_example.log"
TUMOR_BAM="${BAMDIR}/example_tumor_sorted.bam"
NORMAL_BAM="${BAMDIR}/example_BS_GL_sorted.bam"
HLA_FASTA="${ROOT}/tools/lohhla/example-file/correct-example-out/example.patient.hlaFasta.fa"
HLAS="${ROOT}/tools/lohhla/example-file/hlas"
COPYNUM="${ROOT}/tools/lohhla/example-file/solutions.txt"
GATK_DIR="${ROOT}/tools/picard-lohhla"
NOVO_DIR="${POLYSOLVER_HOME:-${ROOT}/tools/polysolver}/binaries"
NOVOINDEX="${POLYSOLVER_HOME:-${ROOT}/tools/polysolver}/scripts/novoindex"

mkdir -p "${WORKDIR}"
export PATH="${NOVO_DIR}:${NOVOINDEX%/*}:${PATH}"
export NOVOALIGN_LICENSE_FILE="${NOVOALIGN_LICENSE_FILE:?ERROR: set NOVOALIGN_LICENSE_FILE=/path/to/novoalign.lic}"
cp "${NOVOALIGN_LICENSE_FILE}" "${NOVO_DIR}/novoalign.lic"

echo "[$(date -Iseconds)] starting LOHHLA example" | tee -a "${LOG}"
Rscript "${ROOT}/tools/lohhla/LOHHLAscript.R" \
  --patientId "${PATIENT}" \
  --outputDir "${WORKDIR}" \
  --normalBAMfile "${NORMAL_BAM}" \
  --tumorBAMfile "${TUMOR_BAM}" \
  --hlaPath "${HLAS}" \
  --HLAfastaLoc "${HLA_FASTA}" \
  --CopyNumLoc "${COPYNUM}" \
  --mappingStep TRUE \
  --fishingStep TRUE \
  --plottingStep FALSE \
  --coverageStep TRUE \
  --minCoverageFilter 10 \
  --cleanUp FALSE \
  --gatkDir "${GATK_DIR}" \
  --novoDir "${NOVO_DIR}" \
  --HLAexonLoc "${ROOT}/tools/lohhla/data/hla.dat" \
  2>&1 | tee -a "${LOG}"

echo "[$(date -Iseconds)] done" | tee -a "${LOG}"
ls -la "${WORKDIR}" | tee -a "${LOG}"
