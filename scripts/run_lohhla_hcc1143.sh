#!/usr/bin/env bash
# LOHHLA smoke test on HCC1143 subset BAMs (tutorial_9183).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${ROOT}/conf/tools.env.sh"

export PATH="${NEOAG_CONDA_BASE}/envs/neoag-fusion/bin:${NEOAG_CONDA_BASE}/envs/neoag-tools/bin:${PATH}"

PATIENT="hcc1143"
WORKDIR="${ROOT}/work/lohhla_hcc1143/out"
BAMDIR="${ROOT}/work/lohhla_hcc1143/bam"
LOG="${ROOT}/work/run_lohhla_hcc1143.log"
TUMOR_SRC="${HCC1143_TUMOR_BAM:-${ROOT}/data/examples/HCC1143/hcc1143_T_subset50K.bam}"
NORMAL_SRC="${HCC1143_NORMAL_BAM:-${ROOT}/data/examples/HCC1143/hcc1143_N_subset50K.bam}"
TUMOR_BAM="${BAMDIR}/tumor.bam"
NORMAL_BAM="${BAMDIR}/normal.bam"
HLA_FASTA="${POLYSOLVER_HLA_FASTA:-${POLYSOLVER_HOME:-${ROOT}/tools/polysolver}/data/abc_complete.fasta}"
HLAS="${ROOT}/work/lohhla_hcc1143/hlas"
GATK_DIR="${ROOT}/tools/picard-lohhla"
NOVO_DIR="${POLYSOLVER_HOME:-${ROOT}/tools/polysolver}/binaries"
NOVOINDEX="${POLYSOLVER_HOME:-${ROOT}/tools/polysolver}/scripts/novoindex"

mkdir -p "${WORKDIR}" "${BAMDIR}"
ln -sf "${TUMOR_SRC}" "${TUMOR_BAM}"
ln -sf "${NORMAL_SRC}" "${NORMAL_BAM}"
export PATH="${NOVO_DIR}:${NOVOINDEX%/*}:${PATH}"
export NOVOALIGN_LICENSE_FILE="${NOVOALIGN_LICENSE_FILE:?ERROR: set NOVOALIGN_LICENSE_FILE=/path/to/novoalign.lic}"
cp "${NOVOALIGN_LICENSE_FILE}" "${NOVO_DIR}/novoalign.lic"

echo "[$(date -Iseconds)] indexing BAMs if needed" | tee -a "${LOG}"
for bam in "${TUMOR_BAM}" "${NORMAL_BAM}"; do
  if [[ ! -f "${bam}.bai" ]] || [[ "${bam}.bai" -ot "${bam}" ]]; then
    samtools index -@ 4 "${bam}"
  fi
done

echo "[$(date -Iseconds)] starting LOHHLA" | tee -a "${LOG}"
Rscript "${ROOT}/tools/lohhla/LOHHLAscript.R" \
  --patientId "${PATIENT}" \
  --outputDir "${WORKDIR}" \
  --normalBAMfile "${NORMAL_BAM}" \
  --tumorBAMfile "${TUMOR_BAM}" \
  --hlaPath "${HLAS}" \
  --HLAfastaLoc "${HLA_FASTA}" \
  --CopyNumLoc FALSE \
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
