#!/usr/bin/env bash
# Run FACETS using a clean biallelic SNP-only pileup, then
# stride-downsample before FACETS fitting. This mode avoids FACETS pseudo-depth
# sites from snp-pileup -P and is useful as a purity robustness check.
#
# Required environment:
#   TUMOR_BAM=/path/tumor.bam
#   NORMAL_BAM=/path/normal.bam
#
# Common overrides:
#   PATIENT_ID=sample_id
#   OUTDIR=results/facets/<sample>/omni2p5_biallelic_snponly_downsample1m
#   FACETS_SNP_VCF=data/facets/reference/1000G_omni2.5.hg38.biallelic.vcf.gz
#   FACETS_TARGET_ROWS=1000000
#   FACETS_CVAL_PRE=25
#   FACETS_CVAL_PROC=25
#   FACETS_MIN_NHET=5
#   FACETS_STEP=all|pileup|downsample|fit|export
#
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
# shellcheck source=/dev/null
source "${ROOT}/conf/tools.env.sh"
export PYTHONPATH="${ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}"
export PATH="${NEOAG_CONDA_BASE}/envs/${NEOAG_FUSION_ENV:-neoag-fusion-r36}/bin:${ROOT}/bin:${PATH}"

PATIENT_ID="${PATIENT_ID:-sample}"
TUMOR_NAME="${TUMOR_NAME:-${PATIENT_ID}}"
NORMAL_NAME="${NORMAL_NAME:-normal}"
TUMOR_BAM="${TUMOR_BAM:?ERROR: set TUMOR_BAM=/path/tumor.bam}"
NORMAL_BAM="${NORMAL_BAM:?ERROR: set NORMAL_BAM=/path/normal.bam}"
OUT="${OUTDIR:-${ROOT}/results/facets_omni2p5_snponly/${PATIENT_ID}}"
LOG="${LOG:-${OUT}/run.log}"
FACETS_STEP="${FACETS_STEP:-all}"

FACETS_SNPSET_NAME="${FACETS_SNPSET_NAME:-omni2p5}"
FACETS_SNP_VCF="${FACETS_SNP_VCF:-${OMNI2P5_VCF:-${ROOT}/data/facets/reference/1000G_omni2.5.hg38.biallelic.vcf.gz}}"
SNP_PILEUP_BIN="${SNP_PILEUP_BIN:-${ROOT}/bin/snp-pileup}"
RSCRIPT="${RSCRIPT:-${NEOAG_CONDA_BASE}/envs/${NEOAG_FUSION_ENV:-neoag-fusion-r36}/bin/Rscript}"
SNP_PILEUP_MIN_READS="${SNP_PILEUP_MIN_READS:-5,5}"
FACETS_TARGET_ROWS="${FACETS_TARGET_ROWS:-1000000}"
FACETS_NDEPTH="${FACETS_NDEPTH:-5}"
FACETS_CVAL_PRE="${FACETS_CVAL_PRE:-25}"
FACETS_CVAL_PROC="${FACETS_CVAL_PROC:-25}"
FACETS_MIN_NHET="${FACETS_MIN_NHET:-5}"

FACETS_HOME="${FACETS_HOME:-${NEOAG_TOOLS_ROOT}/tools/facets}"
FACETS_QUARANTINE="${ROOT}/../neoag_event_pipeline_v03_rc_artifact_quarantine_20260622_091158/tools/facets"
RUN_FACETS_R="${FACETS_HOME}/runFACETS.R"

PILEUP="${OUT}/${PATIENT_ID}.${FACETS_SNPSET_NAME}.snponly.pileup.csv"
DOWNSAMPLED="${OUT}/${PATIENT_ID}.${FACETS_SNPSET_NAME}.snponly.downsample${FACETS_TARGET_ROWS}.csv"
RDS="${OUT}/${PATIENT_ID}.facets.rds"
PURITY_TXT="${OUT}/facets_purity.txt"
CNCF_TSV="${OUT}/facets_cncf.tsv"
PURITY_TSV="${OUT}/purity.tsv"
CNV_TSV="${OUT}/cnv_segments.tsv"
SUMMARY="${OUT}/facets_${FACETS_SNPSET_NAME}_summary.tsv"

mkdir -p "${OUT}"
exec > >(tee -a "${LOG}") 2>&1

step_wanted() {
  [[ "${FACETS_STEP}" == "all" || "${FACETS_STEP}" == "$1" ]]
}

resolve_facets_home() {
  if [[ -f "${RUN_FACETS_R}" ]]; then
    return 0
  fi
  if [[ -f "${FACETS_QUARANTINE}/runFACETS.R" ]]; then
    FACETS_HOME="${FACETS_QUARANTINE}"
    RUN_FACETS_R="${FACETS_HOME}/runFACETS.R"
    return 0
  fi
  echo "ERROR: runFACETS.R not found under ${FACETS_HOME} or ${FACETS_QUARANTINE}" >&2
  exit 1
}

check_prereqs() {
  [[ -x "${SNP_PILEUP_BIN}" ]] || { echo "ERROR: missing snp-pileup: ${SNP_PILEUP_BIN}" >&2; exit 1; }
  [[ -x "${RSCRIPT}" ]] || { echo "ERROR: missing Rscript: ${RSCRIPT}" >&2; exit 1; }
  [[ -s "${FACETS_SNP_VCF}" ]] || { echo "ERROR: missing FACETS_SNP_VCF: ${FACETS_SNP_VCF}" >&2; exit 1; }
  [[ -s "${TUMOR_BAM}" ]] || { echo "ERROR: missing TUMOR_BAM: ${TUMOR_BAM}" >&2; exit 1; }
  [[ -s "${NORMAL_BAM}" ]] || { echo "ERROR: missing NORMAL_BAM: ${NORMAL_BAM}" >&2; exit 1; }
  [[ -s "${TUMOR_BAM}.bai" ]] || { echo "ERROR: missing tumor BAI: ${TUMOR_BAM}.bai" >&2; exit 1; }
  [[ -s "${NORMAL_BAM}.bai" ]] || { echo "ERROR: missing normal BAI: ${NORMAL_BAM}.bai" >&2; exit 1; }
  if ! "${RSCRIPT}" -e 'stopifnot(requireNamespace("facets", quietly=TRUE))' 2>/dev/null; then
    echo "ERROR: facets R package missing in ${RSCRIPT}" >&2
    exit 1
  fi
  resolve_facets_home
}

write_summary() {
  local status="$1"
  cat > "${SUMMARY}" <<EOS
metric	value
sample_id	${PATIENT_ID}
tumor_name	${TUMOR_NAME}
normal_name	${NORMAL_NAME}
tumor_bam	${TUMOR_BAM}
normal_bam	${NORMAL_BAM}
snp_vcf	${FACETS_SNP_VCF}
snp_set_name	${FACETS_SNPSET_NAME}
pileup	${PILEUP}
downsampled_pileup	${DOWNSAMPLED}
rds	${RDS}
purity_tsv	${PURITY_TSV}
cnv_tsv	${CNV_TSV}
mode	${FACETS_SNPSET_NAME}_biallelic_snponly_downsample
pseudo_snps	false
snp_pileup_min_reads	${SNP_PILEUP_MIN_READS}
target_rows	${FACETS_TARGET_ROWS}
facets_ndepth	${FACETS_NDEPTH}
facets_cval_pre	${FACETS_CVAL_PRE}
facets_cval_proc	${FACETS_CVAL_PROC}
facets_min_nhet	${FACETS_MIN_NHET}
status	${status}
EOS
}

run_pileup() {
  echo "==> ${FACETS_SNPSET_NAME} SNP-only snp-pileup $(date -Is)"
  echo "    sample=${PATIENT_ID}"
  echo "    tumor_bam=${TUMOR_BAM}"
  echo "    normal_bam=${NORMAL_BAM}"
  echo "    vcf=${FACETS_SNP_VCF}"
  echo "    output=${PILEUP}"
  echo "    params=-q15 -Q20 -r${SNP_PILEUP_MIN_READS} without -P pseudo-snps"
  "${SNP_PILEUP_BIN}" -q15 -Q20 -r"${SNP_PILEUP_MIN_READS}" \
    "${FACETS_SNP_VCF}" "${PILEUP}" "${NORMAL_BAM}" "${TUMOR_BAM}"
  local rows
  rows="$(tail -n +2 "${PILEUP}" | wc -l | tr -d " ")"
  echo "    pileup_rows=${rows}"
}

run_downsample() {
  [[ -s "${PILEUP}" ]] || { echo "ERROR: missing pileup: ${PILEUP}" >&2; exit 1; }
  local rows stride down_rows
  rows="$(tail -n +2 "${PILEUP}" | wc -l | tr -d " ")"
  stride=$(( rows / FACETS_TARGET_ROWS ))
  [[ "${stride}" -lt 1 ]] && stride=1
  echo "==> stride downsample $(date -Is) rows=${rows} target=${FACETS_TARGET_ROWS} stride=${stride}"
  awk -F, -v stride="${stride}" 'NR==1 || (NR>1 && ((NR-2)%stride)==0)' "${PILEUP}" > "${DOWNSAMPLED}"
  down_rows="$(tail -n +2 "${DOWNSAMPLED}" | wc -l | tr -d " ")"
  cat > "${OUT}/downsample_summary.tsv" <<EOS
metric	value
source	${PILEUP}
source_rows	${rows}
stride	${stride}
downsample_rows	${down_rows}
EOS
  cat "${OUT}/downsample_summary.tsv"
}

run_fit() {
  [[ -s "${DOWNSAMPLED}" ]] || { echo "ERROR: missing downsampled pileup: ${DOWNSAMPLED}" >&2; exit 1; }
  echo "==> FACETS fit $(date -Is)"
  export FACETS_NDEPTH FACETS_CVAL_PRE FACETS_CVAL_PROC FACETS_MIN_NHET
  "${RSCRIPT}" "${ROOT}/scripts/facets_fit_from_pileup.R" "${DOWNSAMPLED}" "${RDS}"
}

run_export() {
  [[ -s "${RDS}" ]] || { echo "ERROR: missing FACETS RDS: ${RDS}" >&2; exit 1; }
  echo "==> FACETS export $(date -Is)"
  if "${RSCRIPT}" "${RUN_FACETS_R}" "${RDS}" "${PURITY_TXT}" "${CNCF_TSV}"; then
    "${ROOT}/bin/neoag-v03" convert-facets \
      --purity-input "${PURITY_TXT}" \
      --purity-output "${PURITY_TSV}" \
      --sample-id "${PATIENT_ID}" \
      --cnv-input "${CNCF_TSV}" \
      --cnv-output "${CNV_TSV}"
  else
    echo "WARNING: FACETS export did not produce purity. Writing indeterminate purity.tsv."
    cat > "${PURITY_TSV}" <<EOS
sample_id	purity	ploidy	source	note
${PATIENT_ID}	NA	NA	facets_${FACETS_SNPSET_NAME}_biallelic_snponly_downsample	FACETS did not export purity; inspect ${RDS} and ${LOG}.
EOS
  fi
}

check_prereqs
write_summary "started"

if step_wanted pileup; then
  run_pileup
fi
if step_wanted downsample; then
  run_downsample
fi
if step_wanted fit; then
  run_fit
fi
if step_wanted export; then
  run_export
fi

write_summary "finished"
echo "==> run_facets_omni2p5_snponly_downsample finished $(date -Is)"
ls -lh "${OUT}"
