#!/usr/bin/env bash
# FACETS smoke test on HCC1143 subset BAMs (tutorial_9183).
#
# Usage:
#   bash scripts/run_facets_hcc1143.sh            # HCC1143 50K-read subset
#   bash scripts/run_facets_hcc1143.sh --example  # bundled stomach WES example (no BAMs)
#
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${ROOT}/conf/tools.env.sh"

RSCRIPT="${NEOAG_CONDA_BASE}/envs/neoag-fusion/bin/Rscript"
export PATH="${NEOAG_CONDA_BASE}/envs/neoag-fusion/bin:${NEOAG_CONDA_BASE}/envs/neoag-tools/bin:${NEOAG_CONDA_BASE}/envs/neoag-gatk/bin:${ROOT}/bin:${PATH}"

MODE="hcc1143"
if [[ "${1:-}" == "--example" ]]; then
  MODE="example"
fi

SAMPLE="hcc1143"
OUTDIR="${ROOT}/work/facets_hcc1143"
LOG="${ROOT}/work/run_facets_hcc1143.log"
DBSNP_VCF="${NEOAG_DBSNP_VCF:?ERROR: set NEOAG_DBSNP_VCF=/path/to/dbsnp.vcf.gz}"
TUMOR_SRC="${HCC1143_TUMOR_BAM:-${ROOT}/data/examples/HCC1143/hcc1143_T_subset50K.bam}"
NORMAL_SRC="${HCC1143_NORMAL_BAM:-${ROOT}/data/examples/HCC1143/hcc1143_N_subset50K.bam}"

mkdir -p "${OUTDIR}"

{
  echo "==> run_facets_hcc1143 $(date -Iseconds)"
  echo "    mode=${MODE}"
  echo "    outdir=${OUTDIR}"
  echo "    dbsnp=${DBSNP_VCF}"
} | tee "${LOG}"

if [[ ! -x "${ROOT}/bin/snp-pileup" ]]; then
  echo "ERROR: snp-pileup not found at ${ROOT}/bin/snp-pileup" | tee -a "${LOG}"
  exit 1
fi
if [[ ! -x "${RSCRIPT}" ]]; then
  echo "ERROR: neoag-fusion Rscript not found: ${RSCRIPT}" | tee -a "${LOG}"
  exit 1
fi
if ! "${RSCRIPT}" -e 'stopifnot(requireNamespace("facets", quietly=TRUE))' 2>/dev/null; then
  echo "ERROR: facets R package missing in neoag-fusion" | tee -a "${LOG}"
  exit 1
fi

RDS="${OUTDIR}/${SAMPLE}.facets.rds"
PILEUP="${OUTDIR}/${SAMPLE}.snp_pileup.csv.gz"
PURITY_TXT="${OUTDIR}/facets_purity.txt"
CNCF_TSV="${OUTDIR}/facets_cncf.tsv"
SUMMARY="${OUTDIR}/facets_summary.txt"

if [[ "${MODE}" == "example" ]]; then
  echo "==> FACETS bundled stomach example (no snp-pileup)" | tee -a "${LOG}"
  SAMPLE="stomach_example"
  RDS="${OUTDIR}/${SAMPLE}.facets.rds"
  PURITY_TXT="${OUTDIR}/facets_purity.txt"
  CNCF_TSV="${OUTDIR}/facets_cncf.tsv"
  export FACETS_OUT_RDS="${RDS}"
  "${RSCRIPT}" - <<'RS' | tee -a "${LOG}"
suppressPackageStartupMessages(library(facets))
set.seed(1234)
datafile <- system.file("extdata", "stomach.csv.gz", package = "facets")
rcmat <- readSnpMatrix(datafile)
xx <- preProcSample(rcmat)
oo <- procSample(xx, cval = 150)
fit <- emcncf(oo)
saveRDS(fit, Sys.getenv("FACETS_OUT_RDS"))
cat("purity:", fit$purity, "ploidy:", fit$ploidy, "\n")
RS
else
  for bam in "${TUMOR_SRC}" "${NORMAL_SRC}"; do
    if [[ ! -f "${bam}" ]]; then
      echo "ERROR: missing BAM: ${bam}" | tee -a "${LOG}"
      exit 1
    fi
  done
  if [[ ! -f "${DBSNP_VCF}" ]]; then
    echo "ERROR: missing dbSNP VCF: ${DBSNP_VCF}" | tee -a "${LOG}"
    exit 1
  fi

  echo "==> Index BAMs if needed" | tee -a "${LOG}"
  for bam in "${NORMAL_SRC}" "${TUMOR_SRC}"; do
    if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" ]]; then
      samtools index -@ 4 "${bam}"
    fi
  done

  # Low-coverage subset: relax depth filters (normal,tumor min reads = 0,0).
  echo "==> snp-pileup (normal then tumor)" | tee -a "${LOG}"
  snp-pileup -g -q15 -Q20 -P100 -r0,0 \
    "${DBSNP_VCF}" "${PILEUP}" "${NORMAL_SRC}" "${TUMOR_SRC}"

  PILEUP_ROWS="$(zcat "${PILEUP}" | tail -n +2 | wc -l | tr -d ' ')"
  echo "    pileup rows: ${PILEUP_ROWS}" | tee -a "${LOG}"
  if [[ "${PILEUP_ROWS}" -lt 10 ]]; then
    echo "WARN: very few SNPs; falling back to bundled stomach example" | tee -a "${LOG}"
    MODE="example_fallback"
    SAMPLE="stomach_example"
    RDS="${OUTDIR}/${SAMPLE}.facets.rds"
    export FACETS_OUT_RDS="${RDS}"
    "${RSCRIPT}" - <<'RS' | tee -a "${LOG}"
suppressPackageStartupMessages(library(facets))
set.seed(1234)
datafile <- system.file("extdata", "stomach.csv.gz", package = "facets")
rcmat <- readSnpMatrix(datafile)
xx <- preProcSample(rcmat)
oo <- procSample(xx, cval = 150)
fit <- emcncf(oo)
saveRDS(fit, Sys.getenv("FACETS_OUT_RDS"))
cat("purity:", fit$purity, "ploidy:", fit$ploidy, "\n")
RS
  else
    echo "==> FACETS fit from pileup" | tee -a "${LOG}"
    export FACETS_NDEPTH=3 FACETS_CVAL_PRE=15 FACETS_CVAL_PROC=15 FACETS_MIN_NHET=3
    "${RSCRIPT}" "${ROOT}/scripts/facets_fit_from_pileup.R" "${PILEUP}" "${RDS}"
  fi
fi

echo "==> Extract purity/CNV tables" | tee -a "${LOG}"
"${RSCRIPT}" "${ROOT}/tools/facets/runFACETS.R" "${RDS}" "${PURITY_TXT}" "${CNCF_TSV}"

{
  echo "sample=${SAMPLE}"
  echo "mode=${MODE}"
  echo "rds=${RDS}"
  echo "pileup=${PILEUP}"
  echo "purity=$(tail -n +2 "${PURITY_TXT}")"
  echo "cncf_segments=$(tail -n +2 "${CNCF_TSV}" | wc -l | tr -d ' ')"
} > "${SUMMARY}"

echo "==> done" | tee -a "${LOG}"
ls -lh "${OUTDIR}" | tee -a "${LOG}"
cat "${SUMMARY}" | tee -a "${LOG}"
