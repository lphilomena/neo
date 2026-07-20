#!/usr/bin/env bash
# ============================================================================
# run_main_all_qc.sh — Launch the main_all_qc Nextflow workflow.
#
# Full pipeline: HLA typing → variant calling → neoantigen scoring
#       + parallel QC: LOH_CHECK (LOHHLA + SpecHLA)
#                     PURITY_CHECK (FACETS + PURPLE)
#
# Prerequisites:
#   source conf/tools.env.sh
#   docker pull sequenza/sequenza:3.99rc1   # if using Docker mode
#
# Usage:
#   # Minimal (conda mode):
#   bash scripts/run_main_all_qc.sh \
#     --normal_bam /path/to/normal.bam \
#     --tumor_bam /path/to/tumor.bam \
#     --sample_id SAMPLE001
#
#   # Full (with dbSNP VCF for FACETS):
#   bash scripts/run_main_all_qc.sh \
#     --normal_bam /path/to/normal.bam \
#     --tumor_bam /path/to/tumor.bam \
#     --sample_id SAMPLE001 \
#     --dbsnp_vcf /path/to/dbsnp_chr.vcf.gz \
#     --outdir /path/to/results
#
#   # Docker mode:
#   NEOAG_RUNNER_MODE=docker \
#   bash scripts/run_main_all_qc.sh \
#     --normal_bam /path/to/normal.bam \
#     --tumor_bam /path/to/tumor.bam \
#     --sample_id SAMPLE001
#
#   # Skip QC subworkflows (main pipeline only):
#   bash scripts/run_main_all_qc.sh --skip_qc ...
#
#   # Dry-run (print commands only):
#   bash scripts/run_main_all_qc.sh --dry-run ...
#
# Environment overrides:
#   NEOAG_RUNNER_MODE          "conda" (default) or "docker"
#   NEOAG_REFERENCE_FASTA      Reference genome FASTA path
#   NEOAG_DBSNP_VCF            dbSNP/common SNP VCF for FACETS snp-pileup
#   NXF_HOME                   Nextflow metadata home (default: work/.nextflow_home)
#   NEOAG_PROFILE               Nextflow -profile value (default: none)
#   NEOAG_RESUME                Set to "1" to add -resume
# ============================================================================

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "${ROOT}"

# ---------------------------------------------------------------------------
# Environment setup
# ---------------------------------------------------------------------------
if [[ -f "${ROOT}/conf/tools.env.sh" ]]; then
  source "${ROOT}/conf/tools.env.sh"
fi

export PYTHONPATH="${ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}"
export PATH="${ROOT}/bin:${PATH}"

# VEP Perl compat (mirrors bin/neoag-nextflow logic)
export PATH="${ROOT}/bin:${NEOAG_CONDA_BASE:-${HOME}/miniforge3}/envs/${NEOAG_VEP_ENV:-neoag-vep}/bin:${PATH}"

# Keep Nextflow metadata out of repository root
export NXF_HOME="${NXF_HOME:-${ROOT}/work/.nextflow_home}"
mkdir -p "${NXF_HOME}"

# ---------------------------------------------------------------------------
# Default parameter values
# ---------------------------------------------------------------------------
SAMPLE_ID="${SAMPLE_ID:-}"
NORMAL_BAM="${NORMAL_BAM:-}"
TUMOR_BAM="${TUMOR_BAM:-}"
REFERENCE_FASTA="${REFERENCE_FASTA:-${NEOAG_REFERENCE_FASTA:-}}"
DBSNP_VCF="${DBSNP_VCF:-${NEOAG_DBSNP_VCF:-}}"
OUTDIR="${OUTDIR:-${ROOT}/results/all_qc}"
PROFILE_NAME="${PROFILE_NAME:-default}"
TUMOR_SAMPLE_NAME="${TUMOR_SAMPLE_NAME:-TUMOR}"
NORMAL_SAMPLE_NAME="${NORMAL_SAMPLE_NAME:-NORMAL}"
STRICT_MODE="${STRICT_MODE:-false}"
UPSTREAM_STUB="${UPSTREAM_STUB:-false}"
SKIP_QC="${SKIP_QC:-false}"
CONFIG_FILE="${CONFIG_FILE:-${ROOT}/conf/main_full.config}"
NEOAG_PROFILE="${NEOAG_PROFILE:-}"
RESUME_FLAG=""
DRY_RUN=false

# ---------------------------------------------------------------------------
# Help
# ---------------------------------------------------------------------------
usage() {
  cat <<'EOF'
Usage: bash scripts/run_main_all_qc.sh --normal_bam <bam> --tumor_bam <bam> --sample_id <id> [options]

Required:
  --normal_bam <path>        Normal blood/sample BAM
  --tumor_bam <path>         Tumor BAM
  --sample_id <id>           Sample identifier

Reference (auto-detected from env if not set):
  --reference_fasta <path>   Reference genome FASTA (default: $NEOAG_REFERENCE_FASTA)

QC (optional — QC subworkflows run by default):
  --dbsnp_vcf <path>         dbSNP/common SNP VCF for FACETS snp-pileup
  --skip_qc                  Skip LOH_CHECK and PURITY_CHECK subworkflows

Pipeline options:
  --outdir <dir>             Output directory (default: results/all_qc)
  --profile_name <name>      Scoring profile (default: default)
  --tumor_sample_name <str>  Tumor sample name in BAM (default: TUMOR)
  --normal_sample_name <str> Normal sample name in BAM (default: NORMAL)
  --strict_mode              Enable strict mode (default: false)
  --upstream_stub            Stub upstream tools (default: false)
  --config <path>            Nextflow config file (default: conf/main_full.config)

Nextflow options:
  --resume                   Resume previous run
  --with-report <path>       Write Nextflow execution report
  --dry-run                  Print the nextflow command without running it

Environment:
  NEOAG_RUNNER_MODE          "conda" (default) or "docker"
  NEOAG_REFERENCE_FASTA      Reference genome FASTA
  NEOAG_DBSNP_VCF            dbSNP/common SNP VCF
  NEOAG_PROFILE               Nextflow -profile value (e.g. "docker,slurm")
EOF
}

# ---------------------------------------------------------------------------
# Parse command-line arguments
# ---------------------------------------------------------------------------
declare -a EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample_id)            SAMPLE_ID="$2";            shift 2 ;;
    --normal_bam)           NORMAL_BAM="$2";           shift 2 ;;
    --tumor_bam)            TUMOR_BAM="$2";            shift 2 ;;
    --reference_fasta)      REFERENCE_FASTA="$2";      shift 2 ;;
    --dbsnp_vcf)            DBSNP_VCF="$2";            shift 2 ;;
    --outdir)               OUTDIR="$2";               shift 2 ;;
    --profile_name)         PROFILE_NAME="$2";         shift 2 ;;
    --tumor_sample_name)    TUMOR_SAMPLE_NAME="$2";    shift 2 ;;
    --normal_sample_name)   NORMAL_SAMPLE_NAME="$2";   shift 2 ;;
    --config)               CONFIG_FILE="$2";          shift 2 ;;
    --strict_mode)          STRICT_MODE=true;          shift ;;
    --upstream_stub)        UPSTREAM_STUB=true;        shift ;;
    --skip_qc)              SKIP_QC=true;              shift ;;
    --resume)               RESUME_FLAG="-resume";     shift ;;
    --with-report)          EXTRA_ARGS+=("-with-report" "$2"); shift 2 ;;
    --dry-run)              DRY_RUN=true;              shift ;;
    -h|--help)              usage; exit 0 ;;
    *)
      echo "ERROR: unknown option: $1" >&2
      usage
      exit 1
      ;;
  esac
done

# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------
fail() { echo "ERROR: $*" >&2; usage >&2; exit 1; }

[[ -n "${SAMPLE_ID}" ]]   || fail "--sample_id is required"
[[ -n "${NORMAL_BAM}" ]]  || fail "--normal_bam is required"
[[ -n "${TUMOR_BAM}" ]]   || fail "--tumor_bam is required"

for f in "${NORMAL_BAM}" "${TUMOR_BAM}"; do
  [[ -f "${f}" ]] || fail "BAM not found: ${f}"
done

if [[ -n "${REFERENCE_FASTA}" ]]; then
  [[ -f "${REFERENCE_FASTA}" ]] || fail "REFERENCE_FASTA not found: ${REFERENCE_FASTA}"
else
  fail "REFERENCE_FASTA not set. Use --reference_fasta or set NEOAG_REFERENCE_FASTA."
fi

if [[ "${SKIP_QC}" != "true" ]] && [[ -n "${DBSNP_VCF}" ]]; then
  [[ -f "${DBSNP_VCF}" ]] || echo "WARNING: dbsnp_vcf not found: ${DBSNP_VCF} (FACETS will produce stub output)"
fi

if [[ ! -f "${CONFIG_FILE}" ]]; then
  fail "Config file not found: ${CONFIG_FILE}"
fi

# ---------------------------------------------------------------------------
# Build Nextflow command
# ---------------------------------------------------------------------------
declare -a NF_ARGS=(
  run
  "${ROOT}/workflows/main_all_qc.nf"
  --sample_id "${SAMPLE_ID}"
  --normal_bam "${NORMAL_BAM}"
  --tumor_bam "${TUMOR_BAM}"
  --reference_fasta "${REFERENCE_FASTA}"
  --outdir "${OUTDIR}"
  --profile_name "${PROFILE_NAME}"
  --tumor_sample_name "${TUMOR_SAMPLE_NAME}"
  --normal_sample_name "${NORMAL_SAMPLE_NAME}"
  -c "${CONFIG_FILE}"
)

[[ -n "${DBSNP_VCF}" ]]   && NF_ARGS+=("--dbsnp_vcf" "${DBSNP_VCF}")
[[ "${STRICT_MODE}"   == "true" ]] && NF_ARGS+=("--strict_mode" "true")
[[ "${UPSTREAM_STUB}" == "true" ]] && NF_ARGS+=("--upstream_stub" "true")
[[ "${SKIP_QC}"       == "true" ]] && NF_ARGS+=("--skip_qc" "true")
[[ -n "${RESUME_FLAG}" ]]          && NF_ARGS+=("${RESUME_FLAG}")

if [[ -n "${NEOAG_PROFILE}" ]]; then
  NF_ARGS+=("-profile" "${NEOAG_PROFILE}")
fi

NF_ARGS+=("${EXTRA_ARGS[@]}")

# ---------------------------------------------------------------------------
# Runtime summary
# ---------------------------------------------------------------------------
RUNNER_MODE="${NEOAG_RUNNER_MODE:-conda}"
echo "============================================================"
echo "  NeoAg main_all_qc — Nextflow Pipeline Launcher"
echo "============================================================"
echo "  sample_id           = ${SAMPLE_ID}"
echo "  normal_bam          = ${NORMAL_BAM}"
echo "  tumor_bam           = ${TUMOR_BAM}"
echo "  reference_fasta     = ${REFERENCE_FASTA}"
echo "  dbsnp_vcf           = ${DBSNP_VCF:-'(not set — FACETS stub)'}"
echo "  outdir              = ${OUTDIR}"
echo "  profile_name        = ${PROFILE_NAME}"
echo "  runner_mode         = ${RUNNER_MODE}"
echo "  skip_qc             = ${SKIP_QC}"
echo "  strict_mode         = ${STRICT_MODE}"
echo "  upstream_stub       = ${UPSTREAM_STUB}"
echo "  config              = ${CONFIG_FILE}"
echo "  nextflow_profile    = ${NEOAG_PROFILE:-(none)}"
echo "  resume              = ${RESUME_FLAG:-(no)}"
echo "============================================================"
echo ""

if [[ "${DRY_RUN}" == "true" ]]; then
  echo "[DRY-RUN] Nextflow command:"
  echo "  cd ${ROOT}"
  echo "  bin/neoag-nextflow ${NF_ARGS[*]}"
  echo ""
  echo "[DRY-RUN] Complete. Remove --dry-run to execute."
  exit 0
fi

# ---------------------------------------------------------------------------
# Launch
# ---------------------------------------------------------------------------
echo "[$(date -Is)] Launching main_all_qc workflow ..."
echo ""

"${ROOT}/bin/neoag-nextflow" "${NF_ARGS[@]}"

RC=$?
echo ""
echo "[$(date -Is)] main_all_qc finished (exit=${RC})"

if [[ ${RC} -eq 0 ]]; then
  echo ""
  echo "Output files:"
  echo "  Scoring:       ${OUTDIR}/"
  echo "  HLA typing:    ${OUTDIR}/hla_typing/"
  if [[ "${SKIP_QC}" != "true" ]]; then
    echo "  QC results:    ${OUTDIR}/qc/"
    echo "    LOHHLA:      ${OUTDIR}/qc/lohhla/"
    echo "    SpecHLA:     ${OUTDIR}/qc/spechla/"
    echo "    FACETS:      ${OUTDIR}/qc/facets/"
    echo "    PURPLE:      ${OUTDIR}/qc/purple/"
  fi
fi

exit ${RC}
