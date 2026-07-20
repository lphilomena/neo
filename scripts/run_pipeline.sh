#!/usr/bin/env bash
# ============================================================================
# run_pipeline.sh — Unified NeoAg Nextflow Pipeline Launcher
#
# Single entry-point for all 5 core neoantigen prediction workflows.
# Handles parameter validation, environment setup, and Nextflow invocation.
#
# Workflows supported:
#   main_fromVCF           TOML config → upstream + scoring
#   main_fromVCF_nohla     TOML + BAM/FASTQ → OptiType HLA → scoring
#   main_all               BAM×2 → OptiType HLA + Mutect2 + scoring
#   main_all_qc            BAM×2 → OptiType HLA + Mutect2 + scoring + QC
#   main_all_nohla         BAM×2 + manual HLA → Mutect2 + scoring
#
# Usage:
#   # TOML-based (pre-computed VCF + HLA)
#   bash scripts/run_pipeline.sh --workflow main_fromVCF \
#     --run_config conf/run.mycase.toml --sample_id SAMPLE001
#
#   # BAM-based with auto HLA + QC
#   bash scripts/run_pipeline.sh --workflow main_all_qc \
#     --normal_bam normal.bam --tumor_bam tumor.bam \
#     --sample_id SAMPLE001 --reference_fasta ref.fa \
#     --dbsnp_vcf dbsnp.vcf.gz
#
#   # BAM-based with manual HLA
#   bash scripts/run_pipeline.sh --workflow main_all_nohla \
#     --normal_bam normal.bam --tumor_bam tumor.bam \
#     --hla_alleles "HLA-A*02:01,HLA-B*07:02,HLA-C*07:02" \
#     --sample_id SAMPLE001 --reference_fasta ref.fa
#
#   # Dry-run (print command only)
#   bash scripts/run_pipeline.sh --dry-run ...
# ============================================================================

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "${ROOT}"

# ---------------------------------------------------------------------------
# Environment setup (mirrors bin/neoag-nextflow logic)
# ---------------------------------------------------------------------------
if [[ -f "${ROOT}/conf/tools.env.sh" ]]; then
  source "${ROOT}/conf/tools.env.sh"
fi

export PYTHONPATH="${ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}"
export PATH="${ROOT}/bin:${PATH}"

# VEP Perl compat
export PATH="${ROOT}/bin:${NEOAG_CONDA_BASE:-${HOME}/miniforge3}/envs/${NEOAG_VEP_ENV:-neoag-vep}/bin:${PATH}"

# Keep Nextflow metadata out of repository root
export NXF_HOME="${NXF_HOME:-${ROOT}/work/.nextflow_home}"
mkdir -p "${NXF_HOME}"

# ---------------------------------------------------------------------------
# Default parameter values
# ---------------------------------------------------------------------------
WORKFLOW="${WORKFLOW:-}"
SAMPLE_ID="${SAMPLE_ID:-}"
NORMAL_BAM="${NORMAL_BAM:-}"
TUMOR_BAM="${TUMOR_BAM:-}"
INPUT_BAM="${INPUT_BAM:-}"
INPUT_FQ1="${INPUT_FQ1:-}"
INPUT_FQ2="${INPUT_FQ2:-}"
RUN_CONFIG="${RUN_CONFIG:-}"
HLA_ALLELES="${HLA_ALLELES:-}"
REFERENCE_FASTA="${REFERENCE_FASTA:-${NEOAG_REFERENCE_FASTA:-}}"
DBSNP_VCF="${DBSNP_VCF:-${NEOAG_DBSNP_VCF:-}}"
OUTDIR="${OUTDIR:-${ROOT}/results/pipeline}"
PROFILE_NAME="${PROFILE_NAME:-default}"
TUMOR_SAMPLE_NAME="${TUMOR_SAMPLE_NAME:-TUMOR}"
NORMAL_SAMPLE_NAME="${NORMAL_SAMPLE_NAME:-NORMAL}"
SEQ_TYPE="${SEQ_TYPE:-dna}"
CONFIG_FILE="${CONFIG_FILE:-${ROOT}/conf/main_full.config}"
NEOAG_PROFILE="${NEOAG_PROFILE:-}"
SKIP_QC="${SKIP_QC:-false}"
SKIP_HLA_TYPING="${SKIP_HLA_TYPING:-false}"
STRICT_MODE="${STRICT_MODE:-false}"
UPSTREAM_STUB="${UPSTREAM_STUB:-false}"
RESUME_FLAG=""
DRY_RUN=false

# ---------------------------------------------------------------------------
# Valid workflow names
# ---------------------------------------------------------------------------
VALID_WORKFLOWS=(
  main_fromVCF
  main_fromVCF_nohla
  main_all
  main_all_qc
  main_all_nohla
)

# ---------------------------------------------------------------------------
# Help
# ---------------------------------------------------------------------------
usage() {
  cat <<'EOF'
Usage: bash scripts/run_pipeline.sh --workflow <name> [options]

Workflows (--workflow):
  main_fromVCF           TOML config with pre-computed VCF + HLA
  main_fromVCF_nohla     TOML config + BAM/FASTQ for auto HLA typing
  main_all               BAM×2 → auto HLA typing + variant calling + scoring
  main_all_qc            BAM×2 → auto HLA + variant calling + scoring + QC
  main_all_nohla         BAM×2 + manual HLA alleles → variant calling + scoring

Required by workflow:

  main_fromVCF:
    --run_config <path>       TOML config file
    --sample_id <id>          Sample identifier

  main_fromVCF_nohla:
    --run_config <path>       TOML config file
    --input_bam <path>        BAM for HLA typing (or --input_fq1/--input_fq2)
    --sample_id <id>          Sample identifier

  main_all / main_all_qc / main_all_nohla:
    --normal_bam <path>       Normal/blood BAM
    --tumor_bam <path>        Tumor BAM
    --sample_id <id>          Sample identifier
    --reference_fasta <path>  Reference genome FASTA

  main_all_nohla (additionally):
    --hla_alleles <str>       Comma-separated HLA, e.g. "HLA-A*02:01,HLA-B*07:02"

  main_all_qc (optional):
    --dbsnp_vcf <path>        dbSNP VCF for FACETS snp-pileup
    --skip_qc                 Skip QC subworkflows (keep main pipeline)

General options:
  --outdir <dir>              Output directory (default: results/pipeline)
  --profile_name <name>       Scoring profile (default: default)
  --tumor_sample_name <str>   Tumor read-group name in BAM (default: TUMOR)
  --normal_sample_name <str>  Normal read-group name in BAM (default: NORMAL)
  --seq_type <str>            Sequence type for HLA typing: dna|rna (default: dna)
  --config <path>             Nextflow config file (default: conf/main_full.config)
  --strict_mode               Enable strict mode
  --upstream_stub             Run upstream in stub mode

Nextflow options:
  --resume                    Resume previous run
  --with-report <path>        Write Nextflow execution report
  --dry-run                   Print the nextflow command without executing

Environment:
  NEOAG_RUNNER_MODE           "conda" (default) or "docker"
  NEOAG_REFERENCE_FASTA       Reference genome FASTA
  NEOAG_DBSNP_VCF             dbSNP/common SNP VCF for FACETS
  NEOAG_PROFILE               Nextflow -profile value (e.g. "docker,slurm")
EOF
}

# ---------------------------------------------------------------------------
# Parse command-line arguments
# ---------------------------------------------------------------------------
declare -a EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --workflow)              WORKFLOW="$2";              shift 2 ;;
    --sample_id)             SAMPLE_ID="$2";             shift 2 ;;
    --normal_bam)            NORMAL_BAM="$2";            shift 2 ;;
    --tumor_bam)             TUMOR_BAM="$2";             shift 2 ;;
    --input_bam)             INPUT_BAM="$2";             shift 2 ;;
    --input_fq1)             INPUT_FQ1="$2";             shift 2 ;;
    --input_fq2)             INPUT_FQ2="$2";             shift 2 ;;
    --run_config)            RUN_CONFIG="$2";            shift 2 ;;
    --hla_alleles)           HLA_ALLELES="$2";           shift 2 ;;
    --reference_fasta)       REFERENCE_FASTA="$2";       shift 2 ;;
    --dbsnp_vcf)             DBSNP_VCF="$2";             shift 2 ;;
    --outdir)                OUTDIR="$2";                shift 2 ;;
    --profile_name)          PROFILE_NAME="$2";          shift 2 ;;
    --tumor_sample_name)     TUMOR_SAMPLE_NAME="$2";     shift 2 ;;
    --normal_sample_name)    NORMAL_SAMPLE_NAME="$2";    shift 2 ;;
    --seq_type)              SEQ_TYPE="$2";              shift 2 ;;
    --config)                CONFIG_FILE="$2";           shift 2 ;;
    --strict_mode)           STRICT_MODE=true;           shift ;;
    --upstream_stub)         UPSTREAM_STUB=true;         shift ;;
    --skip_qc)               SKIP_QC=true;               shift ;;
    --skip_hla_typing)       SKIP_HLA_TYPING=true;       shift ;;
    --resume)                RESUME_FLAG="-resume";      shift ;;
    --with-report)           EXTRA_ARGS+=("-with-report" "$2"); shift 2 ;;
    --dry-run)               DRY_RUN=true;               shift ;;
    -h|--help)               usage; exit 0 ;;
    *)
      echo "ERROR: unknown option: $1" >&2
      usage
      exit 1
      ;;
  esac
done

# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------
fail() { echo "ERROR: $*" >&2; usage >&2; exit 1; }
warn() { echo "WARNING: $*" >&2; }

is_valid_workflow() {
  for w in "${VALID_WORKFLOWS[@]}"; do
    [[ "$w" == "$1" ]] && return 0
  done
  return 1
}

check_file() {
  local desc="$1" path="$2"
  if [[ ! -f "${path}" ]]; then
    fail "${desc} not found: ${path}"
  fi
}

# ---------------------------------------------------------------------------
# Validate workflow selection
# ---------------------------------------------------------------------------
[[ -n "${WORKFLOW}" ]] || fail "--workflow is required. Choose from: ${VALID_WORKFLOWS[*]}"
is_valid_workflow "${WORKFLOW}" || fail "Invalid workflow '${WORKFLOW}'. Choose from: ${VALID_WORKFLOWS[*]}"

# ---------------------------------------------------------------------------
# Workflow-specific parameter validation
# ---------------------------------------------------------------------------
declare -a NF_ARGS=(
  run
  "${ROOT}/workflows/${WORKFLOW}.nf"
)

case "${WORKFLOW}" in
  main_fromVCF)
    [[ -n "${RUN_CONFIG}" ]] || fail "--run_config is required for ${WORKFLOW}"
    [[ -n "${SAMPLE_ID}" ]]  || fail "--sample_id is required for ${WORKFLOW}"
    check_file "Run config" "${RUN_CONFIG}"

    NF_ARGS+=(
      --run_config "${RUN_CONFIG}"
      --sample_id "${SAMPLE_ID}"
      --outdir "${OUTDIR}"
      --profile_name "${PROFILE_NAME}"
      -c "${CONFIG_FILE}"
    )
    ;;

  main_fromVCF_nohla)
    [[ -n "${RUN_CONFIG}" ]] || fail "--run_config is required for ${WORKFLOW}"
    [[ -n "${SAMPLE_ID}" ]]  || fail "--sample_id is required for ${WORKFLOW}"
    check_file "Run config" "${RUN_CONFIG}"

    if [[ "${SKIP_HLA_TYPING}" != "true" ]]; then
      if [[ -z "${INPUT_BAM}" && -z "${INPUT_FQ1}" ]]; then
        fail "HLA typing requires --input_bam or --input_fq1. Use --skip_hla_typing if HLA is pre-configured in TOML."
      fi
      [[ -n "${INPUT_BAM}" ]] && check_file "Input BAM" "${INPUT_BAM}"
      [[ -n "${INPUT_FQ1}" ]] && check_file "Input FASTQ R1" "${INPUT_FQ1}"
      [[ -n "${INPUT_FQ2}" ]] && check_file "Input FASTQ R2" "${INPUT_FQ2}"
    fi

    NF_ARGS+=(
      --run_config "${RUN_CONFIG}"
      --sample_id "${SAMPLE_ID}"
      --outdir "${OUTDIR}"
      --seq_type "${SEQ_TYPE}"
      -c "${CONFIG_FILE}"
    )

    [[ -n "${INPUT_BAM}" ]] && NF_ARGS+=(--input_bam "${INPUT_BAM}")
    [[ -n "${INPUT_FQ1}" ]] && NF_ARGS+=(--input_fq1 "${INPUT_FQ1}")
    [[ -n "${INPUT_FQ2}" ]] && NF_ARGS+=(--input_fq2 "${INPUT_FQ2}")
    [[ "${SKIP_HLA_TYPING}" == "true" ]] && NF_ARGS+=(--skip_hla_typing)
    ;;

  main_all | main_all_qc | main_all_nohla)
    [[ -n "${SAMPLE_ID}" ]]      || fail "--sample_id is required"
    [[ -n "${NORMAL_BAM}" ]]     || fail "--normal_bam is required"
    [[ -n "${TUMOR_BAM}" ]]      || fail "--tumor_bam is required"
    check_file "Normal BAM" "${NORMAL_BAM}"
    check_file "Tumor BAM"  "${TUMOR_BAM}"

    if [[ -n "${REFERENCE_FASTA}" ]]; then
      check_file "Reference FASTA" "${REFERENCE_FASTA}"
      # Check for .fai
      local fai="${REFERENCE_FASTA}.fai"
      [[ -f "${fai}" ]] || warn "Reference index not found: ${fai}"
      # Check for .dict
      local fasta_base="${REFERENCE_FASTA}"
      local dict=""
      if [[ "${fasta_base}" == *.chr.fa ]]; then
        dict="${fasta_base%.chr.fa}.dict"
      elif [[ "${fasta_base}" == *.chr.fasta ]]; then
        dict="${fasta_base%.chr.fasta}.dict"
      elif [[ "${fasta_base}" == *.fa ]]; then
        dict="${fasta_base%.fa}.dict"
      elif [[ "${fasta_base}" == *.fasta ]]; then
        dict="${fasta_base%.fasta}.dict"
      fi
      [[ -n "${dict}" && -f "${dict}" ]] || warn "Reference dict not found: ${dict:-N/A}"
    else
      fail "REFERENCE_FASTA not set. Use --reference_fasta or set NEOAG_REFERENCE_FASTA."
    fi

    NF_ARGS+=(
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

    # Workflow-specific extras
    case "${WORKFLOW}" in
      main_all_nohla)
        [[ -n "${HLA_ALLELES}" ]] || fail "--hla_alleles is required for main_all_nohla (e.g. --hla_alleles \"HLA-A*02:01,HLA-B*07:02,HLA-C*07:02\")"
        NF_ARGS+=(--hla_alleles "${HLA_ALLELES}")
        ;;
      main_all_qc)
        if [[ -n "${DBSNP_VCF}" ]]; then
          check_file "dbSNP VCF" "${DBSNP_VCF}"
          NF_ARGS+=(--dbsnp_vcf "${DBSNP_VCF}")
        else
          warn "dbsnp_vcf not set — FACETS will produce stub output"
        fi
        [[ "${SKIP_QC}" == "true" ]] && NF_ARGS+=(--skip_qc true)
        ;;
    esac
    ;;
esac

# ---------------------------------------------------------------------------
# Common flags
# ---------------------------------------------------------------------------
[[ "${STRICT_MODE}"   == "true" ]] && NF_ARGS+=("--strict_mode" "true")
[[ "${UPSTREAM_STUB}" == "true" ]] && NF_ARGS+=("--upstream_stub" "true")
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
echo "  NeoAg run_pipeline — ${WORKFLOW}"
echo "============================================================"
echo "  sample_id           = ${SAMPLE_ID:-N/A}"
echo "  workflow            = ${WORKFLOW}"
echo "  runner_mode         = ${RUNNER_MODE}"
echo "  outdir              = ${OUTDIR}"
echo "  profile_name        = ${PROFILE_NAME}"
echo "  config              = ${CONFIG_FILE}"
echo "  nextflow_profile    = ${NEOAG_PROFILE:-(none)}"
echo "  resume              = ${RESUME_FLAG:-(no)}"
echo "  dry_run             = ${DRY_RUN}"

case "${WORKFLOW}" in
  main_fromVCF)
    echo "  run_config          = ${RUN_CONFIG}"
    ;;
  main_fromVCF_nohla)
    echo "  run_config          = ${RUN_CONFIG}"
    echo "  input_bam           = ${INPUT_BAM:-(not set)}"
    echo "  input_fq1           = ${INPUT_FQ1:-(not set)}"
    echo "  skip_hla_typing     = ${SKIP_HLA_TYPING}"
    echo "  seq_type            = ${SEQ_TYPE}"
    ;;
  main_all | main_all_qc | main_all_nohla)
    echo "  normal_bam          = ${NORMAL_BAM}"
    echo "  tumor_bam           = ${TUMOR_BAM}"
    echo "  reference_fasta     = ${REFERENCE_FASTA}"
    echo "  tumor_sample_name   = ${TUMOR_SAMPLE_NAME}"
    echo "  normal_sample_name  = ${NORMAL_SAMPLE_NAME}"
    case "${WORKFLOW}" in
      main_all_qc)
        echo "  dbsnp_vcf           = ${DBSNP_VCF:-(not set — FACETS stub)}"
        echo "  skip_qc             = ${SKIP_QC}"
        ;;
      main_all_nohla)
        echo "  hla_alleles         = ${HLA_ALLELES}"
        ;;
    esac
    ;;
esac

echo "  strict_mode         = ${STRICT_MODE}"
echo "  upstream_stub       = ${UPSTREAM_STUB}"
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
echo "[$(date -Is)] Launching ${WORKFLOW} ..."
echo ""

"${ROOT}/bin/neoag-nextflow" "${NF_ARGS[@]}"

RC=$?
echo ""
echo "[$(date -Is)] ${WORKFLOW} finished (exit=${RC})"

if [[ ${RC} -eq 0 ]]; then
  echo ""
  echo "Output files:"
  echo "  Scoring results:  ${OUTDIR}/"
  case "${WORKFLOW}" in
    main_all | main_all_qc)
      echo "  HLA typing:       ${OUTDIR}/hla_typing/"
      ;;
  esac
  case "${WORKFLOW}" in
    main_all_qc)
      if [[ "${SKIP_QC}" != "true" ]]; then
        echo "  QC results:       ${OUTDIR}/qc/"
        echo "    LOHHLA:         ${OUTDIR}/qc/lohhla/"
        echo "    SpecHLA:        ${OUTDIR}/qc/spechla/"
        echo "    FACETS:         ${OUTDIR}/qc/facets/"
        echo "    PURPLE:         ${OUTDIR}/qc/purple/"
      fi
      ;;
  esac
fi

exit ${RC}
