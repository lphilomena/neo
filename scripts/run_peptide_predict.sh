#!/usr/bin/env bash
# Standalone peptide–HLA prediction from flexible CSV/TSV input.
#
# Migrated from neoantigen2 run_peptide_predict.sh, with strict peptide–HLA pair
# deduplication (same peptide + different HLA alleles are kept as separate pairs).
#
# Usage:
#   bash scripts/run_peptide_predict.sh -i peptides.csv [-o results_dir] [options]
#
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "${ROOT}"

if [[ -f "${ROOT}/conf/tools.env.sh" ]]; then
  # shellcheck disable=SC1091
  source "${ROOT}/conf/tools.env.sh"
fi

INPUT=""
OUTDIR="${ROOT}/results/peptide_predict"
SAMPLE_ID="SAMPLE001"
PROFILE="default"
STUB=false
SKIP_NETMHCPAN=false
SKIP_MHCFLURRY=false
SKIP_PRIME=false
SKIP_BIGMHC=false
SKIP_STABPAN=false

usage() {
  cat <<EOF
Usage: $(basename "$0") -i <input.csv|tsv> [options]

Required:
  -i, --input <file>        Input table with peptide + HLA columns

Options:
  -o, --outdir <dir>        Output directory (default: results/peptide_predict)
  --sample-id <id>          Sample identifier (default: SAMPLE001)
  --profile <name>          Scoring profile (default: default)
  --stub                    Fast stub predictors
  --skip-netmhcpan          Skip NetMHCpan
  --skip-mhcflurry          Skip MHCflurry
  --skip-prime              Skip PRIME
  --skip-bigmhc-im          Skip BigMHC_IM
  --skip-stabpan            Skip NetMHCstabpan
  -h, --help                Show help

Supported peptide columns: peptide, peptide_seq, seq, sequence, epitope, mer
Supported HLA columns: hla, hla_allele, allele, mhc, hla_type
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input) INPUT="$2"; shift 2 ;;
    -o|--outdir) OUTDIR="$2"; shift 2 ;;
    --sample-id) SAMPLE_ID="$2"; shift 2 ;;
    --profile) PROFILE="$2"; shift 2 ;;
    --stub) STUB=true; shift ;;
    --skip-netmhcpan) SKIP_NETMHCPAN=true; shift ;;
    --skip-mhcflurry) SKIP_MHCFLURRY=true; shift ;;
    --skip-prime) SKIP_PRIME=true; shift ;;
    --skip-bigmhc-im) SKIP_BIGMHC=true; shift ;;
    --skip-stabpan) SKIP_STABPAN=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "${INPUT}" ]]; then
  echo "ERROR: -i/--input is required" >&2
  usage
  exit 1
fi
if [[ ! -f "${INPUT}" ]]; then
  echo "ERROR: input file not found: ${INPUT}" >&2
  exit 1
fi

ARGS=(
  peptide-predict
  -i "${INPUT}"
  -o "${OUTDIR}"
  --sample-id "${SAMPLE_ID}"
  --profile "${PROFILE}"
)
[[ "${STUB}" == true ]] && ARGS+=(--stub)
[[ "${SKIP_NETMHCPAN}" == true ]] && ARGS+=(--skip-netmhcpan)
[[ "${SKIP_MHCFLURRY}" == true ]] && ARGS+=(--skip-mhcflurry)
[[ "${SKIP_PRIME}" == true ]] && ARGS+=(--skip-prime)
[[ "${SKIP_BIGMHC}" == true ]] && ARGS+=(--skip-bigmhc-im)
[[ "${SKIP_STABPAN}" == true ]] && ARGS+=(--skip-stabpan)

if command -v neoag-v03 >/dev/null 2>&1; then
  neoag-v03 "${ARGS[@]}"
elif [[ -x "${ROOT}/bin/neoag-v03" ]]; then
  "${ROOT}/bin/neoag-v03" "${ARGS[@]}"
else
  python3 -m neoag_v03.cli "${ARGS[@]}"
fi
