#!/usr/bin/env bash
# Standalone DeepImmuno batch prediction on flexible peptide CSV/TSV input.
#
# Usage:
#   bash scripts/run_deepimmuno.sh -i pairs.csv -o results/deepimmuno.tsv
#
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "${ROOT}"

if [[ -f "${ROOT}/conf/tools.env.sh" ]]; then
  # shellcheck disable=SC1091
  source "${ROOT}/conf/tools.env.sh"
fi

INPUT=""
OUT=""
SAMPLE_ID="SAMPLE001"
STUB=false

usage() {
  cat <<EOF
Usage: $(basename "$0") -i <input.csv|tsv> -o <output.tsv> [--sample-id ID] [--stub]
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input) INPUT="$2"; shift 2 ;;
    -o|--output) OUT="$2"; shift 2 ;;
    --sample-id) SAMPLE_ID="$2"; shift 2 ;;
    --stub) STUB=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "${INPUT}" || -z "${OUT}" ]]; then
  echo "ERROR: -i and -o are required" >&2
  usage
  exit 1
fi

WORKDIR="$(dirname "${OUT}")/deepimmuno_work"
mkdir -p "${WORKDIR}"

ARGS=(
  convert-peptide-input
  -i "${INPUT}"
  -o "${WORKDIR}"
  --sample-id "${SAMPLE_ID}"
)
if command -v neoag >/dev/null 2>&1; then
  neoag "${ARGS[@]}"
else
  python3 -m neoag.cli "${ARGS[@]}"
fi

RAW="${WORKDIR}/00_input/raw_peptides.tsv"
RUN_ARGS=(
  run-tool deepimmuno
  --raw-peptides "${RAW}"
  --output "${OUT}"
  --workdir "${WORKDIR}"
  --sample-id "${SAMPLE_ID}"
)
[[ "${STUB}" == true ]] && RUN_ARGS+=(--stub)

if command -v neoag >/dev/null 2>&1; then
  neoag "${RUN_ARGS[@]}"
else
  python3 -m neoag.cli "${RUN_ARGS[@]}"
fi

echo "Wrote DeepImmuno evidence to ${OUT}"
