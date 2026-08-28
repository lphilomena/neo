#!/usr/bin/env bash
# SNAF wrapper for abnormal-splicing neoantigen discovery.
# SNAF workflows are cohort-scale and reference-heavy; this script standardizes
# environment checks and provides a controlled place to run a user-supplied SNAF
# Python workflow script.
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: bash scripts/run_snaf_sample.sh --workflow workflow.py --outdir OUT [--bam-dir DIR] [--hla-file HLA.tsv] [--sample-id ID]

Install first:
  bash scripts/install_splice_tools.sh
USAGE
}

WORKFLOW=""; OUTDIR=""; BAM_DIR=""; HLA_FILE=""; SAMPLE_ID="sample"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --workflow) WORKFLOW="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --bam-dir) BAM_DIR="$2"; shift 2 ;;
    --hla-file) HLA_FILE="$2"; shift 2 ;;
    --sample-id) SAMPLE_ID="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done
[[ -n "$WORKFLOW" && -s "$WORKFLOW" ]] || { echo "ERROR: --workflow required" >&2; exit 2; }
[[ -n "$OUTDIR" ]] || { echo "ERROR: --outdir required" >&2; exit 2; }
SNAF_PY="${SNAF_PYTHON:-}"
if [[ -z "$SNAF_PY" ]]; then
  if [[ -n "${NEOAG_CONDA_BASE:-}" ]]; then
    SNAF_PY="${NEOAG_CONDA_BASE}/bin/conda run -n ${NEOAG_SNAF_ENV:-neoag-snaf} python"
  else
    SNAF_PY="python"
  fi
fi
mkdir -p "$OUTDIR"
export NEOAG_SNAF_OUTDIR="$OUTDIR"
export NEOAG_SNAF_BAM_DIR="$BAM_DIR"
export NEOAG_SNAF_HLA_FILE="$HLA_FILE"
export NEOAG_SNAF_SAMPLE_ID="$SAMPLE_ID"
$SNAF_PY - <<PY
import snaf
print("SNAF import OK")
PY
$SNAF_PY "$WORKFLOW"
echo "SNAF workflow output: $OUTDIR"
