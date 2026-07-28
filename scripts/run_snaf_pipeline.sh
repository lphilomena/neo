#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 --bam BAM --hla-file FILE --sample-id ID --db-dir DIR --outdir DIR [--threads N]" >&2
}

BAM=""; HLA=""; SAMPLE=""; DB=""; OUT=""; THREADS=8
IMAGE="${NEOAG_ALTANALYZE_IMAGE:-neoag-altanalyze:snaf}"
SNAF_PY="${SNAF_PYTHON:-python}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam) BAM="$2"; shift 2;;
    --hla-file) HLA="$2"; shift 2;;
    --sample-id) SAMPLE="$2"; shift 2;;
    --db-dir) DB="$2"; shift 2;;
    --outdir) OUT="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --altanalyze-image) IMAGE="$2"; shift 2;;
    *) usage; exit 2;;
  esac
done
for value in "$BAM" "$HLA" "$SAMPLE" "$DB" "$OUT"; do [[ -n "$value" ]] || { usage; exit 2; }; done
[[ -s "$BAM" && -s "$HLA" ]] || { echo "ERROR: missing BAM or HLA file" >&2; exit 2; }
[[ -s "$DB/controls/GTEx_junction_counts.h5ad" ]] || { echo "ERROR: incomplete SNAF database: $DB" >&2; exit 2; }
command -v docker >/dev/null || { echo "ERROR: docker is required for AltAnalyze" >&2; exit 2; }
docker image inspect "$IMAGE" >/dev/null 2>&1 || { echo "ERROR: missing AltAnalyze image: $IMAGE" >&2; exit 2; }

WORK="$OUT/altanalyze_work"
mkdir -p "$WORK/bam" "$OUT"
mounts=(
  -v "$WORK:/mnt"
  -v "$BAM:/mnt/bam/${SAMPLE}.bam:ro"
  -v "$BAM:/mnt/bam/${SAMPLE}_replicate.bam:ro"
)
if [[ -s "${BAM}.bai" ]]; then
  mounts+=(
    -v "${BAM}.bai:/mnt/bam/${SAMPLE}.bam.bai:ro"
    -v "${BAM}.bai:/mnt/bam/${SAMPLE}_replicate.bam.bai:ro"
  )
fi
docker run --rm "${mounts[@]}" "$IMAGE" identify bam "$THREADS"
MATRIX="$WORK/altanalyze_output/ExpressionInput/counts.original.pruned.txt"
[[ -s "$MATRIX" ]] || { echo "ERROR: AltAnalyze matrix was not created" >&2; exit 1; }

export NEOAG_SNAF_OUTDIR="$OUT"
export NEOAG_SNAF_MATRIX="$MATRIX"
export NEOAG_SNAF_DB="$DB"
export NEOAG_SNAF_HLA_FILE="$HLA"
export NEOAG_SNAF_SAMPLE_ID="$SAMPLE"
export NEOAG_SNAF_CORES="$THREADS"
"$SNAF_PY" "$(dirname "$0")/snaf_sample_workflow.py"
test -s "$OUT/snaf_candidates.tsv"
