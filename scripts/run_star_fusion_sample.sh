#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 --fastq1 R1.fq.gz --fastq2 R2.fq.gz --ctat-genome-lib DIR --sample-id ID --outdir OUT [--threads N]" >&2
}

FQ1=""; FQ2=""; CTAT=""; SAMPLE_ID="sample"; OUTDIR=""; THREADS=16
while [[ $# -gt 0 ]]; do
  case "$1" in
    --fastq1) FQ1="$2"; shift 2 ;;
    --fastq2) FQ2="$2"; shift 2 ;;
    --ctat-genome-lib) CTAT="$2"; shift 2 ;;
    --sample-id) SAMPLE_ID="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done
[[ -s "$FQ1" && -s "$FQ2" ]] || { echo "ERROR: paired FASTQ files are required" >&2; exit 2; }
[[ -d "$CTAT" ]] || { echo "ERROR: CTAT genome library missing: $CTAT" >&2; exit 3; }
STAR_FUSION_BIN="${NEOAG_STAR_FUSION_BIN:-$(command -v STAR-Fusion || true)}"
[[ -n "$STAR_FUSION_BIN" && -x "$STAR_FUSION_BIN" ]] || { echo "ERROR: STAR-Fusion executable not found" >&2; exit 3; }
mkdir -p "$OUTDIR"
"$STAR_FUSION_BIN" --left_fq "$FQ1" --right_fq "$FQ2" \
  --genome_lib_dir "$CTAT" --CPU "$THREADS" --output_dir "$OUTDIR" \
  >"$OUTDIR/star-fusion.log" 2>&1
test -s "$OUTDIR/star-fusion.fusion_predictions.tsv"
