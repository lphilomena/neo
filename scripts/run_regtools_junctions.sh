#!/usr/bin/env bash
# Extract splice junctions from an RNA BAM with RegTools.
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: bash scripts/run_regtools_junctions.sh --bam RNA.bam --out OUT.tsv [--sample-id ID] [--strandness XS|RF|FR] [--junction-anchor N] [--min-intron N] [--max-intron N]

Environment fallback:
  NEOAG_REGTOOLS_BIN   regtools wrapper/executable
USAGE
}

BAM=""; OUT=""; SAMPLE_ID="sample"; STRANDNESS="${REGTOOLS_STRANDNESS:-XS}"; ANCHOR="${REGTOOLS_JUNCTION_ANCHOR:-8}"; MIN_INTRON="${REGTOOLS_MIN_INTRON:-20}"; MAX_INTRON="${REGTOOLS_MAX_INTRON:-500000}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam) BAM="$2"; shift 2 ;;
    --out) OUT="$2"; shift 2 ;;
    --sample-id) SAMPLE_ID="$2"; shift 2 ;;
    --strandness) STRANDNESS="$2"; shift 2 ;;
    --junction-anchor) ANCHOR="$2"; shift 2 ;;
    --min-intron) MIN_INTRON="$2"; shift 2 ;;
    --max-intron) MAX_INTRON="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done
[[ -n "$BAM" && -s "$BAM" ]] || { echo "ERROR: --bam missing or not found: $BAM" >&2; exit 2; }
[[ -n "$OUT" ]] || { echo "ERROR: --out required" >&2; exit 2; }
REGTOOLS_BIN="${NEOAG_REGTOOLS_BIN:-$(command -v regtools-neoag || command -v regtools || true)}"
[[ -n "$REGTOOLS_BIN" && -x "$REGTOOLS_BIN" ]] || { echo "ERROR: regtools not found; run scripts/install_splice_tools.sh" >&2; exit 3; }
mkdir -p "$(dirname "$OUT")"
case "$STRANDNESS" in
  XS|RF|FR) ;;
  *) echo "ERROR: --strandness must be one of XS, RF, FR" >&2; exit 2 ;;
esac
args=(junctions extract -a "$ANCHOR" -m "$MIN_INTRON" -M "$MAX_INTRON" -s "$STRANDNESS" -o "$OUT")
args+=("$BAM")
"$REGTOOLS_BIN" "${args[@]}"
cat > "${OUT%.tsv}.summary.json" <<JSON
{
  "sample_id": "$SAMPLE_ID",
  "method": "regtools junctions extract",
  "strandness": "$STRANDNESS",
  "bam": "$BAM",
  "out": "$OUT"
}
JSON
echo "$OUT"
