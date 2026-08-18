#!/usr/bin/env bash
set -euo pipefail
usage() {
  cat <<USAGE
Usage: bash scripts/run_fusioncatcher_sample.sh --fastq1 R1.fq.gz --fastq2 R2.fq.gz --outdir OUTDIR [--sample-id ID] [--fusioncatcher-ref REF]

Environment fallback:
  FUSIONCATCHER_BIN / NEOAG_FUSIONCATCHER_BIN   FusionCatcher executable
  NEOAG_FUSIONCATCHER_REF                       FusionCatcher data/reference directory
USAGE
}
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FASTQ1=""; FASTQ2=""; OUTDIR=""; SAMPLE_ID="sample"; FC_REF="${NEOAG_FUSIONCATCHER_REF:-}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --fastq1) FASTQ1="$2"; shift 2 ;;
    --fastq2) FASTQ2="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --sample-id) SAMPLE_ID="$2"; shift 2 ;;
    --fusioncatcher-ref|--data-dir) FC_REF="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done
[[ -n "$FASTQ1" && -f "$FASTQ1" ]] || { echo "ERROR: --fastq1 missing or not found: ${FASTQ1:-unset}" >&2; exit 2; }
[[ -n "$FASTQ2" && -f "$FASTQ2" ]] || { echo "ERROR: --fastq2 missing or not found: ${FASTQ2:-unset}" >&2; exit 2; }
[[ -n "$OUTDIR" ]] || { echo "ERROR: --outdir required" >&2; exit 2; }
if [[ -f "$ROOT/conf/tools.env.sh" ]]; then
  # shellcheck disable=SC1091
  source "$ROOT/conf/tools.env.sh" || true
fi
FC_BIN="${FUSIONCATCHER_BIN:-${NEOAG_FUSIONCATCHER_BIN:-}}"
if [[ -z "$FC_BIN" ]]; then
  for candidate in \
    "$ROOT/bin/fusioncatcher-neoag" \
    "$ROOT/../open-neo-deploy/bin/fusioncatcher-neoag" \
    "$ROOT/../open-neo-deploy/env_tool/tools/fusioncatcher/bin/fusioncatcher" \
    "$(command -v fusioncatcher-neoag 2>/dev/null || true)" \
    "$(command -v fusioncatcher 2>/dev/null || true)"; do
    if [[ -n "$candidate" && -x "$candidate" ]]; then
      FC_BIN="$candidate"
      break
    fi
  done
fi
[[ -n "$FC_BIN" && -x "$FC_BIN" ]] || { echo "ERROR: FusionCatcher executable not found; set FUSIONCATCHER_BIN or NEOAG_FUSIONCATCHER_BIN" >&2; exit 3; }
if [[ -z "$FC_REF" && -d "$ROOT/../open-neo-deploy/refs/data/easyfuse/easyfuse_ref_v4/fusioncatcher_index" ]]; then
  FC_REF="$ROOT/../open-neo-deploy/refs/data/easyfuse/easyfuse_ref_v4/fusioncatcher_index"
fi
[[ -n "$FC_REF" && -d "$FC_REF" ]] || { echo "ERROR: FusionCatcher reference/data directory not found; set --fusioncatcher-ref or NEOAG_FUSIONCATCHER_REF: ${FC_REF:-unset}" >&2; exit 3; }
mkdir -p "$OUTDIR"
LOG="$OUTDIR/fusioncatcher.log"
if [[ -x "$ROOT/scripts/patch_easyfuse_fusioncatcher_compat.sh" ]]; then
  NEOAG_FUSIONCATCHER_REF="$FC_REF" bash "$ROOT/scripts/patch_easyfuse_fusioncatcher_compat.sh" >> "$LOG" 2>&1 || true
fi
"$FC_BIN" -d "$FC_REF" -i "$FASTQ1,$FASTQ2" -o "$OUTDIR" >> "$LOG" 2>&1
FINAL=""
for candidate in \
  "$OUTDIR/final-list_candidate-fusion-genes.txt" \
  "$OUTDIR/final-list_candidate-fusion-genes.hg19.txt" \
  "$OUTDIR/final-list_candidate-fusion-genes.hg38.txt" \
  "$OUTDIR/final-list_candidate-fusion-genes.tsv"; do
  if [[ -s "$candidate" ]]; then
    FINAL="$candidate"
    break
  fi
done
[[ -n "$FINAL" ]] || FINAL="$(find "$OUTDIR" -maxdepth 2 -type f -name 'final-list_candidate-fusion-genes*' -size +0c 2>/dev/null | head -1 || true)"
[[ -n "$FINAL" && -s "$FINAL" ]] || { echo "ERROR: FusionCatcher final-list output not found in $OUTDIR" >&2; exit 4; }
ln -sf "$(basename "$FINAL")" "$OUTDIR/fusioncatcher.final-list.txt" 2>/dev/null || cp "$FINAL" "$OUTDIR/fusioncatcher.final-list.txt"
cat > "$OUTDIR/fusioncatcher.summary.json" <<JSON
{
  "sample_id": "$SAMPLE_ID",
  "fastq1": "$FASTQ1",
  "fastq2": "$FASTQ2",
  "fusioncatcher_ref": "$FC_REF",
  "fusioncatcher_bin": "$FC_BIN",
  "final_list": "$FINAL",
  "normalized_final_list": "$OUTDIR/fusioncatcher.final-list.txt",
  "log": "$LOG"
}
JSON
echo "$OUTDIR/fusioncatcher.final-list.txt"
