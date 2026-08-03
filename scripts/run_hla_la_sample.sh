#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
[[ -f "$REPO_ROOT/conf/tools.env.sh" ]] && source "$REPO_ROOT/conf/tools.env.sh"

usage() {
  cat <<USAGE
Usage: $0 --bam FILE --sample-id ID --outdir DIR [options]

Options:
  --graph DIR       Prepared HLA-LA graph directory (default: HLALA_GRAPH)
  --threads N       Maximum threads (default: 8)
  --long-read       Pass --longReads ont2d to HLA-LA
  --force           Replace the completion marker and rerun
USAGE
}

BAM=""
SAMPLE_ID=""
OUTDIR=""
GRAPH=${HLALA_GRAPH:-${HLA_LA_GRAPH:-}}
THREADS=8
LONG_READ=0
FORCE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam) BAM="$2"; shift 2 ;;
    --sample-id) SAMPLE_ID="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --graph) GRAPH="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --long-read) LONG_READ=1; shift ;;
    --force) FORCE=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n "$BAM" && -n "$SAMPLE_ID" && -n "$OUTDIR" ]] || { usage >&2; exit 2; }
[[ -s "$BAM" ]] || { echo "ERROR: BAM/CRAM missing: $BAM" >&2; exit 3; }
[[ -s "$BAM.bai" || -s "${BAM%.bam}.bai" || -s "$BAM.crai" ]] || {
  echo "ERROR: BAM/CRAM index missing for $BAM" >&2; exit 3;
}
[[ -d "$GRAPH" && -s "$GRAPH/serializedGRAPH" ]] || { echo "ERROR: prepared graph missing: $GRAPH" >&2; exit 3; }
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: threads must be a positive integer" >&2; exit 2; }

mkdir -p "$OUTDIR/work"
if [[ -s "$OUTDIR/.complete" && "$FORCE" != 1 ]]; then
  echo "HLA-LA already complete: $OUTDIR"
  exit 0
fi
rm -f "$OUTDIR/.complete"

args=(--BAM "$BAM" --graph "$(basename "$GRAPH")" --customGraphDir "$(dirname "$GRAPH")" --sampleID "$SAMPLE_ID" --maxThreads "$THREADS" --workingDir "$OUTDIR/work")
[[ "$LONG_READ" == 1 ]] && args+=(--longReads ont2d)

"$REPO_ROOT/scripts/run_hla_la_container.sh" "${args[@]}" 2>&1 | tee "$OUTDIR/run.log"

RESULT=$(find "$OUTDIR/work" -type f \( -name 'R1_bestguess_G.txt' -o -name 'R1_bestguess.txt' \) -print -quit)
[[ -n "$RESULT" && -s "$RESULT" ]] || { echo "ERROR: HLA-LA bestguess output missing" >&2; exit 5; }
printf 'key\tvalue\n' > "$OUTDIR/run_metadata.tsv"
printf 'sample_id\t%s\n' "$SAMPLE_ID" >> "$OUTDIR/run_metadata.tsv"
printf 'bam\t%s\n' "$BAM" >> "$OUTDIR/run_metadata.tsv"
printf 'graph\t%s\n' "$GRAPH" >> "$OUTDIR/run_metadata.tsv"
printf 'threads\t%s\n' "$THREADS" >> "$OUTDIR/run_metadata.tsv"
printf 'bestguess\t%s\n' "$RESULT" >> "$OUTDIR/run_metadata.tsv"
date -Is > "$OUTDIR/.complete"
echo "HLA-LA completed: $RESULT"
