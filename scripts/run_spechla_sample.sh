#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
usage() { echo "Usage: $0 --bam NORMAL.GRCh38.bam --sample-id ID --outdir DIR [--threads N] [--force]" >&2; }
BAM=""; SAMPLE_ID=""; OUTDIR=""; THREADS=5; FORCE=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam) BAM="$2"; shift 2 ;;
    --sample-id) SAMPLE_ID="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --force) FORCE=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option $1" >&2; usage; exit 2 ;;
  esac
done
[[ -s "$BAM" && -n "$SAMPLE_ID" && -n "$OUTDIR" ]] || { usage; exit 2; }
[[ -s "$BAM.bai" || -s "${BAM%.bam}.bai" ]] || { echo "ERROR: BAM index missing" >&2; exit 3; }
mkdir -p "$OUTDIR/reads" "$OUTDIR/typing"
if [[ -s "$OUTDIR/.complete" && "$FORCE" != 1 ]]; then echo "SpecHLA already complete: $OUTDIR"; exit 0; fi
rm -f "$OUTDIR/.complete"

SPECHLA_MODE=extract "$ROOT/scripts/run_spechla_container.sh" \
  -s "$SAMPLE_ID" -b "$BAM" -r hg38 -o "$OUTDIR/reads" 2>&1 | tee "$OUTDIR/extract.log"
fq1="$OUTDIR/reads/${SAMPLE_ID}_extract_1.fq.gz"
fq2="$OUTDIR/reads/${SAMPLE_ID}_extract_2.fq.gz"
[[ -s "$fq1" && -s "$fq2" ]] || { echo "ERROR: SpecHLA extracted FASTQs missing" >&2; exit 5; }

"$ROOT/scripts/run_spechla_container.sh" -n "$SAMPLE_ID" -1 "$fq1" -2 "$fq2" \
  -o "$OUTDIR/typing" -u 0 -j "$THREADS" 2>&1 | tee "$OUTDIR/run.log"
result="$(find "$OUTDIR/typing" -type f \( -name 'hla.result.txt' -o -name '*.hla.result.txt' \) -size +0c -print -quit)"
[[ -n "$result" ]] || { echo "ERROR: SpecHLA hla.result.txt missing" >&2; exit 5; }
printf 'key\tvalue\n' > "$OUTDIR/run_metadata.tsv"
printf 'sample_id\t%s\n' "$SAMPLE_ID" >> "$OUTDIR/run_metadata.tsv"
printf 'bam\t%s\n' "$BAM" >> "$OUTDIR/run_metadata.tsv"
printf 'result\t%s\n' "$result" >> "$OUTDIR/run_metadata.tsv"
date -Is > "$OUTDIR/.complete"
echo "SpecHLA completed: $result"
