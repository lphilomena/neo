#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
[[ -f "$ROOT/conf/tools.env.sh" ]] && source "$ROOT/conf/tools.env.sh"

usage() {
  cat <<USAGE
Usage: $0 --bam NORMAL.GRCh38.bam --sample-id ID --outdir DIR [options]

Options:
  --threads N          samtools threads (default: 4)
  --samtools-bin FILE  samtools executable (default: SAMTOOLS_BIN or PATH)
  --optitype-bin FILE  OptiType executable (default: OPTITYPE_BIN/optitype)
  --force              rerun even when .complete exists
USAGE
}

BAM=""; SAMPLE_ID=""; OUTDIR=""; THREADS=4; FORCE=0
OPTITYPE_BIN="${OPTITYPE_BIN:-$(command -v optitype 2>/dev/null || true)}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-$(command -v samtools 2>/dev/null || true)}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam) BAM="$2"; shift 2 ;;
    --sample-id) SAMPLE_ID="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --samtools-bin) SAMTOOLS_BIN="$2"; shift 2 ;;
    --optitype-bin) OPTITYPE_BIN="$2"; shift 2 ;;
    --force) FORCE=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -s "$BAM" && -n "$SAMPLE_ID" && -n "$OUTDIR" ]] || { usage >&2; exit 2; }
[[ -s "$BAM.bai" || -s "${BAM%.bam}.bai" ]] || { echo "ERROR: BAM index missing" >&2; exit 3; }
[[ -x "$OPTITYPE_BIN" ]] || { echo "ERROR: OptiType executable missing" >&2; exit 127; }
export PATH="$(dirname "$OPTITYPE_BIN"):$PATH"
[[ -x "$SAMTOOLS_BIN" ]] || { echo "ERROR: samtools missing" >&2; exit 127; }
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: invalid threads" >&2; exit 2; }

mkdir -p "$OUTDIR/fastq" "$OUTDIR/results"
if [[ -s "$OUTDIR/.complete" && "$FORCE" != 1 ]]; then
  echo "OptiType already complete: $OUTDIR"
  exit 0
fi
rm -f "$OUTDIR/.complete"

contig=chr6
"$SAMTOOLS_BIN" view -H "$BAM" | grep -q $'@SQ\tSN:chr6\t' || contig=6
region="${contig}:28510120-33480577"
region_bam="$OUTDIR/fastq/${SAMPLE_ID}.MHC.bam"
r1="$OUTDIR/fastq/${SAMPLE_ID}.MHC.R1.fastq.gz"
r2="$OUTDIR/fastq/${SAMPLE_ID}.MHC.R2.fastq.gz"
single="$OUTDIR/fastq/${SAMPLE_ID}.MHC.single.fastq.gz"

if [[ ! -s "$r1" || ! -s "$r2" ]]; then
  "$SAMTOOLS_BIN" view -@ "$THREADS" -bh "$BAM" "$region" -o "$region_bam"
  "$SAMTOOLS_BIN" collate -@ "$THREADS" -u -O "$region_bam" | \
    "$SAMTOOLS_BIN" fastq -@ "$THREADS" -1 "$r1" -2 "$r2" -s "$single" -0 /dev/null -n -
else
  echo "Reusing existing MHC FASTQs under $OUTDIR/fastq"
fi

if "$OPTITYPE_BIN" run --help >/dev/null 2>&1; then
  "$OPTITYPE_BIN" run -i "$r1" -i "$r2" --dna --prefix "$SAMPLE_ID" \
    --threads "$THREADS" -o "$OUTDIR/results" 2>&1 | tee "$OUTDIR/run.log"
else
  "$OPTITYPE_BIN" -i "$r1" "$r2" --dna --prefix "$SAMPLE_ID" \
    -o "$OUTDIR/results" 2>&1 | tee "$OUTDIR/run.log"
fi

result="$(find "$OUTDIR/results" -type f -name '*_result.tsv' -size +0c -print -quit)"
[[ -n "$result" ]] || { echo "ERROR: OptiType result missing" >&2; exit 5; }
printf 'key\tvalue\n' > "$OUTDIR/run_metadata.tsv"
printf 'sample_id\t%s\n' "$SAMPLE_ID" >> "$OUTDIR/run_metadata.tsv"
printf 'bam\t%s\n' "$BAM" >> "$OUTDIR/run_metadata.tsv"
printf 'region\t%s\n' "$region" >> "$OUTDIR/run_metadata.tsv"
printf 'result\t%s\n' "$result" >> "$OUTDIR/run_metadata.tsv"
date -Is > "$OUTDIR/.complete"
echo "OptiType completed: $result"
