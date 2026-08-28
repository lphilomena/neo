#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 --fastq1 R1.fq.gz --fastq2 R2.fq.gz --sample-id ID --outdir OUT [--threads N]" >&2
}

FQ1=""; FQ2=""; SAMPLE_ID="sample"; OUTDIR=""; THREADS=4
while [[ $# -gt 0 ]]; do
  case "$1" in
    --fastq1) FQ1="$2"; shift 2 ;;
    --fastq2) FQ2="$2"; shift 2 ;;
    --sample-id) SAMPLE_ID="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done
[[ -s "$FQ1" && -s "$FQ2" ]] || { echo "ERROR: paired FASTQ files are required" >&2; exit 2; }
[[ -n "$OUTDIR" ]] || { echo "ERROR: --outdir is required" >&2; exit 2; }
mkdir -p "$OUTDIR"

tool="basic_fastq_integrity"
check_fastq() {
  local path="$1"
  if [[ "$path" == *.gz ]]; then
    gzip -t "$path"
    local header
    header="$(gzip -cd "$path" | head -n 4 || true)"
    printf '%s\n' "$header" | awk 'NR==1 && substr($0,1,1)!="@"{exit 1} NR==3 && substr($0,1,1)!="+"{exit 1}'
  else
    head -n 4 "$path" | awk 'NR==1 && substr($0,1,1)!="@"{exit 1} NR==3 && substr($0,1,1)!="+"{exit 1}'
  fi
}
check_fastq "$FQ1"
check_fastq "$FQ2"
if command -v fastp >/dev/null 2>&1; then
  tool="fastp"
  fastp -i "$FQ1" -I "$FQ2" --thread "$THREADS" \
    --disable_adapter_trimming --disable_quality_filtering --disable_length_filtering \
    --html "$OUTDIR/fastp.html" --json "$OUTDIR/fastp.json" \
    -o /dev/null -O /dev/null >"$OUTDIR/fastp.log" 2>&1
elif command -v fastqc >/dev/null 2>&1; then
  tool="fastqc"
  fastqc --threads "$THREADS" --outdir "$OUTDIR" "$FQ1" "$FQ2" >"$OUTDIR/fastqc.log" 2>&1
fi

cat >"$OUTDIR/qc.complete.json" <<JSON
{"sample_id":"$SAMPLE_ID","status":"PASS","tool":"$tool","fastq1":"$FQ1","fastq2":"$FQ2"}
JSON
