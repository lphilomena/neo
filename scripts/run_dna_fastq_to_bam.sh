#!/usr/bin/env bash
set -euo pipefail

usage() { echo "Usage: $0 --fastq1 R1 --fastq2 R2 --reference REF --sample-id ID --outdir DIR [--threads N]" >&2; }
FQ1=""; FQ2=""; REF=""; SAMPLE="sample"; OUTDIR=""; THREADS=16
while [[ $# -gt 0 ]]; do
  case "$1" in
    --fastq1) FQ1="$2"; shift 2 ;;
    --fastq2) FQ2="$2"; shift 2 ;;
    --reference) REF="$2"; shift 2 ;;
    --sample-id) SAMPLE="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument $1" >&2; usage; exit 2 ;;
  esac
done
[[ -s "$FQ1" && -s "$FQ2" && -s "$REF" && -n "$OUTDIR" ]] || { usage; exit 2; }
command -v bwa >/dev/null || { echo "ERROR: bwa is unavailable" >&2; exit 127; }
command -v samtools >/dev/null || { echo "ERROR: samtools is unavailable" >&2; exit 127; }
mkdir -p "$OUTDIR"
BAM="$OUTDIR/$SAMPLE.sorted.bam"
if [[ ! -s "$BAM" ]]; then
  bwa mem -t "$THREADS" -R "@RG\\tID:$SAMPLE\\tSM:$SAMPLE\\tPL:ILLUMINA" "$REF" "$FQ1" "$FQ2" \
    | samtools sort -@ "$THREADS" -o "$BAM.tmp" -
  mv "$BAM.tmp" "$BAM"
fi
[[ -s "$BAM.bai" ]] || samtools index -@ "$THREADS" "$BAM"
samtools quickcheck -v "$BAM"
