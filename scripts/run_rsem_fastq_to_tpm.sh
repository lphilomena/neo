#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  run_rsem_fastq_to_tpm.sh --fastq1 R1.fq.gz --fastq2 R2.fq.gz --outdir OUTDIR [options]

Options:
  --sample-id SAMPLE        Sample ID used as output prefix, default: sample
  --rsem-reference PREFIX   RSEM reference prefix; default: $RSEM_REFERENCE
  --threads N              Threads; default: $RSEM_THREADS or 8

Environment:
  RSEM_BIN                 rsem-calculate-expression executable, default: command on PATH
  RSEM_REFERENCE           RSEM reference prefix built from matching transcriptome/GTF
  RSEM_THREADS             thread count
USAGE
}

FASTQ1=""; FASTQ2=""; OUTDIR=""; SAMPLE_ID="sample"; RSEM_REFERENCE="${RSEM_REFERENCE:-}"; THREADS="${RSEM_THREADS:-8}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --fastq1) FASTQ1="$2"; shift 2 ;;
    --fastq2) FASTQ2="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --sample-id) SAMPLE_ID="$2"; shift 2 ;;
    --rsem-reference) RSEM_REFERENCE="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n "$FASTQ1" && -f "$FASTQ1" ]] || { echo "ERROR: --fastq1 missing or not a file: ${FASTQ1:-unset}" >&2; exit 3; }
[[ -n "$FASTQ2" && -f "$FASTQ2" ]] || { echo "ERROR: --fastq2 missing or not a file: ${FASTQ2:-unset}" >&2; exit 3; }
[[ -n "$OUTDIR" ]] || { echo "ERROR: --outdir is required" >&2; exit 3; }
[[ -n "$RSEM_REFERENCE" ]] || { echo "ERROR: RSEM_REFERENCE/--rsem-reference is required" >&2; exit 3; }

RSEM_BIN="${RSEM_BIN:-$(command -v rsem-calculate-expression || true)}"
[[ -n "$RSEM_BIN" && -x "$RSEM_BIN" ]] || { echo "ERROR: rsem-calculate-expression not found; set RSEM_BIN or put it on PATH" >&2; exit 3; }

mkdir -p "$OUTDIR"
LOG="$OUTDIR/rsem_quant.log"
OUT_PREFIX="$OUTDIR/$SAMPLE_ID"

"$RSEM_BIN" --paired-end -p "$THREADS" --estimate-rspd --append-names \
  "$FASTQ1" "$FASTQ2" "$RSEM_REFERENCE" "$OUT_PREFIX" >"$LOG" 2>&1

GENES_RESULTS="$OUT_PREFIX.genes.results"
[[ -f "$GENES_RESULTS" ]] || { echo "ERROR: expected RSEM genes result missing: $GENES_RESULTS" >&2; exit 4; }

awk 'BEGIN{FS=OFS="\t"} NR==1{for(i=1;i<=NF;i++){if($i=="gene_id") gid=i; if($i=="TPM") tpm=i} print "gene_id","tpm"; next} gid && tpm {print $gid,$tpm}' \
  "$GENES_RESULTS" > "$OUTDIR/gene_tpm.tsv"

cat > "$OUTDIR/rna_fastq_to_tpm.summary.json" <<JSON
{
  "method": "rsem",
  "sample_id": "$SAMPLE_ID",
  "fastq1": "$FASTQ1",
  "fastq2": "$FASTQ2",
  "rsem_reference": "$RSEM_REFERENCE",
  "genes_results": "$GENES_RESULTS",
  "gene_tpm": "$OUTDIR/gene_tpm.tsv",
  "log": "$LOG"
}
JSON
