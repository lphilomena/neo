#!/usr/bin/env bash
set -euo pipefail
usage() {
  cat <<USAGE
Usage: bash scripts/run_salmon_fastq_to_tpm.sh --fastq1 R1.fq.gz --fastq2 R2.fq.gz --outdir OUT [--sample-id ID] [--salmon-index IDX] [--tx2gene TSV]

Environment fallback:
  SALMON_BIN       salmon executable, default: salmon on PATH
  SALMON_INDEX     Salmon transcriptome index
  SALMON_TX2GENE   transcript-to-gene TSV/CSV
USAGE
}
FASTQ1=""; FASTQ2=""; OUTDIR=""; SAMPLE_ID="sample"; SALMON_INDEX="${SALMON_INDEX:-}"; TX2GENE="${SALMON_TX2GENE:-}"; THREADS="${SALMON_THREADS:-8}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --fastq1) FASTQ1="$2"; shift 2 ;;
    --fastq2) FASTQ2="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --sample-id) SAMPLE_ID="$2"; shift 2 ;;
    --salmon-index) SALMON_INDEX="$2"; shift 2 ;;
    --tx2gene) TX2GENE="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done
[[ -n "$FASTQ1" && -f "$FASTQ1" ]] || { echo "ERROR: --fastq1 missing or not found: $FASTQ1" >&2; exit 2; }
[[ -n "$FASTQ2" && -f "$FASTQ2" ]] || { echo "ERROR: --fastq2 missing or not found: $FASTQ2" >&2; exit 2; }
[[ -n "$OUTDIR" ]] || { echo "ERROR: --outdir required" >&2; exit 2; }
[[ -n "$SALMON_INDEX" && -d "$SALMON_INDEX" ]] || { echo "ERROR: SALMON_INDEX/--salmon-index missing or not a directory: ${SALMON_INDEX:-unset}" >&2; exit 3; }
[[ -n "$TX2GENE" && -f "$TX2GENE" ]] || { echo "ERROR: SALMON_TX2GENE/--tx2gene missing or not a file: ${TX2GENE:-unset}" >&2; exit 3; }
SALMON_BIN="${SALMON_BIN:-$(command -v salmon || true)}"
[[ -n "$SALMON_BIN" && -x "$SALMON_BIN" ]] || { echo "ERROR: salmon executable not found; set SALMON_BIN or put salmon on PATH" >&2; exit 3; }
mkdir -p "$OUTDIR"
QUANT_DIR="$OUTDIR/salmon_quant"
LOG="$OUTDIR/salmon_quant.log"
"$SALMON_BIN" quant -i "$SALMON_INDEX" -l A -1 "$FASTQ1" -2 "$FASTQ2" -p "$THREADS" --validateMappings -o "$QUANT_DIR" >"$LOG" 2>&1
python "$(dirname "$0")/salmon_quant_to_gene_tpm.py" --quant-sf "$QUANT_DIR/quant.sf" --tx2gene "$TX2GENE" --out "$OUTDIR/gene_tpm.tsv"
cat > "$OUTDIR/rna_fastq_to_tpm.summary.json" <<JSON
{
  "sample_id": "$SAMPLE_ID",
  "method": "salmon",
  "fastq1": "$FASTQ1",
  "fastq2": "$FASTQ2",
  "salmon_index": "$SALMON_INDEX",
  "tx2gene": "$TX2GENE",
  "quant_sf": "$QUANT_DIR/quant.sf",
  "gene_tpm": "$OUTDIR/gene_tpm.tsv",
  "log": "$LOG"
}
JSON
echo "$OUTDIR/gene_tpm.tsv"
