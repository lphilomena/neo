#!/usr/bin/env bash
set -euo pipefail
usage() {
  cat <<USAGE
Usage: bash scripts/run_merge_paired_fastq.sh --outdir OUT --sample-id SAMPLE --fastq1 R1 --fastq2 R2 [--fastq1 R1b --fastq2 R2b ...]

Concatenates one or more paired-end FASTQ(.gz) batches into a single R1/R2 pair.
Input files are copied byte-for-byte in the order supplied; gzip members remain valid when concatenated.
USAGE
}
OUTDIR=""; SAMPLE_ID="sample"; R1S=(); R2S=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2 ;;
    --sample-id) SAMPLE_ID="$2"; shift 2 ;;
    --fastq1) R1S+=("$2"); shift 2 ;;
    --fastq2) R2S+=("$2"); shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done
[[ -n "$OUTDIR" ]] || { echo "ERROR: --outdir required" >&2; exit 2; }
[[ ${#R1S[@]} -gt 0 && ${#R1S[@]} -eq ${#R2S[@]} ]] || { echo "ERROR: matching --fastq1/--fastq2 batches required" >&2; exit 2; }
for f in "${R1S[@]}" "${R2S[@]}"; do
  [[ -f "$f" ]] || { echo "ERROR: FASTQ not found: $f" >&2; exit 2; }
done
mkdir -p "$OUTDIR"
R1_OUT="$OUTDIR/${SAMPLE_ID}_R1.fq.gz"
R2_OUT="$OUTDIR/${SAMPLE_ID}_R2.fq.gz"
: > "$R1_OUT"
: > "$R2_OUT"
append_fastq() {
  local input="$1"
  local output="$2"
  case "$input" in
    *.gz) cat "$input" >> "$output" ;;
    *) gzip -c "$input" >> "$output" ;;
  esac
}
for f in "${R1S[@]}"; do append_fastq "$f" "$R1_OUT"; done
for f in "${R2S[@]}"; do append_fastq "$f" "$R2_OUT"; done
cat > "$OUTDIR/merge_fastq.summary.tsv" <<TSV
sample_id	read	batch_count	output
$SAMPLE_ID	R1	${#R1S[@]}	$R1_OUT
$SAMPLE_ID	R2	${#R2S[@]}	$R2_OUT
TSV
echo "$R1_OUT"
echo "$R2_OUT"
