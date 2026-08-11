#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 --fastq1 R1.fq.gz --fastq2 R2.fq.gz --star-index DIR --gtf FILE --sample-id ID --outdir OUT [--threads N]" >&2
}

FQ1=""; FQ2=""; STAR_INDEX=""; GTF=""; SAMPLE_ID="sample"; OUTDIR=""; THREADS=16
while [[ $# -gt 0 ]]; do
  case "$1" in
    --fastq1) FQ1="$2"; shift 2 ;;
    --fastq2) FQ2="$2"; shift 2 ;;
    --star-index) STAR_INDEX="$2"; shift 2 ;;
    --gtf) GTF="$2"; shift 2 ;;
    --sample-id) SAMPLE_ID="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done
validate_fastq_list() {
  local label="$1" value="$2"
  [[ -n "$value" ]] || { echo "ERROR: ${label} is required" >&2; exit 2; }
  local IFS=,
  local path
  for path in $value; do
    [[ -s "$path" ]] || { echo "ERROR: ${label} FASTQ missing or empty: $path" >&2; exit 2; }
  done
}
validate_fastq_list --fastq1 "$FQ1"
validate_fastq_list --fastq2 "$FQ2"
[[ -d "$STAR_INDEX" ]] || { echo "ERROR: STAR index missing: $STAR_INDEX" >&2; exit 3; }
[[ -s "$GTF" ]] || { echo "ERROR: GTF missing: $GTF" >&2; exit 3; }
[[ -n "$OUTDIR" ]] || { echo "ERROR: --outdir is required" >&2; exit 2; }
STAR_BIN="${NEOAG_STAR_BIN:-$(command -v STAR || true)}"
[[ -n "$STAR_BIN" && -x "$STAR_BIN" ]] || { echo "ERROR: STAR executable not found" >&2; exit 3; }
mkdir -p "$OUTDIR"

read_command=(--readFilesCommand cat)
if [[ "$FQ1,$FQ2" == *.gz* ]]; then
  IFS=, read -r -a _fq_all <<< "$FQ1,$FQ2"
  all_gz=1
  any_gz=0
  for _fq in "${_fq_all[@]}"; do
    [[ "$_fq" == *.gz ]] && any_gz=1 || all_gz=0
  done
  if [[ "$all_gz" == "1" ]]; then
    read_command=(--readFilesCommand zcat)
  elif [[ "$any_gz" == "1" ]]; then
    echo "ERROR: FASTQ mates must use the same compression format" >&2
    exit 2
  fi
fi

"$STAR_BIN" \
  --runThreadN "$THREADS" \
  --genomeDir "$STAR_INDEX" \
  --sjdbGTFfile "$GTF" \
  --readFilesIn "$FQ1" "$FQ2" \
  "${read_command[@]}" \
  --outFileNamePrefix "$OUTDIR/" \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NH HI AS nM NM ch \
  --chimSegmentMin 12 \
  --chimJunctionOverhangMin 12 \
  --chimOutType Junctions WithinBAM \
  --chimMultimapScoreRange 3 \
  --chimScoreJunctionNonGTAG -4 \
  --chimMultimapNmax 20 \
  --alignSJDBoverhangMin 10 \
  --alignMatesGapMax 100000 \
  --alignIntronMax 100000 \
  >"$OUTDIR/star.stdout.log" 2>"$OUTDIR/star.stderr.log"

if command -v samtools >/dev/null 2>&1; then
  samtools index -@ "$THREADS" "$OUTDIR/Aligned.sortedByCoord.out.bam"
fi
test -s "$OUTDIR/Aligned.sortedByCoord.out.bam"
test -s "$OUTDIR/Chimeric.out.junction"
test -s "$OUTDIR/SJ.out.tab"
