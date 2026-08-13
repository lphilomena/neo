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

SAMTOOLS="${SAMTOOLS:-${NEOAG_CONDA_BASE:-/home/na/miniforge3}/envs/${NEOAG_TOOLS_ENV:-neoag-tools}/bin/samtools}"
if [[ ! -x "$SAMTOOLS" ]]; then
  SAMTOOLS="$(command -v samtools || true)"
fi
[[ -x "$SAMTOOLS" ]] || { echo "ERROR: samtools not found for SpecHLA direct MHC extraction" >&2; exit 3; }

fq1="$OUTDIR/reads/${SAMPLE_ID}_extract_1.fq.gz"
fq2="$OUTDIR/reads/${SAMPLE_ID}_extract_2.fq.gz"
single="$OUTDIR/reads/${SAMPLE_ID}_extract_single.fq.gz"
extract_bam="$OUTDIR/reads/${SAMPLE_ID}.tmp.extract.bam"

direct_mhc_extract() {
  local contig
  contig="$(set +o pipefail; "$SAMTOOLS" idxstats "$BAM" | awk -F '\t' '($1 == "chr6" || $1 == "6") && $2 >= 33150000 { print $1; exit }')"
  [[ -n "$contig" ]] || { echo "ERROR: cannot find chr6/6 in BAM header for direct MHC extraction" >&2; return 5; }
  echo "WARN: using direct GRCh38 MHC extraction fallback on ${contig}:29600000-33150000" >&2
  "$SAMTOOLS" view -b "$BAM" "${contig}:29600000-33150000" \
    | "$SAMTOOLS" view -b -F 0x4 - \
    | "$SAMTOOLS" sort -@ "$THREADS" -m 2G -O BAM -o "$extract_bam" -
  "$SAMTOOLS" index "$extract_bam"
}

set +e
SPECHLA_MODE=extract "$ROOT/scripts/run_spechla_container.sh" \
  -s "$SAMPLE_ID" -b "$BAM" -r hg38 -o "$OUTDIR/reads" 2>&1 | tee "$OUTDIR/extract.log"
extract_rc=${PIPESTATUS[0]}
set -e
if [[ ! -s "$fq1" || ! -s "$fq2" ]]; then
  if [[ ! -s "$extract_bam" ]]; then
    echo "WARN: SpecHLA extraction failed (exit=$extract_rc) and extracted BAM is missing; trying direct MHC extraction" >&2
    direct_mhc_extract 2>&1 | tee -a "$OUTDIR/extract.log"
  fi
  [[ -s "$extract_bam" ]] || { echo "ERROR: direct MHC extraction did not produce an extracted BAM" >&2; exit 5; }
  echo "WARN: SpecHLA bamUtil FASTQ conversion failed; using samtools fallback" >&2
  rm -f "$fq1" "$fq2" "$single"
  SPECHLA_CMD="$ROOT/scripts/spechla_bam_to_fastq.sh" "$ROOT/scripts/run_spechla_container.sh" \
    "$extract_bam" "$fq1" "$fq2" "$single" "$THREADS" 2>&1 | tee -a "$OUTDIR/extract.log"
fi
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
