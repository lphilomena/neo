#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 --fastq1 R1.fq.gz --fastq2 R2.fq.gz --ctat-genome-lib DIR --sample-id ID --outdir OUT [--threads N]" >&2
}

FQ1=""; FQ2=""; CTAT=""; SAMPLE_ID="sample"; OUTDIR=""; THREADS=${STAR_FUSION_THREADS:-4}
while [[ $# -gt 0 ]]; do
  case "$1" in
    --fastq1) FQ1="$2"; shift 2 ;;
    --fastq2) FQ2="$2"; shift 2 ;;
    --ctat-genome-lib) CTAT="$2"; shift 2 ;;
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
[[ -d "$CTAT" ]] || { echo "ERROR: CTAT genome library missing: $CTAT" >&2; exit 3; }
STAR_FUSION_BIN="${NEOAG_STAR_FUSION_BIN:-$(command -v STAR-Fusion || true)}"
if [[ -z "$STAR_FUSION_BIN" || ! -x "$STAR_FUSION_BIN" ]]; then
  for candidate in \
    "$(cd "$(dirname "$0")/.." && pwd)/../open-neo-deploy/env_tool/tools/STAR-Fusion/STAR-Fusion" \
    "$(cd "$(dirname "$0")/.." && pwd)/tools/STAR-Fusion/STAR-Fusion"; do
    if [[ -x "$candidate" ]]; then
      STAR_FUSION_BIN="$candidate"
      break
    fi
  done
fi
[[ -n "$STAR_FUSION_BIN" && -x "$STAR_FUSION_BIN" ]] || { echo "ERROR: STAR-Fusion executable not found" >&2; exit 3; }

CONDA_BASE="${NEOAG_CONDA_BASE:-${HOME}/miniforge3}"
if [[ -d "${CONDA_BASE}/envs/neoag-fusion/bin" ]]; then
  export PATH="${CONDA_BASE}/envs/neoag-fusion/bin:${PATH}"
fi
STAR_FUSION_PERL="${STAR_FUSION_PERL:-}"
if [[ -z "$STAR_FUSION_PERL" && -x "${CONDA_BASE}/envs/neoag-vep/bin/perl" ]]; then
  STAR_FUSION_PERL="${CONDA_BASE}/envs/neoag-vep/bin/perl"
fi
if [[ -z "$STAR_FUSION_PERL" ]]; then
  STAR_FUSION_PERL="$(command -v perl || true)"
fi
[[ -n "$STAR_FUSION_PERL" && -x "$STAR_FUSION_PERL" ]] || { echo "ERROR: perl executable not found" >&2; exit 3; }

for perl_dir in \
  "${CONDA_BASE}/envs/neoag-vep/lib/perl5/5.32/site_perl" \
  "${CONDA_BASE}/envs/neoag-vep/lib/perl5/site_perl" \
  "${CONDA_BASE}/envs/neoag-vep/lib/perl5/5.32/vendor_perl" \
  "${CONDA_BASE}/envs/neoag-vep/lib/perl5/vendor_perl" \
  "${CONDA_BASE}/envs/neoag-fusion/lib/perl5/5.32/site_perl" \
  "${CONDA_BASE}/envs/neoag-fusion/lib/perl5/site_perl" \
  "${CONDA_BASE}/envs/neoag-fusion/lib/perl5/5.32/vendor_perl" \
  "${CONDA_BASE}/envs/neoag-fusion/lib/perl5/vendor_perl"; do
  [[ -d "$perl_dir" ]] && export PERL5LIB="${perl_dir}:${PERL5LIB:-}"
done
STAR_FUSION_PERL5LIB="${PERL5LIB:-}"
"$STAR_FUSION_PERL" -MSet::IntervalTree -e 1 >/dev/null 2>&1 || {
  echo "ERROR: Perl module Set::IntervalTree not found; install it or expose neoag-vep Perl libs via PERL5LIB" >&2
  exit 3
}
mkdir -p "$OUTDIR"
PERL5LIB="$STAR_FUSION_PERL5LIB" "$STAR_FUSION_PERL" "$STAR_FUSION_BIN" --left_fq "$FQ1" --right_fq "$FQ2" \
  --genome_lib_dir "$CTAT" --CPU "$THREADS" --output_dir "$OUTDIR" \
  >"$OUTDIR/star-fusion.log" 2>&1
test -s "$OUTDIR/star-fusion.fusion_predictions.tsv"
