#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
usage() {
  echo "Usage: $0 --sample-id ID --tumor-id ID --normal-id ID --tumor-bam BAM --normal-bam BAM --outdir DIR [--threads N]" >&2
}

SAMPLE_ID=""; TUMOR_ID=""; NORMAL_ID=""; TUMOR_BAM=""; NORMAL_BAM=""; OUTDIR=""; THREADS=4
while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample-id) SAMPLE_ID="$2"; shift 2 ;;
    --tumor-id) TUMOR_ID="$2"; shift 2 ;;
    --normal-id) NORMAL_ID="$2"; shift 2 ;;
    --tumor-bam) TUMOR_BAM="$2"; shift 2 ;;
    --normal-bam) NORMAL_BAM="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option $1" >&2; usage; exit 2 ;;
  esac
done
[[ -n "$SAMPLE_ID" && -n "$TUMOR_ID" && -n "$NORMAL_ID" && -n "$TUMOR_BAM" && -n "$NORMAL_BAM" && -n "$OUTDIR" ]] || { usage; exit 2; }

HMF_ENV="${HMF_ENV:-${HMFTOOLS_HOME:-}/.conda}"
REF_ROOT="${HMFTOOLS_REFERENCE_ROOT:-${NEOAG_REF_BUNDLE:-}/data/hmf/purple_reference}"
REF_FASTA="${HMFTOOLS_REFERENCE_FASTA:-${SEQUENZA_FASTA:-${NEOAG_REFERENCE_FASTA:-}}}"
AMBER_LOCI="${HMFTOOLS_AMBER_LOCI:-$REF_ROOT/amber/GermlineHetPon.38.vcf.gz}"
GC_PROFILE="${HMFTOOLS_GC_PROFILE:-$REF_ROOT/cobalt/GC_profile.1000bp.38.cnp}"
ENSEMBL_DIR="${HMFTOOLS_ENSEMBL_DATA_DIR:-$REF_ROOT/ensembl_data_cache_38}"
INPUT_DIR="$OUTDIR/input"
mkdir -p "$INPUT_DIR"
stage_bam_with_index() {
  local source_bam="$1" label="$2"
  if [[ -s "$source_bam.bai" ]]; then
    printf '%s\n' "$source_bam"
    return
  fi
  local alternate_index="${source_bam%.bam}.bai"
  [[ -s "$alternate_index" ]] || { echo "ERROR: missing BAM index for $source_bam" >&2; return 3; }
  local staged_bam="$INPUT_DIR/${label}.bam"
  ln -sfn "$source_bam" "$staged_bam"
  ln -sfn "$alternate_index" "$staged_bam.bai"
  printf '%s\n' "$staged_bam"
}
TUMOR_BAM="$(stage_bam_with_index "$TUMOR_BAM" tumor)"
NORMAL_BAM="$(stage_bam_with_index "$NORMAL_BAM" normal)"

for exe in amber cobalt purple; do [[ -x "$HMF_ENV/bin/$exe" ]] || { echo "ERROR: missing $HMF_ENV/bin/$exe" >&2; exit 127; }; done
for file in "$TUMOR_BAM" "$TUMOR_BAM.bai" "$NORMAL_BAM" "$NORMAL_BAM.bai" "$REF_FASTA" "$REF_FASTA.fai" "$AMBER_LOCI" "$GC_PROFILE"; do
  [[ -s "$file" ]] || { echo "ERROR: missing required input $file" >&2; exit 3; }
done
[[ -d "$ENSEMBL_DIR" ]] || { echo "ERROR: Ensembl cache missing: $ENSEMBL_DIR" >&2; exit 3; }

AMBER_DIR="$OUTDIR/amber"; COBALT_DIR="$OUTDIR/cobalt"; PURPLE_DIR="$OUTDIR/purple"
mkdir -p "$AMBER_DIR" "$COBALT_DIR" "$PURPLE_DIR"
exec > >(tee -a "$OUTDIR/run.log") 2>&1

if [[ ! -s "$AMBER_DIR/.complete" ]]; then
  "$HMF_ENV/bin/amber" -reference "$NORMAL_ID" -reference_bam "$NORMAL_BAM" -tumor "$TUMOR_ID" -tumor_bam "$TUMOR_BAM" \
    -loci "$AMBER_LOCI" -ref_genome "$REF_FASTA" -ref_genome_version 38 -threads "$THREADS" -output_dir "$AMBER_DIR"
  date -Is > "$AMBER_DIR/.complete"
fi
if [[ ! -s "$COBALT_DIR/.complete" ]]; then
  "$HMF_ENV/bin/cobalt" -reference "$NORMAL_ID" -reference_bam "$NORMAL_BAM" -tumor "$TUMOR_ID" -tumor_bam "$TUMOR_BAM" \
    -gc_profile "$GC_PROFILE" -ref_genome "$REF_FASTA" -ref_genome_version 38 -threads "$THREADS" -output_dir "$COBALT_DIR"
  date -Is > "$COBALT_DIR/.complete"
fi
if [[ ! -s "$PURPLE_DIR/.complete" ]]; then
  "$HMF_ENV/bin/purple" -reference "$NORMAL_ID" -tumor "$TUMOR_ID" -amber_dir "$AMBER_DIR" -cobalt_dir "$COBALT_DIR" \
    -ref_genome "$REF_FASTA" -ref_genome_version 38 -gc_profile "$GC_PROFILE" -ensembl_data_dir "$ENSEMBL_DIR" \
    -threads "$THREADS" -no_charts -output_dir "$PURPLE_DIR"
  date -Is > "$PURPLE_DIR/.complete"
fi
purity_file="$(find "$PURPLE_DIR" -maxdepth 1 -type f -name '*.purple.purity.tsv' -size +0c -print -quit)"
cnv_file="$(find "$PURPLE_DIR" -maxdepth 1 -type f -name '*.purple.cnv.somatic.tsv' -size +0c -print -quit)"
[[ -n "$purity_file" && -n "$cnv_file" ]] || { echo "ERROR: PURPLE outputs missing" >&2; exit 5; }
ln -sfn "$purity_file" "$OUTDIR/purple_purity.tsv"
ln -sfn "$cnv_file" "$OUTDIR/purple_cnv_segments.tsv"
date -Is > "$OUTDIR/.complete"
echo "PURPLE suite complete: $OUTDIR"
