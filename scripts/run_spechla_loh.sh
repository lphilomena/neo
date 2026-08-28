#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)

usage() {
  cat <<USAGE
Usage: $0 --sample-id ID --typing-dir DIR --purity FLOAT --ploidy FLOAT --outdir DIR [options]

Run SpecHLA's allele-copy/LOH module from an existing SpecHLA typing result.

Required:
  --sample-id ID          Sample identifier written to SpecHLA output
  --typing-dir DIR        Directory containing hla.result.txt and HLA_*_freq.txt
  --purity FLOAT          Tumor purity in (0, 1]
  --ploidy FLOAT          Tumor ploidy greater than zero
  --outdir DIR            Output directory

Optional:
  --typing-result FILE    Override typing-dir/hla.result.txt
  --het-cutoff FLOAT      Minimum informative heterozygous SNP count (default: 5)
  --lohhla-hla-loh FILE  Build a LOHHLA/SpecHLA cross-check table
  --force                 Replace generated result files in outdir
USAGE
}

SAMPLE_ID=""
TYPING_DIR=""
TYPING_RESULT=""
PURITY=""
PLOIDY=""
OUTDIR=""
HET_CUTOFF=5
LOHHLA_HLA_LOH=""
FORCE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample-id) SAMPLE_ID="$2"; shift 2 ;;
    --typing-dir) TYPING_DIR="$2"; shift 2 ;;
    --typing-result) TYPING_RESULT="$2"; shift 2 ;;
    --purity) PURITY="$2"; shift 2 ;;
    --ploidy) PLOIDY="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --het-cutoff) HET_CUTOFF="$2"; shift 2 ;;
    --lohhla-hla-loh) LOHHLA_HLA_LOH="$2"; shift 2 ;;
    --force) FORCE=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

for value in SAMPLE_ID TYPING_DIR PURITY PLOIDY OUTDIR; do
  [[ -n "${!value}" ]] || { echo "ERROR: --${value,,} is required" >&2; exit 2; }
done
[[ -d "$TYPING_DIR" ]] || { echo "ERROR: typing directory missing: $TYPING_DIR" >&2; exit 3; }
awk -v x="$PURITY" 'BEGIN { exit !(x ~ /^[0-9]+([.][0-9]+)?$/ && x > 0 && x <= 1) }' || {
  echo "ERROR: purity must be numeric and in (0, 1]: $PURITY" >&2; exit 2;
}
awk -v x="$PLOIDY" 'BEGIN { exit !(x ~ /^[0-9]+([.][0-9]+)?$/ && x > 0) }' || {
  echo "ERROR: ploidy must be numeric and greater than zero: $PLOIDY" >&2; exit 2;
}
awk -v x="$HET_CUTOFF" 'BEGIN { exit !(x ~ /^[0-9]+([.][0-9]+)?$/ && x >= 0) }' || {
  echo "ERROR: het-cutoff must be numeric and non-negative: $HET_CUTOFF" >&2; exit 2;
}

if [[ -z "$TYPING_RESULT" ]]; then
  TYPING_RESULT="$TYPING_DIR/hla.result.txt"
fi
[[ -s "$TYPING_RESULT" ]] || { echo "ERROR: SpecHLA typing result missing: $TYPING_RESULT" >&2; exit 3; }
TYPING_RESULT="$(readlink -f "$TYPING_RESULT")"

mapfile -t FREQ_FILES < <(find "$TYPING_DIR" -maxdepth 1 -type f -name 'HLA_*_freq.txt' -print | sort)
(( ${#FREQ_FILES[@]} > 0 )) || { echo "ERROR: no HLA_*_freq.txt files found in $TYPING_DIR" >&2; exit 3; }

mkdir -p "$OUTDIR"
OUTDIR="$(readlink -f "$OUTDIR")"
RAW="$OUTDIR/merge.hla.copy.txt"
NORMALIZED="$OUTDIR/hla_loh.tsv"
DETAIL="$OUTDIR/spechla_loh_evidence.tsv"
FILELIST="$OUTDIR/freq.list"
if [[ "$FORCE" != "1" && ( -e "$RAW" || -e "$NORMALIZED" ) ]]; then
  echo "ERROR: output exists; use --force to replace generated result files: $OUTDIR" >&2
  exit 4
fi

for freq_file in "${FREQ_FILES[@]}"; do
  readlink -f "$freq_file"
done > "$FILELIST"
rm -f "$RAW" "$NORMALIZED" "$DETAIL" "$OUTDIR/hla_loh.crosscheck.tsv" "$OUTDIR/hla_loh.consensus.tsv" "$OUTDIR/.complete"

SPECHLA_MODE=loh "$REPO_ROOT/scripts/run_spechla_container.sh" \
  -S "$SAMPLE_ID" -C "$HET_CUTOFF" -purity "$PURITY" -ploidy "$PLOIDY" \
  -F "$FILELIST" -T "$TYPING_RESULT" -O "$OUTDIR"

[[ -s "$RAW" ]] || { echo "ERROR: SpecHLA did not produce $RAW" >&2; exit 5; }
cp "$RAW" "$DETAIL"
"$REPO_ROOT/bin/neoag" convert-spechla -i "$RAW" -o "$NORMALIZED" --min-het "$HET_CUTOFF"
[[ -s "$NORMALIZED" ]] || { echo "ERROR: normalized SpecHLA LOH result missing" >&2; exit 5; }

if [[ -n "$LOHHLA_HLA_LOH" ]]; then
  [[ -s "$LOHHLA_HLA_LOH" ]] || { echo "ERROR: LOHHLA table missing: $LOHHLA_HLA_LOH" >&2; exit 3; }
  "$REPO_ROOT/bin/neoag" crosscheck-hla-loh \
    --lohhla-hla-loh "$LOHHLA_HLA_LOH" \
    --spechla-hla-loh "$NORMALIZED" \
    --out "$OUTDIR/hla_loh.crosscheck.tsv" \
    --consensus-out "$OUTDIR/hla_loh.consensus.tsv"
fi

{
  printf 'key\tvalue\n'
  printf 'sample_id\t%s\n' "$SAMPLE_ID"
  printf 'purity\t%s\n' "$PURITY"
  printf 'ploidy\t%s\n' "$PLOIDY"
  printf 'het_cutoff\t%s\n' "$HET_CUTOFF"
  printf 'typing_result\t%s\n' "$TYPING_RESULT"
  printf 'frequency_file_count\t%s\n' "${#FREQ_FILES[@]}"
} > "$OUTDIR/run_metadata.tsv"
date -Is > "$OUTDIR/.complete"
echo "SpecHLA LOH completed: $NORMALIZED"
