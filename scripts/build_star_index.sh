#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 --reference-fasta FASTA --gtf GTF --star-index DIR [--threads N] [--sjdb-overhang N]" >&2
}

FASTA=""; GTF=""; STAR_INDEX=""; THREADS=16; SJDB_OVERHANG=149
while [[ $# -gt 0 ]]; do
  case "$1" in
    --reference-fasta) FASTA="$2"; shift 2 ;;
    --gtf) GTF="$2"; shift 2 ;;
    --star-index) STAR_INDEX="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --sjdb-overhang) SJDB_OVERHANG="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done

[[ -s "$FASTA" ]] || { echo "ERROR: reference FASTA missing: $FASTA" >&2; exit 2; }
[[ -s "$GTF" ]] || { echo "ERROR: GTF missing: $GTF" >&2; exit 2; }
[[ -n "$STAR_INDEX" ]] || { echo "ERROR: --star-index is required" >&2; exit 2; }
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: --threads must be positive" >&2; exit 2; }
[[ "$SJDB_OVERHANG" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: --sjdb-overhang must be positive" >&2; exit 2; }

required=(Genome SA SAindex genomeParameters.txt)
complete=1
for name in "${required[@]}"; do
  [[ -s "$STAR_INDEX/$name" ]] || complete=0
done
if [[ "$complete" == 1 ]]; then
  echo "STAR index already complete: $STAR_INDEX"
  exit 0
fi

STAR_BIN="${NEOAG_STAR_BIN:-$(command -v STAR || true)}"
[[ -n "$STAR_BIN" && -x "$STAR_BIN" ]] || { echo "ERROR: STAR executable not found" >&2; exit 3; }

parent="$(dirname "$STAR_INDEX")"
mkdir -p "$parent"
building="${STAR_INDEX}.building.$$"
mkdir -p "$building"
cleanup() { rm -rf "$building"; }
trap cleanup EXIT

"$STAR_BIN" \
  --runMode genomeGenerate \
  --runThreadN "$THREADS" \
  --genomeDir "$building" \
  --genomeFastaFiles "$FASTA" \
  --sjdbGTFfile "$GTF" \
  --sjdbOverhang "$SJDB_OVERHANG"

for name in "${required[@]}"; do
  [[ -s "$building/$name" ]] || { echo "ERROR: STAR did not create $name" >&2; exit 4; }
done

if [[ -e "$STAR_INDEX" ]]; then
  backup="${STAR_INDEX}.incomplete.$(date +%Y%m%d_%H%M%S)"
  mv "$STAR_INDEX" "$backup"
  echo "Moved incomplete index to: $backup"
fi
mv "$building" "$STAR_INDEX"
trap - EXIT
echo "STAR index complete: $STAR_INDEX"
