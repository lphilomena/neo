#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 4 || $# -gt 5 ]]; then
  echo "Usage: $0 EXTRACTED.bam READ1.fq.gz READ2.fq.gz SINGLE.fq.gz [THREADS]" >&2
  exit 2
fi

bam=$1
fq1=$2
fq2=$3
single=$4
threads=${5:-4}
[[ -s "$bam" ]] || { echo "ERROR: extracted BAM missing: $bam" >&2; exit 3; }
samtools quickcheck "$bam"
samtools collate -@ "$threads" -u -O "$bam" |
  samtools fastq -@ "$threads" -n -1 "$fq1" -2 "$fq2" -s "$single" -0 /dev/null -
[[ -s "$fq1" && -s "$fq2" ]] || { echo "ERROR: samtools FASTQ fallback produced empty mate files" >&2; exit 4; }
