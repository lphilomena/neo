#!/usr/bin/env bash
set -euo pipefail

usage() { echo "Usage: $0 --tumor-bam T --normal-bam N --reference REF --sample-id ID --outdir DIR [--germline-resource VCF] [--pon VCF]" >&2; }
TUMOR=""; NORMAL=""; REF=""; SAMPLE="sample"; OUTDIR=""; GERMLINE=""; PON=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --tumor-bam) TUMOR="$2"; shift 2 ;;
    --normal-bam) NORMAL="$2"; shift 2 ;;
    --reference) REF="$2"; shift 2 ;;
    --sample-id) SAMPLE="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --germline-resource) GERMLINE="$2"; shift 2 ;;
    --pon) PON="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument $1" >&2; usage; exit 2 ;;
  esac
done
[[ -s "$TUMOR" && -s "$NORMAL" && -s "$REF" && -n "$OUTDIR" ]] || { usage; exit 2; }
command -v gatk >/dev/null || { echo "ERROR: gatk is unavailable" >&2; exit 127; }
command -v samtools >/dev/null || { echo "ERROR: samtools is unavailable" >&2; exit 127; }
mkdir -p "$OUTDIR"
[[ -s "$TUMOR.bai" ]] || samtools index "$TUMOR"
[[ -s "$NORMAL.bai" ]] || samtools index "$NORMAL"
sample_name() {
  local bam="$1" sample=""
  sample="$(samtools samples "$bam" 2>/dev/null | head -1 || true)"
  if [[ -z "$sample" ]]; then
    sample="$(samtools view -H "$bam" | awk -F '\t' '$1 == "@RG" {for (i=2; i<=NF; i++) if ($i ~ /^SM:/) {sub(/^SM:/, "", $i); print $i}}' | sort -u | head -1)"
  fi
  [[ -n "$sample" ]] || { echo "ERROR: BAM has no readable @RG SM sample name: $bam" >&2; return 2; }
  printf '%s\n' "$sample"
}
TUMOR_SAMPLE="$(sample_name "$TUMOR")"
NORMAL_SAMPLE="$(sample_name "$NORMAL")"
[[ "$TUMOR_SAMPLE" != "$NORMAL_SAMPLE" ]] || { echo "ERROR: tumor and normal BAM use the same SM value: $TUMOR_SAMPLE" >&2; exit 2; }
EXTRA=()
[[ -n "$GERMLINE" ]] && EXTRA+=(--germline-resource "$GERMLINE")
[[ -n "$PON" ]] && EXTRA+=(--panel-of-normals "$PON")
RAW="$OUTDIR/somatic.raw.vcf.gz"
PASS="$OUTDIR/somatic.pass.vcf.gz"
gatk Mutect2 -R "$REF" -I "$TUMOR" -tumor "$TUMOR_SAMPLE" -I "$NORMAL" -normal "$NORMAL_SAMPLE" "${EXTRA[@]}" -O "$RAW"
gatk FilterMutectCalls -R "$REF" -V "$RAW" -O "$PASS"
gatk IndexFeatureFile -I "$PASS"
