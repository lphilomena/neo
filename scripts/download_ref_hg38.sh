#!/usr/bin/env bash
# Download Broad GRCh38 reference + gnomAD + PoN for GATK Mutect2 (WES/WGS).
#
# Usage:
#   bash scripts/download_ref_hg38.sh
#   bash scripts/download_ref_hg38.sh /path/to/ref_dir
#
# Background (survives terminal close):
#   nohup bash scripts/download_ref_hg38.sh > data/ref/hg38/download.log 2>&1 &
#   tail -f data/ref/hg38/download.log
#
# Output layout (default: data/ref/hg38/):
#   Homo_sapiens_assembly38.fasta[.fai|.dict]
#   GRCh38.fa -> Homo_sapiens_assembly38.fasta   (symlink for docs/CLI)
#   GRCh38.dict -> Homo_sapiens_assembly38.fasta.dict
#   af-only-gnomad.hg38.vcf.gz[.tbi]
#   1000g_pon.hg38.vcf.gz[.tbi]
#   mutect2-exome-panel.vcf.gz -> 1000g_pon.hg38.vcf.gz
#
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
REF_DIR="${1:-${ROOT}/data/ref/hg38}"
ARIA2="${ARIA2:-/usr/bin/aria2c}"
WGET="${WGET:-/usr/bin/wget}"

BROAD_REF="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0"
SOMATIC="https://storage.googleapis.com/gatk-best-practices/somatic-hg38"

mkdir -p "${REF_DIR}"
cd "${REF_DIR}"

fetch() {
  local url="$1"
  local out="$2"
  local min_bytes="${3:-1}"
  if [[ -f "${out}" ]] && [[ "$(stat -c%s "${out}" 2>/dev/null || echo 0)" -ge "${min_bytes}" ]]; then
    echo "==> skip (exists): ${out}"
    return 0
  fi
  rm -f "${out}"
  echo "==> download: ${out}"
  if [[ -x "${ARIA2}" ]]; then
    "${ARIA2}" -c -x 8 -s 8 -o "${out}" "${url}"
  else
    "${WGET}" -c -O "${out}" "${url}"
  fi
}

echo "==> Reference directory: ${REF_DIR}"

fetch "${BROAD_REF}/Homo_sapiens_assembly38.fasta" "Homo_sapiens_assembly38.fasta" 3249912778
fetch "${BROAD_REF}/Homo_sapiens_assembly38.fasta.fai" "Homo_sapiens_assembly38.fasta.fai" 1000
fetch "${BROAD_REF}/Homo_sapiens_assembly38.dict" "Homo_sapiens_assembly38.fasta.dict" 1000

fetch "${SOMATIC}/af-only-gnomad.hg38.vcf.gz" "af-only-gnomad.hg38.vcf.gz" 3184275189
fetch "${SOMATIC}/af-only-gnomad.hg38.vcf.gz.tbi" "af-only-gnomad.hg38.vcf.gz.tbi" 1000

fetch "${SOMATIC}/1000g_pon.hg38.vcf.gz" "1000g_pon.hg38.vcf.gz" 17273497
fetch "${SOMATIC}/1000g_pon.hg38.vcf.gz.tbi" "1000g_pon.hg38.vcf.gz.tbi" 1000

ln -sfn Homo_sapiens_assembly38.fasta GRCh38.fa
ln -sfn Homo_sapiens_assembly38.fasta.fai GRCh38.fa.fai
ln -sfn Homo_sapiens_assembly38.fasta.dict GRCh38.dict
ln -sfn af-only-gnomad.hg38.vcf.gz gnomad.vcf.gz
ln -sfn af-only-gnomad.hg38.vcf.gz.tbi gnomad.vcf.gz.tbi
ln -sfn 1000g_pon.hg38.vcf.gz mutect2-exome-panel.vcf.gz
ln -sfn 1000g_pon.hg38.vcf.gz.tbi mutect2-exome-panel.vcf.gz.tbi

echo "==> Done. Example CLI paths:"
echo "    --reference-fasta ${REF_DIR}/GRCh38.fa"
echo "    --reference-dict  ${REF_DIR}/GRCh38.dict"
echo "    --gnomad-vcf      ${REF_DIR}/af-only-gnomad.hg38.vcf.gz"
echo "    --panel-of-normals ${REF_DIR}/1000g_pon.hg38.vcf.gz"
