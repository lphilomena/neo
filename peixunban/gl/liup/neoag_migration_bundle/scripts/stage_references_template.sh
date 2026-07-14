#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "${SCRIPT_DIR}/common.sh"
REF_ROOT="${NEOAG_REF_ROOT:-${TOOLS_ROOT}/refs}"
mkdir -p \
  "${REF_ROOT}/vep_cache/homo_sapiens/115_GRCh38" \
  "${REF_ROOT}/facets" \
  "${REF_ROOT}/hla" \
  "${REF_ROOT}/lohhla/hla_reference" \
  "${REF_ROOT}/ascat" \
  "${REF_ROOT}/arriba" \
  "${REF_ROOT}/ctat" \
  "${REF_ROOT}/easyfuse"
cat > "${BUNDLE_ROOT}/refs/reference_manifest.local_template.yaml" <<MANIFEST
genome_build: GRCh38
reference_fasta: ${REF_ROOT}/GRCh38.fa
gencode_gtf: ${REF_ROOT}/gencode.v44.annotation.gtf
vep_cache: ${REF_ROOT}/vep_cache/homo_sapiens/115_GRCh38
protein_fasta: ${REF_ROOT}/protein.fa
normal_proteome: ${REF_ROOT}/normal_proteome.fa
normal_ligandome: ${REF_ROOT}/normal_ligandome.tsv
hla_reference: ${REF_ROOT}/hla
facets_snp_vcf: ${REF_ROOT}/facets/common_snp.vcf.gz
ascat_loci: ${REF_ROOT}/ascat/G1000_loci_hg38.txt
ascat_alleles: ${REF_ROOT}/ascat/G1000_alleles_hg38.txt
arriba_reference: ${REF_ROOT}/arriba
star_fusion_reference: ${REF_ROOT}/ctat/current
lohhla_reference: ${REF_ROOT}/lohhla/hla_reference
MANIFEST
cat <<MSG
Reference directory skeleton created at:
  ${REF_ROOT}

Template manifest written to:
  ${BUNDLE_ROOT}/refs/reference_manifest.local_template.yaml

Place real production reference files in that tree, then run Doctor with this manifest.
The tiny refs/test_refs files are only for migration self-test, not production.
MSG
