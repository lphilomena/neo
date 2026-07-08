#!/usr/bin/env bash
# Verify a portable neoag reference bundle layout.
set -euo pipefail

BUNDLE="${1:-${NEOAG_REF_BUNDLE:-${NEOAG_TOOLS_ROOT:-}/}}"
if [[ -z "$BUNDLE" || ! -d "$BUNDLE" ]]; then
  echo "ERROR: reference bundle directory not found. Usage: bash scripts/verify_reference_bundle.sh /path/to/neodata4git" >&2
  exit 2
fi
BUNDLE="$(cd "$BUNDLE" && pwd)"
FAILED=0
check_file() { [[ -f "$1" || -L "$1" ]] && echo "[OK] file $1" || { echo "[MISS] file $1"; FAILED=1; }; }
check_dir() { [[ -d "$1" || -L "$1" ]] && echo "[OK] dir  $1" || { echo "[MISS] dir  $1"; FAILED=1; }; }
optional_file() { [[ -f "$1" || -L "$1" ]] && echo "[OK] optional file $1" || echo "[WARN] optional file missing $1"; }
optional_dir() { [[ -d "$1" || -L "$1" ]] && echo "[OK] optional dir  $1" || echo "[WARN] optional dir missing $1"; }

echo "==> Reference bundle: $BUNDLE"
check_file "$BUNDLE/data/ref/hg38/Homo_sapiens_assembly38.fasta"
check_file "$BUNDLE/data/ref/hg38/Homo_sapiens_assembly38.fasta.fai"
check_file "$BUNDLE/data/ref/hg38/Homo_sapiens_assembly38.dict"
check_file "$BUNDLE/data/ref/hg38/gencode.gtf"
optional_file "$BUNDLE/data/ref/hg38/capture.bed"
check_dir "$BUNDLE/data/vep/homo_sapiens/105_GRCh38"
check_file "$BUNDLE/data/vep/homo_sapiens/105_GRCh38/info.txt"
check_file "$BUNDLE/data/facets/reference/common_snp.hg38.vcf.gz"
check_file "$BUNDLE/data/facets/reference/common_snp.hg38.vcf.gz.tbi"
check_dir "$BUNDLE/work/vep_plugins"
check_file "$BUNDLE/work/vep_plugins/Wildtype.pm"
check_file "$BUNDLE/work/vep_plugins/Frameshift.pm"
check_dir "$BUNDLE/work/nextflow_cache"

optional_dir "$BUNDLE/data/ref/ctat/current"
optional_dir "$BUNDLE/data/easyfuse/current"
optional_dir "$BUNDLE/data/ascat/reference/WGS_hg38"
optional_file "$BUNDLE/data/sequenza/reference/GRCh38.primary_assembly.chr.fa"
optional_file "$BUNDLE/data/sequenza/reference/gc.wig.gz"
optional_dir "$BUNDLE/data/hla/spechla_db"
optional_dir "$BUNDLE/data/hla/PRG_MHC_GRCh38_withIMGT"
optional_dir "$BUNDLE/data/hla/optitype_reference"
optional_dir "$BUNDLE/data/lohhla/polysolver"
optional_file "$BUNDLE/data/lohhla/novoalign.lic"
optional_dir "$BUNDLE/data/predictors/netMHCpan"
optional_dir "$BUNDLE/data/predictors/prime"
optional_dir "$BUNDLE/data/predictors/mixMHCpred_install"
optional_dir "$BUNDLE/data/hmf/purple_reference"
optional_file "$BUNDLE/data/normal/proteome/Homo_sapiens.GRCh38.pep.all.fa"

if [[ -f "$BUNDLE/neodata4git.env.sh" ]]; then
  echo "[OK] env  $BUNDLE/neodata4git.env.sh"
else
  echo "[WARN] env file missing $BUNDLE/neodata4git.env.sh"
fi
if [[ -f "$BUNDLE/reference_manifest.tsv" ]]; then
  echo "[OK] manifest $BUNDLE/reference_manifest.tsv"
fi
if [[ -f "$BUNDLE/tool_reference_manifest.tsv" ]]; then
  echo "[OK] manifest $BUNDLE/tool_reference_manifest.tsv"
fi

if [[ "$FAILED" != "0" ]]; then
  echo "==> Required reference bundle checks failed" >&2
  exit 1
fi
echo "==> Reference bundle verification passed"
