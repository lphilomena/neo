#!/usr/bin/env bash
# Run pVACsplice for one sample from RegTools junctions and VEP-annotated VCF.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
[[ -f "${ROOT}/conf/tools.env.sh" ]] && source "${ROOT}/conf/tools.env.sh"

usage() {
  cat <<'USAGE'
Usage: bash scripts/run_pvacsplice_sample.sh --junctions J.tsv --annotated-vcf VEP.vcf.gz --sample-id ID --hla HLA-A*02:01,HLA-B*07:02 --outdir OUT [--ref-fasta REF.fa] [--gtf genes.gtf] [--algorithm MHCflurry]

This is the production pVACsplice driver. For a lightweight adapter smoke, use:
  bash scripts/run_splice_tool_smoke.sh
USAGE
}

JUNCTIONS=""; VCF=""; SAMPLE_ID=""; HLA=""; OUTDIR=""; REF_FASTA="${NEOAG_REFERENCE_FASTA:-${NEOAG_REF_BUNDLE:-}/data/ref/hg38/Homo_sapiens_assembly38.fasta}"; GTF="${NEOAG_GENCODE_GTF:-${NEOAG_REF_BUNDLE:-}/data/ref/hg38/gencode.gtf}"; ALGORITHM="${PVACSPLICE_ALGORITHM:-MHCflurry}"; THREADS="${PVACSPLICE_THREADS:-4}"; PASS_ONLY=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --junctions) JUNCTIONS="$2"; shift 2 ;;
    --annotated-vcf) VCF="$2"; shift 2 ;;
    --sample-id) SAMPLE_ID="$2"; shift 2 ;;
    --hla|--hla-alleles) HLA="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --ref-fasta) REF_FASTA="$2"; shift 2 ;;
    --gtf) GTF="$2"; shift 2 ;;
    --algorithm) ALGORITHM="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --pass-only) PASS_ONLY=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done
[[ -n "$JUNCTIONS" && -s "$JUNCTIONS" ]] || { echo "ERROR: --junctions missing or not found: $JUNCTIONS" >&2; exit 2; }
[[ -n "$VCF" && -s "$VCF" ]] || { echo "ERROR: --annotated-vcf missing or not found: $VCF" >&2; exit 2; }
[[ -n "$SAMPLE_ID" ]] || { echo "ERROR: --sample-id required" >&2; exit 2; }
[[ -n "$HLA" ]] || { echo "ERROR: --hla required" >&2; exit 2; }
[[ -n "$OUTDIR" ]] || { echo "ERROR: --outdir required" >&2; exit 2; }
[[ -s "$REF_FASTA" ]] || { echo "ERROR: reference FASTA missing: $REF_FASTA" >&2; exit 3; }
[[ -s "$GTF" ]] || { echo "ERROR: GTF missing: $GTF" >&2; exit 3; }
PVACSPLICE_BIN="${NEOAG_PVACSPLICE_BIN:-$(command -v pvacsplice-neoag || command -v pvacsplice || true)}"
[[ -n "$PVACSPLICE_BIN" && -x "$PVACSPLICE_BIN" ]] || { echo "ERROR: pvacsplice not found; run scripts/setup_tools_env.sh or scripts/install_splice_tools.sh" >&2; exit 3; }
mkdir -p "$OUTDIR"
args=(run "$JUNCTIONS" "$SAMPLE_ID" "$HLA" "$ALGORITHM" "$OUTDIR" "$VCF" "$REF_FASTA" "$GTF" -t "$THREADS")
[[ "$PASS_ONLY" == "1" ]] && args+=(--pass-only)
"$PVACSPLICE_BIN" "${args[@]}"
echo "pVACsplice output: $OUTDIR"
