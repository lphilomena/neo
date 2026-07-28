#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
usage() { echo "Usage: $0 --bam1 NORMAL.bam --bam2 TUMOR.bam --reference GRCh38.fa --loci hg38.identity.vcf --outdir DIR [--depth 15]" >&2; }
BAM1=""; BAM2=""; REF=""; LOCI=""; OUTDIR=""; DEPTH=15
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam1) BAM1="$2"; shift 2 ;;
    --bam2) BAM2="$2"; shift 2 ;;
    --reference) REF="$2"; shift 2 ;;
    --loci) LOCI="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --depth) DEPTH="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument $1" >&2; usage; exit 2 ;;
  esac
done
for value in "$BAM1" "$BAM2" "$REF" "$LOCI"; do [[ -s "$value" ]] || { echo "ERROR: missing input $value" >&2; exit 2; }; done
[[ -n "$OUTDIR" ]] || { usage; exit 2; }
command -v bam-matcher >/dev/null 2>&1 || { echo "ERROR: bam-matcher wrapper unavailable" >&2; exit 127; }
ENV_NAME="${NEOAG_BAM_MATCHER_ENV:-neoag-bam-matcher}"
ENV_BIN="${NEOAG_CONDA_BASE:-}/envs/$ENV_NAME/bin"
SAMTOOLS_BIN="${BAM_MATCHER_SAMTOOLS:-$(command -v samtools 2>/dev/null || true)}"
FREEBAYES_BIN="${BAM_MATCHER_FREEBAYES:-$(command -v freebayes 2>/dev/null || true)}"
[[ -x "$SAMTOOLS_BIN" ]] || SAMTOOLS_BIN="$ENV_BIN/samtools"
[[ -x "$FREEBAYES_BIN" ]] || FREEBAYES_BIN="$ENV_BIN/freebayes"
[[ -x "$SAMTOOLS_BIN" ]] || { echo "ERROR: samtools unavailable" >&2; exit 127; }
[[ -x "$FREEBAYES_BIN" ]] || { echo "ERROR: freebayes unavailable" >&2; exit 127; }
[[ -s "$REF.fai" ]] || { echo "ERROR: reference FASTA index missing: $REF.fai" >&2; exit 2; }
mkdir -p "$OUTDIR/cache"
CONFIG="$OUTDIR/bam-matcher.conf"
RAW="$OUTDIR/bam_matcher.short.tsv"
cat > "$CONFIG" <<EOF
[VariantCallers]
caller: freebayes
GATK:
freebayes: $FREEBAYES_BIN
samtools: $SAMTOOLS_BIN
varscan:
java: $(command -v java || true)
[ScriptOptions]
DP_threshold: $DEPTH
number_of_SNPs:
fast_freebayes: True
VCF_file: $LOCI
[VariantCallerParameters]
GATK_MEM: 4
GATK_nt: 1
VARSCAN_MEM: 4
[GenomeReference]
REFERENCE: $REF
REF_ALTERNATE:
CHROM_MAP:
[BatchOperations]
CACHE_DIR: $OUTDIR/cache
[Miscellaneous]
EOF
bam-matcher --bam1 "$BAM1" --bam2 "$BAM2" --config "$CONFIG" --output "$RAW" --short-output
PYTHONPATH="$ROOT/src${PYTHONPATH:+:$PYTHONPATH}" python "$ROOT/scripts/parse_bam_matcher.py" \
  --input "$RAW" --output "$OUTDIR/sample_identity.tsv" --fail-on-mismatch
