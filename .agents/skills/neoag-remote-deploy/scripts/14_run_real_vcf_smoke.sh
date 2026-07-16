#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
TOOLS_ROOT="/root/neo/env_tool"
LICENSED_ROOT="/root/neo/licensed_tools"
CONDA_BASE=""
OUTDIR=""
SAMPLE_ID="M1ML150017383_L01_438"
TOP_N=50
RAW_VCF="/mnt/zjl-bgi-zzb/peixunban/gl/data/chenxiaoliang_data/data/liver_0520_WGS_shortReads/somatic/M1ML150017383_L01_438.align.somatic.pass.vcf.gz"
ANNOTATED_VCF="/mnt/zzb/peixunban/gl/data/chenxiaoliang_data/work/neoag_sliding/M1ML150017383_L01_438_agent_20260708/upstream/tools/M1ML150017383_L01_438.vep.annotated.vcf"
HLA_ALLELES="HLA-A*02:06,HLA-A*30:01,HLA-B*13:02,HLA-B*48:01,HLA-C*06:02,HLA-C*08:01"
HLA_FILE=""
SKIP_MHCFLURRY=0
SKIP_STABPAN=1
SKIP_EXTRACTION=0

usage() {
  cat <<USAGE
Usage: 14_run_real_vcf_smoke.sh [options]

Run a real-data smoke test from the default M1ML150017383 somatic VCF. The
script uses a VEP-annotated VCF for peptide extraction, builds a top-N
peptide-HLA table, then runs peptide-predict with NetMHCpan, PRIME, BigMHC_IM
and optional predictors.

Options:
  --project-root DIR       Project checkout (default: auto-detected)
  --tools-root DIR         Tool/env root (default: /root/neo/env_tool)
  --licensed-root DIR      Licensed tool root (default: /root/neo/licensed_tools)
  --conda-base DIR         Conda base (default: tools-root/miniforge3)
  --outdir DIR             Output directory (default: work/agent_deploy/real_vcf_smoke_<timestamp>)
  --raw-vcf FILE           Raw somatic VCF used as the test anchor
  --annotated-vcf FILE     VEP-annotated VCF for peptide extraction
  --sample-id ID           Sample id (default: M1ML150017383_L01_438)
  --hla-alleles LIST       Comma-separated HLA alleles
  --hla-file FILE          Optional file containing HLA alleles
  --top-n N                Number of unique peptides to predict (default: 50)
  --skip-extraction        Reuse outdir/variant_peptides.annotated.tsv
  --skip-mhcflurry         Skip MHCflurry
  --run-stabpan            Do not skip NetMHCstabpan
  -h, --help               Show help

Defaults run MHCflurry and skip NetMHCstabpan: NetMHCstabpan IEDB smoke tests are slow.
Use --skip-mhcflurry if this host has TensorFlow compatibility issues.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --project-root) PROJECT_ROOT="$2"; shift 2 ;;
    --tools-root) TOOLS_ROOT="$2"; shift 2 ;;
    --licensed-root) LICENSED_ROOT="$2"; shift 2 ;;
    --conda-base) CONDA_BASE="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --raw-vcf) RAW_VCF="$2"; shift 2 ;;
    --annotated-vcf) ANNOTATED_VCF="$2"; shift 2 ;;
    --sample-id) SAMPLE_ID="$2"; shift 2 ;;
    --hla-alleles) HLA_ALLELES="$2"; shift 2 ;;
    --hla-file) HLA_FILE="$2"; shift 2 ;;
    --top-n) TOP_N="$2"; shift 2 ;;
    --skip-extraction) SKIP_EXTRACTION=1; shift ;;
    --skip-mhcflurry) SKIP_MHCFLURRY=1; shift ;;
    --run-stabpan) SKIP_STABPAN=0; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

cd "$PROJECT_ROOT"
if [[ -z "$OUTDIR" ]]; then
  OUTDIR="work/agent_deploy/real_vcf_smoke_$(date +%Y%m%d_%H%M%S)"
fi
mkdir -p "$OUTDIR"
LOG="$OUTDIR/real_vcf_smoke.log"
: > "$LOG"
log() { printf "%s\n" "$*" | tee -a "$LOG"; }

CONDA_BASE="${CONDA_BASE:-$TOOLS_ROOT/miniforge3}"
PY="$CONDA_BASE/envs/neoag-tools/bin/python"
[[ -x "$PY" ]] || PY="$(command -v python3)"

if [[ -f "conf/tools.env.sh" ]]; then
  # shellcheck source=/dev/null
  source conf/tools.env.sh || true
fi

export PYTHONPATH="$PROJECT_ROOT/src:${PYTHONPATH:-}"
export NEOAG_PROJECT_ROOT="$PROJECT_ROOT"
export NEOAG_TOOLS_ROOT="$TOOLS_ROOT"
export NEOAG_CONDA_BASE="$CONDA_BASE"
export NEOAG_NETMHCPAN_BIN="${NEOAG_NETMHCPAN_BIN:-$TOOLS_ROOT/bin/netMHCpan}"
export NEOAG_NETMHCPAN_HOME="${NEOAG_NETMHCPAN_HOME:-$LICENSED_ROOT/netMHCpan}"
export NETMHCPAN_HOME="${NETMHCPAN_HOME:-$LICENSED_ROOT/netMHCpan}"
export NETMHCpan="${NETMHCpan:-$LICENSED_ROOT/netMHCpan}"
export NETMHCSTABPAN_HOME="${NETMHCSTABPAN_HOME:-$LICENSED_ROOT/netMHCstabpan}"
export PRIME_HOME="${PRIME_HOME:-$TOOLS_ROOT/tools/prime}"
export MIXMHCPRED_HOME="${MIXMHCPRED_HOME:-$TOOLS_ROOT/tools/mixMHCpred_install}"
export BIGMHC_DIR="${BIGMHC_DIR:-$TOOLS_ROOT/tools/bigmhc}"
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-}"
export NEOAG_FORCE_CPU="${NEOAG_FORCE_CPU:-1}"
export TF_USE_LEGACY_KERAS="${TF_USE_LEGACY_KERAS:-1}"
export PATH="$TOOLS_ROOT/bin:$LICENSED_ROOT/netMHCstabpan:$PRIME_HOME:$MIXMHCPRED_HOME:$TOOLS_ROOT/tools/bin:$PROJECT_ROOT/bin:$CONDA_BASE/envs/neoag-tools/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:${PATH:-}"

if [[ -n "$HLA_FILE" && -f "$HLA_FILE" ]]; then
  HLA_ALLELES="$($PY - "$HLA_FILE" <<'PYHLA'
import re, sys
text = open(sys.argv[1], errors="ignore").read()
vals = []
for m in re.finditer(r"HLA-[ABC]\*\d{2}:\d{2}(?::\d{2})?", text):
    if m.group(0) not in vals:
        vals.append(m.group(0))
print(",".join(vals))
PYHLA
)"
fi
[[ -n "$HLA_ALLELES" ]] || { echo "ERROR: no HLA alleles resolved" >&2; exit 44; }

log "raw_vcf=$RAW_VCF"
log "annotated_vcf=$ANNOTATED_VCF"
log "hla_alleles=$HLA_ALLELES"
log "outdir=$OUTDIR"

[[ -f "$RAW_VCF" ]] || { echo "REAL_VCF_MISSING: $RAW_VCF" >&2; exit 41; }

VARIANT_TSV="$OUTDIR/variant_peptides.annotated.tsv"
if [[ "$SKIP_EXTRACTION" != "1" ]]; then
  [[ -f "$ANNOTATED_VCF" ]] || { echo "ANNOTATED_VCF_MISSING: $ANNOTATED_VCF" >&2; exit 42; }
  log "==> extract variant peptides"
  "$PY" -m neoag_v03.cli extract-variant-peptides \
    --input-vcf "$ANNOTATED_VCF" \
    --output "$VARIANT_TSV" \
    --sample-id "$SAMPLE_ID" \
    --tumor-sample-name "$SAMPLE_ID" \
    --hla-alleles "$HLA_ALLELES" 2>&1 | tee "$OUTDIR/extract_variant_peptides.log"
else
  [[ -s "$VARIANT_TSV" ]] || { echo "VARIANT_TSV_MISSING_FOR_SKIP_EXTRACTION: $VARIANT_TSV" >&2; exit 43; }
fi

PAIR_TSV="$OUTDIR/peptide_hla_top${TOP_N}.tsv"
log "==> build top-${TOP_N} peptide-HLA input"
"$PY" - "$VARIANT_TSV" "$PAIR_TSV" "$TOP_N" "$SAMPLE_ID" "$HLA_ALLELES" <<'PYPAIRS'
import csv, sys
inp, outp, topn_s, sample, hla_csv = sys.argv[1:6]
topn = int(topn_s)
hlas = [x.strip() for x in hla_csv.split(',') if x.strip()]
seen = []
seen_set = set()
with open(inp, newline='') as fh:
    reader = csv.DictReader(fh, delimiter='\t')
    for row in reader:
        pep = (row.get('mutant_peptide') or row.get('peptide') or '').strip()
        if not pep or pep in seen_set:
            continue
        seen.append(pep)
        seen_set.add(pep)
        if len(seen) >= topn:
            break
with open(outp, 'w', newline='') as out:
    w = csv.writer(out, delimiter='\t')
    w.writerow(['sample_id', 'peptide', 'hla_allele'])
    for pep in seen:
        for hla in hlas:
            w.writerow([sample, pep, hla])
print(f"peptides={len(seen)} hla={len(hlas)} pairs={len(seen) * len(hlas)}")
PYPAIRS

PREDICT_OUT="$OUTDIR/peptide_predict_top${TOP_N}"
args=(peptide-predict -i "$PAIR_TSV" -o "$PREDICT_OUT" --sample-id "$SAMPLE_ID")
[[ "$SKIP_MHCFLURRY" == "1" ]] && args+=(--skip-mhcflurry)
[[ "$SKIP_STABPAN" == "1" ]] && args+=(--skip-stabpan)

log "==> peptide-predict ${args[*]}"
NEOAG_PRIME_JOBS="${NEOAG_PRIME_JOBS:-2}" "$PY" -m neoag_v03.cli "${args[@]}" 2>&1 | tee "$OUTDIR/peptide_predict_top${TOP_N}.log"

log "==> summarize outputs"
find "$OUTDIR" -maxdepth 3 -type f -printf "%p\t%s bytes\n" | sort | tee "$OUTDIR/output_files.tsv"
log "real_vcf_smoke_done=$OUTDIR"
