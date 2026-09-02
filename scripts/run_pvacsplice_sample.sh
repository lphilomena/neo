#!/usr/bin/env bash
# Run pVACsplice and emit an auditable list of all_epitopes outputs.
set -euo pipefail
PYTHON_BIN="${NEOAG_PYTHON:-$(command -v python3 || command -v python || true)}"
[[ -n "$PYTHON_BIN" ]] || { echo "ERROR: python3 or python is required" >&2; exit 3; }
usage() {
  cat <<'USAGE'
Usage:
  bash scripts/run_pvacsplice_sample.sh \
    --junctions regtools.cis_splice_effects.tsv \
    --annotated-vcf tumor.vep.vcf.gz \
    --sample-id SAMPLE001 \
    (--hla HLA-A*02:01,HLA-B*07:02 | --hla-file hla.txt) \
    --algorithms MHCflurry,MHCflurryEL \
    --outdir results/pvacsplice \
    --ref-fasta GRCh38.fa --gtf gencode.gtf [options]

Options:
  --threads N                 Default: 4
  --junction-score N          Optional pVACsplice junction cutoff
  --variant-distance N        Optional pVACsplice cis distance
  --pass-only
  --extra-arg ARG             Repeatable reviewed pVACsplice argument
  --pvacsplice-bin PATH

Outputs:
  OUTDIR/pvacsplice_all_epitopes.list
  OUTDIR/pvacsplice_run_manifest.json
USAGE
}
JUNCTIONS=""; VCF=""; SAMPLE=""; HLA=""; HLA_FILE=""; OUTDIR=""; REF=""; GTF=""
ALGORITHMS="${NEOAG_PVACSPLICE_ALGORITHMS:-MHCflurry}"; THREADS="4"; SCORE=""; DIST=""; PASS=0; EXTRA=()
BIN="${NEOAG_PVACSPLICE_BIN:-$(command -v pvacsplice || true)}"
while [[ $# -gt 0 ]]; do
 case "$1" in
  --junctions) JUNCTIONS="$2"; shift 2;; --annotated-vcf) VCF="$2"; shift 2;;
  --sample-id) SAMPLE="$2"; shift 2;; --hla|--hla-alleles) HLA="$2"; shift 2;;
  --hla-file) HLA_FILE="$2"; shift 2;; --algorithms|--algorithm) ALGORITHMS="$2"; shift 2;;
  --outdir) OUTDIR="$2"; shift 2;; --ref-fasta) REF="$2"; shift 2;; --gtf) GTF="$2"; shift 2;;
  --threads) THREADS="$2"; shift 2;; --junction-score) SCORE="$2"; shift 2;;
  --variant-distance) DIST="$2"; shift 2;; --pass-only) PASS=1; shift;;
  --extra-arg) EXTRA+=("$2"); shift 2;; --pvacsplice-bin) BIN="$2"; shift 2;;
  -h|--help) usage; exit 0;; *) echo "ERROR: unknown argument $1" >&2; usage >&2; exit 2;;
 esac
done
for pair in JUNCTIONS VCF REF GTF; do value="${!pair}"; [[ -n "$value" && -s "$value" ]] || { echo "ERROR: $pair missing/empty: $value" >&2; exit 2; }; done
[[ -n "$SAMPLE" && -n "$OUTDIR" ]] || { echo "ERROR: --sample-id and --outdir required" >&2; exit 2; }
if [[ -n "$HLA_FILE" ]]; then [[ -s "$HLA_FILE" ]] || { echo "ERROR: HLA file missing" >&2; exit 2; }; HLA="$(tr '\n; ' ',,,' < "$HLA_FILE" | sed -E 's/,+/,/g;s/^,//;s/,$//')"; fi
[[ -n "$HLA" ]] || { echo "ERROR: HLA alleles required" >&2; exit 2; }
[[ -n "$BIN" ]] || { echo "ERROR: pvacsplice not found" >&2; exit 3; }
mkdir -p "$OUTDIR"
# pVACsplice expects algorithm names as separate positional args (not a single CSV token).
read -r -a ALG_ARR <<< "$(echo "$ALGORITHMS" | tr ',;' '  ')"
[[ ${#ALG_ARR[@]} -ge 1 ]] || { echo "ERROR: empty --algorithms" >&2; exit 2; }
cmd=("$BIN" run "$JUNCTIONS" "$SAMPLE" "$HLA" "${ALG_ARR[@]}" "$OUTDIR" "$VCF" "$REF" "$GTF" -t "$THREADS")
[[ -n "$SCORE" ]] && cmd+=(--junction-score "$SCORE")
[[ -n "$DIST" ]] && cmd+=(--variant-distance "$DIST")
[[ "$PASS" == 1 ]] && cmd+=(--pass-only)
cmd+=("${EXTRA[@]}")
printf 'pVACsplice command:' >&2; printf ' %q' "${cmd[@]}" >&2; printf '\n' >&2
"${cmd[@]}"
LIST="$OUTDIR/pvacsplice_all_epitopes.list"
find "$OUTDIR" -type f -name '*all_epitopes.tsv' -size +0c -print | sort -u > "$LIST"
[[ -s "$LIST" ]] || { echo "ERROR: pVACsplice completed but no all_epitopes.tsv was found" >&2; exit 4; }
"$PYTHON_BIN" - "$JUNCTIONS" "$VCF" "$SAMPLE" "$HLA" "$ALGORITHMS" "$OUTDIR" "$REF" "$GTF" "$BIN" "$LIST" "${cmd[*]}" <<'PY'
from __future__ import annotations
import hashlib,json,subprocess,sys
from datetime import datetime,timezone
from pathlib import Path
junctions,vcf,sample,hla,algorithms,outdir,ref,gtf,exe,list_path,command=sys.argv[1:]
def sha(p):
 p=Path(p); h=hashlib.sha256()
 with p.open('rb') as f:
  for c in iter(lambda:f.read(1024*1024),b''): h.update(c)
 return h.hexdigest()
try:
 cp=subprocess.run([exe,'--version'],text=True,capture_output=True,check=False)
 version=(cp.stdout or cp.stderr).strip().splitlines()[0]
except Exception: version='UNASSESSED'
outs=[Path(x) for x in Path(list_path).read_text().splitlines() if x.strip()]
payload={'schema_version':'neoag-pvacsplice-run-v2','created_at':datetime.now(timezone.utc).isoformat(),
 'sample_id':sample,'hla_alleles':[x for x in hla.split(',') if x], 'algorithms':[x for x in algorithms.split(',') if x],
 'tool':'pVACsplice','version':version or 'UNASSESSED','executable':exe,'command':command,
 'junctions':{'path':str(Path(junctions).resolve()),'sha256':sha(junctions)},
 'annotated_vcf':{'path':str(Path(vcf).resolve()),'sha256':sha(vcf)},
 'reference_fasta':{'path':str(Path(ref).resolve()),'sha256':sha(ref)},
 'gtf':{'path':str(Path(gtf).resolve()),'sha256':sha(gtf)},
 'all_epitopes':[{'path':str(p.resolve()),'sha256':sha(p)} for p in outs]}
Path(outdir,'pvacsplice_run_manifest.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
PY
printf '%s\n' "$LIST"
