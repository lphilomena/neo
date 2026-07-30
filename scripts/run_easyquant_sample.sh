#!/usr/bin/env bash
# Run TRON EasyQuant/bp-quant against a NeoAg-generated exact query table.
set -euo pipefail
PYTHON_BIN="${NEOAG_PYTHON:-$(command -v python3 || command -v python || true)}"
[[ -n "$PYTHON_BIN" ]] || { echo "ERROR: python3 or python is required" >&2; exit 3; }

usage() {
  cat <<'USAGE'
Usage:
  bash scripts/run_easyquant_sample.sh \
    --query-table splice_easyquant_input.tsv \
    (--bam tumor.rna.bam | --fq1 tumor_R1.fastq.gz --fq2 tumor_R2.fastq.gz) \
    --outdir results/easyquant [options]

Required:
  --query-table PATH          NeoAg-generated table: name,sequence,position
  --bam PATH                  Coordinate-independent RNA BAM input
    or
  --fq1 PATH --fq2 PATH       Paired RNA FASTQ input
  --outdir DIR

Options:
  --bp-distance N             Breakpoint distance (default: 10)
  --mapper star|bowtie2       Mapping backend (default: star)
  --threads N                 Threads (default: 4)
  --allow-mismatches          Pass EasyQuant --allow_mismatches
  --interval-mode             Pass EasyQuant --interval_mode
  --skip-singleton            Pass EasyQuant --skip_singleton
  --no-stringent              Do not pass --stringent_params
  --keep-all                  Pass EasyQuant --keep_all
  --bp-quant-bin PATH         Override bp_quant executable
  --dry-run                   Print the command and exit

Outputs:
  OUTDIR/easyquant_quantification.list
  OUTDIR/easyquant_run_manifest.json
USAGE
}

QUERY=""; BAM=""; FQ1=""; FQ2=""; OUTDIR=""
DISTANCE="10"; MAPPER="star"; THREADS="4"; STRINGENT=1; DRY=0
ALLOW=0; INTERVAL=0; SKIP_SINGLETON=0; KEEP_ALL=0
BP_BIN="${NEOAG_BP_QUANT_BIN:-$(command -v bp_quant || true)}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --query-table) QUERY="$2"; shift 2 ;;
    --bam) BAM="$2"; shift 2 ;;
    --fq1) FQ1="$2"; shift 2 ;;
    --fq2) FQ2="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --bp-distance) DISTANCE="$2"; shift 2 ;;
    --mapper) MAPPER="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --allow-mismatches) ALLOW=1; shift ;;
    --interval-mode) INTERVAL=1; shift ;;
    --skip-singleton) SKIP_SINGLETON=1; shift ;;
    --no-stringent) STRINGENT=0; shift ;;
    --keep-all) KEEP_ALL=1; shift ;;
    --bp-quant-bin) BP_BIN="$2"; shift 2 ;;
    --dry-run) DRY=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n "$QUERY" && -s "$QUERY" ]] || { echo "ERROR: --query-table missing/empty: $QUERY" >&2; exit 2; }
[[ -n "$OUTDIR" ]] || { echo "ERROR: --outdir required" >&2; exit 2; }
[[ "$MAPPER" == "star" || "$MAPPER" == "bowtie2" ]] || { echo "ERROR: --mapper must be star or bowtie2" >&2; exit 2; }
if [[ -n "$BAM" ]]; then
  [[ -z "$FQ1" && -z "$FQ2" ]] || { echo "ERROR: choose BAM or FASTQ, not both" >&2; exit 2; }
  [[ -s "$BAM" ]] || { echo "ERROR: BAM missing/empty: $BAM" >&2; exit 2; }
else
  [[ -n "$FQ1" && -n "$FQ2" && -s "$FQ1" && -s "$FQ2" ]] || { echo "ERROR: provide --bam or paired --fq1/--fq2" >&2; exit 2; }
fi
[[ -n "$BP_BIN" ]] || { echo "ERROR: bp_quant not found" >&2; exit 3; }
if [[ "$DRY" != 1 ]]; then [[ -x "$BP_BIN" || "$(command -v "$BP_BIN" 2>/dev/null || true)" ]] || { echo "ERROR: bp_quant executable not usable: $BP_BIN" >&2; exit 3; }; fi
mkdir -p "$OUTDIR"

cmd=("$BP_BIN" pipeline -s "$QUERY" -o "$OUTDIR" -d "$DISTANCE" -m "$MAPPER" -t "$THREADS")
if [[ -n "$BAM" ]]; then cmd+=(-b "$BAM"); else cmd+=(-1 "$FQ1" -2 "$FQ2"); fi
[[ "$STRINGENT" == 1 ]] && cmd+=(--stringent_params)
[[ "$ALLOW" == 1 ]] && cmd+=(--allow_mismatches)
[[ "$INTERVAL" == 1 ]] && cmd+=(--interval_mode)
[[ "$SKIP_SINGLETON" == 1 ]] && cmd+=(--skip_singleton)
[[ "$KEEP_ALL" == 1 ]] && cmd+=(--keep_all)

printf 'EasyQuant command:' >&2; printf ' %q' "${cmd[@]}" >&2; printf '\n' >&2
[[ "$DRY" == 1 ]] && exit 0
"${cmd[@]}"

LIST="$OUTDIR/easyquant_quantification.list"
find "$OUTDIR" -type f \( -name 'quantification.tsv' -o -name '*quantification*.tsv' \) -size +0c -print | sort -u > "$LIST"
[[ -s "$LIST" ]] || { echo "ERROR: EasyQuant completed but no non-empty quantification TSV was found" >&2; exit 4; }

"$PYTHON_BIN" - "$QUERY" "$BAM" "$FQ1" "$FQ2" "$OUTDIR" "$BP_BIN" "$LIST" "${cmd[*]}" <<'PY'
from __future__ import annotations
import hashlib, json, subprocess, sys
from datetime import datetime, timezone
from pathlib import Path
query, bam, fq1, fq2, outdir, exe, list_path, command = sys.argv[1:]
def sha(p: Path) -> str:
    h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024), b''): h.update(chunk)
    return h.hexdigest()
def meta(value: str):
    if not value: return None
    p=Path(value)
    return {'path':str(p.resolve()), 'sha256':sha(p)}
try:
    cp=subprocess.run([exe,'--version'],text=True,capture_output=True,check=False)
    version=(cp.stdout or cp.stderr).strip().splitlines()[0]
except Exception:
    version='UNASSESSED'
outputs=[Path(x) for x in Path(list_path).read_text().splitlines() if x.strip()]
payload={
 'schema_version':'neoag-easyquant-run-v1','created_at':datetime.now(timezone.utc).isoformat(),
 'tool':'EasyQuant/bp-quant','version':version or 'UNASSESSED','executable':exe,'command':command,
 'query_table':meta(query),'bam':meta(bam),'fq1':meta(fq1),'fq2':meta(fq2),
 'quantification':[{'path':str(p.resolve()),'sha256':sha(p)} for p in outputs],
}
Path(outdir,'easyquant_run_manifest.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
PY
printf '%s\n' "$LIST"
