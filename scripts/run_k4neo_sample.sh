#!/usr/bin/env bash
# Run k4neo sequence screening with explicit upstream-license acknowledgement.
set -euo pipefail
PYTHON_BIN="${NEOAG_PYTHON:-$(command -v python3 || command -v python || true)}"
[[ -n "$PYTHON_BIN" ]] || { echo "ERROR: python3 or python is required" >&2; exit 3; }

usage() {
  cat <<'USAGE'
Usage:
  bash scripts/run_k4neo_sample.sh \
    --query-table splice_k4neo_input.tsv \
    --database /path/k4neo_db \
    --index /path/index_manifest.tsv \
    --outdir results/k4neo \
    --license-accepted [options]

Required:
  --query-table PATH          NeoAg-generated cts_id/cts_seq query table
  --database PATH             k4neo database directory
  --index PATH                k4neo index manifest
  --outdir DIR
  --license-accepted          Confirms the operator reviewed the upstream license

Options:
  --prefix NAME               Output prefix (default: neoag)
  --k4neo-bin PATH            Override k4neo-annotator
  --extra-arg ARG             Repeatable reviewed raw option
  --dry-run

Outputs:
  OUTDIR/k4neo_healthy_sample_rate.list
  OUTDIR/k4neo_annotated.list
  OUTDIR/k4neo_uniqueness.list
  OUTDIR/k4neo_run_manifest.json
USAGE
}
QUERY=""; DB=""; INDEX=""; OUTDIR=""; PREFIX="neoag"; ACCEPT=0; DRY=0
K4_BIN="${NEOAG_K4NEO_BIN:-$(command -v k4neo-annotator || true)}"; EXTRA=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --query-table) QUERY="$2"; shift 2 ;;
    --database) DB="$2"; shift 2 ;;
    --index) INDEX="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --prefix) PREFIX="$2"; shift 2 ;;
    --k4neo-bin) K4_BIN="$2"; shift 2 ;;
    --extra-arg) EXTRA+=("$2"); shift 2 ;;
    --license-accepted) ACCEPT=1; shift ;;
    --dry-run) DRY=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done
[[ "$ACCEPT" == 1 ]] || { echo "ERROR: k4neo requires --license-accepted after reviewing the upstream license" >&2; exit 2; }
[[ -n "$QUERY" && -s "$QUERY" ]] || { echo "ERROR: query table missing/empty: $QUERY" >&2; exit 2; }
[[ -n "$DB" && -e "$DB" ]] || { echo "ERROR: database not found: $DB" >&2; exit 2; }
[[ -n "$INDEX" && -s "$INDEX" ]] || { echo "ERROR: index manifest missing/empty: $INDEX" >&2; exit 2; }
[[ -n "$OUTDIR" ]] || { echo "ERROR: --outdir required" >&2; exit 2; }
[[ -n "$K4_BIN" ]] || { echo "ERROR: k4neo-annotator not found" >&2; exit 3; }
mkdir -p "$OUTDIR/work"
OUT_PREFIX="$OUTDIR/$PREFIX"
cmd=("$K4_BIN" --database "$DB" --index "$INDEX" --queries "$QUERY" --working-dir "$OUTDIR/work" --output "$OUT_PREFIX")
cmd+=("${EXTRA[@]}")
printf 'k4neo command:' >&2; printf ' %q' "${cmd[@]}" >&2; printf '\n' >&2
[[ "$DRY" == 1 ]] && exit 0
"${cmd[@]}"

HEALTHY="$OUTDIR/k4neo_healthy_sample_rate.list"; ANNOTATED="$OUTDIR/k4neo_annotated.list"; UNIQUE="$OUTDIR/k4neo_uniqueness.list"
find "$OUTDIR" -type f \( -iname '*healthy*sample*rate*.tsv' -o -iname '*healthy*sample*rate*.tsv.gz' \) -size +0c -print | sort -u > "$HEALTHY"
find "$OUTDIR" -type f \( -iname '*annotated*.tsv' -o -iname '*annotated*.tsv.gz' \) -size +0c -print | sort -u > "$ANNOTATED"
find "$OUTDIR" -type f \( -iname '*unique*.tsv' -o -iname '*uniqueness*.tsv' -o -iname '*unique*.tsv.gz' -o -iname '*uniqueness*.tsv.gz' \) -size +0c -print | sort -u > "$UNIQUE"
if [[ ! -s "$HEALTHY" && ! -s "$ANNOTATED" && ! -s "$UNIQUE" ]]; then
  echo "ERROR: k4neo completed but no recognized result TSV was found" >&2; exit 4
fi
"$PYTHON_BIN" - "$QUERY" "$DB" "$INDEX" "$OUTDIR" "$K4_BIN" "$HEALTHY" "$ANNOTATED" "$UNIQUE" "${cmd[*]}" <<'PY'
from __future__ import annotations
import hashlib,json,subprocess,sys
from datetime import datetime,timezone
from pathlib import Path
query,db,index,outdir,exe,healthy,annotated,unique,command=sys.argv[1:]
def sha(p):
 p=Path(p); h=hashlib.sha256()
 with p.open('rb') as f:
  for c in iter(lambda:f.read(1024*1024),b''): h.update(c)
 return h.hexdigest()
def files(list_path):
 return [{'path':str(Path(x).resolve()),'sha256':sha(x)} for x in Path(list_path).read_text().splitlines() if x.strip()]
try:
 cp=subprocess.run([exe,'--version'],text=True,capture_output=True,check=False)
 version=(cp.stdout or cp.stderr).strip().splitlines()[0]
except Exception: version='UNASSESSED'
payload={'schema_version':'neoag-k4neo-run-v1','created_at':datetime.now(timezone.utc).isoformat(),
 'license_acknowledged':True,'tool':'k4neo','version':version or 'UNASSESSED','executable':exe,'command':command,
 'query_table':{'path':str(Path(query).resolve()),'sha256':sha(query)},
 'database':str(Path(db).resolve()),'index':{'path':str(Path(index).resolve()),'sha256':sha(index)},
 'healthy_sample_rate':files(healthy),'annotated':files(annotated),'uniqueness':files(unique)}
Path(outdir,'k4neo_run_manifest.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
PY
printf '%s\n%s\n%s\n' "$HEALTHY" "$ANNOTATED" "$UNIQUE"
