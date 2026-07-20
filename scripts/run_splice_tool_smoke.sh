#!/usr/bin/env bash
# Lightweight splice tool smoke: checks installed CLIs and exercises the
# repository splice-junction adapter using fixture junctions.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
[[ -f conf/tools.env.sh ]] && source conf/tools.env.sh

PY="${NEOAG_TOOLS_PYTHON:-${NEOAG_CONDA_BASE:-}/envs/neoag-tools/bin/python}"
[[ -x "$PY" ]] || PY="$(command -v python3)"
OUTDIR="${OUTDIR:-work/splice_tool_smoke}"
mkdir -p "$OUTDIR"

echo "==> pVACsplice CLI"
PVACSPLICE_BIN="${NEOAG_PVACSPLICE_BIN:-$(command -v pvacsplice-neoag || command -v pvacsplice || true)}"
if [[ -n "$PVACSPLICE_BIN" && -x "$PVACSPLICE_BIN" ]]; then
  "$PVACSPLICE_BIN" --help | head -8
else
  echo "WARN: pvacsplice missing"
fi

echo "==> RegTools CLI"
REGTOOLS_BIN="${NEOAG_REGTOOLS_BIN:-$(command -v regtools-neoag || command -v regtools || true)}"
if [[ -n "$REGTOOLS_BIN" && -x "$REGTOOLS_BIN" ]]; then
  "$REGTOOLS_BIN" junctions extract -h | head -8 || true
else
  echo "WARN: regtools missing"
fi

echo "==> Splice-junction adapter fixture"
PYTHONPATH=src "$PY" -m neoag_v03.cli build-intermediates   --entry-mode splice_junction   --splice-junction-tsv data/fixtures/regtools_splice_junctions.tsv   --sample-id SPLICE_SMOKE   --outdir "$OUTDIR/intermediates"

test -s "$OUTDIR/intermediates/parsed/raw_events.tsv"
echo "PASS: splice smoke wrote $OUTDIR/intermediates/parsed/raw_events.tsv"
