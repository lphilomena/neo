#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
HOST="${NEOAG_AGENT_WEB_HOST:-0.0.0.0}"
PORT="${NEOAG_AGENT_WEB_PORT:-8787}"
CONDA_BASE="${NEOAG_CONDA_BASE:-$(conda info --base 2>/dev/null || true)}"
PYTHON="${NEOAG_LLM_PYTHON:-}"
if [[ -z "$PYTHON" && -n "$CONDA_BASE" && -x "$CONDA_BASE/envs/neoag-vllm/bin/python" ]]; then
  PYTHON="$CONDA_BASE/envs/neoag-vllm/bin/python"
elif [[ -z "$PYTHON" ]]; then
  PYTHON="$(command -v python3 || command -v python || true)"
fi
if [[ -z "$PYTHON" || ! -x "$PYTHON" ]]; then
  echo "ERROR: set NEOAG_LLM_PYTHON=/path/to/python with FastAPI/uvicorn installed" >&2
  exit 127
fi
export NEOAG_PROJECT_ROOT="$ROOT"
export PYTHONPATH="$ROOT/src:${PYTHONPATH:-}"
exec "$PYTHON" -m uvicorn neoag.web.agent_app:app --host "$HOST" --port "$PORT"
