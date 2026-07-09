#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
HOST="${NEOAG_AGENT_WEB_HOST:-0.0.0.0}"
PORT="${NEOAG_AGENT_WEB_PORT:-8787}"
PYTHON="${NEOAG_LLM_PYTHON:-/home/na/miniforge3/envs/neoag-vllm/bin/python}"
if [[ ! -x "$PYTHON" ]]; then
  echo "ERROR: Python not executable: $PYTHON" >&2
  exit 127
fi
export NEOAG_PROJECT_ROOT="$ROOT"
export PYTHONPATH="$ROOT/src:${PYTHONPATH:-}"
exec "$PYTHON" -m uvicorn neoag_v03.web.agent_app:app --host "$HOST" --port "$PORT"
