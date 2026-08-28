#!/usr/bin/env bash
set -euo pipefail

CONFIG="${CONFIG:-configs/llm/litellm_config.example.yaml}"
PORT="${PORT:-4000}"

if ! command -v litellm >/dev/null 2>&1; then
  echo "litellm is not installed. Install optional agent deps first:" >&2
  echo "  python -m pip install -e '.[agent-llm]'" >&2
  exit 2
fi

litellm --config "${CONFIG}" --port "${PORT}"
