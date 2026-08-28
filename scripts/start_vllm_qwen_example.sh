#!/usr/bin/env bash
set -euo pipefail

# Example only. Adjust model path/name, tensor parallelism, GPU memory, and port
# for your local/HPC environment. Do not run this blindly on shared servers.
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CONDA_BASE="${CONDA_BASE:-${NEOAG_CONDA_BASE:-$(conda info --base 2>/dev/null || true)}}"
VLLM_ENV="${VLLM_ENV:-neoag-vllm}"
VLLM_PYTHON="${VLLM_PYTHON:-}"

if [[ -z "${VLLM_PYTHON}" ]]; then
  if [[ -x "${CONDA_BASE}/envs/${VLLM_ENV}/bin/python" ]]; then
    VLLM_PYTHON="${CONDA_BASE}/envs/${VLLM_ENV}/bin/python"
  else
    VLLM_PYTHON="$(command -v python3 || command -v python || true)"
  fi
fi

if [[ -z "${VLLM_PYTHON}" ]]; then
  echo "ERROR: no python interpreter found; set VLLM_PYTHON=/path/to/python" >&2
  exit 127
fi

if ! "${VLLM_PYTHON}" - <<PY >/dev/null 2>&1
import vllm
PY
then
  cat >&2 <<MSG
ERROR: vLLM is not installed in: ${VLLM_PYTHON}

Recommended setup on this server:
  ${CONDA_BASE}/bin/conda create -y -n ${VLLM_ENV} python=3.11
  ${CONDA_BASE}/envs/${VLLM_ENV}/bin/python -m pip install -U pip
  ${CONDA_BASE}/envs/${VLLM_ENV}/bin/python -m pip install vllm

Then run:
  VLLM_PYTHON=${CONDA_BASE}/envs/${VLLM_ENV}/bin/python MODEL_PATH=/path/to/Qwen3-32B bash ${ROOT}/scripts/start_vllm_qwen_example.sh
MSG
  exit 127
fi

MODEL_PATH="${MODEL_PATH:-/models/Qwen3-32B}"
SERVED_MODEL_NAME="${SERVED_MODEL_NAME:-qwen3-32b}"
HOST="${HOST:-0.0.0.0}"
PORT="${PORT:-8000}"
TP="${TP:-1}"
MAX_LEN="${MAX_LEN:-32768}"
API_KEY="${LOCAL_VLLM_API_KEY:-local-dev-key}"

if [[ "${MODEL_PATH}" == /* && ! -e "${MODEL_PATH}" ]]; then
  echo "ERROR: local MODEL_PATH does not exist: ${MODEL_PATH}" >&2
  echo "       Set MODEL_PATH=/path/to/local/model, or use a HuggingFace model id such as Qwen/Qwen3-32B." >&2
  exit 2
elif [[ ! -e "${MODEL_PATH}" ]]; then
  echo "WARN: MODEL_PATH is not a local path; vLLM will treat it as a HuggingFace model id: ${MODEL_PATH}" >&2
fi

exec "${VLLM_PYTHON}" -m vllm.entrypoints.openai.api_server \
  --model "${MODEL_PATH}" \
  --served-model-name "${SERVED_MODEL_NAME}" \
  --host "${HOST}" \
  --port "${PORT}" \
  --tensor-parallel-size "${TP}" \
  --max-model-len "${MAX_LEN}" \
  --api-key "${API_KEY}"
