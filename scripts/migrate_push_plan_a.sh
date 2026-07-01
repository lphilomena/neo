#!/usr/bin/env bash
# Push project to remote via rsync (Plan A). Requires SSH access to target.
#
# Usage:
#   REMOTE=na@10.200.50.134 REMOTE_DIR=/home/na/project bash scripts/migrate_push_plan_a.sh
#
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
REMOTE="${REMOTE:-na@10.200.50.134}"
REMOTE_DIR="${REMOTE_DIR:-/home/na/project}"
DEST="${REMOTE}:${REMOTE_DIR}/neoag_event_pipeline_v03_rc/"

echo "==> Rsync Plan A to ${DEST}"
ssh -o ConnectTimeout=15 "${REMOTE}" "mkdir -p ${REMOTE_DIR}"

rsync -avz --progress \
  --exclude work/ \
  --exclude results/ \
  --exclude .venv/ \
  --exclude __pycache__/ \
  --exclude 'tools/*.tmp/' \
  --exclude 'tools/EasyFuse.tmp/' \
  --exclude 'tools/STAR-Fusion.tmp/' \
  --exclude 'tools/fusioncatcher.tmp/' \
  "${ROOT}/" "${DEST}"

echo "==> Remote setup:"
echo "  ssh ${REMOTE} 'bash ${REMOTE_DIR}/neoag_event_pipeline_v03_rc/scripts/migrate_setup_plan_a.sh'"
