#!/usr/bin/env bash
# Plan A export: code + data/ref + tools + optional conda_packs (exclude work/results/venv).
#
# Usage:
#   bash scripts/migrate_export_plan_a.sh
#   bash scripts/migrate_export_plan_a.sh /path/to/neoag_v03_plan_a.tar.gz
#
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="${1:-$(dirname "${ROOT}")/neoag_v03_plan_a.tar.gz}"
OUT="$(readlink -f "${OUT}")"
OUT_NAME="$(basename "${OUT}")"

echo "==> Exporting Plan A bundle to ${OUT}"
echo "    Source: ${ROOT}"

tar -czf "${OUT}" \
  --exclude='./work' \
  --exclude='./results' \
  --exclude='./.venv' \
  --exclude='./__pycache__' \
  --exclude="./${OUT_NAME}" \
  --exclude='./migrate_export.log' \
  --exclude='./tools/EasyFuse.tmp' \
  --exclude='./tools/STAR-Fusion.tmp' \
  --exclude='./tools/fusioncatcher.tmp' \
  --exclude='./*.pyc' \
  -C "$(dirname "${ROOT}")" \
  "$(basename "${ROOT}")"

ls -lh "${OUT}"
echo "==> Done. Transfer with:"
echo "  scp ${OUT} na@10.200.50.134:/home/na/project/"
