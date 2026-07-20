#!/usr/bin/env bash
# Install DeepImmuno-CNN under tools/DeepImmuno (optional immunogenicity source).
#
# Usage:
#   bash scripts/install_deepimmuno.sh
#   bash scripts/install_deepimmuno.sh /path/to/existing/DeepImmuno
#
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TARGET="${1:-${ROOT}/tools/DeepImmuno}"
REPO="${DEEPIMMUNO_GIT_URL:-https://github.com/frankligy/DeepImmuno.git}"

if [[ -d "${TARGET}" && -f "${TARGET}/deepimmuno-cnn.py" ]]; then
  echo "==> DeepImmuno already present at ${TARGET}"
else
  mkdir -p "$(dirname "${TARGET}")"
  echo "==> Cloning DeepImmuno into ${TARGET} ..."
  git clone --depth 1 "${REPO}" "${TARGET}"
fi

required=(
  "deepimmuno-cnn.py"
  "data/after_pca.txt"
  "data/hla2paratopeTable_aligned.txt"
  "models/cnn_model_331_3_7"
)
for rel in "${required[@]}"; do
  if [[ ! -e "${TARGET}/${rel}" ]]; then
    echo "ERROR: missing ${TARGET}/${rel}" >&2
    exit 1
  fi
done

TOOLS_ENV="${ROOT}/conf/tools.env.sh"
line="export DEEPIMMUNO_DIR=\"${TARGET}\""
if [[ -f "${TOOLS_ENV}" ]]; then
  if grep -q 'DEEPIMMUNO_DIR' "${TOOLS_ENV}"; then
    sed -i "s|^export DEEPIMMUNO_DIR=.*|${line}|" "${TOOLS_ENV}"
  else
    printf '\n# DeepImmuno-CNN (optional immunogenicity)\n%s\n' "${line}" >> "${TOOLS_ENV}"
  fi
else
  mkdir -p "${ROOT}/conf"
  printf '%s\n' "${line}" > "${TOOLS_ENV}"
fi

echo "==> Done. Test:"
echo "  source ${TOOLS_ENV}"
echo "  neoag check-tools | grep deepimmuno"
