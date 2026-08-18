#!/usr/bin/env bash
# Patch EasyFuse requantification so STAR_CUSTOM removes STAR-version-specific
# transient-index parameters before reading the index.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
# shellcheck source=/dev/null
source "${ROOT}/conf/tools.env.sh" 2>/dev/null || true

if [[ ! -f "${NEOAG_EASYFUSE_HOME:-}/modules/05_requantification.nf" && -f "${ROOT}/../open-neo-deploy/env_tool/tools/EasyFuse/modules/05_requantification.nf" ]]; then
  NEOAG_EASYFUSE_HOME="${ROOT}/../open-neo-deploy/env_tool/tools/EasyFuse"
fi

module="${NEOAG_EASYFUSE_HOME:?ERROR: set NEOAG_EASYFUSE_HOME}/modules/05_requantification.nf"
[[ -f "${module}" ]] || { echo "ERROR: EasyFuse module not found: ${module}" >&2; exit 1; }

python3 - "${module}" <<'PY'
from pathlib import Path
import re
import sys

p = Path(sys.argv[1])
text = p.read_text()
match = re.search(r'(process STAR_CUSTOM \{.*?script:\s*\n\s*"""\n)(.*?)(\n\s*"""\n\})', text, re.S)
if not match:
    raise SystemExit(f"STAR_CUSTOM block not found in {p}")

prefix, body, suffix = match.groups()
cleanup = """\n    if [ -f ${star_index}/genomeParameters.txt ]; then
        sed -i.bak -e '/^genomeType[[:space:]]/d' -e '/^genomeTransform/d' ${star_index}/genomeParameters.txt
    fi
"""

if "genomeTransform" in body and "genomeParameters.txt" in body:
    print(f"already patched {p}")
else:
    star_pos = body.find("\n    STAR --genomeDir ${star_index}")
    if star_pos < 0:
        raise SystemExit(f"STAR_CUSTOM STAR command not found in {p}")
    body = body[:star_pos] + cleanup + body[star_pos:]
    p.write_text(text[:match.start()] + prefix + body + suffix + text[match.end():])
    print(f"patched {p}")
PY
