#!/usr/bin/env bash
# Patch easy-fuse entrypoints to use the conda env's Python (not neoag-tools via env python).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CONDA_CACHE="${ROOT}/work/.nextflow_conda"

fix_env_prefix() {
  local prefix="$1"
  [[ -d "${prefix}/bin" ]] || return 0
  local py="${prefix}/bin/python"
  local rscript="${prefix}/bin/Rscript"
  [[ -x "${py}" ]] || return 0

  for script in "${prefix}/bin/easy-fuse"; do
    [[ -f "${script}" ]] || continue
    if head -1 "${script}" | grep -q '^#!/usr/bin/env python'; then
      sed -i "1s|^#!/usr/bin/env python.*|#!${py}|" "${script}"
      echo "    patched ${script}"
    fi
  done

  if [[ -x "${rscript}" ]]; then
    shopt -s nullglob
    for script in "${prefix}"/lib/python*/site-packages/pyeasyfuse/resources/R/R_model_prediction.R; do
      [[ -f "${script}" ]] || continue
      if head -1 "${script}" | grep -q '^#!/usr/bin/env Rscript'; then
        sed -i "1s|^#!/usr/bin/env Rscript.*|#!${rscript}|" "${script}"
        echo "    patched ${script}"
      fi
    done
    shopt -u nullglob
  fi

  if ! "${prefix}/bin/easy-fuse" --help >/dev/null 2>&1; then
    echo "ERROR: easy-fuse still broken in ${prefix}" >&2
    return 1
  fi
  patch_fusionannotation_gff_name "${prefix}"
  echo "    verified easy-fuse in ${prefix}"
}

patch_fusionannotation_gff_name() {
  local prefix="$1"
  shopt -s nullglob
  for fa in "${prefix}"/lib/python*/site-packages/pyeasyfuse/fusionannotation.py; do
    [[ -f "${fa}" ]] || continue
    if grep -q 'NEOAG_PATCH: gff Name fallback' "${fa}" 2>/dev/null; then
      echo "    already patched ${fa}"
      continue
    fi
    python3 - "${fa}" <<'PY'
import sys
from pathlib import Path
path = Path(sys.argv[1])
text = path.read_text(encoding="utf-8")
old = """            elif parent.id.startswith(\"gene:\"):
                gene_id = parent.id
                gene_name = parent.attributes[\"Name\"][0]
                gene_biotype = parent.attributes[\"biotype\"][0]"""
new = """            elif parent.id.startswith(\"gene:\"):
                gene_id = parent.id
                # NEOAG_PATCH: gff Name fallback (Ensembl 110: ~1.6k genes lack Name)
                gene_name = (
                    parent.attributes.get(\"Name\", [None])[0]
                    or parent.attributes.get(\"gene_name\", [None])[0]
                    or parent.attributes.get(\"gene_id\", [None])[0]
                    or parent.id.replace(\"gene:\", \"\")
                )
                gene_biotype = parent.attributes[\"biotype\"][0]"""
if old not in text:
    sys.exit(f"fusionannotation.py pattern not found in {path}")
path.write_text(text.replace(old, new, 1), encoding="utf-8")
print(f"    patched {path}")
PY
  done
  shopt -u nullglob
}

echo "==> fix_easyfuse_pyeasyfuse_env $(date -Is)"

# shellcheck source=/dev/null
source "${ROOT}/conf/tools.env.sh"

shopt -s nullglob
for prefix in "${CONDA_CACHE}"/env-*; do
  if [[ -f "${prefix}/bin/easy-fuse" ]]; then
    echo "==> ${prefix}"
    fix_env_prefix "${prefix}"
  fi
done

echo "==> fix_easyfuse_pyeasyfuse_env done"
