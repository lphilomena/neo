#!/usr/bin/env bash
# EasyFuse star_index is built with STAR-avx2; bioconda alignment env ships a monolithic STAR
# that can segfault on large libraries. Install the SIMD dispatch wrapper + STAR-* binaries.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
# shellcheck source=/dev/null
source "${ROOT}/conf/tools.env.sh" 2>/dev/null || true

CONDA_CACHE="${ROOT}/work/.nextflow_conda"
STAR_SRC="${NEOAG_STAR_FUSION_HOME:-${ROOT}/tools/STAR-Fusion}"
if [[ -n "${NEOAG_CONDA_BASE:-}" && -x "${NEOAG_CONDA_BASE}/envs/${NEOAG_FUSION_ENV:-neoag-fusion}/bin/STAR-avx2" ]]; then
  STAR_SRC="${NEOAG_CONDA_BASE}/envs/${NEOAG_FUSION_ENV:-neoag-fusion}/bin"
fi

if [[ -n "${NEOAG_CONDA_BASE:-}" && ! -x "${STAR_SRC}/STAR-avx2" && -x "${NEOAG_CONDA_BASE}/envs/neoag-fusion/bin/STAR-avx2" ]]; then
  STAR_SRC="${NEOAG_CONDA_BASE}/envs/neoag-fusion/bin"
fi

[[ -x "${STAR_SRC}/STAR-avx2" ]] || {
  echo "ERROR: STAR-avx2 source not found (tried ${STAR_SRC})" >&2
  exit 1
}

patch_star_runtime_libs() {
  local prefix="$1"
  local lib
  [[ -d "${prefix}/lib" ]] || mkdir -p "${prefix}/lib"
  for lib in libhts.so libhts.so.3 libhts.so.1.21 libdeflate.so.0 liblzma.so.5 libbz2.so.1.0; do
    if [[ -e "${STAR_SRC}/../lib/${lib}" && ! -e "${prefix}/lib/${lib}" ]]; then
      cp -a "${STAR_SRC}/../lib/${lib}" "${prefix}/lib/${lib}"
    fi
  done
  if [[ -x "${prefix}/bin/STAR-avx2" ]]; then
    if ldd "${prefix}/bin/STAR-avx2" | grep -q 'not found'; then
      echo "ERROR: STAR-avx2 still has unresolved libraries in ${prefix}" >&2
      ldd "${prefix}/bin/STAR-avx2" >&2
      return 1
    fi
  fi
}


patch_starfusion_version_check() {
  local prefix="$1"
  local sf="${prefix}/lib/STAR-Fusion/STAR-Fusion"
  [[ -f "${sf}" ]] || return 0
  if grep -q '^(?:STAR_)?' "${sf}"; then
    echo "    STAR-Fusion version check already patched ${prefix}"
    return 0
  fi
  python3 - "${sf}" <<'PY'
from pathlib import Path
import sys
p = Path(sys.argv[1])
text = p.read_text(errors="ignore")
old = r'if ($star_version_info =~ /^STAR_(\d+)\.(\d+)\./) {'
new = r'if ($star_version_info =~ /^(?:STAR_)?(\d+)\.(\d+)(?:\.|[A-Za-z]|\s|$)/) {'
if old not in text:
    print(f"STAR-Fusion version check pattern not found in {p}; leaving unchanged")
else:
    p.write_text(text.replace(old, new))
    print(f"patched STAR-Fusion version check {p}")
PY
  perl -c "${sf}" >/dev/null
  echo "    verified STAR-Fusion version check ${prefix}"
}

patch_prefix() {
  local prefix="$1"
  [[ -d "${prefix}/bin" ]] || return 0
  [[ -x "${prefix}/bin/STAR" ]] || return 0

  # FusionCatcher 1.33 hard-requires STAR 2.7.2b; do not replace with neoag-fusion 2.7.11b.
  if [[ -x "${prefix}/bin/fusioncatcher" || -f "${prefix}/etc/configuration.cfg" ]]; then
    echo "    skip fusioncatcher env ${prefix}"
    return 0
  fi

  if [[ -f "${prefix}/bin/STAR.orig-bioconda" ]]; then
    echo "    already patched ${prefix}"
    patch_star_runtime_libs "${prefix}"
    patch_starfusion_version_check "${prefix}"
    return 0
  fi

  echo "    patching ${prefix}"
  cp -a "${prefix}/bin/STAR" "${prefix}/bin/STAR.orig-bioconda"
  [[ -x "${prefix}/bin/STARlong" ]] && cp -a "${prefix}/bin/STARlong" "${prefix}/bin/STARlong.orig-bioconda"

  for bin in STAR STAR-avx2 STAR-avx STAR-sse4.1 STAR-ssse3 STAR-sse3 STAR-sse2 STAR-sse STAR-plain \
             STARlong STARlong-avx2 STARlong-avx STARlong-sse4.1 STARlong-ssse3 STARlong-sse3 STARlong-plain; do
    [[ -f "${STAR_SRC}/${bin}" ]] || continue
    cp -a "${STAR_SRC}/${bin}" "${prefix}/bin/${bin}"
    chmod +x "${prefix}/bin/${bin}"
  done

  patch_star_runtime_libs "${prefix}"
  patch_starfusion_version_check "${prefix}"
  "${prefix}/bin/STAR" --version >/dev/null
  echo "    verified STAR wrapper -> $("${prefix}/bin/STAR" --version 2>&1 | head -1)"
}

echo "==> patch_easyfuse_star_avx2 $(date -Is)"
echo "    source=${STAR_SRC}"

shopt -s nullglob
for prefix in "${CONDA_CACHE}"/env-*; do
  [[ -x "${prefix}/bin/STAR" ]] || continue
  patch_prefix "${prefix}"
done

echo "==> patch_easyfuse_star_avx2 done"
