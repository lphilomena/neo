#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

echo "== Git summary =="
git status --short | awk '{count[$1]++} END {for (key in count) print key, count[key]}' | sort || true

echo
echo "== Ignored artifact directories present =="
for path in work results tools dist conda_packs .nextflow .nextflow_user; do
  if [[ -e "$path" ]]; then
    printf '%-18s present
' "$path"
  fi
done

echo
echo "== Root-owned files within release tree =="
find . -maxdepth 3 \( -user root -o -group root \) -printf '%u:%g %p\n' | sort | sed -n '1,120p'

echo
echo "== Nextflow home =="
NXF_HOME="${NXF_HOME:-$ROOT/work/.nextflow_home}"
mkdir -p "$NXF_HOME"
if [[ -w "$NXF_HOME" ]]; then
  echo "writable: $NXF_HOME"
else
  echo "not writable: $NXF_HOME"
  exit 1
fi

echo
echo "== Tool environment hint =="
echo "source $ROOT/conf/tools.env.sh"
echo "neoag-v03 check-tools"
