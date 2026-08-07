#!/usr/bin/env bash
# Run a conda transaction with a pre-fetched Bioconductor data package.
# The patched package cache is isolated and restored before this script exits.
set -euo pipefail

CONDA_BASE="${NEOAG_CONDA_BASE:-}"
CACHE_ROOT="${NEOAG_TOOLS_ROOT:-}/install_cache"
PACKAGE_KEY=""

usage() {
  cat <<'USAGE'
Usage: with_bioc_data_cache.sh --conda-base DIR --cache-root DIR --package-key KEY -- COMMAND [ARG...]

Example package keys: genomeinfodbdata-1.2.9, genomeinfodbdata-1.2.13.
The command must be a conda create/update operation. The helper downloads the
data tarball once, verifies its MD5, and uses an isolated, temporary package
cache patch so failed Bioconductor mirrors do not break the transaction.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --conda-base) CONDA_BASE="$2"; shift 2 ;;
    --cache-root) CACHE_ROOT="$2"; shift 2 ;;
    --package-key) PACKAGE_KEY="$2"; shift 2 ;;
    --) shift; break ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n "$CONDA_BASE" && -x "$CONDA_BASE/bin/conda" ]] || { echo "ERROR: valid --conda-base is required" >&2; exit 2; }
[[ -n "$CACHE_ROOT" ]] || { echo "ERROR: --cache-root is required" >&2; exit 2; }
[[ -n "$PACKAGE_KEY" ]] || { echo "ERROR: --package-key is required" >&2; exit 2; }
[[ $# -gt 0 ]] || { echo "ERROR: command is required after --" >&2; exit 2; }

PKG_CACHE="$CACHE_ROOT/conda_pkgs_bioc"
DATA_CACHE="$CACHE_ROOT/bioconductor"
mkdir -p "$PKG_CACHE" "$DATA_CACHE"

curl_supports_retry_all_errors() {
  command -v curl >/dev/null 2>&1 || return 1
  { curl --help all 2>/dev/null || curl --help 2>/dev/null; } | grep -q -- '--retry-all-errors'
}

download_file() {
  local url="$1" destination="$2"
  local -a curl_args=(-fL --retry 5 --connect-timeout 30)
  if curl_supports_retry_all_errors; then
    curl_args+=(--retry-all-errors)
  else
    echo "INFO: curl lacks --retry-all-errors; using portable retry options." >&2
  fi
  curl "${curl_args[@]}" -o "$destination" "$url"
}

META_DIR="$(find "$PKG_CACHE" -maxdepth 1 -type d -name 'bioconductor-data-packages-*' | sort | tail -1)"
if [[ -z "$META_DIR" ]]; then
  echo "==> Populate isolated Bioconductor metadata cache"
  mapfile -t CONDA_META < <("$CONDA_BASE/bin/conda" search --json --override-channels \
    -c bioconda bioconductor-data-packages | "$CONDA_BASE/bin/python" -c '
import json, sys
records = json.load(sys.stdin)["bioconductor-data-packages"]
record = sorted(records, key=lambda x: (x["version"], x.get("build_number", 0)))[-1]
print(record["fn"])
print(record.get("url") or record["channel"].rstrip("/") + "/" + record["fn"])
')
  META_PACKAGE="$PKG_CACHE/${CONDA_META[0]}"
  if [[ ! -s "$META_PACKAGE" ]]; then
    download_file "${CONDA_META[1]}" "$META_PACKAGE"
  fi
  META_DIR="$PKG_CACHE/${CONDA_META[0]%.conda}"
  rm -rf "$META_DIR"
  "$CONDA_BASE/bin/python" - "$META_PACKAGE" "$META_DIR" <<'PY'
import sys
from conda_package_handling.api import extract
extract(sys.argv[1], dest_dir=sys.argv[2])
PY
fi
[[ -n "$META_DIR" ]] || { echo "ERROR: bioconductor-data-packages metadata was not extracted" >&2; exit 3; }
JSON="$META_DIR/share/bioconductor-data-packages/dataURLs.json"
HOOK="$META_DIR/bin/installBiocDataPackage.sh"
PATHS="$META_DIR/info/paths.json"
[[ -f "$JSON" && -f "$HOOK" && -f "$PATHS" ]] || { echo "ERROR: incomplete Bioconductor metadata package: $META_DIR" >&2; exit 3; }

mapfile -t META < <(python3 - "$JSON" "$PACKAGE_KEY" <<'PY'
import json, sys
entry = json.load(open(sys.argv[1], encoding="utf-8")).get(sys.argv[2])
if not entry:
    raise SystemExit(f"package key not found: {sys.argv[2]}")
print(entry["fn"])
print(entry["md5"])
print(*entry["urls"], sep="\n")
PY
)
FN="${META[0]}"
MD5="${META[1]}"
TARBALL="$DATA_CACHE/$FN"

if ! echo "$MD5  $TARBALL" | md5sum -c - >/dev/null 2>&1; then
  rm -f "$TARBALL"
  if echo "$MD5  $CACHE_ROOT/$FN" | md5sum -c - >/dev/null 2>&1; then
    echo "==> Import existing verified cache: $CACHE_ROOT/$FN"
    cp "$CACHE_ROOT/$FN" "$TARBALL"
  fi
fi
if ! echo "$MD5  $TARBALL" | md5sum -c - >/dev/null 2>&1; then
  for url in "${META[@]:2}"; do
    echo "==> Download $FN from $url"
    if command -v aria2c >/dev/null 2>&1; then
      aria2c -x 1 -s 1 --file-allocation=none --allow-overwrite=true \
        --dir "$DATA_CACHE" --out "$FN" "$url" || true
    else
      download_file "$url" "$TARBALL" || true
    fi
    echo "$MD5  $TARBALL" | md5sum -c - >/dev/null 2>&1 && break
    rm -f "$TARBALL"
  done
fi
echo "$MD5  $TARBALL" | md5sum -c - >/dev/null 2>&1 || {
  echo "ERROR: unable to download an MD5-valid $FN" >&2
  exit 4
}

BACKUP_DIR="$(mktemp -d "$CACHE_ROOT/bioc-patch.XXXXXX")"
cp -a "$HOOK" "$BACKUP_DIR/hook"
cp -a "$PATHS" "$BACKUP_DIR/paths.json"
restore_cache() {
  [[ -d "$BACKUP_DIR" ]] || return 0
  cp -a "$BACKUP_DIR/hook" "$HOOK"
  cp -a "$BACKUP_DIR/paths.json" "$PATHS"
  rm -rf "$BACKUP_DIR"
}
trap restore_cache EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

cat > "$HOOK" <<'HOOK'
#!/bin/bash
set -euo pipefail
PACKAGE_KEY="$1"
SCRIPT_DIR="$(dirname -- "${BASH_SOURCE[0]}")/../share/bioconductor-data-packages"
JSON="$SCRIPT_DIR/dataURLs.json"
FN="$(yq ".\"$PACKAGE_KEY\".fn" "$JSON" | tr -d '"')"
MD5="$(yq ".\"$PACKAGE_KEY\".md5" "$JSON" | tr -d '"')"
CACHED="${NEOAG_BIOC_DATA_CACHE_DIR:?Bioconductor data cache is not configured}/$FN"
echo "$MD5  $CACHED" | md5sum -c -
STAGING="$PREFIX/share/$PACKAGE_KEY"
mkdir -p "$STAGING"
cp "$CACHED" "$STAGING/$FN"
R CMD INSTALL --library="$PREFIX/lib/R/library" "$STAGING/$FN"
rm -f "$STAGING/$FN"
rmdir "$STAGING"
HOOK
chmod +x "$HOOK"

python3 - "$PATHS" "$HOOK" <<'PY'
import hashlib, json, os, sys
paths_file, hook = sys.argv[1:]
data = json.load(open(paths_file, encoding="utf-8"))
payload = open(hook, "rb").read()
for entry in data["paths"]:
    if entry["_path"] == "bin/installBiocDataPackage.sh":
        entry["sha256"] = hashlib.sha256(payload).hexdigest()
        entry["size_in_bytes"] = len(payload)
        break
else:
    raise SystemExit("installBiocDataPackage.sh missing from info/paths.json")
with open(paths_file, "w", encoding="utf-8") as handle:
    json.dump(data, handle, separators=(",", ":"))
PY

echo "==> Run conda transaction with cached $FN"
CONDA_PKGS_DIRS="$PKG_CACHE" NEOAG_BIOC_DATA_CACHE_DIR="$DATA_CACHE" "$@"
