#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(pwd)"
MANIFEST="configs/assets/production_assets.tsv"
SOURCE_HOST=""
OUTDIR="work/remote_deploy/assets"
EXECUTE=0

usage() {
  cat <<'USAGE'
Usage: 15_sync_asset_manifest.sh [options]

Synchronize large deployment assets listed in a TSV manifest. Default mode is
dry-run; add --execute to copy. The manifest is intentionally data-only so large
models/references stay out of Git. Source symlinks are dereferenced so manifests
can point at stable /mnt links while targets receive real files/directories.

Options:
  --project-root DIR       Project checkout (default: current directory)
  --asset-manifest FILE    TSV manifest (default: configs/assets/production_assets.tsv)
  --asset-source-host HOST Default source host for relative/local source paths
  --outdir DIR             Report/log directory (default: work/remote_deploy/assets)
  --execute                Actually copy assets
  -h, --help               Show help

Manifest columns:
  asset_name     Required stable id, e.g. bigmhc_models
  source_path    Source path or user@host:/path. Directories should omit a trailing /*
  target_path    Target file/directory path on this machine
  kind           dir or file (default: dir)
  required       1/0, yes/no, true/false (default: 1)
  sha256         Optional checksum for files, or sha256sum manifest file for dirs
  marker         Optional path inside target directory, or target file marker

Lines beginning with # are ignored. The first non-comment line must be the
header. Unknown columns are ignored.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --project-root) PROJECT_ROOT="$2"; shift 2 ;;
    --asset-manifest) MANIFEST="$2"; shift 2 ;;
    --asset-source-host) SOURCE_HOST="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --execute) EXECUTE=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

cd "$PROJECT_ROOT"
[[ -f "$MANIFEST" ]] || { echo "ASSET_MANIFEST_MISSING: $MANIFEST" >&2; exit 50; }
mkdir -p "$OUTDIR"
LOG="$OUTDIR/asset_sync.log"
REPORT="$OUTDIR/asset_sync_report.tsv"
: > "$LOG"
MODE="DRY_RUN"
[[ "$EXECUTE" == "1" ]] && MODE="EXECUTE"

log() { printf '%s\n' "$*" | tee -a "$LOG"; }
csv_get() {
  local idx="$1"; shift
  local -a fields=("$@")
  if [[ "$idx" =~ ^[0-9]+$ && "$idx" -ge 0 && "$idx" -lt "${#fields[@]}" ]]; then
    printf '%s' "${fields[$idx]}"
  fi
}
is_truthy() {
  case "${1,,}" in
    ""|1|yes|true|required) return 0 ;;
    0|no|false|optional) return 1 ;;
    *) return 0 ;;
  esac
}
source_spec() {
  local src="$1"
  if [[ "$src" == *:* || -z "$SOURCE_HOST" ]]; then
    printf '%s' "$src"
  else
    printf '%s:%s' "$SOURCE_HOST" "$src"
  fi
}
target_has_marker() {
  local kind="$1" dst="$2" marker="$3"
  [[ "$marker" == "-" ]] && marker=""
  if [[ -n "$marker" ]]; then
    [[ -e "$dst/$marker" || -e "$marker" ]]
  elif [[ "$kind" == "file" ]]; then
    [[ -s "$dst" ]]
  else
    [[ -d "$dst" ]]
  fi
}
verify_sha256() {
  local kind="$1" dst="$2" sha="$3"
  [[ -n "$sha" && "$sha" != "-" ]] || return 0
  if [[ "$kind" == "file" ]]; then
    echo "$sha  $dst" | sha256sum -c -
  else
    [[ -f "$dst/$sha" ]] && (cd "$dst" && sha256sum -c "$sha")
  fi
}

header_seen=0
declare -A col=()
{
  echo -e "asset_name\tkind\trequired\tsource\ttarget\tstatus\tdetail"
} > "$REPORT"

while IFS= read -r line || [[ -n "$line" ]]; do
  [[ -z "$line" || "$line" =~ ^[[:space:]]*# ]] && continue
  IFS=$'\t' read -r -a fields <<< "$line"
  if [[ "$header_seen" == "0" ]]; then
    for i in "${!fields[@]}"; do col["${fields[$i]}"]="$i"; done
    for required_col in asset_name source_path target_path; do
      [[ -n "${col[$required_col]+x}" ]] || { echo "ASSET_MANIFEST_BAD_HEADER: missing $required_col" >&2; exit 51; }
    done
    header_seen=1
    continue
  fi

  name="$(csv_get "${col[asset_name]}" "${fields[@]}")"
  src="$(csv_get "${col[source_path]}" "${fields[@]}")"
  dst="$(csv_get "${col[target_path]}" "${fields[@]}")"
  kind="$(csv_get "${col[kind]:--1}" "${fields[@]}")"
  required="$(csv_get "${col[required]:--1}" "${fields[@]}")"
  sha="$(csv_get "${col[sha256]:--1}" "${fields[@]}")"
  marker="$(csv_get "${col[marker]:--1}" "${fields[@]}")"
  kind="${kind:-dir}"
  required="${required:-1}"
  [[ -n "$name" && -n "$src" && -n "$dst" ]] || continue

  if target_has_marker "$kind" "$dst" "$marker"; then
    if verify_sha256 "$kind" "$dst" "$sha" >/dev/null 2>&1; then
      status="present"
      detail="already available"
    else
      status="checksum_failed"
      detail="existing target failed sha256 verification"
      is_truthy "$required" && { echo -e "$name\t$kind\t$required\t$src\t$dst\t$status\t$detail" >> "$REPORT"; echo "ASSET_CHECKSUM_FAILED: $name" >&2; exit 52; }
    fi
    echo -e "$name\t$kind\t$required\t$src\t$dst\t$status\t$detail" >> "$REPORT"
    log "$name: $status ($detail)"
    continue
  fi

  spec="$(source_spec "$src")"
  if [[ "$kind" == "file" ]]; then
    cmd="mkdir -p '$(dirname "$dst")' && rsync -aL '$spec' '$dst'"
  else
    cmd="mkdir -p '$dst' && rsync -aL '$spec/' '$dst/'"
  fi
  log ""
  log "==> [$MODE] sync asset $name"
  log "+ $cmd"
  if [[ "$EXECUTE" == "1" ]]; then
    if bash -lc "$cmd" 2>&1 | tee -a "$LOG"; then
      if target_has_marker "$kind" "$dst" "$marker" && verify_sha256 "$kind" "$dst" "$sha" >/dev/null 2>&1; then
        status="synced"
        detail="copy and verification completed"
      else
        status="verify_failed"
        detail="target marker/checksum missing after copy"
        echo -e "$name\t$kind\t$required\t$src\t$dst\t$status\t$detail" >> "$REPORT"
        is_truthy "$required" && { echo "ASSET_VERIFY_FAILED: $name" >&2; exit 53; }
      fi
    else
      status="sync_failed"
      detail="rsync failed"
      echo -e "$name\t$kind\t$required\t$src\t$dst\t$status\t$detail" >> "$REPORT"
      is_truthy "$required" && { echo "ASSET_SYNC_FAILED: $name" >&2; exit 54; }
    fi
  else
    status="planned"
    detail="dry-run only"
  fi
  echo -e "$name\t$kind\t$required\t$src\t$dst\t$status\t$detail" >> "$REPORT"
done < "$MANIFEST"

[[ "$header_seen" == "1" ]] || { echo "ASSET_MANIFEST_EMPTY: $MANIFEST" >&2; exit 55; }
log ""
log "asset_sync_report=$REPORT"
