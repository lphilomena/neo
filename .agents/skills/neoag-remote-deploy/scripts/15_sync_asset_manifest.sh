#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(pwd)"
MANIFEST="configs/assets/production_assets.tsv"
SOURCE_HOST=""
SSH_KEY="${NEOAG_ASSET_SSH_KEY:-}"
SHARED_ASSET_ROOT="${NEOAG_SHARED_ASSET_ROOT:-}"
TOOLS_ROOT="${NEOAG_TOOLS_ROOT:-/opt/neoag/env_tool}"
REFERENCE_ROOT="${NEOAG_REFERENCE_ROOT:-/opt/neoag/refs}"
LICENSED_ROOT="${NEOAG_LICENSED_ROOT:-/opt/neoag/licensed_tools}"
OUTDIR="work/remote_deploy/assets"
EXECUTE=0

usage() {
  cat <<'USAGE'
Usage: 15_sync_asset_manifest.sh [options]

Synchronize large deployment assets listed in a TSV manifest. Default mode is
dry-run; add --execute to copy. The manifest is intentionally data-only so large
models/references stay out of Git. With --shared-asset-root, local manifest
sources are linked into the install tree to avoid duplicating large assets.

Options:
  --project-root DIR       Project checkout (default: current directory)
  --asset-manifest FILE    TSV manifest (default: configs/assets/production_assets.tsv)
  --asset-source-host HOST Default source host for relative/local source paths
  --asset-ssh-key FILE    SSH private key used by rsync for remote asset paths
  --shared-asset-root DIR  Link assets from a locally mounted shared root
  --tools-root DIR         Resolve /srv/neoag-tools targets below this root
  --reference-root DIR     Resolve reference targets below this root
  --licensed-root DIR      Resolve licensed-tool targets below this root
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
    --asset-ssh-key) SSH_KEY="$2"; shift 2 ;;
    --shared-asset-root) SHARED_ASSET_ROOT="$2"; shift 2 ;;
    --tools-root) TOOLS_ROOT="$2"; shift 2 ;;
    --reference-root) REFERENCE_ROOT="$2"; shift 2 ;;
    --licensed-root) LICENSED_ROOT="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --execute) EXECUTE=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ -n "$SHARED_ASSET_ROOT" && -n "$SOURCE_HOST" ]]; then
  echo "ERROR: --shared-asset-root cannot be combined with --asset-source-host" >&2
  exit 2
fi

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
resolve_shared_source() {
  local src="$1"
  if [[ -n "$SHARED_ASSET_ROOT" && "$src" == /srv/neoag-assets/source/* ]]; then
    printf '%s/%s' "${SHARED_ASSET_ROOT%/}" "${src#/srv/neoag-assets/source/}"
  else
    printf '%s' "$src"
  fi
}
resolve_target() {
  local dst="$1"
  case "$dst" in
    /srv/neoag-assets/install/*) printf '%s/%s' "${REFERENCE_ROOT%/}" "${dst#/srv/neoag-assets/install/}" ;;
    /srv/neoag-tools/*) printf '%s/%s' "${TOOLS_ROOT%/}" "${dst#/srv/neoag-tools/}" ;;
    /srv/neoag-licensed/*) printf '%s/%s' "${LICENSED_ROOT%/}" "${dst#/srv/neoag-licensed/}" ;;
    *) printf '%s' "$dst" ;;
  esac
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
rsync_transport() {
  if [[ -n "$SSH_KEY" ]]; then
    printf '%s' "-e 'ssh -i $(printf '%q' "$SSH_KEY") -o BatchMode=yes -o ConnectTimeout=15'"
  else
    printf '%s' "-e 'ssh -o BatchMode=yes -o ConnectTimeout=15'"
  fi
}

RSYNC_TRANSPORT="$(rsync_transport)"

rsync_flags() {
  local name="$1" src="$2" dst="$3"
  case "$dst:$name:$src" in
    "$LICENSED_ROOT"/*:*|*polysolver*|*licensed*)
      # Licensed tool trees may contain vendor binaries or links that are not
      # readable through their source targets. Preserve those links instead of
      # dereferencing them with -L.
      printf '%s' "-a"
      ;;
    *)
      # Reference assets commonly contain absolute symlinks that only resolve on
      # the source host, so keep dereferencing for non-licensed assets.
      printf '%s' "-aL"
      ;;
  esac
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
  src="$(resolve_shared_source "$src")"
  dst="$(resolve_target "$dst")"

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
  if [[ -n "$SHARED_ASSET_ROOT" ]]; then
    if [[ ! -e "$src" ]]; then
      status="source_missing"
      detail="shared source does not exist"
      echo -e "$name\t$kind\t$required\t$src\t$dst\t$status\t$detail" >> "$REPORT"
      log "$name: $status ($src)"
      is_truthy "$required" && { echo "ASSET_SOURCE_MISSING: $name: $src" >&2; exit 54; }
      continue
    fi
    if [[ -e "$dst" || -L "$dst" ]]; then
      status="target_conflict"
      detail="refusing to replace existing target"
      echo -e "$name\t$kind\t$required\t$src\t$dst\t$status\t$detail" >> "$REPORT"
      echo "ASSET_TARGET_CONFLICT: $name: $dst" >&2
      is_truthy "$required" && exit 56
      continue
    fi
    cmd="mkdir -p '$(dirname "$dst")' && ln -s '$src' '$dst'"
  elif [[ "$kind" == "file" ]]; then
    rsync_opts="$(rsync_flags "$name" "$src" "$dst")"
    cmd="mkdir -p '$(dirname "$dst")' && rsync $rsync_opts $RSYNC_TRANSPORT '$spec' '$dst'"
  else
    rsync_opts="$(rsync_flags "$name" "$src" "$dst")"
    cmd="mkdir -p '$dst' && rsync $rsync_opts $RSYNC_TRANSPORT '$spec/' '$dst/'"
  fi
  log ""
  log "==> [$MODE] sync asset $name"
  log "+ $cmd"
  if [[ "$EXECUTE" == "1" ]]; then
    if bash -lc "$cmd" 2>&1 | tee -a "$LOG"; then
      if target_has_marker "$kind" "$dst" "$marker" && verify_sha256 "$kind" "$dst" "$sha" >/dev/null 2>&1; then
        status="synced"
        if [[ -n "$SHARED_ASSET_ROOT" ]]; then detail="link and verification completed"; else detail="copy and verification completed"; fi
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
