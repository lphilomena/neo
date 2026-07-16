#!/usr/bin/env bash
set -euo pipefail
OUTDIR="work/remote_deploy"
LICENSED_ROOT="/root/neo/licensed_tools"
NETMHCPAN_TAR=""
NETMHCPAN_DIR=""
NETMHCPAN_URL=""
MIXMHCPRED_DIR=""
MIXMHCPRED_ARCHIVE=""
MIXMHCPRED_URL=""
NETMHCSTABPAN_DIR=""
NETMHCSTABPAN_ARCHIVE=""
NETMHCSTABPAN_URL=""
EXECUTE=0
ALLOW_DOWNLOAD=0

usage() {
  cat <<'USAGE'
Usage: 12_install_local_licensed_tools.sh [options]

Install licensed/restricted tools from source files, directories, or approved
URLs visible to the target machine. The installer copies/extracts assets into
--licensed-root and never creates target symlinks to /mnt, /home, or old-machine
paths.

Options:
  --licensed-root DIR          Target licensed-tool root (default: /root/neo/licensed_tools)
  --netmhcpan-tar FILE         Local NetMHCpan .tar.gz/.tgz/.tar archive
  --netmhcpan-dir DIR          Existing NetMHCpan install directory to copy
  --netmhcpan-url URL          Approved NetMHCpan archive URL to download
  --mixmhcpred-dir DIR         Existing MixMHCpred install directory to copy
  --mixmhcpred-archive FILE    Local MixMHCpred .zip/.tar.gz/.tgz/.tar archive
  --mixmhcpred-url URL         Approved MixMHCpred archive URL to download
  --netmhcstabpan-dir DIR      Existing NetMHCstabpan install directory to copy
  --netmhcstabpan-archive FILE Local NetMHCstabpan archive
  --netmhcstabpan-url URL      Approved NetMHCstabpan archive URL to download
  --allow-download             Permit network download from supplied URLs
  --outdir DIR                 Output directory for logs/downloads (default: work/remote_deploy)
  --execute                    Actually write files; default is dry-run
  -h, --help                   Show this help

Examples:
  bash 12_install_local_licensed_tools.sh \
    --netmhcpan-tar /mnt/.../netMHCpan-4.2c.Linux.tar.gz \
    --mixmhcpred-dir /mnt/.../tools/mixMHCpred_install \
    --licensed-root /root/neo/licensed_tools \
    --execute

  bash 12_install_local_licensed_tools.sh \
    --mixmhcpred-url https://approved.example/MixMHCpred.zip \
    --allow-download --execute

For licensed tools, use only official or user-approved URLs and verify the local
license permits installation on this machine. This script does not bypass login,
registration, click-through licenses, or institutional access controls.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --licensed-root) LICENSED_ROOT="$2"; shift 2 ;;
    --netmhcpan-tar) NETMHCPAN_TAR="$2"; shift 2 ;;
    --netmhcpan-dir) NETMHCPAN_DIR="$2"; shift 2 ;;
    --netmhcpan-url) NETMHCPAN_URL="$2"; shift 2 ;;
    --mixmhcpred-dir) MIXMHCPRED_DIR="$2"; shift 2 ;;
    --mixmhcpred-archive) MIXMHCPRED_ARCHIVE="$2"; shift 2 ;;
    --mixmhcpred-url) MIXMHCPRED_URL="$2"; shift 2 ;;
    --netmhcstabpan-dir) NETMHCSTABPAN_DIR="$2"; shift 2 ;;
    --netmhcstabpan-archive) NETMHCSTABPAN_ARCHIVE="$2"; shift 2 ;;
    --netmhcstabpan-url) NETMHCSTABPAN_URL="$2"; shift 2 ;;
    --allow-download) ALLOW_DOWNLOAD=1; shift ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --execute) EXECUTE=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

mkdir -p "$OUTDIR"
LOG="$OUTDIR/local_licensed_tool_install.log"
REPORT="$OUTDIR/local_licensed_tool_install_report.md"
DOWNLOAD_DIR="$OUTDIR/downloads"
: > "$LOG"
MODE="DRY_RUN"
[[ "$EXECUTE" == "1" ]] && MODE="EXECUTE"

log() { printf '%s\n' "$*" | tee -a "$LOG"; }
backup_target() {
  local dst="$1"
  if [[ -e "$dst" || -L "$dst" ]]; then
    local bak="${dst}.bak_$(date +%Y%m%d_%H%M%S)"
    log "backup $dst -> $bak"
    [[ "$EXECUTE" == "1" ]] && mv "$dst" "$bak"
  fi
}
copy_dir() {
  local name="$1" src="$2" dst="$3"
  [[ -n "$src" ]] || return 0
  [[ -d "$src" ]] || { echo "ASSET_SOURCE_MISSING: $name source is not a directory: $src" >&2; exit 20; }
  log "[$MODE] install $name directory: $src -> $dst"
  if [[ "$EXECUTE" == "1" ]]; then
    mkdir -p "$(dirname "$dst")"
    backup_target "$dst"
    mkdir -p "$dst"
    rsync -aH --copy-links --delete "$src"/ "$dst"/ | tee -a "$LOG"
  fi
}
archive_ext() {
  case "$1" in
    *.tar.gz|*.tgz) echo tar.gz ;;
    *.tar) echo tar ;;
    *.zip) echo zip ;;
    *) echo unknown ;;
  esac
}
extract_archive() {
  local name="$1" src="$2" dst="$3"
  [[ -n "$src" ]] || return 0
  [[ -f "$src" ]] || { echo "ASSET_SOURCE_MISSING: $name archive is not a file: $src" >&2; exit 20; }
  local kind; kind="$(archive_ext "$src")"
  [[ "$kind" != "unknown" ]] || { echo "UNSUPPORTED_ARCHIVE: $src" >&2; exit 22; }
  log "[$MODE] install $name archive: $src -> $dst"
  if [[ "$EXECUTE" == "1" ]]; then
    mkdir -p "$(dirname "$dst")"
    backup_target "$dst"
    mkdir -p "$dst"
    case "$kind" in
      tar.gz) tar -xzf "$src" -C "$dst" --strip-components=1 ;;
      tar) tar -xf "$src" -C "$dst" --strip-components=1 ;;
      zip)
        local tmp
        tmp="$(mktemp -d)"
        python3 - "$src" "$tmp" <<'PY'
import sys, zipfile
zipfile.ZipFile(sys.argv[1]).extractall(sys.argv[2])
PY
        local entries
        entries="$(find "$tmp" -mindepth 1 -maxdepth 1 | wc -l)"
        if [[ "$entries" == "1" && -d "$(find "$tmp" -mindepth 1 -maxdepth 1)" ]]; then
          rsync -aH --copy-links --delete "$(find "$tmp" -mindepth 1 -maxdepth 1)"/ "$dst"/ | tee -a "$LOG"
        else
          rsync -aH --copy-links --delete "$tmp"/ "$dst"/ | tee -a "$LOG"
        fi
        rm -rf "$tmp"
        ;;
    esac
  fi
}
download_url() {
  local name="$1" url="$2"
  [[ -n "$url" ]] || return 0
  [[ "$ALLOW_DOWNLOAD" == "1" ]] || { echo "DOWNLOAD_NOT_APPROVED: $name URL provided but --allow-download is missing" >&2; exit 23; }
  mkdir -p "$DOWNLOAD_DIR"
  local base dest
  base="$(basename "${url%%\?*}")"
  [[ -n "$base" && "$base" != "/" ]] || base="$name.download"
  dest="$DOWNLOAD_DIR/$base"
  log "[$MODE] download $name: $url -> $dest"
  if [[ "$EXECUTE" == "1" ]]; then
    if command -v curl >/dev/null 2>&1; then
      curl -fL --retry 3 --connect-timeout 20 -o "$dest" "$url"
    elif command -v wget >/dev/null 2>&1; then
      wget -O "$dest" "$url"
    else
      echo "DOWNLOAD_TOOL_MISSING: curl or wget is required" >&2
      exit 24
    fi
  fi
  printf '%s\n' "$dest"
}

if [[ -z "$NETMHCPAN_TAR" && -z "$NETMHCPAN_DIR" && -z "$NETMHCPAN_URL" && -z "$MIXMHCPRED_DIR" && -z "$MIXMHCPRED_ARCHIVE" && -z "$MIXMHCPRED_URL" && -z "$NETMHCSTABPAN_DIR" && -z "$NETMHCSTABPAN_ARCHIVE" && -z "$NETMHCSTABPAN_URL" ]]; then
  echo "ASSET_SOURCE_MISSING: provide a local source or approved URL. If the package is not local, the agent must retrieve an official/user-approved URL first." >&2
  exit 20
fi

mkdir -p "$LICENSED_ROOT"
[[ -n "$NETMHCPAN_URL" ]] && NETMHCPAN_TAR="$(download_url netMHCpan "$NETMHCPAN_URL")"
[[ -n "$MIXMHCPRED_URL" ]] && MIXMHCPRED_ARCHIVE="$(download_url MixMHCpred "$MIXMHCPRED_URL")"
[[ -n "$NETMHCSTABPAN_URL" ]] && NETMHCSTABPAN_ARCHIVE="$(download_url NetMHCstabpan "$NETMHCSTABPAN_URL")"

extract_archive "netMHCpan" "$NETMHCPAN_TAR" "$LICENSED_ROOT/netMHCpan"
copy_dir "netMHCpan" "$NETMHCPAN_DIR" "$LICENSED_ROOT/netMHCpan"
extract_archive "MixMHCpred" "$MIXMHCPRED_ARCHIVE" "$LICENSED_ROOT/mixMHCpred_install"
copy_dir "MixMHCpred" "$MIXMHCPRED_DIR" "$LICENSED_ROOT/mixMHCpred_install"
extract_archive "NetMHCstabpan" "$NETMHCSTABPAN_ARCHIVE" "$LICENSED_ROOT/netMHCstabpan"
copy_dir "NetMHCstabpan" "$NETMHCSTABPAN_DIR" "$LICENSED_ROOT/netMHCstabpan"

{
  echo "# Local licensed tool install report"
  echo
  echo "Mode: \`$MODE\`"
  echo "Licensed root: \`$LICENSED_ROOT\`"
  echo "Log: \`$LOG\`"
  echo "Downloads: \`$DOWNLOAD_DIR\`"
  echo
  echo "Installed sources:"
  [[ -n "$NETMHCPAN_TAR" ]] && echo "- NetMHCpan archive: \`$NETMHCPAN_TAR\`"
  [[ -n "$NETMHCPAN_DIR" ]] && echo "- NetMHCpan directory: \`$NETMHCPAN_DIR\`"
  [[ -n "$NETMHCPAN_URL" ]] && echo "- NetMHCpan URL: \`$NETMHCPAN_URL\`"
  [[ -n "$MIXMHCPRED_ARCHIVE" ]] && echo "- MixMHCpred archive: \`$MIXMHCPRED_ARCHIVE\`"
  [[ -n "$MIXMHCPRED_DIR" ]] && echo "- MixMHCpred directory: \`$MIXMHCPRED_DIR\`"
  [[ -n "$MIXMHCPRED_URL" ]] && echo "- MixMHCpred URL: \`$MIXMHCPRED_URL\`"
  [[ -n "$NETMHCSTABPAN_ARCHIVE" ]] && echo "- NetMHCstabpan archive: \`$NETMHCSTABPAN_ARCHIVE\`"
  [[ -n "$NETMHCSTABPAN_DIR" ]] && echo "- NetMHCstabpan directory: \`$NETMHCSTABPAN_DIR\`"
  [[ -n "$NETMHCSTABPAN_URL" ]] && echo "- NetMHCstabpan URL: \`$NETMHCSTABPAN_URL\`"
  echo
  if [[ "$EXECUTE" != "1" ]]; then
    echo "Dry run only. Re-run with \`--execute\` and \`--allow-download\` when applicable after license/transfer approval."
  else
    echo "Files were installed without target symlinks. Next run activation rewrite and runtime validation."
  fi
} > "$REPORT"
echo "local_licensed_tool_install_report=$REPORT"
