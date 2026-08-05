#!/usr/bin/env bash
set -euo pipefail

OUTDIR="work/remote_deploy/claude_code"
INSTALLER_URL="https://claude.ai/install.sh"
CHANNEL="stable"
EXECUTE=0
ALLOW_DOWNLOAD=0
FORCE=0

usage() {
  cat <<'USAGE'
Usage: 17_install_claude_code.sh [options]

Install Claude Code with Anthropic's official native installer. The default is
the stable release channel. Installation is per-user and does not authenticate,
store an API key, or modify project settings.

Options:
  --outdir DIR              Log/report directory
  --channel VALUE           stable, latest, or an exact version (default: stable)
  --installer-url URL       Official/user-approved installer URL
  --allow-download          Approve the installer and binary download
  --force                   Re-run the installer when claude is already available
  --execute                 Perform installation; otherwise print a dry run
  -h, --help                Show this help
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2 ;;
    --channel) CHANNEL="$2"; shift 2 ;;
    --installer-url) INSTALLER_URL="$2"; shift 2 ;;
    --allow-download) ALLOW_DOWNLOAD=1; shift ;;
    --force) FORCE=1; shift ;;
    --execute) EXECUTE=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ ! "$CHANNEL" =~ ^(stable|latest|[0-9]+\.[0-9]+\.[0-9]+)$ ]]; then
  echo "CLAUDE_CODE_CHANNEL_INVALID: use stable, latest, or X.Y.Z" >&2
  exit 2
fi
if [[ "$INSTALLER_URL" != https://claude.ai/install.sh ]]; then
  echo "CLAUDE_CODE_INSTALLER_URL_REJECTED: expected https://claude.ai/install.sh" >&2
  exit 2
fi

mkdir -p "$OUTDIR"
LOG="$OUTDIR/claude_code_install.log"
REPORT="$OUTDIR/claude_code_install_report.md"
: > "$LOG"

log() { printf '%s\n' "$*" | tee -a "$LOG"; }
find_claude() {
  if command -v claude >/dev/null 2>&1; then
    command -v claude
  elif [[ -x "$HOME/.local/bin/claude" ]]; then
    printf '%s\n' "$HOME/.local/bin/claude"
  else
    return 1
  fi
}

MODE="DRY_RUN"
[[ "$EXECUTE" == "1" ]] && MODE="EXECUTE"
log "mode=$MODE"
log "channel=$CHANNEL"
log "installer_url=$INSTALLER_URL"

CLAUDE_BIN="$(find_claude 2>/dev/null || true)"
SKIP_INSTALL=0
if [[ -n "$CLAUDE_BIN" && "$FORCE" != "1" && "$CHANNEL" =~ ^[0-9] ]]; then
  current_version="$($CLAUDE_BIN --version 2>/dev/null || true)"
  [[ "$current_version" == "$CHANNEL "* ]] && SKIP_INSTALL=1
fi

if [[ "$SKIP_INSTALL" == "1" ]]; then
  log "Requested Claude Code version is already installed: $CLAUDE_BIN"
elif [[ "$EXECUTE" != "1" ]]; then
  log "[DRY_RUN] download the official installer and install Claude Code channel/version: $CHANNEL"
else
  [[ "$ALLOW_DOWNLOAD" == "1" ]] || {
    echo "DOWNLOAD_NOT_APPROVED: Claude Code requires --allow-download" >&2
    exit 23
  }
  command -v curl >/dev/null 2>&1 || { echo "CLAUDE_CODE_PREREQUISITE_MISSING: curl" >&2; exit 31; }
  command -v bash >/dev/null 2>&1 || { echo "CLAUDE_CODE_PREREQUISITE_MISSING: bash" >&2; exit 31; }
  case "$(uname -s)" in
    Linux|Darwin) ;;
    *) echo "CLAUDE_CODE_PLATFORM_UNSUPPORTED: $(uname -s)" >&2; exit 31 ;;
  esac

  INSTALLER="$OUTDIR/claude-code-install.sh"
  log "Downloading Anthropic installer to $INSTALLER"
  curl -fsSL --retry 3 "$INSTALLER_URL" -o "$INSTALLER" 2>&1 | tee -a "$LOG"
  chmod 0700 "$INSTALLER"
  log "Running official native installer for channel/version: $CHANNEL"
  bash "$INSTALLER" "$CHANNEL" 2>&1 | tee -a "$LOG"
  CLAUDE_BIN="$(find_claude 2>/dev/null || true)"
fi

VERSION="UNASSESSED"
if [[ -n "$CLAUDE_BIN" && -x "$CLAUDE_BIN" ]]; then
  VERSION="$($CLAUDE_BIN --version 2>&1)"
  log "claude_bin=$CLAUDE_BIN"
  log "claude_version=$VERSION"
elif [[ "$EXECUTE" == "1" ]]; then
  echo "CLAUDE_CODE_VERIFY_FAILED: claude binary not found after installation" >&2
  exit 42
fi

{
  echo "# Claude Code install report"
  echo
  echo "Mode: \`$MODE\`"
  echo "Release channel/version: \`$CHANNEL\`"
  echo "Installer: \`$INSTALLER_URL\`"
  echo "Binary: \`${CLAUDE_BIN:-not installed in dry-run}\`"
  echo "Version: \`$VERSION\`"
  echo "Log: \`$LOG\`"
  echo
  echo "Authentication was not attempted. Run \`claude\` interactively or configure an approved enterprise provider after deployment."
} > "$REPORT"

log "claude_code_install_report=$REPORT"
