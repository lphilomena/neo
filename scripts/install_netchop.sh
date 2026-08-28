#!/usr/bin/env bash
set -euo pipefail

ARCHIVE="${1:-}"
LICENSED_ROOT="${NEOAG_LICENSED_ROOT:-/opt/neoag/licensed_tools}"
TOOLS_ROOT="${NEOAG_TOOLS_ROOT:-/opt/neoag/env_tool}"

[[ -n "$ARCHIVE" && -f "$ARCHIVE" ]] || {
  echo "Usage: install_netchop.sh /path/to/netchop-3.1d.Linux.tar.gz" >&2
  exit 2
}

home="$LICENSED_ROOT/netchop"
mkdir -p "$home" "$TOOLS_ROOT/bin"
tar -xzf "$ARCHIVE" -C "$home"
binary="$(find "$home" -type f -path '*/Linux_x86_64/bin/netChop' -print -quit)"
[[ -n "$binary" && -x "$binary" ]] || {
  echo "NetChop binary missing after extracting $ARCHIVE" >&2
  exit 3
}
ln -sfn "$binary" "$TOOLS_ROOT/bin/netChop"
if ! "$TOOLS_ROOT/bin/netChop" -h >/dev/null 2>&1; then
  test_fasta="$(find "$home" -type f -path '*/test/test.fsa' -print -quit)"
  [[ -n "$test_fasta" ]] && "$TOOLS_ROOT/bin/netChop" "$test_fasta" >/dev/null
fi
echo "NetChop 3.1d installed: $TOOLS_ROOT/bin/netChop"
