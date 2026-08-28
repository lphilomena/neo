#!/usr/bin/env bash
set -euo pipefail
RELEASE=""
SHA256=""
DEST="."
while [[ $# -gt 0 ]]; do
  case "$1" in
    --release) RELEASE="$2"; shift 2 ;;
    --sha256) SHA256="$2"; shift 2 ;;
    --dest) DEST="$2"; shift 2 ;;
    -h|--help) echo "Usage: $0 --release file.tar.gz [--sha256 file.sha256] [--dest DIR]"; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done
[[ -n "$RELEASE" && -f "$RELEASE" ]] || { echo "CHECKSUM_FAILED: release missing" >&2; exit 10; }
if [[ -n "$SHA256" ]]; then
  [[ -f "$SHA256" ]] || { echo "CHECKSUM_FAILED: sha256 file missing" >&2; exit 10; }
  sha256sum -c "$SHA256" || { echo "CHECKSUM_FAILED" >&2; exit 10; }
fi
mkdir -p "$DEST"
tar -xzf "$RELEASE" -C "$DEST"
echo "unpacked_to=$DEST"
