#!/usr/bin/env bash
set -euo pipefail

SOURCE="${NEOAG_SNAF_SOURCE:-}"
IMAGE="${NEOAG_ALTANALYZE_IMAGE:-neoag-altanalyze:snaf}"
DEBIAN_ARCHIVE="${NEOAG_DEBIAN_ARCHIVE:-http://archive.debian.org/debian}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --source) SOURCE="$2"; shift 2;;
    --image) IMAGE="$2"; shift 2;;
    *) echo "ERROR: unknown option $1" >&2; exit 2;;
  esac
done
[[ -d "$SOURCE/AltAnalyze" ]] && SOURCE="$SOURCE/AltAnalyze"
for file in Dockerfile AltAnalyze.sh prune.py Hs.bed requirements_slim_test.txt; do
  [[ -s "$SOURCE/$file" ]] || { echo "ERROR: missing $SOURCE/$file" >&2; exit 2; }
done
command -v docker >/dev/null || { echo "ERROR: docker is required" >&2; exit 2; }

context="$(mktemp -d)"
cleanup() { rm -rf "$context"; }
trap cleanup EXIT
cp -a "$SOURCE/." "$context/"

# The upstream Python 2 image is Debian Buster. Its live mirrors were retired,
# so use the immutable Debian archive and disable expired metadata checks.
awk '/^FROM[[:space:]]/ && !done {
  print
  print "RUN printf \047deb __DEBIAN_ARCHIVE__ buster main\\n\047 > /etc/apt/sources.list && rm -f /etc/apt/sources.list.d/* && printf \047Acquire::Check-Valid-Until \\\"false\\\";\\nAcquire::AllowInsecureRepositories \\\"true\\\";\\n\047 > /etc/apt/apt.conf.d/99archive"
  done = 1
  next
}
{ print }' "$SOURCE/Dockerfile" > "$context/Dockerfile"
sed -i.bak 's/apt-get install -y parallel/apt-get install -y --no-install-recommends parallel/' "$context/Dockerfile"
sed -i.bak "s|__DEBIAN_ARCHIVE__|$DEBIAN_ARCHIVE|g" "$context/Dockerfile"

docker build --pull=false -t "$IMAGE" "$context"
docker image inspect "$IMAGE" --format 'AltAnalyze image ready: {{.Id}} {{.Size}} bytes'
