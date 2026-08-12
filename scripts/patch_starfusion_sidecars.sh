#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 --star-fusion-home DIR [--source DIR]" >&2
}

STAR_FUSION_HOME=""
SOURCE_ROOT="${NEOAG_STAR_FUSION_SIDECAR_SRC:-}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --star-fusion-home) STAR_FUSION_HOME="$2"; shift 2 ;;
    --source) SOURCE_ROOT="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done

[[ -n "$STAR_FUSION_HOME" && -d "$STAR_FUSION_HOME" ]] || {
  echo "ERROR: STAR-Fusion home missing: ${STAR_FUSION_HOME:-<unset>}" >&2
  exit 2
}

PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
if [[ -z "$SOURCE_ROOT" ]]; then
  shopt -s nullglob
  for candidate in \
    "$PROJECT_ROOT/../open-neo-deploy/env_tool/conda_pkgs"/star-fusion-*/lib/STAR-Fusion \
    "$PROJECT_ROOT/../open-neo-deploy/env_tool/conda_pkgs"/star-fusion-*/lib/STAR-Fusion/FusionInspector; do
    if [[ -x "$candidate/FusionFilter/blast_and_promiscuity_filter.pl" && -x "$candidate/FusionAnnotator/FusionAnnotator" ]]; then
      SOURCE_ROOT="$candidate"
      break
    fi
  done
  shopt -u nullglob
fi

need_filter=0
need_annotator=0
[[ -x "$STAR_FUSION_HOME/FusionFilter/blast_and_promiscuity_filter.pl" ]] || need_filter=1
[[ -x "$STAR_FUSION_HOME/FusionAnnotator/FusionAnnotator" ]] || need_annotator=1

if [[ "$need_filter" -eq 1 || "$need_annotator" -eq 1 ]]; then
  [[ -n "$SOURCE_ROOT" && -d "$SOURCE_ROOT" ]] || {
    echo "ERROR: STAR-Fusion sidecar source not found; set NEOAG_STAR_FUSION_SIDECAR_SRC" >&2
    exit 3
  }
  [[ "$need_filter" -eq 0 ]] || cp -a "$SOURCE_ROOT/FusionFilter/." "$STAR_FUSION_HOME/FusionFilter/"
  [[ "$need_annotator" -eq 0 ]] || cp -a "$SOURCE_ROOT/FusionAnnotator/." "$STAR_FUSION_HOME/FusionAnnotator/"
fi

chmod +x "$STAR_FUSION_HOME/FusionFilter/"*.pl "$STAR_FUSION_HOME/FusionFilter/util/"*.pl 2>/dev/null || true
chmod +x "$STAR_FUSION_HOME/FusionAnnotator/FusionAnnotator" "$STAR_FUSION_HOME/FusionAnnotator/util/"*.pl 2>/dev/null || true

filter_script="$STAR_FUSION_HOME/FusionFilter/util/filter_by_annotation_rules.pl"
[[ -s "$filter_script" ]] || { echo "ERROR: missing $filter_script" >&2; exit 3; }
if ! grep -q 'annot_filter.pass' "$filter_script"; then
  perl -0pi -e 's/my \$pass_file = "\$predictions_file\.pass";/my \$pass_file = "\$predictions_file.annot_filter.pass";/' "$filter_script"
fi
perl -c "$filter_script" >/dev/null

test -x "$STAR_FUSION_HOME/FusionFilter/blast_and_promiscuity_filter.pl"
test -x "$STAR_FUSION_HOME/FusionAnnotator/FusionAnnotator"
