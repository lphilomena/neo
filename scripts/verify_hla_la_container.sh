#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
[[ -f "$REPO_ROOT/conf/tools.env.sh" ]] && source "$REPO_ROOT/conf/tools.env.sh"

HLALA_HOME=${HLALA_HOME:-${HLA_LA_HOME:-${NEOAG_TOOLS_ROOT:-$REPO_ROOT}/tools/HLA-LA}}
HLALA_ENV_PREFIX=${HLALA_ENV_PREFIX:-${HLA_LA_ENV_PREFIX:-$HLALA_HOME/.conda}}
HLALA_BIN=${HLALA_BIN:-${HLA_LA_BIN:-$HLALA_ENV_PREFIX/bin/HLA-LA.pl}}
GRAPH=${HLALA_GRAPH:-${HLA_LA_GRAPH:-$HLALA_HOME/graphs/PRG_MHC_GRCh38_withIMGT}}

[[ -x "$HLALA_BIN" ]] || { echo "ERROR: real HLA-LA executable missing: $HLALA_BIN" >&2; exit 1; }
[[ -x "$HLALA_ENV_PREFIX/bin/perl" ]] || { echo "ERROR: HLA-LA Perl runtime missing" >&2; exit 1; }
[[ -x "$HLALA_ENV_PREFIX/bin/samtools" ]] || { echo "ERROR: HLA-LA samtools runtime missing" >&2; exit 1; }
[[ -s "$GRAPH/serializedGRAPH" ]] || { echo "ERROR: prepared HLA-LA graph marker missing: $GRAPH/serializedGRAPH" >&2; exit 1; }
[[ -s "$GRAPH/PRG/graph.txt" ]] || { echo "ERROR: HLA-LA PRG graph missing: $GRAPH/PRG/graph.txt" >&2; exit 1; }

VERIFY_WORK=${TMPDIR:-/tmp}/neoag_hlala_verify
mkdir -p "$VERIFY_WORK"
NEOAG_HLALA_BACKEND=native "$REPO_ROOT/scripts/run_hla_la_container.sh" \
  --testing 1 --workingDir "$VERIFY_WORK" --sampleID VERIFY >/dev/null
echo "PASS: HLA-LA executable starts: $HLALA_BIN"
echo "PASS: HLA-LA environment: $HLALA_ENV_PREFIX"
echo "PASS: prepared graph: $GRAPH"
