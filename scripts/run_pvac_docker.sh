#!/usr/bin/env bash
# Run pVACseq/pvacfuse via official Docker image when conda install is not used.
# Example:
#   bash scripts/run_pvac_docker.sh pvacseq SAMPLE normal.vcf tumor.vcf "HLA-A*02:01" /outdir
#
set -euo pipefail

IMAGE="${NEOAG_PVAC_DOCKER:-griffithlab/pvactools:6.1.0}"
TOOL="${1:?tool: pvacseq|pvacfuse}"
shift

WORKDIR="${NEOAG_PVAC_WORKDIR:-$(pwd)/work/pvac_docker}"
mkdir -p "${WORKDIR}"

echo "==> Pulling ${IMAGE} (if needed) ..."
docker pull "${IMAGE}"

docker run --rm \
  -v "${WORKDIR}:/work" \
  -v "$(pwd):/data:ro" \
  "${IMAGE}" \
  "${TOOL}" "$@"
