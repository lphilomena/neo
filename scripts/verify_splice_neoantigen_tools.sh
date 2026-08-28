#!/usr/bin/env bash
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONDA_BASE="${NEOAG_CONDA_BASE:-$(conda info --base 2>/dev/null || true)}"
ENV_ROOT="${NEOAG_ENV_ROOT:-${CONDA_BASE}/envs}"
TOOLS_ROOT="${NEOAG_SPLICE_TOOLS_ROOT:-${NEOAG_TOOLS_ROOT:-${ROOT}}/tools/splice}"
ARCHIVE_ROOT="${NEOAG_CONTAINER_ARCHIVE_ROOT:-${ROOT}/work/container_images}"

printf 'tool\tstatus\tdetail\n'
if [[ -x "${ENV_ROOT}/neoag-splice/bin/regtools" ]]; then
  printf 'RegTools\tREADY\t%s\n' "${ENV_ROOT}/neoag-splice/bin/regtools"
else printf 'RegTools\tBLOCKED\tmissing executable\n'; fi
pvac=""
for candidate in \
  "${NEOAG_PVACSPLICE_BIN:-}" \
  "$(command -v pvacsplice 2>/dev/null || true)" \
  "${CONDA_BASE}/envs/neoag-tools/bin/pvacsplice" \
  "${HOME}/miniforge3/envs/neoag-tools/bin/pvacsplice" \
  "${HOME}/miniconda3/envs/neoag-tools/bin/pvacsplice"
do
  if [[ -n "${candidate}" && -x "${candidate}" ]]; then
    pvac="${candidate}"
    break
  fi
done
if [[ -n "${pvac}" ]]; then printf 'pVACsplice\tREADY\t%s\n' "${pvac}"; else printf 'pVACsplice\tBLOCKED\tmissing executable\n'; fi
if "${ENV_ROOT}/neoag-snaf/bin/python" -c 'import snaf' >/dev/null 2>&1; then printf 'SNAF\tREADY\t%s\n' "${ENV_ROOT}/neoag-snaf"; else printf 'SNAF\tPARTIAL\tenvironment/source incomplete\n'; fi
if [[ -f "${TOOLS_ROOT}/ASNEO/ASNEO.py" ]]; then printf 'ASNEO\tPARTIAL\tGRCh37-only source staged\n'; else printf 'ASNEO\tBLOCKED\tsource missing\n'; fi
if [[ -f "${TOOLS_ROOT}/NeoSplice-0.0.3/README.md" ]]; then printf 'NeoSplice\tPARTIAL\tmatched tumor/normal RNA required\n'; else printf 'NeoSplice\tBLOCKED\tsource missing or incomplete\n'; fi
if docker image inspect ghcr.io/tron-bioinformatics/splice2neo:v0.6.14 >/dev/null 2>&1; then
  printf 'splice2neo\tREADY\tcontainer image present\n'
elif [[ -s "${ARCHIVE_ROOT}/splice2neo_v0.6.14.tar" ]]; then printf 'splice2neo\tPARTIAL\timage archive present but not loaded\n'
else printf 'splice2neo\tBLOCKED\timage missing\n'; fi
splicemutr="${TOOLS_ROOT}/SpliceMutr-main"
if [[ -f "${splicemutr}/Rscripts/form_transcripts.R" \
   && -x "${ENV_ROOT}/neoag-splicemutr/bin/Rscript" \
   && -x "${ENV_ROOT}/neoag-splicemutr-leafcutter/bin/Rscript" ]]; then
  printf 'SpliceMutr\tREADY\t%s\n' "${splicemutr}"
elif [[ -f "${splicemutr}/Rscripts/form_transcripts.R" ]]; then
  printf 'SpliceMutr\tPARTIAL\tsource staged; environments incomplete\n'
else
  printf 'SpliceMutr\tBLOCKED\tsource missing\n'
fi
