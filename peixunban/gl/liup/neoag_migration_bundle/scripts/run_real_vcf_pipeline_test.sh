#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=common.sh
source "${SCRIPT_DIR}/common.sh"

PRODUCTION_ACTIVATE="${NEOAG_TOOLS_ROOT:-${TOOLS_ROOT}}/activate_neoag_production_refs.sh"
if [[ -f "${PRODUCTION_ACTIVATE}" ]]; then
  # shellcheck source=/dev/null
  source "${PRODUCTION_ACTIVATE}"
else
  activate_paths
fi

REAL_VCF="${REAL_VCF:-/mnt/zjl-bgi-zzb/peixunban/gl/data/chenxiaoliang_data/data/liver_0520_WGS_shortReads/somatic/M1ML150017383_L01_438.align.somatic.pass.vcf.gz}"
TUMOR_SAMPLE_NAME="${TUMOR_SAMPLE_NAME:-M1ML150017383_L01_438}"
SAMPLE_ID="${SAMPLE_ID:-M1ML150017383_L01_438}"
REAL_VCF_OUTDIR="${REAL_VCF_OUTDIR:-${PROJECT_ROOT}/results/${SAMPLE_ID}_real_vcf_release_test}"
REAL_VCF_CONFIG="${REAL_VCF_CONFIG:-${PROJECT_ROOT}/conf/run.real_vcf.release_test.toml}"
NORMAL_EXPRESSION="${NORMAL_EXPRESSION:-resources/normal_expression.example.tsv}"
NORMAL_HLA_LIGANDS="${NORMAL_HLA_LIGANDS:-resources/normal_hla_ligands.example.tsv}"
HLA_ALLELES="${HLA_ALLELES:-HLA-A*02:06,HLA-A*30:01,HLA-B*13:02,HLA-C*06:02,HLA-C*08:01}"

if [[ ! -f "${REAL_VCF}" ]]; then
  cat >&2 <<MSG
ERROR: REAL_VCF not found: ${REAL_VCF}

Provide the real VCF test input or override it, for example:
  REAL_VCF=/path/to/sample.somatic.pass.vcf.gz \\
  TUMOR_SAMPLE_NAME=${TUMOR_SAMPLE_NAME} \\
  bash ${BASH_SOURCE[0]}

The default path is the validated internal test VCF used during the migration release check.
MSG
  exit 2
fi

if [[ ! -f "${REAL_VCF}.tbi" && ! -f "${REAL_VCF}.csi" ]]; then
  echo "WARN: no tabix index found next to ${REAL_VCF}; VEP/VCF tools may require .tbi or .csi" >&2
fi

python - <<PY
from pathlib import Path
hlas = [x.strip() for x in "${HLA_ALLELES}".split(",") if x.strip()]
config = Path("${REAL_VCF_CONFIG}")
config.parent.mkdir(parents=True, exist_ok=True)
quoted_hlas = ", ".join('"%s"' % h for h in hlas)
config.write_text(f"""[sample]
id = "${SAMPLE_ID}"
profile = "default"

[tools]
stub = false
enabled = ["netmhcpan", "mhcflurry"]
immunogenicity_stub = false

[inputs]
entry_mode = "snv_indel"
variant_peptide_extraction = true
variants_vcf = "${REAL_VCF}"
tumor_sample_name = "${TUMOR_SAMPLE_NAME}"
hla_alleles = [{quoted_hlas}]
extract_appm_from_vcf = false
normal_expression = "${NORMAL_EXPRESSION}"
normal_hla_ligands = "${NORMAL_HLA_LIGANDS}"
""")
print(config)
PY

cd "${PROJECT_ROOT}"
log "running real VCF release test"
log "config: ${REAL_VCF_CONFIG}"
log "outdir: ${REAL_VCF_OUTDIR}"
rm -rf "${REAL_VCF_OUTDIR}"
neoag-v03 run-full --config "${REAL_VCF_CONFIG}" --outdir "${REAL_VCF_OUTDIR}"

for path in \
  "${REAL_VCF_OUTDIR}/scoring/ranked_peptides.v03.tsv" \
  "${REAL_VCF_OUTDIR}/scoring/ranked_events.v03.tsv" \
  "${REAL_VCF_OUTDIR}/presentation/mhcflurry_evidence.tsv" \
  "${REAL_VCF_OUTDIR}/presentation/prime_evidence.tsv" \
  "${REAL_VCF_OUTDIR}/presentation/bigmhc_im_evidence.tsv" \
  "${REAL_VCF_OUTDIR}/reports/evidence_report.v03.html"; do
  [[ -s "${path}" ]] || { echo "ERROR: expected non-empty output missing: ${path}" >&2; exit 3; }
done

wc -l \
  "${REAL_VCF_OUTDIR}/scoring/ranked_peptides.v03.tsv" \
  "${REAL_VCF_OUTDIR}/scoring/ranked_events.v03.tsv" \
  "${REAL_VCF_OUTDIR}/presentation/presentation_evidence.tsv"
log "real VCF release test OK"
