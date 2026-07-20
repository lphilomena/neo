// @Deprecated: Use APPM_2 (modules/appm_2/main.nf) instead — APPM_2 provides
// superset functionality including CNV/purity integration and expanded outputs.
// This module is kept for backward compatibility and will be removed in a
// future release.
process APPM_LITE {
  tag "$sample_id"
  publishDir "${params.outdir}/appm", mode: 'copy'

  input:
    val sample_id
    val profile_name
    path vep_appm_file
    path expression_file
    path hla_loh_file

  output:
    path "appm_lite.tsv", emit: appm_lite
    path "appm_summary.tsv", emit: appm_summary
    path "versions.yml", emit: versions

  script:
  """
  neoag-v03 appm-lite \
    --sample-id '${sample_id}' \
    --profile '${profile_name}' \
    --vep-tsv '${vep_appm_file}' \
    --expression '${expression_file}' \
    --hla-loh '${hla_loh_file}' \
    --outdir .

  echo "APPM_LITE:" > versions.yml
  echo "  neoag-v03: \$(python -c 'import neoag_v03; print(neoag_v03.__version__)')" >> versions.yml
  """
}
