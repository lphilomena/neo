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
  neoag appm-lite \
    --sample-id '${sample_id}' \
    --profile '${profile_name}' \
    --vep-tsv '${vep_appm_file}' \
    --expression '${expression_file}' \
    --hla-loh '${hla_loh_file}' \
    --outdir .

  echo "APPM_LITE:" > versions.yml
  echo "  neoag: \$(python -c 'import neoag; print(neoag.__version__)')" >> versions.yml
  """
}
