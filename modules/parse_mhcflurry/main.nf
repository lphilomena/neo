process PARSE_MHCFLURRY {
  tag "$sample_id"
  publishDir "${params.outdir}/presentation", mode: 'copy'

  input:
    val sample_id
    path mhcflurry_file

  output:
    path "mhcflurry_evidence.tsv", emit: mhcflurry_evidence
    path "versions.yml", emit: versions

  script:
  """
  neoag parse-mhcflurry \
    --sample-id '${sample_id}' \
    --input '${mhcflurry_file}' \
    --out mhcflurry_evidence.tsv

  echo "PARSE_MHCFLURRY:" > versions.yml
  echo "  neoag: \$(python -c 'import neoag; print(neoag.__version__)')" >> versions.yml
  """
}
