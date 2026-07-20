process PARSE_NETMHCPAN {
  tag "$sample_id"
  publishDir "${params.outdir}/presentation", mode: 'copy'

  input:
    val sample_id
    path netmhcpan_file

  output:
    path "netmhcpan_evidence.tsv", emit: netmhcpan_evidence
    path "versions.yml", emit: versions

  script:
  """
  neoag parse-netmhcpan \
    --sample-id '${sample_id}' \
    --input '${netmhcpan_file}' \
    --out netmhcpan_evidence.tsv

  echo "PARSE_NETMHCPAN:" > versions.yml
  echo "  neoag: \$(python -c 'import neoag; print(neoag.__version__)')" >> versions.yml
  """
}
