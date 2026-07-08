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
  neoag-v03 parse-netmhcpan \
    --sample-id '${sample_id}' \
    --input '${netmhcpan_file}' \
    --out netmhcpan_evidence.tsv

  echo "PARSE_NETMHCPAN:" > versions.yml
  echo "  neoag-v03: \$(python -c 'import neoag_v03; print(neoag_v03.__version__)')" >> versions.yml
  """
}
