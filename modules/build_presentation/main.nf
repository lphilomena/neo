process BUILD_PRESENTATION {
  tag "presentation"
  publishDir "${params.outdir}/presentation", mode: 'copy'

  input:
    path raw_peptides
    path netmhcpan_evidence
    path mhcflurry_evidence
    val profile_name

  output:
    path "presentation_evidence.tsv", emit: presentation_evidence
    path "versions.yml", emit: versions

  script:
  """
  neoag build-presentation-evidence \
    --raw-peptides '${raw_peptides}' \
    --netmhcpan '${netmhcpan_evidence}' \
    --mhcflurry '${mhcflurry_evidence}' \
    --profile '${profile_name}' \
    --out presentation_evidence.tsv

  echo "BUILD_PRESENTATION:" > versions.yml
  echo "  neoag: \$(python -c 'import neoag; print(neoag.__version__)')" >> versions.yml
  """
}
