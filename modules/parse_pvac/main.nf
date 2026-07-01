process PARSE_PVAC {
  tag "$sample_id"
  publishDir "${params.outdir}/parsed", mode: 'copy'

  input:
    val sample_id
    val profile_name
    path pvac_files

  output:
    path "raw_events.tsv", emit: raw_events
    path "raw_peptides.tsv", emit: raw_peptides
    path "versions.yml", emit: versions

  script:
  def files = pvac_files.collect{ it.toString() }.join(' ')
  """
  neoag-v03 parse-pvac \
    --sample-id '${sample_id}' \
    --profile '${profile_name}' \
    --pvac ${files} \
    --events-out raw_events.tsv \
    --peptides-out raw_peptides.tsv

  echo "PARSE_PVAC:" > versions.yml
  echo "  neoag-v03: \$(python -c 'import neoag_v03; print(neoag_v03.__version__)')" >> versions.yml
  """
}
