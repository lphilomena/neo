process CCF_LITE {
  tag "ccf-lite"
  publishDir "${params.outdir}/clonality", mode: 'copy'

  input:
    path raw_events
    val profile_name
    path purity_file
    path cnv_file

  output:
    path "ccf_lite.tsv", emit: ccf_lite
    path "versions.yml", emit: versions

  script:
  """
  neoag-v03 ccf-lite \
    --events '${raw_events}' \
    --profile '${profile_name}' \
    --purity '${purity_file}' \
    --cnv '${cnv_file}' \
    --out ccf_lite.tsv

  echo "CCF_LITE:" > versions.yml
  echo "  neoag-v03: \$(python -c 'import neoag_v03; print(neoag_v03.__version__)')" >> versions.yml
  """
}
