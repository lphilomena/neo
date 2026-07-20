process CCF_2 {
  tag "ccf-2"
  label 'medium'
  publishDir "${params.outdir}/clonality", mode: 'copy'

  input:
    path raw_events
    val profile_name
    path purity_file
    path cnv_file

  output:
    path "ccf_2.tsv", emit: ccf_2
    path "ccf_lite.tsv", emit: ccf_lite
    path "ccf_input_qc.tsv", emit: ccf_input_qc
    path "ccf_conflicts.tsv", emit: ccf_conflicts
    path "ccf_cluster.tsv", emit: ccf_cluster
    path "versions.yml", emit: versions

  script:
  """
  neoag-v03 ccf-2 \
    --events '${raw_events}' \
    --profile '${profile_name}' \
    --purity '${purity_file}' \
    --cnv '${cnv_file}' \
    --out ccf_2.tsv

  cp ccf_2.tsv ccf_lite.tsv

  echo "CCF_2:" > versions.yml
  echo "  neoag-v03: \$(python -c 'import neoag_v03; print(neoag_v03.__version__)')" >> versions.yml
  """
}
