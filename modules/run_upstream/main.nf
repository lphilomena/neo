process RUN_UPSTREAM {
  tag "$sample_id"
  label 'large'
  time '24.h'
  publishDir "${params.outdir}", mode: 'copy'

  input:
    val sample_id
    path run_config

  output:
    path "upstream", emit: upstream_dir
    path "upstream/tools/*", emit: tool_outputs, optional: true

  script:
  """
  neoag-v03 run-upstream --config '${run_config}' --outdir "\$PWD/upstream"
  """
}
