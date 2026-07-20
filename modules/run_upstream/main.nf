process RUN_UPSTREAM {
  tag "$sample_id"
  publishDir "${params.outdir}", mode: 'copy'

  input:
    val sample_id
    path run_config

  output:
    path "upstream/**", emit: upstream_dir
    path "upstream/tools/*", emit: tool_outputs, optional: true

  script:
  """
  neoag run-upstream --config '${run_config}' --outdir upstream
  """
}
