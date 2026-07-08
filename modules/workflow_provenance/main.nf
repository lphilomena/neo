process WORKFLOW_PROVENANCE {
  tag "workflow-provenance"
  publishDir "${params.outdir}/provenance", mode: 'copy'

  input:
    path version_files, stageAs: 'versions_*.yml'

  output:
    path "workflow_provenance.yml", emit: workflow_provenance

  script:
  """
  echo "workflow_provenance:" > workflow_provenance.yml
  echo "  generated_by: Nextflow" >> workflow_provenance.yml
  echo "  modules:" >> workflow_provenance.yml
  for f in ${version_files}; do
    echo "    - file: \$(basename "\$f")" >> workflow_provenance.yml
    sed 's/^/      /' "\$f" >> workflow_provenance.yml
  done
  """
}
