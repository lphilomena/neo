process VALIDATION_PLAN {
  tag "validation-plan"
  publishDir "${params.outdir}/scoring", mode: 'copy'

  input:
    path ranked_peptides

  output:
    path "validation_plan.tsv", emit: validation_plan
    path "versions.yml", emit: versions

  script:
  """
  neoag validation-plan \
    --ranked-peptides '${ranked_peptides}' \
    --outdir '${params.outdir}' \
    --out validation_plan.tsv

  echo "VALIDATION_PLAN:" > versions.yml
  echo "  neoag: \$(python -c 'import neoag; print(neoag.__version__)')" >> versions.yml
  """
}
