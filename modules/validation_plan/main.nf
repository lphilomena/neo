process VALIDATION_PLAN_V03 {
  tag "validation-plan"
  publishDir "${params.outdir}/scoring", mode: 'copy'

  input:
    path ranked_peptides

  output:
    path "validation_plan.v03.tsv", emit: validation_plan
    path "versions.yml", emit: versions

  script:
  """
  neoag-v03 validation-plan-v03 \
    --ranked-peptides '${ranked_peptides}' \
    --outdir '${params.outdir}' \
    --out validation_plan.v03.tsv

  echo "VALIDATION_PLAN_V03:" > versions.yml
  echo "  neoag-v03: \$(python -c 'import neoag_v03; print(neoag_v03.__version__)')" >> versions.yml
  """
}
