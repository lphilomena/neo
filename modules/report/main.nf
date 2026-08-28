process REPORT {
  tag "report"
  publishDir "${params.outdir}/reports", mode: 'copy'

  input:
    path ranked_events
    path ranked_peptides
    path appm_summary
    path validation_plan
    val profile_name

  output:
    path "evidence_report.html", emit: report
    path "evidence_report.patient.html", emit: patient_report
    path "evidence_report.technical.html", emit: technical_report
    path "versions.yml", emit: versions

  script:
  """
  neoag report \
    --profile '${profile_name}' \
    --ranked-events '${ranked_events}' \
    --ranked-peptides '${ranked_peptides}' \
    --appm-summary '${appm_summary}' \
    --validation-plan '${validation_plan}' \
    --outdir '${params.outdir}' \
    --audience both \
    --out evidence_report.html

  echo "REPORT:" > versions.yml
  echo "  neoag: \$(python -c 'import neoag; print(neoag.__version__)')" >> versions.yml
  """
}
