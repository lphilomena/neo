process REPORT_V03 {
  tag "report-v03"
  publishDir "${params.outdir}/reports", mode: 'copy'

  input:
    path ranked_events
    path ranked_peptides
    path appm_summary
    path validation_plan
    val profile_name

  output:
    path "evidence_report.v03.html", emit: report
    path "evidence_report.patient.html", emit: patient_report
    path "evidence_report.technical.html", emit: technical_report
    path "versions.yml", emit: versions

  script:
  """
  neoag-v03 report-v03 \
    --profile '${profile_name}' \
    --ranked-events '${ranked_events}' \
    --ranked-peptides '${ranked_peptides}' \
    --appm-summary '${appm_summary}' \
    --validation-plan '${validation_plan}' \
    --outdir '${params.outdir}' \
    --audience both \
    --out evidence_report.v03.html

  echo "REPORT_V03:" > versions.yml
  echo "  neoag-v03: \$(python -c 'import neoag_v03; print(neoag_v03.__version__)')" >> versions.yml
  """
}
