process REPORT_V041 {
  tag "report-v041"
  label 'small'
  container params.neoag_container ?: null
  publishDir "${params.outdir}/reports", mode: 'copy'

  input:
    path ranked_events
    path ranked_peptides
    path appm_summary
    path validation_plan
    path appm_gene_status
    path appm_module_scores
    path appm_submodule_scores
    path appm_conflicts
    path appm_peptide_modifiers
    path immune_escape_summary
    path peptide_escape_flags
    path peptide_safety
    path ccf_lite
    val profile_name

  output:
    path "evidence_report.v041.html", emit: report
    path "versions.yml", emit: versions

  script:
  """
  neoag-v03 report-v041 \
    --profile '${profile_name}' \
    --ranked-events '${ranked_events}' \
    --ranked-peptides '${ranked_peptides}' \
    --appm-summary '${appm_summary}' \
    --appm-gene-status '${appm_gene_status}' \
    --appm-module-scores '${appm_module_scores}' \
    --appm-submodule-scores '${appm_submodule_scores}' \
    --appm-conflicts '${appm_conflicts}' \
    --appm-peptide-modifiers '${appm_peptide_modifiers}' \
    --immune-escape-summary '${immune_escape_summary}' \
    --peptide-escape-flags '${peptide_escape_flags}' \
    --peptide-safety '${peptide_safety}' \
    --ccf '${ccf_lite}' \
    --validation-plan '${validation_plan}' \
    --out evidence_report.v041.html

  echo "REPORT_V041:" > versions.yml
  echo "  neoag-v03: \$(python -c 'import neoag_v03; print(neoag_v03.__version__)')" >> versions.yml
  """
}
