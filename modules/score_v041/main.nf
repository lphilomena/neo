process SCORE_V041 {
  tag "score-v041"
  publishDir "${params.outdir}/scoring", mode: 'copy'

  input:
    path raw_events
    path raw_peptides
    path presentation_evidence
    path appm_summary
    path ccf_lite
    path normal_expression
    path normal_hla_ligands
    path peptide_safety
    path peptide_escape_flags
    path appm_peptide_modifiers
    val profile_name

  output:
    path "ranked_events.tsv", emit: ranked_events
    path "ranked_peptides.tsv", emit: ranked_peptides
    path "versions.yml", emit: versions

  script:
  """
  neoag score \
    --raw-events '${raw_events}' \
    --raw-peptides '${raw_peptides}' \
    --presentation '${presentation_evidence}' \
    --appm-summary '${appm_summary}' \
    --ccf '${ccf_lite}' \
    --normal-expression '${normal_expression}' \
    --normal-hla-ligands '${normal_hla_ligands}' \
    --peptide-safety '${peptide_safety}' \
    --peptide-escape-flags '${peptide_escape_flags}' \
    --appm-peptide-modifiers '${appm_peptide_modifiers}' \
    --profile '${profile_name}' \
    --out-events ranked_events.tsv \
    --out-peptides ranked_peptides.tsv

  echo "SCORE_V041:" > versions.yml
  echo "  neoag: \$(python -c 'import neoag; print(neoag.__version__)')" >> versions.yml
  """
}
