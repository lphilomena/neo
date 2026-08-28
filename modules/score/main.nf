process SCORE {
  tag "score"
  publishDir "${params.outdir}/scoring", mode: 'copy'

  input:
    path raw_events
    path raw_peptides
    path presentation_evidence
    path appm_summary
    path ccf_lite
    path normal_expression
    path normal_hla_ligands
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
    --profile '${profile_name}' \
    --out-events ranked_events.tsv \
    --out-peptides ranked_peptides.tsv

  echo "SCORE:" > versions.yml
  echo "  neoag: \$(python -c 'import neoag; print(neoag.__version__)')" >> versions.yml
  """
}
