process SCORE_V03 {
  tag "score-v03"
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
    path "ranked_events.v03.tsv", emit: ranked_events
    path "ranked_peptides.v03.tsv", emit: ranked_peptides
    path "versions.yml", emit: versions

  script:
  """
  neoag-v03 score-v03 \
    --raw-events '${raw_events}' \
    --raw-peptides '${raw_peptides}' \
    --presentation '${presentation_evidence}' \
    --appm-summary '${appm_summary}' \
    --ccf '${ccf_lite}' \
    --normal-expression '${normal_expression}' \
    --normal-hla-ligands '${normal_hla_ligands}' \
    --profile '${profile_name}' \
    --out-events ranked_events.v03.tsv \
    --out-peptides ranked_peptides.v03.tsv

  echo "SCORE_V03:" > versions.yml
  echo "  neoag-v03: \$(python -c 'import neoag_v03; print(neoag_v03.__version__)')" >> versions.yml
  """
}
