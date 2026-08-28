process PEPTIDE_SAFETY {
  tag "peptide-safety"
  publishDir "${params.outdir}/safety", mode: 'copy'

  input:
    path raw_events
    path raw_peptides
    val profile_name
    path normal_expression
    path normal_hla_ligands

  output:
    path "peptide_safety.tsv", emit: peptide_safety
    path "event_safety.tsv", emit: event_safety
    path "versions.yml", emit: versions

  script:
  """
  neoag peptide-safety \
    --raw-events '${raw_events}' \
    --raw-peptides '${raw_peptides}' \
    --profile '${profile_name}' \
    --normal-expression '${normal_expression}' \
    --normal-hla-ligands '${normal_hla_ligands}' \
    --out peptide_safety.tsv \
    --event-out event_safety.tsv

  echo "PEPTIDE_SAFETY:" > versions.yml
  echo "  neoag: \$(python -c 'import neoag; print(neoag.__version__)')" >> versions.yml
  """
}
