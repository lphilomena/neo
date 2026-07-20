process APPM_2 {
  tag "$sample_id"
  publishDir "${params.outdir}/appm", mode: 'copy'

  input:
    val sample_id
    val profile_name
    path vep_appm_file
    path expression_file
    path hla_loh_file
    path cnv_file
    path raw_peptides
    path purity_file

  output:
    path "appm_gene_status.tsv", emit: appm_gene_status
    path "appm_pathway_status.tsv", emit: appm_pathway_status
    path "appm_module_scores.tsv", emit: appm_module_scores
    path "appm_submodule_scores.tsv", emit: appm_submodule_scores
    path "appm_immune_context.tsv", emit: appm_immune_context
    path "appm_evidence_completeness.tsv", emit: appm_evidence_completeness
    path "appm_input_status.tsv", emit: appm_input_status
    path "appm_conflicts.tsv", emit: appm_conflicts
    path "peptide_appm_flags.tsv", emit: peptide_appm_flags
    path "appm_peptide_modifiers.tsv", emit: appm_peptide_modifiers
    path "appm_summary.tsv", emit: appm_summary
    path "versions.yml", emit: versions

  script:
  """
  neoag appm-2 \
    --sample-id '${sample_id}' \
    --profile '${profile_name}' \
    --vep-tsv '${vep_appm_file}' \
    --expression '${expression_file}' \
    --hla-loh '${hla_loh_file}' \
    --cnv '${cnv_file}' \
    --raw-peptides '${raw_peptides}' \
    --tumor-purity '${purity_file}' \
    --outdir .

  echo "APPM_2:" > versions.yml
  echo "  neoag: \$(python -c 'import neoag; print(neoag.__version__)')" >> versions.yml
  """
}
