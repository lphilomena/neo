process IMMUNE_ESCAPE {
  tag "$sample_id"
  publishDir "${params.outdir}/immune_escape", mode: 'copy'

  input:
    val sample_id
    path raw_peptides
    val profile_name
    path vep_appm_file
    path cnv_file
    path expression_file
    path hla_loh_file
    path appm_gene_status
    path appm_pathway_status
    path ccf_lite

  output:
    path "immune_escape_summary.tsv", emit: immune_escape_summary
    path "peptide_escape_flags.tsv", emit: peptide_escape_flags
    path "versions.yml", emit: versions

  script:
  """
  neoag-v03 immune-escape \
    --sample-id '${sample_id}' \
    --raw-peptides '${raw_peptides}' \
    --profile '${profile_name}' \
    --vep-tsv '${vep_appm_file}' \
    --cnv '${cnv_file}' \
    --expression '${expression_file}' \
    --hla-loh '${hla_loh_file}' \
    --appm-gene-status '${appm_gene_status}' \
    --appm-pathway-status '${appm_pathway_status}' \
    --ccf '${ccf_lite}' \
    --outdir .

  echo "IMMUNE_ESCAPE:" > versions.yml
  echo "  neoag-v03: \$(python -c 'import neoag_v03; print(neoag_v03.__version__)')" >> versions.yml
  """
}
