process SV_BUILD_RAW {
    tag "${sample_id}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    val sample_id
    path sv_list
    path reference_fasta
    path gencode_gtf

    output:
    path "parsed/raw_events.tsv", emit: raw_events
    path "parsed/raw_peptides.tsv", emit: raw_peptides
    path "sv/sv_events.full.tsv", emit: sv_events_full
    path "sv/sv_event_to_peptide.tsv", emit: sv_event_to_peptide
    path "sv/sv_mutant_proteins.fa", emit: sv_mutant_proteins

    script:
    def profile = params.profile_name ?: params.profile ?: 'sv_wgs_phase1'
    def hla = params.hla ?: "${projectDir}/../data/fixtures_sv/hla.txt"
    def expr = params.expression_tsv ? "--expression ${params.expression_tsv}" : ""
    def rna = params.rna_junction_tsv ? "--rna-junctions ${params.rna_junction_tsv}" : ""
    def nexpr = params.normal_expression ? "--normal-expression ${params.normal_expression}" : ""
    def nlig = params.normal_hla_ligands ? "--normal-hla-ligands ${params.normal_hla_ligands}" : ""
    def callers = params.sv_callers ? "--callers ${params.sv_callers.join(' ')}" : ""
    def wes = params.wes_mode ? 'sv-build-raw-wes' : 'sv-build-raw'
    def samples = params.tumor_sample_name && params.normal_sample_name ? "--tumor-sample-name '${params.tumor_sample_name}' --normal-sample-name '${params.normal_sample_name}'" : ""
    def build = params.genome_build ?: 'GRCh38'
    def products = params.expressed_products_tsv ? "--expressed-products '${params.expressed_products_tsv}'" : ""
    def capture = params.wes_mode ? "--capture-bed '${params.capture_bed}' --capture-near-bp ${params.capture_near_bp ?: 250} --capture-slop-bp ${params.capture_slop_bp ?: 1000}" : ""
    """
    if [[ "${sv_list}" == *.list ]]; then
      SV_VCFS=\$(grep -v '^#' ${sv_list} | tr '\n' ' ')
    else
      SV_VCFS=${sv_list}
    fi
    neoag ${wes} \
      --sample-id ${sample_id} \
      --profile ${profile} \
      --sv-vcf \${SV_VCFS} \
      ${callers} \
      --reference-fasta ${reference_fasta} \
      --gencode-gtf ${gencode_gtf} \
      --hla '${hla}' \
      --outdir . \
      --merge-distance-bp ${params.merge_distance_bp ?: 200} \
      --genome-build '${build}' \
      ${samples} ${products} ${capture} \
      ${expr} ${rna} ${nexpr} ${nlig}
    """
}
