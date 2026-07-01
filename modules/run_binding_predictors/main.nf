process RUN_BINDING_PREDICTORS {
    tag "${sample_id}"
    publishDir "${params.outdir}/tools", mode: 'copy'

    input:
    val sample_id
    path raw_peptides
    val binding_stub

    output:
    path "netmhcpan.xls", emit: netmhcpan
    path "mhcflurry.csv", emit: mhcflurry
    path "versions.yml", emit: versions

    script:
    def stub = binding_stub ? '--stub' : ''
    def strict = params.strict_mode ? 'export NEOAG_STRICT_MODE=1' : ''
    """
    ${strict}
    mkdir -p tools
    neoag-v03 run-tool netmhcpan \
      --sample-id ${sample_id} \
      --raw-peptides '${raw_peptides}' \
      --output netmhcpan.xls \
      --workdir . \
      ${stub}

    neoag-v03 run-tool mhcflurry \
      --sample-id ${sample_id} \
      --raw-peptides '${raw_peptides}' \
      --output mhcflurry.csv \
      --workdir . \
      ${stub}

    echo "RUN_BINDING_PREDICTORS:" > versions.yml
    echo "  neoag-v03: \$(neoag-v03 run-demo --help >/dev/null 2>&1 && echo ok || echo unknown)" >> versions.yml
    """
}
