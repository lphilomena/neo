/*
 * FACETS — Tumor purity, ploidy and allele-specific copy number.
 *
 * Runs snp-pileup + runFACETS.R on paired tumor-normal BAM to estimate
 * tumor purity and clonal/subclonal copy-number segments.
 *
 * Docker image: neoag-purple-suite:ubuntu22.04 (自建镜像)
 */
process FACETS {
  tag "$sample_id"
  label 'qc_check'
  time '8.h'
  publishDir "${params.outdir}/qc/facets", mode: 'copy'

  input:
    val sample_id
    path tumor_bam
    path normal_bam
    path ref_fasta
    path dbsnp_vcf            // dbSNP/common SNP VCF for snp-pileup

  output:
    path "purity.tsv",           emit: purity_tsv
    path "facets_output/**",    emit: raw_output,   optional: true
    path "versions.yml",        emit: versions

  script:
  def runner_mode = System.getenv('NEOAG_RUNNER_MODE') ?: 'conda'
  def pileup_bin = System.getenv('NEOAG_SNP_PILEUP_BIN') ?: 'snp-pileup'
  """
  mkdir -p facets_output

  if [ "${runner_mode}" = "docker" ]; then
    # --- Docker mode: run FACETS via self-built image ---
    docker run --rm \\
      -v "\$PWD":"\$PWD" -w "\$PWD" \\
      -v ${ref_fasta.getParent()}:/ref:ro \\
      neoag-purple-suite:ubuntu22.04 \\
      bash -c "
        mkdir -p /work/facets_output
        cd /work
        # Step 1: snp-pileup
        snp-pileup \\
          --germline=false \\
          --pseudo-snps=100 \\
          --min-map-quality=30 \\
          --min-base-quality=20 \\
          --max-depth=10000 \\
          ${dbsnp_vcf} \\
          facets_output/${sample_id}.pileup.gz \\
          ${normal_bam} \\
          ${tumor_bam}

        # Step 2: runFACETS.R
        gunzip -c facets_output/${sample_id}.pileup.gz > facets_output/${sample_id}.pileup.txt
        runFACETS.R \\
          facets_output/${sample_id}.pileup.txt \\
          facets_output/${sample_id}_purity.txt \\
          facets_output/${sample_id}_cncf.tsv
      " || true
  else
    # --- Conda mode: FACETS from installed tools ---
    if command -v snp-pileup >/dev/null 2>&1; then
      snp-pileup \\
        --germline=false \\
        --pseudo-snps=100 \\
        --min-map-quality=30 \\
        --min-base-quality=20 \\
        --max-depth=10000 \\
        ${dbsnp_vcf} \\
        "facets_output/${sample_id}.pileup.gz" \\
        "${normal_bam}" \\
        "${tumor_bam}" || true

      if [ -f "facets_output/${sample_id}.pileup.gz" ]; then
        gunzip -c "facets_output/${sample_id}.pileup.gz" > "facets_output/${sample_id}.pileup.txt"
      fi
    fi
    if [ -f "facets_output/${sample_id}.pileup.txt" ]; then
      runFACETS.R \\
        "facets_output/${sample_id}.pileup.txt" \\
        "facets_output/${sample_id}_purity.txt" \\
        "facets_output/${sample_id}_cncf.tsv" || true
    fi
  fi

  # --- Convert FACETS purity output to standardised purity.tsv ---
  PURITY_FILE="\$(find facets_output -name '*_purity.txt' 2>/dev/null | head -1)"
  if [ -n "\${PURITY_FILE}" ] && [ -f "\${PURITY_FILE}" ]; then
    CNCF_FILE="\$(find facets_output -name '*_cncf.tsv' 2>/dev/null | head -1)"
    CNV_ARG=""
    if [ -n "\${CNCF_FILE}" ] && [ -f "\${CNCF_FILE}" ]; then
      CNV_ARG="--cnv-input \${CNCF_FILE} --cnv-output facets_output/cnv_segments.tsv"
    fi
    neoag-v03 convert-facets \\
      --purity-input "\${PURITY_FILE}" \\
      --purity-output purity.tsv \\
      --sample-id "${sample_id}" \\
      \${CNV_ARG} || true
  else
    # Emit empty-but-valid TSV if FACETS produced no output.
    echo -e "sample_id\\tpurity\\tploidy\\tmethod\\tstatus" > purity.tsv
    echo -e "${sample_id}\\tNA\\tNA\\tfacets\\tno_data" >> purity.tsv
  fi

  echo "FACETS:" > versions.yml
  echo "  runner: \${runner_mode}" >> versions.yml
  """
}
