/*
 * LOHHLA — HLA allele-specific loss-of-heterozygosity detection.
 *
 * Runs LOHHLA on paired tumor-normal BAM to detect HLA LOH events.
 * Outputs a standardized hla_loh.tsv for downstream immune escape analysis.
 *
 * Docker image: quay.io/biocontainers/lohhla:20171108--hdfd78af_3
 *   or via neoag-v03 run-tool lohhla when NEOAG_RUNNER_MODE=docker
 */
process LOHHLA {
  tag "$sample_id"
  label 'qc_check'
  time '8.h'
  publishDir "${params.outdir}/qc/lohhla", mode: 'copy'

  input:
    val sample_id
    path tumor_bam
    path normal_bam
    val hla_alleles           // list of HLA alleles from OptiType
    path ref_fasta

  output:
    path "hla_loh.tsv",         emit: hla_loh_tsv
    path "LOHHLA_output/**",    emit: raw_output,    optional: true
    path "versions.yml",        emit: versions

  script:
  def runner_mode = System.getenv('NEOAG_RUNNER_MODE') ?: 'conda'
  def hla_str = hla_alleles instanceof List ? hla_alleles.join(';') : hla_alleles.toString()
  """
  mkdir -p LOHHLA_output

  if [ "${runner_mode}" = "docker" ]; then
    # --- Docker mode: run LOHHLA via biocontainers image ---
    docker run --rm \\
      -v "\$PWD":"\$PWD" -w "\$PWD" \\
      -v ${ref_fasta.getParent()}:/ref:ro \\
      quay.io/biocontainers/lohhla:20171108--hdfd78af_3 \\
      Rscript /usr/local/bin/LOHHLAscript.R \\
        --patientId "${sample_id}" \\
        --outputDir "\$PWD/LOHHLA_output" \\
        --normalBAMfile "${normal_bam}" \\
        --tumorBAMfile "${tumor_bam}" \\
        --HLAfastaLoc /usr/local/bin/hla.dat \\
        --HLAexonLoc /usr/local/bin/hla.dat \\
        --CopyNumLoc /usr/local/bin/hlas \\\\
        2>&1 | tee LOHHLA_output/run.log || true
  else
    # --- Conda mode: LOHHLA via neoag-tools conda env ---
    export LOHHLA_HOME="\${LOHHLA_HOME:-\${NEOAG_TOOLS_ROOT}/tools/lohhla}"
    if [ -f "\${LOHHLA_HOME}/LOHHLAscript.R" ]; then
      Rscript "\${LOHHLA_HOME}/LOHHLAscript.R" \\
        --patientId "${sample_id}" \\
        --outputDir "\$PWD/LOHHLA_output" \\
        --normalBAMfile "${normal_bam}" \\
        --tumorBAMfile "${tumor_bam}" \\
        --HLAfastaLoc "\${LOHHLA_HOME}/hla.dat" \\
        --HLAexonLoc "\${LOHHLA_HOME}/hla.dat" \\
        --CopyNumLoc "\${LOHHLA_HOME}/hlas" \\\\
        2>&1 | tee LOHHLA_output/run.log || true
    else
      echo "WARNING: LOHHLA_HOME not set or LOHHHLAscript.R not found. Creating stub output."
    fi
  fi

  # --- Convert raw LOHHLA output to standardised hla_loh.tsv ---
  PRED_FILE=\$(find LOHHLA_output -name "*HLAlossPrediction_CI*" 2>/dev/null | head -1)
  if [ -n "\${PRED_FILE}" ] && [ -f "\${PRED_FILE}" ]; then
    neoag-v03 convert-lohhla -i "\${PRED_FILE}" -o hla_loh.tsv
  else
    # LOHHLA may not produce output if no HLA LOH detected or tool failed.
    # Emit an empty-but-valid hla_loh.tsv so downstream doesn't break.
    echo -e "hla_allele\\tloh_status" > hla_loh.tsv
    for a in ${hla_str//;/ }; do
      echo -e "\${a}\\tno_data" >> hla_loh.tsv
    done
  fi

  echo "LOHHLA:" > versions.yml
  echo "  runner: \${runner_mode}" >> versions.yml
  """
}
