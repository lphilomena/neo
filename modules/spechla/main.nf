/*
 * SpecHLA — HLA typing with allele-specific copy number.
 *
 * Runs SpecHLA on paired tumor-normal BAM to obtain HLA typing and
 * allele-level copy number, then extracts LOH status from merge.hla.copy.txt.
 *
 * Docker image: neoag-spechla:ubuntu22.04 (自建镜像)
 */
process SPECHLA {
  tag "$sample_id"
  label 'qc_check'
  time '12.h'
  publishDir "${params.outdir}/qc/spechla", mode: 'copy'

  input:
    val sample_id
    path tumor_bam
    path normal_bam
    path ref_fasta

  output:
    path "hla_loh.tsv",           emit: hla_loh_tsv
    path "spechla_output/**",    emit: raw_output,    optional: true
    path "versions.yml",         emit: versions

  script:
  def runner_mode = System.getenv('NEOAG_RUNNER_MODE') ?: 'conda'
  """
  mkdir -p spechla_output

  if [ "${runner_mode}" = "docker" ]; then
    # --- Docker mode: run SpecHLA via self-built image ---
    docker run --rm \\
      -v "\$PWD":"\$PWD" -w "\$PWD" \\
      -v ${ref_fasta.getParent()}:/ref:ro \\
      neoag-spechla:ubuntu22.04 \\
      python3 /opt/spechla/main.py \\
        --normal "${normal_bam}" \\
        --tumor "${tumor_bam}" \\
        --reference "${ref_fasta}" \\
        --sample_id "${sample_id}" \\
        --outdir "\$PWD/spechla_output" \\\\
        2>&1 | tee spechla_output/run.log || true
  else
    # --- Conda mode: SpecHLA from conda env ---
    SPECHLA_BIN="\${NEOAG_CONDA_BASE:-\${HOME}/miniconda3}/envs/neoag-spechla/bin/spechla"
    if command -v spechla >/dev/null 2>&1; then
      spechla \\
        --normal "${normal_bam}" \\
        --tumor "${tumor_bam}" \\
        --reference "${ref_fasta}" \\
        --sample_id "${sample_id}" \\
        --outdir "\$PWD/spechla_output" \\\\
        2>&1 | tee spechla_output/run.log || true
    else
      echo "WARNING: specHLA not found in PATH. Creating stub output."
    fi
  fi

  # --- Convert SpecHLA merge.hla.copy.txt to standardised hla_loh.tsv ---
  MERGE_FILE=\$(find spechla_output -name "merge.hla.copy.txt" 2>/dev/null | head -1)
  if [ -n "\${MERGE_FILE}" ] && [ -f "\${MERGE_FILE}" ]; then
    neoag-v03 convert-spechla -i "\${MERGE_FILE}" -o hla_loh.tsv
  else
    # Emit empty-but-valid TSV if SpecHLA produced no output.
    echo -e "hla_allele\\tloh_status" > hla_loh.tsv
    echo -e "HLA-A*00:00\\tno_data" >> hla_loh.tsv
  fi

  echo "SPECHLA:" > versions.yml
  echo "  runner: \${runner_mode}" >> versions.yml
  """
}
