/*
 * PURPLE — HMF PURPLE tumor purity / ploidy / copy number estimator.
 *
 * Runs PURPLE on paired tumor-normal BAM to produce purity/ploidy and
 * somatic/germline CNV calls.  Typically used alongside AMBER + COBALT
 * as part of the HMF pipeline, but can also accept pre-computed inputs.
 *
 * Docker image: neoag-purple-suite:ubuntu22.04 (自建镜像, shared with FACETS)
 */
process PURPLE {
  tag "$sample_id"
  label 'qc_check'
  time '12.h'
  publishDir "${params.outdir}/qc/purple", mode: 'copy'

  input:
    val sample_id
    path tumor_bam
    path normal_bam
    path ref_fasta

  output:
    path "purity.tsv",           emit: purity_tsv
    path "purple_output/**",    emit: raw_output,   optional: true
    path "versions.yml",        emit: versions

  script:
  def runner_mode = System.getenv('NEOAG_RUNNER_MODE') ?: 'conda'
  """
  mkdir -p purple_output

  if [ "${runner_mode}" = "docker" ]; then
    # --- Docker mode: run PURPLE via self-built image ---
    docker run --rm \\
      -v "\$PWD":"\$PWD" -w "\$PWD" \\
      -v ${ref_fasta.getParent()}:/ref:ro \\
      neoag-purple-suite:ubuntu22.04 \\
      bash -c "
        # PURPLE expects tumor + normal BAM + reference + optional AMBER/COBALT dirs.
        # If AMBER/COBALT outputs are not available, PURPLE runs in reference-only
        # mode which provides purity/ploidy from read-depth + BAF alone.
        java -Xmx16G -jar /opt/hmftools/purple.jar \\
          -reference '\$(basename ${ref_fasta})' \\
          -ref_genome /ref/\$(basename ${ref_fasta}) \\
          -tumor ${tumor_bam} \\
          -reference ${normal_bam} \\
          -output_dir \$PWD/purple_output \\
          -threads ${task.cpus} \\\\
          2>&1 | tee purple_output/run.log
      " || true
  else
    # --- Conda mode: PURPLE from installed tools ---
    PURPLE_JAR="\${PURPLE_HOME:-/opt/hmftools}/purple.jar"
    if [ -f "\${PURPLE_JAR}" ]; then
      java -Xmx16G -jar "\${PURPLE_JAR}" \\
        -reference "\$(basename ${ref_fasta})" \\
        -ref_genome "${ref_fasta}" \\
        -tumor "${tumor_bam}" \\
        -reference "${normal_bam}" \\
        -output_dir "\$PWD/purple_output" \\
        -threads ${task.cpus} \\\\
        2>&1 | tee purple_output/run.log || true
    else
      echo "WARNING: PURPLE jar not found at \${PURPLE_JAR}. Creating stub output."
    fi
  fi

  # --- Parse PURPLE purity output to standardised purity.tsv ---
  PURPLE_TSV="\$(find purple_output -name '*.purple.purity.tsv' -o -name '*.purity' 2>/dev/null | head -1)"
  if [ -n "\${PURPLE_TSV}" ] && [ -f "\${PURPLE_TSV}" ]; then
    python3 -c "
import csv, sys
from pathlib import Path
purple_tsv = Path('\${PURPLE_TSV}')
reader = csv.DictReader(purple_tsv.open(), delimiter='\\t') if purple_tsv.suffix == '.tsv' else None
rows = list(reader) if reader else []
purity, ploidy = 'NA', 'NA'
if rows and len(rows) >= 1:
    r = rows[0]
    purity = r.get('purity', r.get('Purity', 'NA'))
    ploidy = r.get('ploidy', r.get('Ploidy', 'NA'))
else:
    for line in purple_tsv.read_text().split('\\n'):
        parts = line.strip().split('\\t')
        if len(parts) >= 2:
            try:
                float(parts[0])
                purity = parts[0]
                ploidy = parts[1]
                break
            except ValueError:
                continue
Path('purity.tsv').write_text(f'sample_id\\tpurity\\tploidy\\tmethod\\tstatus\\n${sample_id}\\t{purity}\\t{ploidy}\\tpurple\\tok\\n')
" || echo -e "sample_id\\tpurity\\tploidy\\tmethod\\tstatus\\n${sample_id}\\tNA\\tNA\\tpurple\\tparse_error" > purity.tsv
  else
    # Emit empty-but-valid TSV if PURPLE produced no output.
    echo -e "sample_id\\tpurity\\tploidy\\tmethod\\tstatus" > purity.tsv
    echo -e "${sample_id}\\tNA\\tNA\\tpurple\\tno_data" >> purity.tsv
  fi

  echo "PURPLE:" > versions.yml
  echo "  runner: \${runner_mode}" >> versions.yml
  """
}
