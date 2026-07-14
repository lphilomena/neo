/*
 * OPTITYPE — Precision 4-digit HLA typing from NGS data.
 *
 * Wraps the `optitype run` command.  Accepts BAM (preferred) or paired FASTQ.
 * Outputs a plain-text hla_alleles.txt file (one normalised allele per line)
 * suitable for injection into the run config TOML.
 */
process OPTITYPE {
  tag "$sample_id"
  label 'large'
  container params.neoag_container ?: null
  time '12.h'
  publishDir "${params.outdir}/hla_typing", mode: 'copy'

  input:
    val sample_id
    path input_bam
    val seq_type               // 'dna' or 'rna'

  output:
    path "hla_alleles.txt", emit: hla_alleles
    path "optitype_result.csv", emit: typing_csv, optional: true
    path "versions.yml", emit: versions

  script:
  def stype = seq_type == 'rna' ? '--rna' : '--dna'
  // neoag-tools/bin/optitype uses #!/usr/bin/env python3.11 which may not
  // exist; neoag-optitype/bin/optitype uses an embedded correct Python path.
  def opti_bin = "${System.getenv('NEOAG_CONDA_BASE') ?: '/home/na/miniforge3'}/envs/neoag-optitype/bin/optitype"
  """
  # Index the input BAM if needed (shared mount BAMs may lack .bai).
  if [ ! -f '${input_bam}.bai' ] && [ ! -f '${input_bam}.csi' ]; then
    echo "Indexing input BAM for HLA typing..."
    samtools index '${input_bam}'
  fi

  # Extract chr6 HLA region and downsample to avoid OptiType uint16 overflow
  # (WGS BAMs can have >65535 reads per HLA position).
  echo "Extracting chr6 HLA region and downsampling..."
  samtools view -b '${input_bam}' chr6:29800000-33200000 | \
    samtools view -s 42.2 -b -o optitype_input.bam 2>/dev/null
  if [ ! -s optitype_input.bam ]; then
    echo "HLA region extraction failed, using full BAM with downsampling"
    samtools view -s 42.01 -b '${input_bam}' -o optitype_input.bam
  fi
  samtools index optitype_input.bam

  ${opti_bin} run \\
    --input optitype_input.bam \\
    ${stype} \\
    --outdir optitype_out \\
    --prefix '${sample_id}' \\
    --threads ${task.cpus}

  # Extract HLA alleles from OptiType result
  python3 -c "
import csv, sys
from pathlib import Path

results = list(Path('optitype_out').glob('*_result.tsv'))
if not results:
    sys.exit('OptiType produced no _result.tsv')

with open(results[0], newline='') as f:
    for row in csv.DictReader(f):
        alleles = []
        for col in ('A1','A2','B1','B2','C1','C2'):
            v = (row.get(col) or '').strip()
            if v and v not in alleles:
                alleles.append(v)
        Path('hla_alleles.txt').write_text('\\n'.join(alleles) + '\\n')
        Path('optitype_result.csv').write_text(open(results[0]).read())
        break
  "

  echo "OPTITYPE:" > versions.yml
  echo "  optitype: \$(optitype --version 2>&1 | head -1)" >> versions.yml
  """
}
