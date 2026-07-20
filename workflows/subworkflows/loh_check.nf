/*
 * LOH_CHECK — HLA Loss-of-Heterozygosity quality-control subworkflow.
 *
 * Runs LOHHLA and SpecHLA in parallel, collecting their per-allele LOH calls.
 * Merges results into a single quality_control.tsv with cross-tool consensus.
 *
 * Included by: workflows/main_all_qc.nf
 */

include { LOHHLA } from '../../modules/lohhla/main.nf'
include { SPECHLA } from '../../modules/spechla/main.nf'

workflow LOH_CHECK {
  take:
    sample_id
    tumor_bam
    normal_bam
    hla_alleles             // list of HLA alleles (from OptiType)
    ref_fasta

  main:
    // --- Run LOHHLA and SpecHLA in parallel -------------------------------
    LOHHLA(
      sample_id,
      tumor_bam,
      normal_bam,
      hla_alleles,
      ref_fasta,
    )

    SPECHLA(
      sample_id,
      tumor_bam,
      normal_bam,
      ref_fasta,
    )

    // --- Merge LOH results into quality_control.tsv -----------------------
    // Collect both TSVs and cross-check with the neoag-v03 CLI.
    // If either tool fails to produce output, the module emits a stub TSV
    // with status "no_data", so this merge always succeeds.
    quality_control = LOHHLA.out.hla_loh_tsv
      .combine(SPECHLA.out.hla_loh_tsv, by: 0)
      .map { lohhla_tsv, spechla_tsv ->
        def outdir = file("qc_loh_check")
        outdir.mkdirs()

        """
        mkdir -p ${outdir}

        # Cross-check LOHHLA and SpecHLA results
        neoag-v03 crosscheck-hla-loh \\
          --lohhla-hla-loh ${lohhla_tsv} \\
          --spechla-hla-loh ${spechla_tsv} \\
          --out ${outdir}/hla_loh_crosscheck.tsv \\
          --consensus-out ${outdir}/hla_loh_consensus.tsv \\\\
          2>&1 || echo "crosscheck-hla-loh failed, using individual results"

        # Build quality_control.tsv
        python3 -c "
import csv
from pathlib import Path

out = Path('${outdir}')
qc_path = out / 'quality_control.tsv'

# Collect LOH status per allele from both tools
loh_data = {}  # allele -> {lohhla: status, spechla: status, consensus: status}

# Read LOHHLA
lh = out / '..' / '..'
lh_file = Path('${lohhla_tsv}')
if lh_file.exists():
    with lh_file.open() as f:
        for row in csv.DictReader(f, delimiter='\\t'):
            allele = row.get('hla_allele', '')
            status = row.get('loh_status', '')
            if allele:
                loh_data.setdefault(allele, {})['lohhla'] = status

# Read SpecHLA
sh_file = Path('${spechla_tsv}')
if sh_file.exists():
    with sh_file.open() as f:
        for row in csv.DictReader(f, delimiter='\\t'):
            allele = row.get('hla_allele', '')
            status = row.get('loh_status', '')
            if allele:
                loh_data.setdefault(allele, {})['spechla'] = status

# Read consensus if available
cc = out / 'hla_loh_consensus.tsv'
if cc.exists():
    with cc.open() as f:
        for row in csv.DictReader(f, delimiter='\\t'):
            allele = row.get('hla_allele', '')
            status = row.get('loh_status', '')
            if allele:
                loh_data.setdefault(allele, {})['consensus'] = status

# Write quality_control.tsv
with qc_path.open('w') as f:
    f.write('check_type\\tcheck_name\\tallele\\tlohhla_status\\tspechla_status\\tconsensus_status\\n')
    if not loh_data:
        f.write('loh\\tno_data\\t.\\tno_data\\tno_data\\tno_data\\n')
    for allele in sorted(loh_data):
        d = loh_data[allele]
        f.write(f'loh\\thla_loh\\t{allele}\\t{d.get(\"lohhla\", \"NA\")}\\t{d.get(\"spechla\", \"NA\")}\\t{d.get(\"consensus\", \"NA\")}\\n')

print(f'Wrote {qc_path}')
"
        """
      }
      .collectFile(
        name: 'quality_control.tsv',
        newLine: true,
        sort: true,
        storeDir: "${params.outdir}/qc"
      )

  emit:
    quality_control_tsv = quality_control
    lohhla_hla_loh       = LOHHLA.out.hla_loh_tsv
    spechla_hla_loh      = SPECHLA.out.hla_loh_tsv
    lohhla_raw           = LOHHLA.out.raw_output
    spechla_raw          = SPECHLA.out.raw_output
}
