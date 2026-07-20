/*
 * PURITY_CHECK — Tumor purity / ploidy quality-control subworkflow.
 *
 * Runs FACETS and PURPLE in parallel, collecting tumor purity and ploidy
 * estimates.  Merges results into a single quality_control.tsv.
 *
 * Included by: workflows/main_all_qc.nf
 */

include { FACETS } from '../../modules/facets/main.nf'
include { PURPLE } from '../../modules/purple/main.nf'

workflow PURITY_CHECK {
  take:
    sample_id
    tumor_bam
    normal_bam
    ref_fasta
    dbsnp_vcf               // dbSNP/common SNP VCF (required for FACETS snp-pileup)

  main:
    // --- Run FACETS and PURPLE in parallel ---------------------------------
    FACETS(
      sample_id,
      tumor_bam,
      normal_bam,
      ref_fasta,
      dbsnp_vcf,
    )

    PURPLE(
      sample_id,
      tumor_bam,
      normal_bam,
      ref_fasta,
    )

    // --- Merge purity results into quality_control.tsv ---------------------
    quality_control = FACETS.out.purity_tsv
      .combine(PURPLE.out.purity_tsv, by: 0)
      .map { facets_tsv, purple_tsv ->
        def outdir = file("qc_purity_check")
        outdir.mkdirs()

        """
        mkdir -p ${outdir}

        python3 -c "
import csv
from pathlib import Path

out = Path('${outdir}')
qc_path = out / 'quality_control.tsv'

def read_purity(path_str):
    '''Read a purity.tsv file and return (purity, ploidy, method, status).'''
    p = Path(path_str)
    if not p.exists():
        return ('NA', 'NA', 'NA', 'missing')
    with p.open() as f:
        for row in csv.DictReader(f, delimiter='\\t'):
            purity = row.get('purity', row.get('Purity', 'NA'))
            ploidy = row.get('ploidy', row.get('Ploidy', 'NA'))
            method = row.get('method', 'NA')
            status = row.get('status', 'ok')
            return (purity, ploidy, method, status)
    return ('NA', 'NA', 'NA', 'parse_error')

facets_purity, facets_ploidy, facets_method, facets_status = read_purity('${facets_tsv}')
purple_purity, purple_ploidy, purple_method, purple_status = read_purity('${purple_tsv}')

with qc_path.open('w') as f:
    f.write('check_type\\tcheck_name\\tfacets_purity\\tfacets_ploidy\\tpurple_purity\\tpurple_ploidy\\tconsensus_purity\\tconsensus_ploidy\\n')

    # Calculate consensus if both tools returned numeric values
    cons_purity, cons_ploidy = 'NA', 'NA'
    try:
        fp, pp = float(facets_purity), float(purple_purity)
        fpl, ppl = float(facets_ploidy), float(purple_ploidy)
        cons_purity = f'{((fp + pp) / 2.0):.4f}'
        cons_ploidy = f'{((fpl + ppl) / 2.0):.4f}'
    except (ValueError, TypeError):
        pass

    f.write(f'purity\\tpurity_check\\t{facets_purity}\\t{facets_ploidy}\\t{purple_purity}\\t{purple_ploidy}\\t{cons_purity}\\t{cons_ploidy}\\n')

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
    facets_purity        = FACETS.out.purity_tsv
    purple_purity        = PURPLE.out.purity_tsv
    facets_raw           = FACETS.out.raw_output
    purple_raw           = PURPLE.out.raw_output
}
