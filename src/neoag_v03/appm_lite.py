from __future__ import annotations
from pathlib import Path
from typing import Any, Mapping

from .appm_v2 import build_appm_2, legacy_appm_rows_from_v2, MHC_I_GENES, MHC_II_GENES, IFNG_GENES, DAMAGING_TERMS

MHC_I = list(dict.fromkeys(MHC_I_GENES + IFNG_GENES))
MHC_II = MHC_II_GENES
DAMAGING = DAMAGING_TERMS
from .utils import read_tsv, write_tsv


def build_appm_lite(
    sample_id: str,
    vep_tsv,
    expression_tsv,
    hla_loh_tsv,
    profile: Mapping[str, Any],
    outdir,
    cnv_tsv=None,
    raw_peptides=None,
    hla_typing_tsv=None,
    tumor_purity_tsv=None,
    proteomics_tsv=None,
    phosphoproteomics_tsv=None,
    hla_ligandome_tsv=None,
):
    """Backward-compatible APPM-lite entry point backed by APPM 2.0.

    Existing callers still receive (legacy_rows, summary) and find
    appm_lite.tsv/appm_summary.tsv. New users also get APPM 2.0 sidecars:
    appm_gene_status.tsv, appm_pathway_status.tsv, peptide_appm_flags.tsv.
    """
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    paths = build_appm_2(
        sample_id=sample_id,
        vep_tsv=vep_tsv,
        expression_tsv=expression_tsv,
        hla_loh_tsv=hla_loh_tsv,
        cnv_tsv=cnv_tsv,
        raw_peptides=raw_peptides,
        hla_typing_tsv=hla_typing_tsv,
        tumor_purity_tsv=tumor_purity_tsv,
        proteomics_tsv=proteomics_tsv,
        phosphoproteomics_tsv=phosphoproteomics_tsv,
        hla_ligandome_tsv=hla_ligandome_tsv,
        profile=profile,
        outdir=out,
    )
    gene_rows = read_tsv(paths["appm_gene_status"])
    summary = read_tsv(paths["appm_summary"])[0]
    expression_assessed = summary.get("expression_assessment_status") == "assessed"
    legacy_rows = legacy_appm_rows_from_v2(gene_rows, expression_assessed)
    write_tsv(out / "appm_lite.tsv", legacy_rows)
    return legacy_rows, summary
