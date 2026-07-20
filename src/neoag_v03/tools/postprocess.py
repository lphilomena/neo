from __future__ import annotations
from pathlib import Path
from ..utils import read_tsv, write_tsv, first
from ..adapters.facets import parse_facets_purity, parse_facets_cncf, write_purity_evidence, write_cnv_segments
from ..adapters.lohhla import parse_lohhla_prediction, write_hla_loh_evidence
from ..adapters.spechla import parse_spechla_loh_merge, write_spechla_hla_loh_evidence
from ..adapters.optitype import parse_optitype_result
from ..appm_lite import MHC_I, MHC_II

APPM_GENES = set(MHC_I + MHC_II)


def vep_to_appm_tsv(vep_path: str | Path, out_path: str | Path) -> None:
    """Normalize VEP tab output to gene\\tconsequence for APPM-lite."""
    genes: dict[str, list[str]] = {}
    for row in read_tsv(vep_path):
        gene = first(row, ["SYMBOL", "Gene", "gene", "Uploaded_variation"], "")
        if gene in APPM_GENES or first(row, ["Gene", "gene"], "") in APPM_GENES:
            g = gene if gene in APPM_GENES else first(row, ["Gene", "gene"], gene)
            cons = first(row, ["Consequence", "consequence", "CSQ"], "")
            if g:
                genes.setdefault(g, [])
                if cons and cons not in genes[g]:
                    genes[g].append(cons.split(",")[0])
    rows = [{"gene": g, "consequence": cons[0] if cons else ""} for g, cons in sorted(genes.items())]
    write_tsv(out_path, rows, ["gene", "consequence"])


def lohhla_to_hla_loh_tsv(lohhla_path: str | Path, out_path: str | Path) -> None:
    rows = parse_lohhla_prediction(lohhla_path)
    if not rows:
        for row in read_tsv(lohhla_path):
            allele = first(row, ["hla_allele", "allele", "HLA", "gene"], "")
            status = first(row, ["loh_status", "LOH", "status", "LOH_status"], "no")
            if allele:
                rows.append({"hla_allele": allele, "loh_status": status})
    if not rows:
        rows = [{"hla_allele": "HLA-A*02:01", "loh_status": "no"}]
    write_hla_loh_evidence(out_path, rows)


def spechla_to_hla_loh_tsv(spechla_path: str | Path, out_path: str | Path) -> None:
    rows = parse_spechla_loh_merge(spechla_path)
    if not rows:
        raise ValueError(f"No HLA LOH rows parsed from SpecHLA output: {spechla_path}")
    write_spechla_hla_loh_evidence(out_path, rows, source_path=spechla_path)


def facets_to_purity_tsv(facets_path: str | Path, sample_id: str, out_path: str | Path) -> None:
    write_purity_evidence(out_path, parse_facets_purity(facets_path, sample_id))


def facets_to_cnv_tsv(cncf_path: str | Path, out_path: str | Path) -> None:
    write_cnv_segments(out_path, parse_facets_cncf(cncf_path))


def optitype_to_hla_alleles(result_csv: str | Path) -> list[str]:
    """Parse OptiType result CSV and return normalised HLA alleles list."""
    return parse_optitype_result(result_csv)
