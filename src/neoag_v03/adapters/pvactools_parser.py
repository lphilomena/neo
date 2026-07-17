from __future__ import annotations
from pathlib import Path
from ..model_layers import enrich_peptide_layers
from ..utils import read_tsv, write_tsv, first, safe_id, to_float
from ..schemas import EVENT_FIELDS, PEPTIDE_FIELDS
from ..driver_gene_db import lookup_driver_relevance
from ..tumor_specificity import compute_tumor_specificity

def discover_tsvs(paths: list[str | Path]) -> list[Path]:
    found = []
    for raw in paths:
        p = Path(raw)
        if not p.exists():
            raise FileNotFoundError(f"Missing pVACtools path: {p}")
        if p.is_file():
            found.append(p)
        else:
            for pat in ["*aggregated*.tsv", "*filtered*.tsv", "*all_epitopes*.tsv", "*.tsv"]:
                found.extend(p.rglob(pat))
    uniq, seen = [], set()
    for x in found:
        if x not in seen:
            uniq.append(x); seen.add(x)
    return uniq

def infer_tool(path: Path) -> str:
    s = str(path).lower()
    if "fuse" in s:
        return "pVACfuse"
    if "splice" in s:
        return "pVACsplice"
    if "bind" in s:
        return "pVACbind"
    return "pVACseq"

def infer_event_type(row: dict[str, str], tool: str) -> str:
    et = first(row, ["Event Type", "event_type", "Variant Type", "variant_type"], "")
    if et:
        return et
    if tool == "pVACfuse":
        return "Fusion"
    if tool == "pVACsplice":
        return "Splice"
    return "SNV"

def event_from_row(row: dict[str, str], sample_id: str, profile_name: str, tool: str, profile: dict | None = None) -> dict[str, str]:
    gene = first(row, ["Gene", "Gene Name", "Hugo_Symbol", "gene"], "UNKNOWN")
    event_type = infer_event_type(row, tool)
    event_name = first(row, ["Mutation", "Variant", "Protein Change", "HGVSp", "Fusion", "fusion_gene", "Event"], gene)
    event_id = first(row, ["event_id", "Index", "Mutation", "Variant"], "")
    if not event_id:
        event_id = safe_id(f"{sample_id}_{event_type}_{gene}_{event_name}")
    expression = to_float(first(row, ["Transcript Expression", "Gene Expression", "Expression", "TPM", "RNA Expr", "Allele Expr"], "0"), 0.0)
    return {
        "event_id": event_id,
        "sample_id": sample_id,
        "disease_profile": profile_name,
        "event_type": event_type,
        "gene": gene,
        "event_name": event_name,
        "chrom": first(row, ["Chromosome", "chrom", "CHROM"], ""),
        "pos": first(row, ["Start", "Position", "POS", "pos"], ""),
        "ref": first(row, ["Reference", "REF", "ref"], ""),
        "alt": first(row, ["Variant", "ALT", "alt"], ""),
        "transcript_id": first(row, ["Transcript", "Feature", "transcript_id"], ""),
        "consequence": first(row, ["Consequence", "consequence"], ""),
        "event_confidence": first(row, ["event_confidence", "Variant Confidence"], "0.7"),
        "event_expression": str(expression),
        "driver_relevance": first(row, ["driver_relevance", "Driver Relevance", "driver"], "")
            or str(lookup_driver_relevance(gene, profile)),
        "tumor_vaf": str(to_float(first(row, ["DNA VAF", "dna_vaf", "VAF", "tumor_vaf"], "0"), 0.0)),
        "tumor_depth": first(row, ["Tumor DNA Depth", "tumor_depth", "DP"], ""),
        "tumor_alt_count": first(row, ["Tumor DNA Alt Count", "tumor_alt_count", "AD_ALT"], ""),
        "rna_vaf": first(row, ["RNA VAF", "rna_vaf", "Tumor RNA VAF"], ""),
        "rna_alt_reads": first(row, ["RNA Alt Count", "Tumor RNA Alt Count", "rna_alt_reads", "rna_alt_count"], ""),
        "rna_depth": first(row, ["RNA Depth", "Tumor RNA Depth", "rna_depth"], ""),
        "clonality": first(row, ["clonality", "CCF", "ccf"], "0.5"),
        "persistence": first(row, ["persistence", "MRD Persistence", "relapse_retained"], "0.5"),
        "tumor_specificity": first(row, ["tumor_specificity", "Tumor Specificity"], "")
            or str(compute_tumor_specificity(gene, expression, profile)),
        "source": tool,
    }

def peptide_from_row(row: dict[str, str], sample_id: str, event: dict[str, str], tool: str) -> dict[str, str]:
    peptide = first(row, ["MT Epitope Seq", "MT Epitope", "Peptide", "peptide", "epitope", "Best Peptide"], "")
    wt = first(row, ["WT Epitope Seq", "WT Epitope", "wildtype_peptide"], "")
    hla = first(row, ["HLA Allele", "Allele", "hla_allele"], "")
    mhc = first(row, ["MHC Class", "mhc_class"], "")
    if not mhc:
        mhc = "II" if any(x in hla for x in ["DR", "DQ", "DP"]) else "I"
    pid = first(row, ["peptide_id", "Index"], "")
    if not pid:
        pid = safe_id(f"{event['event_id']}_{hla}_{peptide}")
    base = {
        "peptide_id": pid,
        "event_id": event["event_id"],
        "sample_id": sample_id,
        "event_type": event["event_type"],
        "mutation_source": event.get("mutation_source", ""),
        "peptide_consequence": event.get("peptide_consequence", ""),
        "gene": event["gene"],
        "peptide": peptide,
        "wildtype_peptide": wt,
        "crosses_junction": first(row, ["crosses_junction", "Crosses Junction"], ""),
        "contains_novel_aa": first(row, ["contains_novel_aa", "Contains Novel AA"], ""),
        "rna_junction_reads": first(row, ["rna_junction_reads", "RNA Junction Reads"], event.get("rna_junction_reads", "")),
        "hla_allele": hla,
        "mhc_class": mhc,
        "source_tool": tool,
        "binding_rank": str(to_float(first(row, ["Best MT Score", "Median MT Score", "MT %Rank", "%ile MT", "binding_rank"], "99"), 99.0)),
        "el_rank": str(to_float(first(row, ["EL Rank", "el_rank", "Best MT EL Score", "%ile MT"], "99"), 99.0)),
        "presentation_score": first(row, ["presentation_score", "Presentation Score"], "0.0"),
        "immunogenicity_score": first(row, ["immunogenicity_score", "Immunogenicity Score"], "0.5"),
        "wildtype_binding_rank": str(
            to_float(
                first(
                    row,
                    [
                        "netmhcpan_wt_rank_ba",
                        "Best WT Score",
                        "WT %Rank",
                        "%ile WT",
                        "wildtype_binding_rank",
                    ],
                    "99",
                ),
                99.0,
            )
        ),
        "self_similarity_score": first(row, ["self_similarity_score", "Self Similarity"], "0.0"),
        "normal_hla_ligand_overlap": first(row, ["normal_hla_ligand_overlap"], "no"),
    }
    for field in (
        "netmhcpan_mt_ic50",
        "netmhcpan_mt_rank_ba",
        "netmhcpan_mt_rank_el",
        "netmhcpan_wt_ic50",
        "netmhcpan_wt_rank_ba",
        "netmhcpan_wt_rank_el",
    ):
        value = first(row, [field], "")
        if value:
            base[field] = value
    return enrich_peptide_layers(base, event)

def parse_pvactools_outputs(paths, sample_id, profile_name, out_events=None, out_peptides=None, profile=None):
    events = {}
    peptides = []
    for tsv in discover_tsvs(paths):
        tool = infer_tool(tsv)
        for row in read_tsv(tsv):
            ev = event_from_row(row, sample_id, profile_name, tool, profile=profile)
            if ev["event_id"] not in events:
                events[ev["event_id"]] = ev
            peptides.append(peptide_from_row(row, sample_id, events[ev["event_id"]], tool))
    evs = list(events.values())
    if out_events:
        write_tsv(out_events, evs, EVENT_FIELDS)
    if out_peptides:
        write_tsv(out_peptides, peptides, PEPTIDE_FIELDS)
    return evs, peptides
