from __future__ import annotations

from pathlib import Path

from ..model_layers import enrich_event_layers, enrich_peptide_layers
from ..provenance import merge_rows_preserving_provenance
from ..schemas import EVENT_FIELDS, PEPTIDE_FIELDS
from ..utils import first, read_tsv, safe_id, to_float, write_tsv


def discover_tsvs(paths: list[str | Path]) -> list[Path]:
    found: list[Path] = []
    for raw in paths:
        path = Path(raw)
        if not path.exists():
            raise FileNotFoundError(f"Missing pVACtools path: {path}")
        if path.is_file():
            found.append(path)
        else:
            for pattern in (
                "*aggregated*.tsv",
                "*filtered*.tsv",
                "*all_epitopes*.tsv",
                "*.tsv",
            ):
                found.extend(path.rglob(pattern))
    unique: list[Path] = []
    seen: set[Path] = set()
    for path in found:
        if path not in seen:
            unique.append(path)
            seen.add(path)
    return unique


def infer_tool(path: Path) -> str:
    value = str(path).lower()
    if "fuse" in value:
        return "pVACfuse"
    if "splice" in value:
        return "pVACsplice"
    if "bind" in value:
        return "pVACbind"
    return "pVACseq"


def infer_event_type(row: dict[str, str], tool: str) -> str:
    event_type = first(
        row, ["Event Type", "event_type", "Variant Type", "variant_type"], ""
    )
    if event_type:
        return event_type
    if tool == "pVACfuse":
        return "Fusion"
    if tool == "pVACsplice":
        return "Splice"
    variant_text = " ".join(
        first(row, [key], "") for key in ("Index", "AA Change", "Consequence")
    ).lower()
    if any(
        token in variant_text
        for token in ("frameshift", ".fs.", "inframe_del", "inframe_ins")
    ):
        return "InDel"
    return "SNV"


def event_from_row(
    row: dict[str, str],
    sample_id: str,
    profile_name: str,
    tool: str,
    *,
    source_file: str = "",
    source_row_number: int = 0,
) -> dict[str, str]:
    gene = first(row, ["Gene", "Gene Name", "Hugo_Symbol", "gene"], "UNKNOWN")
    event_type = infer_event_type(row, tool)
    if tool == "pVACsplice":
        event_name = first(
            row,
            [
                "Splice Junction",
                "splice_junction",
                "junction_id",
                "Junction Name",
                "Junction",
                "event_id",
                "Event",
            ],
            gene,
        )
        event_id = first(
            row,
            ["event_id", "Splice Junction", "splice_junction", "junction_id", "Index"],
            "",
        )
    else:
        event_name = first(
            row,
            [
                "Mutation",
                "Variant",
                "Protein Change",
                "HGVSp",
                "Fusion",
                "fusion_gene",
                "Event",
            ],
            gene,
        )
        event_id = first(row, ["event_id", "Index", "Mutation", "Variant"], "")
    if not event_id:
        event_id = safe_id(f"{sample_id}_{event_type}_{gene}_{event_name}")

    source_record_id = (
        first(row, ["source_record_id", "Index", "ID", "id"], "")
        or f"{tool}:{source_row_number or event_id}"
    )
    provided_reads = first(
        row, ["RNA Junction Reads", "rna_junction_reads", "Read Support"], ""
    )
    canonical = event_name if str(event_name).startswith("SJ|") else ""
    base = {
        "event_id": event_id,
        "sample_id": sample_id,
        "disease_profile": profile_name,
        "event_type": event_type,
        "gene": gene,
        "event_name": event_name,
        "canonical_junction_id": canonical if tool == "pVACsplice" else "",
        "source_junction_id": event_name if tool == "pVACsplice" else "",
        "junction_resolution_status": (
            "RESOLVED_FROM_CANONICAL_ID" if canonical else "AWAITING_JUNCTION_REGISTRY"
        )
        if tool == "pVACsplice"
        else "",
        "junction_resolution_reason": (
            "pVACsplice source junction requires exact registry resolution"
            if tool == "pVACsplice" and not canonical
            else ""
        ),
        "junction_support_status": (
            "PROVIDED_UNVERIFIED" if tool == "pVACsplice" and provided_reads else ""
        ),
        "provided_rna_junction_reads": (
            provided_reads if tool in {"pVACsplice", "pVACfuse"} else ""
        ),
        # Unverified pVACsplice counts do not enter scoring until the exact
        # RegTools junction adapter resolves the source junction.
        "rna_junction_reads": (
            "0"
            if tool == "pVACsplice"
            else provided_reads
            if tool == "pVACfuse"
            else ""
        ),
        "source_file": source_file,
        "source_row_number": str(source_row_number or ""),
        "source_record_id": source_record_id,
        "source_tools": tool,
        "source_records": source_record_id,
        "provenance_record_count": "1",
        "evidence_conflict_status": "NONE",
        "chrom": first(row, ["Chromosome", "chrom", "CHROM"], ""),
        "pos": first(row, ["Start", "Position", "POS", "pos"], ""),
        "ref": first(row, ["Reference", "REF", "ref"], ""),
        "alt": first(row, ["Variant", "ALT", "alt"], ""),
        "transcript_id": first(
            row, ["Transcript", "Best Transcript", "Feature", "transcript_id"], ""
        ),
        "consequence": first(row, ["Consequence", "consequence"], ""),
        "event_confidence": first(row, ["event_confidence", "Variant Confidence"], "0.7"),
        "event_expression": str(
            to_float(
                first(
                    row,
                    [
                        "Transcript Expression",
                        "Gene Expression",
                        "Expression",
                        "TPM",
                        "RNA Expr",
                        "Allele Expr",
                        "Expr",
                    ],
                    "0",
                ),
                0.0,
            )
        ),
        "driver_relevance": first(row, ["driver_relevance", "Driver Relevance", "driver"], "0.0"),
        "tumor_vaf": str(
            to_float(first(row, ["DNA VAF", "dna_vaf", "VAF", "tumor_vaf"], "0"), 0.0)
        ),
        "tumor_depth": first(row, ["Tumor DNA Depth", "tumor_depth", "DP"], ""),
        "tumor_alt_count": first(
            row, ["Tumor DNA Alt Count", "tumor_alt_count", "AD_ALT"], ""
        ),
        "rna_vaf": first(row, ["RNA VAF", "rna_vaf", "Tumor RNA VAF"], ""),
        "rna_alt_reads": first(
            row,
            ["RNA Alt Count", "Tumor RNA Alt Count", "rna_alt_reads", "rna_alt_count"],
            "",
        ),
        "rna_depth": first(row, ["RNA Depth", "Tumor RNA Depth", "rna_depth"], ""),
        "clonality": first(row, ["clonality", "CCF", "ccf"], "0.5"),
        "persistence": first(row, ["persistence", "MRD Persistence", "relapse_retained"], "0.5"),
        "tumor_specificity": first(
            row, ["tumor_specificity", "Tumor Specificity"], "0.7"
        ),
        "source": tool,
    }
    return enrich_event_layers(base)


def peptide_from_row(
    row: dict[str, str],
    sample_id: str,
    event: dict[str, str],
    tool: str,
    *,
    source_file: str = "",
    source_row_number: int = 0,
) -> dict[str, str]:
    peptide = first(
        row,
        ["MT Epitope Seq", "MT Epitope", "Peptide", "peptide", "epitope", "Best Peptide"],
        "",
    )
    wildtype = first(row, ["WT Epitope Seq", "WT Epitope", "wildtype_peptide"], "")
    hla = first(row, ["HLA Allele", "Allele", "hla_allele"], "")
    mhc_class = first(row, ["MHC Class", "mhc_class"], "")
    if not mhc_class:
        mhc_class = "II" if any(token in hla for token in ("DR", "DQ", "DP")) else "I"
    peptide_id = first(row, ["peptide_id", "Index"], "")
    if not peptide_id:
        peptide_id = safe_id(f"{event['event_id']}_{hla}_{peptide}")
    source_record_id = (
        first(row, ["source_record_id", "Index", "ID", "id"], "")
        or f"{tool}:{source_row_number or peptide_id}"
    )
    provided_reads = first(
        row,
        ["rna_junction_reads", "RNA Junction Reads", "Read Support"],
        event.get("provided_rna_junction_reads", ""),
    )

    base = {
        "peptide_id": peptide_id,
        "event_id": event["event_id"],
        "sample_id": sample_id,
        "event_type": event["event_type"],
        "mutation_source": event.get("mutation_source", ""),
        "peptide_consequence": event.get("peptide_consequence", ""),
        "gene": event["gene"],
        "peptide": peptide,
        "wildtype_peptide": wildtype,
        "canonical_junction_id": event.get("canonical_junction_id", ""),
        "source_junction_id": event.get("source_junction_id", ""),
        "crosses_junction": first(row, ["crosses_junction", "Crosses Junction"], ""),
        "contains_novel_aa": first(row, ["contains_novel_aa", "Contains Novel AA"], ""),
        "provided_rna_junction_reads": provided_reads if tool == "pVACsplice" else "",
        "rna_junction_reads": "0" if tool == "pVACsplice" else provided_reads,
        "junction_support_status": (
            "PROVIDED_UNVERIFIED" if tool == "pVACsplice" and provided_reads else ""
        ),
        "hla_allele": hla,
        "mhc_class": mhc_class,
        "source_tool": tool,
        "source_file": source_file,
        "source_row_number": str(source_row_number or ""),
        "source_record_id": source_record_id,
        "source_tools": tool,
        "source_records": source_record_id,
        "provenance_record_count": "1",
        "evidence_conflict_status": "NONE",
        "generation_status": "provided_by_pvactools",
        "binding_rank": str(
            to_float(
                first(
                    row,
                    ["Best MT Score", "Median MT Score", "MT %Rank", "%ile MT", "binding_rank"],
                    "99",
                ),
                99.0,
            )
        ),
        "el_rank": str(
            to_float(
                first(row, ["EL Rank", "el_rank", "Best MT EL Score", "%ile MT"], "99"),
                99.0,
            )
        ),
        "presentation_score": first(row, ["presentation_score", "Presentation Score"], "0.0"),
        "immunogenicity_score": first(
            row, ["immunogenicity_score", "Immunogenicity Score"], "0.5"
        ),
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
        "self_similarity_score": first(
            row, ["self_similarity_score", "Self Similarity"], "0.0"
        ),
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


def parse_pvactools_outputs(
    paths,
    sample_id,
    profile_name,
    out_events=None,
    out_peptides=None,
):
    event_rows: list[dict[str, str]] = []
    peptide_rows: list[dict[str, str]] = []
    for tsv in discover_tsvs(paths):
        tool = infer_tool(tsv)
        for row_number, row in enumerate(read_tsv(tsv), 1):
            event = event_from_row(
                row,
                sample_id,
                profile_name,
                tool,
                source_file=str(tsv),
                source_row_number=row_number,
            )
            event_rows.append(event)
            peptide_rows.append(
                peptide_from_row(
                    row,
                    sample_id,
                    event,
                    tool,
                    source_file=str(tsv),
                    source_row_number=row_number,
                )
            )

    events, _, _ = merge_rows_preserving_provenance(
        event_rows,
        EVENT_FIELDS,
        ("event_id",),
        entity_type="pvactools_event",
    )
    peptides, _, _ = merge_rows_preserving_provenance(
        peptide_rows,
        PEPTIDE_FIELDS,
        ("event_id", "peptide", "hla_allele"),
        entity_type="pvactools_peptide",
    )
    if out_events:
        write_tsv(out_events, events, EVENT_FIELDS)
    if out_peptides:
        write_tsv(out_peptides, peptides, PEPTIDE_FIELDS)
    return events, peptides
