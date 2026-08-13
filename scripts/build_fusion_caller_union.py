#!/usr/bin/env python3
"""Build a provenance-preserving fusion event/peptide union from completed callers."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

from neoag.input_router import build_raw_intermediates
from neoag.model_layers import enrich_event_layers, enrich_peptide_layers, infer_mutation_source, infer_peptide_consequence
from neoag.provenance import merge_rows_preserving_provenance
from neoag.schemas import EVENT_FIELDS, PEPTIDE_FIELDS
from neoag.utils import first, safe_id, to_float, write_tsv

HLA_RE = re.compile(r"(?:HLA-)?(?:A|B|C)\*[0-9]{2,3}(?::[0-9A-Z]{2,3}){1,4}", re.I)
AA_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")


def read_hla(path: Path) -> list[str]:
    values = HLA_RE.findall(path.read_text(encoding="utf-8", errors="replace"))
    return list(dict.fromkeys(value.upper() if value.upper().startswith("HLA-") else "HLA-" + value.upper() for value in values))


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        header = handle.readline()
        handle.seek(0)
        delimiter = "\t" if "\t" in header else ","
        return [{str(key): str(value or "") for key, value in row.items()} for row in csv.DictReader(handle, delimiter=delimiter)]


def clean_gene(value: str) -> str:
    return str(value or "").split("^", 1)[0].strip()


def gene_pair(row: dict[str, str]) -> tuple[str, str]:
    combined = first(row, ["FusionName", "#FusionName", "fusion", "fusion_name", "Fusion_Gene"], "")
    for sep in ("::", "--", "_"):
        if sep in combined:
            left, right = combined.split(sep, 1)
            return clean_gene(left), clean_gene(right)
    return (
        clean_gene(first(row, ["LeftGene", "left_gene", "gene1", "Gene1", "gene5"], "")),
        clean_gene(first(row, ["RightGene", "right_gene", "gene2", "Gene2", "gene3"], "")),
    )


def generic_caller_rows(path: Path, tool: str, sample_id: str, profile: str, hla: list[str]):
    events, peptides, audit = [], [], []
    for index, row in enumerate(read_table(path), 1):
        left_gene, right_gene = gene_pair(row)
        if not left_gene or not right_gene:
            continue
        left_bp = first(row, ["LeftBreakpoint", "breakpoint1", "left_breakpoint", "breakpoint_1"], "")
        right_bp = first(row, ["RightBreakpoint", "breakpoint2", "right_breakpoint", "breakpoint_2"], "")
        direction = first(row, ["direction", "strand", "Strand1", "strand1(gene/fusion)"], "") + "/" + first(row, ["Strand2", "strand2(gene/fusion)"], "")
        pair = f"{left_gene}::{right_gene}"
        event_id = safe_id(f"FUSION|{pair}|{left_bp}|{right_bp}|{direction}")
        reads = first(row, ["JunctionReadCount", "junction_reads", "split_reads", "supporting_reads", "split_reads1"], "")
        frame = first(row, ["frame", "reading_frame", "in_frame", "reading_frame_status"], "")
        base = {field: "" for field in EVENT_FIELDS}
        base.update({
            "event_id": event_id, "sample_id": sample_id, "disease_profile": profile,
            "event_type": "Fusion", "gene": pair, "event_name": pair,
            "consequence": frame or "fusion_orf_unassessed", "rna_junction_reads": reads,
            "event_confidence": "0.7", "event_expression": "0.0", "driver_relevance": "0.0",
            "clonality": "0.5", "persistence": "0.5", "tumor_specificity": "0.7",
            "source": f"{tool}:{path}", "source_file": str(path), "source_row_number": str(index),
            "source_tools": tool, "mutation_source": infer_mutation_source(event_type="Fusion", tool=tool),
            "peptide_consequence": infer_peptide_consequence(event_type="Fusion", consequence="fusion", tool=tool),
        })
        events.append(enrich_event_layers(base))
        sequence = first(row, ["junction_peptide", "neo_peptide_sequence", "fusion_peptide", "mutant_peptide", "peptide"], "").upper().replace("*", "")
        explicit_hla = first(row, ["hla", "hla_allele", "allele"], "")
        row_hla = read_hla_text(explicit_hla) or hla
        windows: list[str] = []
        if AA_RE.fullmatch(sequence or ""):
            if 8 <= len(sequence) <= 12:
                windows = [sequence]
            else:
                windows = list(dict.fromkeys(sequence[start:start + length] for length in (8, 9, 10, 11) for start in range(max(0, len(sequence) - length + 1))))
        for peptide in windows:
            for allele in row_hla:
                pbase = {field: "" for field in PEPTIDE_FIELDS}
                pbase.update({
                    "peptide_id": safe_id(f"{event_id}|{allele}|{peptide}"), "event_id": event_id,
                    "sample_id": sample_id, "event_type": "Fusion", "gene": pair,
                    "peptide": peptide, "hla_allele": allele, "mhc_class": "I",
                    "source_tool": tool, "source_file": str(path), "crosses_junction": "true",
                    "contains_novel_aa": "true", "rna_junction_reads": reads,
                    "mutation_source": base["mutation_source"], "peptide_consequence": base["peptide_consequence"],
                    "binding_rank": "99", "el_rank": "99", "presentation_score": "0.0",
                    "immunogenicity_score": "0.5", "wildtype_binding_rank": "99", "self_similarity_score": "0.0",
                })
                peptides.append(enrich_peptide_layers(pbase))
        audit.append({"event_id": event_id, "gene_pair": pair, "left_breakpoint": left_bp, "right_breakpoint": right_bp, "direction": direction, "source_tool": tool, "source_file": str(path), "source_row": str(index), "peptide_status": "PROVIDED_ORF_PEPTIDE" if windows else "ORF_PEPTIDE_UNAVAILABLE_REVIEW_ONLY"})
    return events, peptides, audit


def read_hla_text(value: str) -> list[str]:
    return list(dict.fromkeys(match.upper() if match.upper().startswith("HLA-") else "HLA-" + match.upper() for match in HLA_RE.findall(value or "")))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True); parser.add_argument("--profile", default="default")
    parser.add_argument("--hla-file", required=True, type=Path); parser.add_argument("--easyfuse", type=Path)
    parser.add_argument("--star-fusion", type=Path); parser.add_argument("--arriba", type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    hla = read_hla(args.hla_file)
    if not hla:
        raise SystemExit("HLA consensus has no class-I alleles")
    args.outdir.mkdir(parents=True, exist_ok=True)
    events, peptides, audit = [], [], []
    if args.easyfuse:
        cfg = {"sample": {"id": args.sample_id, "profile": args.profile}, "inputs": {"entry_mode": "fusion", "easyfuse_tsv": str(args.easyfuse.resolve()), "hla_alleles": hla}}
        easyfuse_out = args.outdir / "easyfuse"
        build_raw_intermediates(cfg, easyfuse_out, root=Path.cwd())
        events.extend(read_table(easyfuse_out / "parsed/raw_events.tsv"))
        peptides.extend(read_table(easyfuse_out / "parsed/raw_peptides.tsv"))
        for event in events:
            audit.append({"event_id": event.get("event_id", ""), "gene_pair": event.get("gene", ""), "left_breakpoint": "", "right_breakpoint": "", "direction": "", "source_tool": "EasyFuse", "source_file": str(args.easyfuse), "source_row": "", "peptide_status": "GENERATED_FROM_EASYFUSE_ORF"})
    for path, tool in ((args.star_fusion, "STAR-Fusion"), (args.arriba, "Arriba")):
        if path:
            e, p, a = generic_caller_rows(path, tool, args.sample_id, args.profile, hla)
            events.extend(e); peptides.extend(p); audit.extend(a)
    if not events:
        raise SystemExit("No fusion events were parsed from supplied caller outputs")
    merged_events, _, _ = merge_rows_preserving_provenance(events, EVENT_FIELDS, ("event_id",), entity_type="fusion_union_event")
    merged_peptides, _, _ = merge_rows_preserving_provenance(peptides, PEPTIDE_FIELDS, ("event_id", "peptide", "hla_allele"), entity_type="fusion_union_peptide")
    write_tsv(args.outdir / "raw_events.tsv", merged_events, EVENT_FIELDS)
    write_tsv(args.outdir / "raw_peptides.tsv", merged_peptides, PEPTIDE_FIELDS)
    write_tsv(args.outdir / "fusion_caller_union.tsv", audit, ["event_id", "gene_pair", "left_breakpoint", "right_breakpoint", "direction", "source_tool", "source_file", "source_row", "peptide_status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
