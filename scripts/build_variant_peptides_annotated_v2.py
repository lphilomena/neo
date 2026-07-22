#!/usr/bin/env python3
"""Build a release-ready annotated peptide catalog from a phased NeoAg run.

The current ranked peptide table is authoritative. A legacy annotated catalog
is used only to recover matching VEP and minigene design fields.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


LEGACY_FIELDS = [
    "peptide_id", "gene", "ensembl_gene_id", "transcript_id", "hgvsc", "hgvsp",
    "chrom", "pos", "ref", "alt", "vaf", "tumor_depth", "tumor_alt_count",
    "rna_vaf", "rna_alt_reads", "rna_depth", "consequence", "protein_position",
    "amino_acids", "multi_aa_flag", "peptide_length", "peptide_start_aa",
    "peptide_end_aa", "mutation_position_in_peptide", "mutant_peptide",
    "wildtype_peptide", "minigene", "minigene_nt", "in_normal_proteome",
    "sample_hla_alleles", "hla_allele", "netmhcpan_mt_ic50",
    "netmhcpan_mt_rank_ba", "netmhcpan_mt_rank_el", "netmhcpan_wt_ic50",
    "netmhcpan_wt_rank_ba", "netmhcpan_wt_rank_el", "mhcflurry_mt_affinity",
    "mhcflurry_mt_affinity_percentile", "mhcflurry_mt_processing_score",
    "mhcflurry_mt_presentation_score", "mhcflurry_wt_affinity",
    "mhcflurry_wt_affinity_percentile", "mhcflurry_wt_processing_score",
    "mhcflurry_wt_presentation_score", "variant_key", "peptide_label",
    "peptide_source", "generation_method", "crosses_junction", "contains_novel_aa",
    "fusion_window_type", "fusion_breakpoint_position_raw", "fusion_generation_method",
    "netmhcstabpan_score", "netmhcstabpan_rank", "netmhcstabpan_wt_score",
    "netmhcstabpan_wt_rank", "prime_score", "prime_rank", "prime_wt_score",
    "prime_wt_rank", "bigmhc_im_score", "bigmhc_im_wt_score",
    "iedb_immunogenicity_score", "iedb_immunogenicity_wt_score",
]

BUILD_FIELDS = [
    "annotation_source", "legacy_annotation_match", "minigene_design_status",
]


def read_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader), list(reader.fieldnames or [])


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def first_nonempty(*values: str) -> str:
    for value in values:
        if str(value or "").strip():
            return str(value)
    return ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ranked-peptides", type=Path, required=True)
    parser.add_argument("--ranked-events", type=Path, required=True)
    parser.add_argument("--legacy-annotated", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--qc-out", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    ranked, ranked_fields = read_tsv(args.ranked_peptides)
    events, _ = read_tsv(args.ranked_events)
    legacy, _ = read_tsv(args.legacy_annotated)
    event_map = {row.get("event_id", ""): row for row in events}

    exact: dict[tuple[str, str, str, str], dict[str, str]] = {}
    event_peptide: dict[tuple[str, str], dict[str, str]] = {}
    gene_peptide_rows: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in legacy:
        event_id = row.get("variant_key", "")
        peptide = row.get("mutant_peptide", "").upper()
        gene = row.get("gene", "").upper()
        hla = row.get("hla_allele", "").upper()
        exact[(event_id, gene, peptide, hla)] = row
        event_peptide.setdefault((event_id, peptide), row)
        gene_peptide_rows[(gene, peptide)].append(row)

    hla_background = ",".join(sorted({row.get("hla_allele", "") for row in ranked if row.get("hla_allele")}))
    match_counts: Counter[str] = Counter()
    output: list[dict[str, str]] = []
    for current in ranked:
        event_id = current.get("event_id", "")
        gene = current.get("gene", "")
        peptide = current.get("peptide", "")
        hla = current.get("hla_allele", "")
        match = exact.get((event_id, gene.upper(), peptide.upper(), hla.upper()))
        match_type = "exact_event_peptide_hla" if match else ""
        if not match:
            match = event_peptide.get((event_id, peptide.upper()))
            match_type = "event_peptide" if match else ""
        if not match:
            candidates = gene_peptide_rows.get((gene.upper(), peptide.upper()), [])
            if candidates:
                match = candidates[0]
                match_type = "gene_peptide"
        match_counts[match_type or "none"] += 1

        row = {field: "" for field in LEGACY_FIELDS}
        if match:
            row.update({field: match.get(field, "") for field in LEGACY_FIELDS})
        event = event_map.get(event_id, {})

        row.update({
            "peptide_id": current.get("peptide_id", ""),
            "gene": gene,
            "transcript_id": first_nonempty(event.get("transcript_id", ""), row.get("transcript_id", "")),
            "chrom": first_nonempty(event.get("chrom", ""), row.get("chrom", "")),
            "pos": first_nonempty(event.get("pos", ""), row.get("pos", "")),
            "ref": first_nonempty(event.get("ref", ""), row.get("ref", "")),
            "alt": first_nonempty(event.get("alt", ""), row.get("alt", "")),
            "vaf": first_nonempty(event.get("tumor_vaf", ""), row.get("vaf", "")),
            "tumor_depth": first_nonempty(event.get("tumor_depth", ""), row.get("tumor_depth", "")),
            "tumor_alt_count": first_nonempty(event.get("tumor_alt_count", ""), row.get("tumor_alt_count", "")),
            "rna_vaf": current.get("rna_vaf", ""),
            "rna_alt_reads": current.get("rna_alt_reads", ""),
            "rna_depth": current.get("rna_depth", ""),
            "consequence": first_nonempty(event.get("consequence", ""), current.get("peptide_consequence", ""), row.get("consequence", "")),
            "peptide_length": str(len(peptide)) if peptide else "",
            "mutation_position_in_peptide": first_nonempty(current.get("mutation_positions_in_peptide", ""), row.get("mutation_position_in_peptide", "")),
            "mutant_peptide": peptide,
            "wildtype_peptide": current.get("wildtype_peptide", ""),
            "in_normal_proteome": {"yes": "yes", "no": "no"}.get(current.get("reference_proteome_exact_match", ""), "not_assessed"),
            "sample_hla_alleles": hla_background,
            "hla_allele": hla,
            "variant_key": event_id,
            "peptide_source": first_nonempty(current.get("mutation_source", ""), current.get("event_type", "")).lower(),
            "generation_method": first_nonempty(current.get("source_tool", ""), row.get("generation_method", "")),
            "crosses_junction": current.get("crosses_junction", ""),
            "contains_novel_aa": current.get("contains_novel_aa", ""),
            "netmhcpan_mt_ic50": current.get("netmhcpan_mt_ic50", ""),
            "netmhcpan_mt_rank_ba": current.get("netmhcpan_mt_rank_ba", ""),
            "netmhcpan_mt_rank_el": current.get("netmhcpan_mt_rank_el", ""),
            "netmhcpan_wt_ic50": current.get("netmhcpan_wt_ic50", ""),
            "netmhcpan_wt_rank_ba": current.get("netmhcpan_wt_rank_ba", ""),
            "netmhcpan_wt_rank_el": current.get("netmhcpan_wt_rank_el", ""),
            "mhcflurry_mt_affinity_percentile": current.get("mhcflurry_affinity_percentile", ""),
            "mhcflurry_mt_processing_score": current.get("mhcflurry_processing_score", ""),
            "mhcflurry_mt_presentation_score": current.get("mhcflurry_presentation_score", ""),
            "mhcflurry_wt_affinity_percentile": current.get("mhcflurry_wt_affinity_percentile", ""),
            "mhcflurry_wt_processing_score": current.get("mhcflurry_wt_processing_score", ""),
            "mhcflurry_wt_presentation_score": current.get("mhcflurry_wt_presentation_score", ""),
            "netmhcstabpan_score": current.get("netmhcstabpan_score", ""),
            "netmhcstabpan_rank": current.get("netmhcstabpan_rank", ""),
            "prime_score": current.get("prime_score", ""),
            "prime_rank": current.get("prime_rank", ""),
            "prime_wt_score": current.get("prime_wt_score", ""),
            "prime_wt_rank": current.get("prime_wt_rank", ""),
            "bigmhc_im_score": current.get("bigmhc_im_score", ""),
            "bigmhc_im_wt_score": current.get("bigmhc_im_wt_score", ""),
            "iedb_immunogenicity_score": current.get("iedb_immunogenicity_score", ""),
        })

        phased = current.get("haplotype_status") == "PHASED_CIS_COMBINED"
        if phased:
            minigene_status = "REQUIRES_PHASED_MINIGENE_REDESIGN"
            row["minigene"] = ""
            row["minigene_nt"] = ""
        elif row.get("minigene"):
            minigene_status = "REUSED_MATCHED_LEGACY_DESIGN"
        else:
            minigene_status = "UNAVAILABLE_REQUIRES_DESIGN"
        row["annotation_source"] = "phased_run" + (f"+legacy_{match_type}" if match_type else "")
        row["legacy_annotation_match"] = "yes" if match else "no"
        row["minigene_design_status"] = minigene_status
        for field in ranked_fields:
            row[field] = current.get(field, "")
        output.append(row)

    output_fields = LEGACY_FIELDS + BUILD_FIELDS + [field for field in ranked_fields if field not in LEGACY_FIELDS]
    write_tsv(args.out, output, output_fields)

    qc_rows = [
        {"metric": "rows", "value": str(len(output))},
        {"metric": "unique_peptide_ids", "value": str(len({row.get('peptide_id', '') for row in output}))},
        {"metric": "phased_tbr1_rows", "value": str(sum(row.get("gene") == "TBR1" and row.get("haplotype_status") == "PHASED_CIS_COMBINED" for row in output))},
        {"metric": "unphased_tbr1_rows", "value": str(sum(row.get("gene") == "TBR1" and row.get("haplotype_status") != "PHASED_CIS_COMBINED" for row in output))},
        {"metric": "rna_alt_reads_assessed", "value": str(sum(bool(row.get("rna_alt_reads", "").strip()) for row in output))},
        {"metric": "safety_status_assessed", "value": str(sum(bool(row.get("safety_status", "").strip()) for row in output))},
        {"metric": "final_priority_assessed", "value": str(sum(bool(row.get("final_priority", "").strip()) for row in output))},
    ]
    qc_rows.extend({"metric": f"legacy_match_{key}", "value": str(value)} for key, value in sorted(match_counts.items()))
    write_tsv(args.qc_out, qc_rows, ["metric", "value"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
