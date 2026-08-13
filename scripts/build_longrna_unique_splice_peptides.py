#!/usr/bin/env python3
"""Build junction-spanning peptides for SQANTI3-only abnormal isoforms."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path

ABNORMAL = {
    "novel_in_catalog",
    "novel_not_in_catalog",
    "antisense",
    "intergenic",
    "fusion",
    "genic",
    "genic_intron",
}
AA = set("ACDEFGHIKLMNPQRSTVWY")


def number(value: str | None, default: float = 0.0) -> float:
    try:
        return float(value or default)
    except (TypeError, ValueError):
        return default


def rows(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def load_hla(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text().splitlines() if line.strip() and not line.startswith("#")]


def load_proteins(path: Path) -> dict[str, str]:
    proteins: dict[str, str] = {}
    isoform = ""
    sequence: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if isoform:
                    proteins.setdefault(isoform, "".join(sequence).rstrip("*"))
                token = line[1:].split()[0]
                isoform = re.sub(r"\.p\d+$", "", token)
                sequence = []
            else:
                sequence.append(line.strip())
    if isoform:
        proteins.setdefault(isoform, "".join(sequence).rstrip("*"))
    return proteins


def load_exons(path: Path, wanted: set[str]) -> dict[str, list[tuple[int, int, str]]]:
    exons = defaultdict(list)
    pattern = re.compile(r'transcript_id "([^"]+)"')
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            match = pattern.search(fields[8])
            if match and match.group(1) in wanted:
                exons[match.group(1)].append((int(fields[3]), int(fields[4]), fields[6]))
    for isoform, items in exons.items():
        items.sort(key=lambda item: item[0], reverse=items[0][2] == "-")
    return exons


def junction_tx_boundaries(exons: list[tuple[int, int, str]]) -> dict[tuple[int, int], int]:
    result = {}
    consumed = 0
    for left, right in zip(exons, exons[1:]):
        consumed += left[1] - left[0] + 1
        genomic = (min(left[1], right[1]) + 1, max(left[0], right[0]) - 1)
        result[genomic] = consumed
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sqanti-classification", required=True, type=Path)
    parser.add_argument("--sqanti-junctions", required=True, type=Path)
    parser.add_argument("--corrected-gtf", required=True, type=Path)
    parser.add_argument("--td2-peptides", required=True, type=Path)
    parser.add_argument("--shortread-comparisons", required=True, type=Path)
    parser.add_argument("--hla", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--peptide-lengths", default="", help="Explicit comma-separated MHC-I lengths")
    parser.add_argument("--high-recall-12mer", action="store_true", help="Use 8,9,10,11,12 unless explicit lengths are provided")
    args = parser.parse_args()
    if args.peptide_lengths:
        peptide_lengths = tuple(sorted({int(value) for value in args.peptide_lengths.split(",") if value.strip()}))
    else:
        peptide_lengths = (8, 9, 10, 11, 12) if args.high_recall_12mer else (8, 9, 10, 11)
    if not peptide_lengths or any(length < 8 or length > 14 for length in peptide_lengths):
        parser.error("--peptide-lengths must contain comma-separated values from 8 through 14")

    shortread_keys = {row.get("junction_key", "") for row in rows(args.shortread_comparisons)}
    classes = {
        row["isoform"]: row
        for row in rows(args.sqanti_classification)
        if row.get("structural_category") in ABNORMAL
    }
    selected = []
    for junction in rows(args.sqanti_junctions):
        isoform = junction["isoform"]
        if isoform not in classes:
            continue
        chrom = junction["chrom"]
        chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
        exon_key = (
            chrom,
            str(int(junction["genomic_start_coord"]) - 1),
            str(int(junction["genomic_end_coord"]) + 1),
            junction["strand"],
        )
        key = ":".join(exon_key)
        if key in shortread_keys:
            continue
        structural_category = classes[isoform]["structural_category"]
        unique_coverage = int(float(junction.get("total_coverage_unique") or 0))
        if structural_category != "fusion" and junction.get("junction_category") != "novel":
            continue
        if unique_coverage < 3 or junction.get("canonical", "").lower() != "canonical":
            continue
        item = dict(junction)
        item["junction_key"] = key
        item["structural_category"] = structural_category
        item["associated_gene"] = classes[isoform].get("associated_gene", "")
        selected.append(item)

    wanted = {item["isoform"] for item in selected}
    exons = load_exons(args.corrected_gtf, wanted)
    proteins = load_proteins(args.td2_peptides)
    hla = load_hla(args.hla)
    boundaries = {isoform: junction_tx_boundaries(items) for isoform, items in exons.items()}

    selected.sort(
        key=lambda item: (
            item["junction_key"],
            -int(float(item.get("total_coverage_unique") or 0)),
            -number(classes[item["isoform"]].get("psauron_score")),
        )
    )
    peptide_rows = []
    prediction_pairs: dict[tuple[str, str], dict[str, str]] = {}
    seen_event_pairs: set[tuple[str, str, str]] = set()
    reason_counts = defaultdict(int)
    for item in selected:
        isoform = item["isoform"]
        cls = classes[isoform]
        if item["structural_category"] == "fusion":
            reason_counts["fusion_requires_jaffal_breakpoint_translation"] += 1
            continue
        if cls.get("CDS_type") != "complete":
            reason_counts["incomplete_cds"] += 1
            continue
        if cls.get("predicted_NMD", "").upper() == "TRUE":
            reason_counts["predicted_nmd"] += 1
            continue
        try:
            if number(cls.get("psauron_score")) < 0.8:
                reason_counts["low_psauron_score"] += 1
                continue
        except ValueError:
            reason_counts["low_psauron_score"] += 1
            continue
        protein = proteins.get(isoform, "")
        try:
            cds_start = int(float(cls.get("CDS_start", "")))
            cds_end = int(float(cls.get("CDS_end", "")))
        except ValueError:
            reason_counts["non_coding_or_missing_cds"] += 1
            continue
        genomic = (int(item["genomic_start_coord"]), int(item["genomic_end_coord"]))
        tx_boundary = boundaries.get(isoform, {}).get(genomic)
        if tx_boundary is None:
            reason_counts["junction_not_mapped_to_gtf"] += 1
            continue
        if not (cds_start <= tx_boundary < cds_end):
            reason_counts["junction_outside_cds"] += 1
            continue
        if not protein or not set(protein).issubset(AA):
            reason_counts["missing_or_nonstandard_protein"] += 1
            continue

        coding_nt_before = tx_boundary - cds_start + 1
        aa_cut = coding_nt_before // 3
        starts_by_length = defaultdict(set)
        for length in peptide_lengths:
            for start in range(max(0, aa_cut - length + 1), min(aa_cut + 1, len(protein) - length + 1)):
                end = start + length
                if start < aa_cut + 1 and end > max(0, aa_cut - 1):
                    starts_by_length[length].add(start)
        for length, starts in starts_by_length.items():
            for start in sorted(starts):
                peptide = protein[start : start + length]
                if len(peptide) != length or not set(peptide).issubset(AA):
                    continue
                for allele in hla:
                    event_id = f"LR_SPLICE_{item['junction_key']}"
                    event_pair = (event_id, peptide, allele)
                    if event_pair in seen_event_pairs:
                        continue
                    seen_event_pairs.add(event_pair)
                    output_row = {
                            "sample_id": "M1ML150017383",
                            "event_id": event_id,
                            "peptide_id": f"{isoform}_{item['junction_number']}_{length}_{start + 1}",
                            "peptide": peptide,
                            "hla_allele": allele,
                            "mhc_class": "I",
                            "event_type": "FUSION" if item["structural_category"] == "fusion" else "SPLICE",
                            "mutation_source": "long_read_rna",
                            "source_tool": "IsoQuant+SQANTI3+TD2",
                            "generation_method": f"neoag_junction_spanning_mhc1_{','.join(map(str, peptide_lengths))}mer",
                            "isoform": isoform,
                            "gene": item["associated_gene"],
                            "structural_category": item["structural_category"],
                            "junction_key": item["junction_key"],
                            "junction_category": item["junction_category"],
                            "canonical": item["canonical"],
                            "longrna_unique_coverage": item.get("total_coverage_unique", ""),
                            "longrna_multi_coverage": item.get("total_coverage_multi", ""),
                            "junction_tx_boundary": str(tx_boundary),
                            "cds_start": str(cds_start),
                            "cds_end": str(cds_end),
                            "protein_length": str(len(protein)),
                        }
                    peptide_rows.append(output_row)
                    prediction_pairs.setdefault(
                        (peptide, allele),
                        {
                            "sample_id": "M1ML150017383",
                            "peptide": peptide,
                            "hla_allele": allele,
                            "source_tool": "IsoQuant+SQANTI3+TD2",
                        },
                    )

    args.outdir.mkdir(parents=True, exist_ok=True)
    junction_fields = list(selected[0]) if selected else ["isoform", "junction_key"]
    with (args.outdir / "longrna_unique_abnormal_junctions.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=junction_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(selected)
    peptide_fields = list(peptide_rows[0]) if peptide_rows else ["sample_id", "peptide", "hla_allele"]
    with (args.outdir / "longrna_unique_junction_peptides.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=peptide_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(peptide_rows)
    prediction_fields = ["sample_id", "peptide", "hla_allele", "source_tool"]
    with (args.outdir / "longrna_unique_prediction_input.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=prediction_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(prediction_pairs.values())

    summary = {
        "status": "PASS",
        "abnormal_categories": sorted(ABNORMAL),
        "shortread_junction_keys": len(shortread_keys),
        "longrna_unique_abnormal_junction_rows": len(selected),
        "longrna_unique_abnormal_isoforms": len(wanted),
        "junction_peptide_hla_rows": len(peptide_rows),
        "unique_peptides": len({row["peptide"] for row in peptide_rows}),
        "unique_peptide_hla_pairs": len(prediction_pairs),
        "mhc1_peptide_lengths": list(peptide_lengths),
        "high_recall_12mer": 12 in peptide_lengths,
        "selection_policy": "novel junction (fusion excepted), unique coverage >=3, canonical, complete CDS, PSAURON >=0.8, non-NMD",
        "excluded_reasons": dict(reason_counts),
    }
    (args.outdir / "longrna_unique_candidate_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
