#!/usr/bin/env python3
"""Add STAR, BAM, usage, and normal-panel QC to SNAF splice candidates."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path

import pysam


STRAND = {0: ".", 1: "+", 2: "-"}


def chrom(value: str) -> str:
    value = value.strip()
    return value if value.startswith("chr") else f"chr{value}"


def snaf_key(row: dict[str, str]) -> tuple[str, int, int]:
    # SNAF exports outer exon-boundary coordinates. STAR SJ uses the closed
    # intron interval, hence +1/-1.
    return chrom(row["chrom"]), int(row["start"]) + 1, int(row["end"]) - 1


def median_hist(hist: Counter[int]) -> str:
    total = sum(hist.values())
    if not total:
        return ""
    target = (total - 1) // 2
    seen = 0
    for value in sorted(hist):
        seen += hist[value]
        if seen > target:
            return str(value)
    return ""


def read_candidates(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        return rows, list(reader.fieldnames or [])


def load_star(path: Path, targets: set[tuple[str, int, int]]):
    hits: dict[tuple[str, int, int], list[dict[str, str]]] = defaultdict(list)
    donor_total: Counter[tuple[str, int, str]] = Counter()
    acceptor_total: Counter[tuple[str, int, str]] = Counter()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            cells = line.rstrip("\n").split("\t")
            if len(cells) < 9:
                continue
            key = (chrom(cells[0]), int(cells[1]), int(cells[2]))
            strand = STRAND.get(int(cells[3]), ".")
            unique = int(cells[6])
            multi = int(cells[7])
            total = unique + multi
            donor_total[(key[0], key[1], strand)] += total
            acceptor_total[(key[0], key[2], strand)] += total
            if key in targets:
                hits[key].append({
                    "junction_strand": strand,
                    "splice_motif_code": cells[4],
                    "star_annotation_code": cells[5],
                    "star_unique_reads": str(unique),
                    "star_multi_reads": str(multi),
                    "star_total_reads": str(total),
                    "max_overhang": cells[8],
                })
    return hits, donor_total, acceptor_total


def bam_qc(path: Path, targets: set[tuple[str, int, int]]):
    stats = defaultdict(lambda: {
        "support": 0, "nondup": 0, "duplicates": 0, "starts": set(),
        "mapq": Counter(), "max_anchor": 0,
    })
    with pysam.AlignmentFile(str(path), "rb") as bam:
        for read in bam.fetch(until_eof=False):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            ref = read.reference_start
            cigars = read.cigartuples or []
            for index, (operation, length) in enumerate(cigars):
                if operation == 3:
                    key = (chrom(read.reference_name), ref + 1, ref + length)
                    if key in targets:
                        left = cigars[index - 1][1] if index and cigars[index - 1][0] in {0, 7, 8} else 0
                        right = cigars[index + 1][1] if index + 1 < len(cigars) and cigars[index + 1][0] in {0, 7, 8} else 0
                        item = stats[key]
                        item["support"] += 1
                        item["duplicates"] += int(read.is_duplicate)
                        if not read.is_duplicate:
                            item["nondup"] += 1
                            item["starts"].add((read.reference_start, read.is_reverse))
                            item["mapq"][read.mapping_quality] += 1
                            item["max_anchor"] = max(item["max_anchor"], min(left, right))
                if operation in {0, 2, 3, 7, 8}:
                    ref += length
    return stats


def load_normal(path: Path, targets: set[tuple[str, int, int]]):
    result: dict[tuple[str, int, int], dict[str, str]] = {}
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            try:
                key = (chrom(row.get("chromosome") or row.get("chrom") or ""), int(row["start"]), int(row["end"]))
            except (KeyError, TypeError, ValueError):
                continue
            if key in targets:
                result[key] = row
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--snaf", required=True, type=Path)
    parser.add_argument("--star-sj", required=True, type=Path)
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--normal-junctions", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--star-output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()

    rows, original_fields = read_candidates(args.snaf)
    targets = {snaf_key(row) for row in rows}
    star_hits, donor_total, acceptor_total = load_star(args.star_sj, targets)
    bam_stats = bam_qc(args.bam, targets)
    normal = load_normal(args.normal_junctions, targets)

    extra = [
        "canonical_chrom", "canonical_intron_start", "canonical_intron_end", "strand", "junction_strand",
        "splice_motif_code", "star_annotation_code", "provided_snaf_junction_reads",
        "star_unique_reads", "star_multi_reads", "star_total_reads", "max_overhang",
        "bam_junction_support_reads", "bam_nonduplicate_support_reads", "bam_duplicate_support_reads",
        "unique_start_count", "anchor_size", "junction_mapq", "junction_mapq_min", "junction_mapq_max",
        "psi", "psi_donor", "psi_acceptor", "psi_method",
        "normal_panel_detection_status", "normal_panel_samples", "normal_panel_reads",
        "normal_panel_tissues", "normal_junction_coverage_status", "junction_qc_status",
    ]
    fields = original_fields + [field for field in extra if field not in original_fields]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    star_rows: dict[tuple[str, int, int], dict[str, str]] = {}
    exact = 0
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            key = snaf_key(row)
            row["provided_snaf_junction_reads"] = row.get("junction_reads", "")
            row["canonical_chrom"], row["canonical_intron_start"], row["canonical_intron_end"] = map(str, key)
            matches = star_hits.get(key, [])
            if len(matches) == 1:
                exact += 1
                row.update(matches[0])
                row["strand"] = row["junction_strand"]
                row["junction_qc_status"] = "STAR_EXACT_MATCH"
            elif len(matches) > 1:
                row["junction_qc_status"] = "STAR_STRAND_CONFLICT"
            else:
                row["junction_qc_status"] = "STAR_NOT_DETECTED"
            strand = row.get("junction_strand", ".")
            star_total = int(row.get("star_total_reads") or 0)
            donor = donor_total[(key[0], key[1], strand)]
            acceptor = acceptor_total[(key[0], key[2], strand)]
            psi_d = star_total / donor if donor else 0.0
            psi_a = star_total / acceptor if acceptor else 0.0
            row["psi_donor"] = f"{psi_d:.6g}" if donor else ""
            row["psi_acceptor"] = f"{psi_a:.6g}" if acceptor else ""
            row["psi"] = f"{min(psi_d, psi_a):.6g}" if donor and acceptor else ""
            row["psi_method"] = "MIN_STAR_DONOR_ACCEPTOR_USAGE" if row["psi"] else "UNASSESSED"
            bq = bam_stats.get(key)
            if bq:
                row["bam_junction_support_reads"] = str(bq["support"])
                row["bam_nonduplicate_support_reads"] = str(bq["nondup"])
                row["bam_duplicate_support_reads"] = str(bq["duplicates"])
                row["unique_start_count"] = str(len(bq["starts"]))
                row["anchor_size"] = str(bq["max_anchor"])
                row["junction_mapq"] = median_hist(bq["mapq"])
                row["junction_mapq_min"] = str(min(bq["mapq"])) if bq["mapq"] else ""
                row["junction_mapq_max"] = str(max(bq["mapq"])) if bq["mapq"] else ""
            nrow = normal.get(key)
            if nrow:
                row["normal_panel_detection_status"] = "DETECTED_BROAD_NORMAL"
                row["normal_panel_samples"] = nrow.get("normal_samples", "")
                row["normal_panel_reads"] = nrow.get("normal_reads", "")
                row["normal_panel_tissues"] = nrow.get("normal_tissues", "")
                row["normal_junction_coverage_status"] = "DETECTED_COVERAGE_PRESENT"
            else:
                row["normal_panel_detection_status"] = "NOT_DETECTED"
                row["normal_junction_coverage_status"] = "NOT_DETECTED_COVERAGE_UNASSESSED"
            writer.writerow(row)
            if len(matches) == 1:
                star_rows[key] = row

    star_fields = [
        "chrom", "start", "end", "strand", "junction_reads", "unique_split_reads",
        "multi_split_reads", "total_split_reads", "max_overhang", "unique_start_count",
        "anchor_size", "junction_mapq", "junction_mapq_min", "junction_mapq_max",
        "psi", "psi_donor", "psi_acceptor", "psi_method", "normal_panel_detection_status",
        "normal_panel_samples", "normal_panel_reads", "normal_panel_tissues",
        "normal_junction_coverage_status", "source_coordinate_system",
    ]
    with args.star_output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, star_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for key, row in sorted(star_rows.items()):
            writer.writerow({
                "chrom": key[0], "start": key[1], "end": key[2],
                "strand": row.get("junction_strand", "."),
                "junction_reads": row.get("star_total_reads", "0"),
                "unique_split_reads": row.get("star_unique_reads", "0"),
                "multi_split_reads": row.get("star_multi_reads", "0"),
                "total_split_reads": row.get("star_total_reads", "0"),
                "max_overhang": row.get("max_overhang", ""),
                "unique_start_count": row.get("unique_start_count", ""),
                "anchor_size": row.get("anchor_size", ""),
                "junction_mapq": row.get("junction_mapq", ""),
                "junction_mapq_min": row.get("junction_mapq_min", ""),
                "junction_mapq_max": row.get("junction_mapq_max", ""),
                "psi": row.get("psi", ""), "psi_donor": row.get("psi_donor", ""),
                "psi_acceptor": row.get("psi_acceptor", ""), "psi_method": row.get("psi_method", ""),
                "normal_panel_detection_status": row.get("normal_panel_detection_status", ""),
                "normal_panel_samples": row.get("normal_panel_samples", ""),
                "normal_panel_reads": row.get("normal_panel_reads", ""),
                "normal_panel_tissues": row.get("normal_panel_tissues", ""),
                "normal_junction_coverage_status": row.get("normal_junction_coverage_status", ""),
                "source_coordinate_system": "intron_1based_closed",
            })

    summary = {
        "candidate_rows": len(rows), "unique_candidate_junctions": len(targets),
        "star_exact_candidate_junctions": len(star_rows), "star_exact_candidate_rows": exact,
        "bam_qc_candidate_junctions": len(bam_stats), "normal_detected_candidate_junctions": len(normal),
        "psi_definition": "minimum of STAR junction usage among junctions sharing the same donor or acceptor",
        "normal_absence_interpretation": "NOT_DETECTED_COVERAGE_UNASSESSED when no exact panel row exists",
    }
    args.summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
