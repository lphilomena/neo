#!/usr/bin/env python3
"""Enrich splice candidates with exact STAR/BAM junction QC evidence."""

from __future__ import annotations

import argparse
import csv
import re
import sqlite3
import statistics
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

COORD_RE = re.compile(r"^(chr[^:]+):(\d+)-(\d+):([+\-?])$")
CANONICAL_RE = re.compile(r"^SJ\|[^|]+\|(chr[^|]+)\|(\d+)\|(\d+)\|([+\-?])$")
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def read_rows(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader), list(reader.fieldnames or [])


def write_rows(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def add_fields(fields: list[str], additions: list[str]) -> list[str]:
    return list(dict.fromkeys([*fields, *additions]))


def norm_chrom(value: str) -> str:
    return value[3:] if value.lower().startswith("chr") else value


def event_junction(row: dict[str, str]) -> tuple[str, int, int, str] | None:
    """Resolve an exact junction from canonical, structured, or display fields."""
    for field in ("canonical_junction_id", "source_junction_id", "event_name"):
        value = str(row.get(field) or "").strip()
        match = CANONICAL_RE.match(value) or COORD_RE.match(value)
        if match:
            chrom, start, end, strand = match.groups()
            if strand in {"+", "-"}:
                return norm_chrom(chrom), int(start), int(end), strand
    chrom = str(row.get("junction_chrom") or row.get("chrom") or "").strip()
    start = str(row.get("junction_start") or "").strip()
    end = str(row.get("junction_end") or "").strip()
    strand = str(row.get("junction_strand") or "").strip()
    if chrom and start.isdigit() and end.isdigit() and strand in {"+", "-"}:
        return norm_chrom(chrom), int(start), int(end), strand
    return None


def load_star(path: Path):
    records = {}
    donor_totals = defaultdict(int)
    acceptor_totals = defaultdict(int)
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, start, end = norm_chrom(fields[0]), int(fields[1]), int(fields[2])
            strand = {"1": "+", "2": "-"}.get(fields[3], "?")
            unique_reads, multi_reads = int(fields[6]), int(fields[7])
            key = (chrom, start, end, strand)
            records[key] = {
                "star_motif": fields[4],
                "star_annotated": fields[5],
                "star_unique_reads": unique_reads,
                "star_multi_reads": multi_reads,
                "star_max_overhang": int(fields[8]),
            }
            donor_totals[(chrom, start, strand)] += unique_reads
            acceptor_totals[(chrom, end, strand)] += unique_reads
    return records, donor_totals, acceptor_totals


def load_crossvalidated_keys(path: Path) -> tuple[set[str], set[str]]:
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        snaf: set[str] = set()
        splicemutr: set[str] = set()
        for row in reader:
            tools = {
                item.strip().lower()
                for item in re.split(r"[;,]", str(row.get("support_tools") or row.get("source_tools") or ""))
                if item.strip()
            }
            status = str(row.get("status") or row.get("evidence_status") or "").upper()
            if status and not any(token in status for token in ("CONFIRMED", "COMPLETED", "PASS", "SUPPORTED")):
                continue
            for field in ("uid", "event_id", "canonical_junction_id", "source_junction_id"):
                value = str(row.get(field) or "").strip()
                if value:
                    if "snaf" in tools:
                        snaf.add(value)
                    if "splicemutr" in tools:
                        splicemutr.add(value)
        return snaf, splicemutr


def splicemutr_origin_events(rows: list[dict[str, str]]) -> set[str]:
    """Events with an exact peptide -> SpliceMutr translated-origin backlink."""
    result: set[str] = set()
    for row in rows:
        source = ";".join(
            str(row.get(field) or "")
            for field in ("source_tool", "source_tools", "source_generator", "source_records")
        ).lower()
        if str(row.get("origin_peptide_id") or "").strip() or "splicemutr" in source:
            event_id = str(row.get("event_id") or "").strip()
            if event_id:
                result.add(event_id)
    return result


def bam_contigs(samtools: str, bam_path: Path) -> dict[str, str]:
    process = subprocess.run(
        [samtools, "view", "-H", str(bam_path)],
        check=True,
        capture_output=True,
        text=True,
    )
    result: dict[str, str] = {}
    for line in process.stdout.splitlines():
        if not line.startswith("@SQ\t"):
            continue
        for field in line.split("\t"):
            if field.startswith("SN:"):
                contig = field[3:]
                result[norm_chrom(contig)] = contig
                break
    return result


def load_bam_metrics(samtools: str, bam_path: Path, keys: set[tuple[str, int, int, str]], workdir: Path):
    workdir.mkdir(parents=True, exist_ok=True)
    by_coord = defaultdict(list)
    for key in keys:
        by_coord[key[:3]].append(key)
    bed = workdir / "candidate_junction_windows.bed"
    contigs = bam_contigs(samtools, bam_path)
    with bed.open("w", encoding="utf-8") as handle:
        for chrom, start, end, _strand in sorted(keys):
            bam_chrom = contigs.get(norm_chrom(chrom), chrom)
            handle.write(f"{bam_chrom}\t{max(0, start - 12)}\t{start + 10}\n")
            handle.write(f"{bam_chrom}\t{max(0, end - 11)}\t{end + 10}\n")
    reads_by_key: dict[tuple[str, int, int, str], dict[str, tuple[int, int, int, int]]] = defaultdict(dict)
    command = [samtools, "view", "-M", "-L", str(bed), str(bam_path)]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    assert process.stdout is not None
    for line in process.stdout:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 11:
            continue
        qname, flag_text, chrom, pos_text, mapq_text, cigar = fields[:6]
        flag = int(flag_text)
        if flag & (0x4 | 0x100 | 0x400 | 0x800):
            continue
        pos, mapq = int(pos_text), int(mapq_text)
        pnext = int(fields[7]) if fields[7].isdigit() else 0
        tags = {item.split(":", 2)[0]: item.rsplit(":", 1)[-1] for item in fields[11:] if ":" in item}
        nh = int(tags.get("NH", "1"))
        operations = [(int(length), op) for length, op in CIGAR_RE.findall(cigar)]
        ref_pos = pos
        for index, (length, op) in enumerate(operations):
            if op == "N":
                intron_start, intron_end = ref_pos, ref_pos + length - 1
                matching = by_coord.get((norm_chrom(chrom), intron_start, intron_end), [])
                if matching:
                    left = operations[index - 1][0] if index > 0 and operations[index - 1][1] in {"M", "=", "X"} else 0
                    right = operations[index + 1][0] if index + 1 < len(operations) and operations[index + 1][1] in {"M", "=", "X"} else 0
                    fragment_start = min(value for value in (pos - 1, pnext - 1) if value >= 0)
                    current = (mapq, nh, min(left, right), fragment_start)
                    for key in matching:
                        previous = reads_by_key[key].get(qname)
                        if previous is None or current[0] > previous[0]:
                            reads_by_key[key][qname] = current
            if op in {"M", "D", "N", "=", "X"}:
                ref_pos += length
    stderr = process.stderr.read() if process.stderr is not None else ""
    if process.wait() != 0:
        raise RuntimeError(f"samtools view failed: {stderr.strip()}")
    result = {}
    for key in keys:
        values = list(reads_by_key.get(key, {}).values())
        unique = [value for value in values if value[1] <= 1]
        multi = [value for value in values if value[1] > 1]
        result[key] = {
            "bam_unique_split_reads": str(len(unique)),
            "bam_multi_split_reads": str(len(multi)),
            "bam_total_split_reads": str(len(values)),
            "unique_fragment_starts": str(len({value[3] for value in unique})),
            "median_mapq": f"{statistics.median(value[0] for value in values):g}" if values else "",
            "bam_max_overhang": str(max((value[2] for value in values), default=0)),
            "bam_multimapping_fraction": f"{len(multi) / len(values):.6g}" if values else "",
        }
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--events", required=True)
    parser.add_argument("--peptides", required=True)
    parser.add_argument("--star-sj", required=True)
    parser.add_argument("--rna-bam", required=True)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--crossvalidated-snaf-splicemutr")
    parser.add_argument("--splice-consensus", help="Consensus TSV containing exact SNAF + SpliceMutr support")
    parser.add_argument("--normal-junction-sqlite")
    parser.add_argument("--matched-normal-star-sj")
    parser.add_argument("--matched-normal-rna-bam")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--genome-build", default="GRCh38")
    parser.add_argument("--min-unique-reads", type=int, default=3)
    parser.add_argument("--min-unique-fragment-starts", type=int, default=2)
    parser.add_argument("--min-overhang", type=int, default=10)
    parser.add_argument("--min-mapq", type=float, default=20.0)
    parser.add_argument("--max-multimapping-fraction", type=float, default=0.20)
    parser.add_argument("--min-tumor-psi", type=float, default=0.05)
    args = parser.parse_args()

    validation_source = args.splice_consensus or args.crossvalidated_snaf_splicemutr
    if not validation_source:
        parser.error("one of --splice-consensus or --crossvalidated-snaf-splicemutr is required")
    if bool(args.matched_normal_star_sj) != bool(args.matched_normal_rna_bam):
        parser.error("--matched-normal-star-sj and --matched-normal-rna-bam must be supplied together")

    events, event_fields = read_rows(Path(args.events))
    peptides, peptide_fields = read_rows(Path(args.peptides))
    star, donor_totals, acceptor_totals = load_star(Path(args.star_sj))
    snaf_validated, splicemutr_validated = load_crossvalidated_keys(Path(validation_source))
    splicemutr_validated.update(splicemutr_origin_events(peptides))
    normal_db = sqlite3.connect(args.normal_junction_sqlite) if args.normal_junction_sqlite else None
    matched_normal_star = matched_normal_donors = matched_normal_acceptors = None
    if args.matched_normal_star_sj:
        matched_normal_star, matched_normal_donors, matched_normal_acceptors = load_star(Path(args.matched_normal_star_sj))
    metrics_by_event: dict[str, dict[str, str]] = {}
    qc_rows: list[dict[str, str]] = []

    splice_events = [row for row in events if str(row.get("event_type") or "").lower() == "splice"]
    junction_keys = set()
    for row in splice_events:
        junction = event_junction(row)
        if junction:
            junction_keys.add(junction)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    bam_by_key = load_bam_metrics(args.samtools, Path(args.rna_bam), junction_keys, outdir)
    matched_normal_bam_by_key = (
        load_bam_metrics(
            args.samtools,
            Path(args.matched_normal_rna_bam),
            junction_keys,
            outdir / "matched_normal",
        )
        if args.matched_normal_rna_bam
        else {}
    )
    for index, row in enumerate(splice_events, 1):
        event_id = str(row.get("event_id") or "")
        junction = event_junction(row)
        metrics: dict[str, str] = {"event_id": event_id}
        if not junction:
            metrics.update({"splice_alignment_qc_status": "FAIL", "junction_read_qc_status": "FAIL", "caller_filter_status": "FAIL", "junction_support_reason": "no exact stranded genomic junction in canonical or structured fields"})
            metrics_by_event[event_id] = metrics
            qc_rows.append(metrics)
            continue
        chrom, start, end, strand = junction
        key = (chrom, start, end, strand)
        star_row = star.get(key)
        uid = event_id.split("|", 1)[1] if "|" in event_id else event_id
        canonical = f"SJ|{args.genome_build}|chr{chrom}|{start}|{end}|{strand}"
        normal_id = f"chr{chrom}:{start}-{end}:{strand}"
        normal_record = None
        normal_has_frequency = False
        if normal_db:
            columns = {item[1] for item in normal_db.execute("pragma table_info(junction_ids)")}
            normal_has_frequency = "normal_samples" in columns
            if normal_has_frequency:
                normal_record = normal_db.execute(
                    "select normal_samples,normal_reads,normal_total_reads,normal_tissues,tissue,source,dataset "
                    "from junction_ids where junction_id=? limit 1", (normal_id,),
                ).fetchone()
            else:
                normal_record = normal_db.execute(
                    "select 1 from junction_ids where junction_id=? limit 1", (normal_id,)
                ).fetchone()
        normal_hit = bool(normal_record)
        normal_samples = int(normal_record[0] or 0) if normal_record and normal_has_frequency else 0
        normal_reads = int(normal_record[1] or 0) if normal_record and normal_has_frequency else 0
        normal_tissues = int(normal_record[3] or 0) if normal_record and normal_has_frequency else 0
        normal_status = (
            "DETECTED_BROAD_NORMAL" if normal_hit and normal_has_frequency and (normal_samples >= 2 or normal_tissues >= 2)
            else "LOW_LEVEL_NONCRITICAL_NORMAL" if normal_hit and normal_has_frequency
            else "DETECTED_NORMAL_CATALOG_FREQUENCY_UNAVAILABLE" if normal_hit
            else "NOT_LISTED_IN_NORMAL_CATALOG" if normal_db
            else "UNASSESSED_NO_NORMAL_COHORT"
        )
        metrics.update({
            "junction_chrom": f"chr{chrom}", "junction_start": str(start), "junction_end": str(end),
            "junction_strand": strand, "canonical_junction_id": canonical, "source_junction_id": normal_id,
            "junction_coordinate_system": "STAR_SJ_1BASED_INTRON_INCLUSIVE",
            "normal_cohort_junction_status": normal_status,
            "normal_cohort_normal_samples": str(normal_samples) if normal_has_frequency and normal_hit else "",
            "normal_cohort_normal_reads": str(normal_reads) if normal_has_frequency and normal_hit else "",
            "normal_cohort_normal_tissues": str(normal_tissues) if normal_has_frequency and normal_hit else "",
            "normal_junction_assessment_status": (
                "DETECTED" if normal_hit else "UNASSESSED_COVERAGE" if normal_db else "UNASSESSED"
            ),
        })
        if matched_normal_star is None:
            metrics.update({
                "matched_normal_junction_status": "UNASSESSED_NO_MATCHED_NORMAL_RNA",
                "matched_normal_junction_reads": "",
                "matched_normal_junction_coverage": "",
            })
        else:
            normal_star_row = matched_normal_star.get(key)
            normal_bam_row = matched_normal_bam_by_key.get(key, {})
            normal_unique = int(normal_bam_row.get("bam_unique_split_reads") or (normal_star_row or {}).get("star_unique_reads") or 0)
            normal_multi = int(normal_bam_row.get("bam_multi_split_reads") or (normal_star_row or {}).get("star_multi_reads") or 0)
            normal_coverage = max(
                (matched_normal_donors or {}).get((chrom, start, strand), 0),
                (matched_normal_acceptors or {}).get((chrom, end, strand), 0),
                normal_unique + normal_multi,
            )
            metrics.update({
                "matched_normal_junction_status": (
                    "DETECTED_MATCHED_NORMAL" if normal_unique + normal_multi > 0
                    else "NOT_DETECTED_ADEQUATE_COVERAGE" if normal_coverage >= 10
                    else "NOT_DETECTED_LOW_COVERAGE"
                ),
                "matched_normal_junction_reads": str(normal_unique + normal_multi),
                "matched_normal_junction_coverage": str(normal_coverage),
            })
        if star_row is None:
            metrics.update({"splice_alignment_qc_status": "FAIL", "junction_read_qc_status": "FAIL", "caller_filter_status": "FAIL", "junction_match_status": "UNRESOLVED", "junction_support_reason": "not found exactly in STAR SJ.out.tab"})
        else:
            bam_row = bam_by_key.get(key, {})
            unique_reads = int(bam_row.get("bam_unique_split_reads") or star_row["star_unique_reads"])
            multi_reads = int(bam_row.get("bam_multi_split_reads") or star_row["star_multi_reads"])
            total_reads = unique_reads + multi_reads
            denominator = max(donor_totals[(chrom, start, strand)], acceptor_totals[(chrom, end, strand)], unique_reads)
            psi = unique_reads / denominator if denominator else 0.0
            aliases = (uid, event_id, canonical, normal_id)
            caller_ok = (
                any(value in snaf_validated for value in aliases)
                and any(value in splicemutr_validated for value in aliases)
            )
            unique_starts = int(bam_row.get("unique_fragment_starts") or 0)
            overhang = max(star_row["star_max_overhang"], int(bam_row.get("bam_max_overhang") or 0))
            median_mapq = float(bam_row.get("median_mapq") or 0)
            multimap_fraction = float(bam_row.get("bam_multimapping_fraction") or (multi_reads / total_reads if total_reads else 1.0))
            qc_checks = {
                "unique_reads": unique_reads >= args.min_unique_reads,
                "unique_fragment_starts": unique_starts >= args.min_unique_fragment_starts,
                "overhang": overhang >= args.min_overhang,
                "mapq": median_mapq >= args.min_mapq,
                "multimapping_fraction": multimap_fraction <= args.max_multimapping_fraction,
                "tumor_psi": psi >= args.min_tumor_psi,
                "caller_filter": caller_ok,
            }
            read_qc_pass = all(qc_checks.values())
            failed_checks = [name for name, passed in qc_checks.items() if not passed]
            metrics.update({
                "splice_alignment_qc_status": "PASS", "junction_resolution_status": "RESOLVED",
                "junction_match_status": "EXACT", "junction_match_method": "STAR_SJ_EXACT_COORDINATES",
                "junction_support_status": "PASS" if read_qc_pass else "FAIL",
                "junction_read_qc_status": "PASS" if read_qc_pass else "FAIL",
                "junction_read_qc_failed_checks": ";".join(failed_checks),
                "unique_junction_reads": str(unique_reads),
                "unique_split_reads": str(unique_reads), "multi_split_reads": str(multi_reads),
                "junction_total_coverage": str(total_reads), "total_split_reads": str(total_reads),
                "unique_fragment_starts": bam_row.get("unique_fragment_starts", ""),
                "max_overhang": str(overhang),
                "median_mapq": bam_row.get("median_mapq", ""), "mapping_quality": bam_row.get("median_mapq", ""),
                "multimapping_fraction": bam_row.get("bam_multimapping_fraction", f"{multi_reads / total_reads:.6g}" if total_reads else ""),
                "tumor_psi": f"{psi:.6g}", "splice_psi": f"{psi:.6g}",
                "tumor_psi_method": "junction_unique_reads/max(shared_donor_unique_reads,shared_acceptor_unique_reads)",
                "caller_filter_status": "PASS" if caller_ok else "FAIL",
                "known_junction": "TRUE" if star_row["star_annotated"] == "1" else "FALSE",
                "splice_annotation_status": "ANNOTATED_NORMAL" if star_row["star_annotated"] == "1" else "UNANNOTATED",
                "annotated_normal_isoform_status": "KNOWN_NORMAL" if star_row["star_annotated"] == "1" else "NOVEL",
                "junction_source_assay_id": "tumor_short_rna_STAR_BAM",
                **bam_row,
            })
        metrics_by_event[event_id] = metrics
        qc_rows.append(metrics)
        if index % 100 == 0:
            print(f"processed {index}/{len(splice_events)} splice events", file=sys.stderr, flush=True)

    if normal_db:
        normal_db.close()
    additions = sorted({key for row in qc_rows for key in row if key != "event_id"})
    enriched_events = [{**row, **metrics_by_event.get(str(row.get("event_id") or ""), {})} for row in events]
    enriched_peptides = [{**row, **metrics_by_event.get(str(row.get("event_id") or ""), {})} for row in peptides]
    write_rows(outdir / "raw_events.enriched.tsv", enriched_events, add_fields(event_fields, additions))
    write_rows(outdir / "raw_peptides.enriched.tsv", enriched_peptides, add_fields(peptide_fields, additions))
    qc_fields = ["event_id", *additions]
    write_rows(outdir / "splice_junction_qc.enriched.tsv", qc_rows, qc_fields)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
