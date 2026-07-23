#!/usr/bin/env python3
"""Build reproducible NeoAg ligand references from an IEDB full TSV export ZIP."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import re
import sqlite3
import sys
import zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


BUILDER_VERSION = "1.0.0"
SCHEMA_VERSION = "1.0"
CANONICAL_PEPTIDE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")
DETAIL_FIELDS = (
    "peptide",
    "host",
    "disease",
    "assay_method",
    "culture_condition",
    "mhc_class",
    "mhc_restriction_name",
)
SIMPLE_FIELDS = ("peptide",)
NORMAL_DISEASE_VALUES = {"", "healthy"}
DIRECT_EX_VIVO = "Direct Ex Vivo"

REQUIRED_COLUMNS = {
    "assay_id": "Assay ID::IEDB IRI",
    "epitope_iri": "Epitope::Epitope IRI",
    "object_type": "Epitope::Object Type",
    "peptide": "Epitope::Name",
    "host": "Host::Name",
    "disease": "in vivo Process::Disease",
    "assay_method": "Assay::Method",
    "qualitative": "Assay::Qualitative Measurement",
    "culture_condition": "Antigen Presenting Cell::Culture Condition",
    "mhc_restriction_name": "MHC Restriction::Name",
    "mhc_class": "MHC Restriction::Class",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def qualified_headers(groups: list[str], fields: list[str]) -> list[str]:
    if len(groups) != len(fields):
        raise ValueError("IEDB two-row header has different group and field column counts")
    return [f"{group.strip()}::{field.strip()}" for group, field in zip(groups, fields)]


def resolve_columns(headers: list[str]) -> dict[str, int]:
    positions: dict[str, list[int]] = {}
    for pos, header in enumerate(headers):
        positions.setdefault(header, []).append(pos)
    resolved: dict[str, int] = {}
    missing: list[str] = []
    ambiguous: list[str] = []
    for logical, header in REQUIRED_COLUMNS.items():
        matches = positions.get(header, [])
        if not matches:
            missing.append(header)
        elif len(matches) > 1:
            ambiguous.append(header)
        else:
            resolved[logical] = matches[0]
    if missing or ambiguous:
        raise ValueError(f"IEDB columns invalid; missing={missing}; ambiguous={ambiguous}")
    return resolved


def choose_member(archive: zipfile.ZipFile, requested: str | None) -> zipfile.ZipInfo:
    files = [item for item in archive.infolist() if not item.is_dir()]
    if requested:
        try:
            return archive.getinfo(requested)
        except KeyError as exc:
            raise ValueError(f"ZIP member not found: {requested}") from exc
    tsv_files = [item for item in files if item.filename.lower().endswith(".tsv")]
    if len(tsv_files) != 1:
        raise ValueError("Archive must contain exactly one TSV or --member must be supplied")
    return tsv_files[0]


def normalize_row(row: list[str], index: dict[str, int]) -> dict[str, str]:
    return {name: row[pos].strip() if pos < len(row) else "" for name, pos in index.items()}


def rejection_reason(record: dict[str, str]) -> str | None:
    if "homo sapiens" not in record["host"].lower():
        return "non_human_host"
    if record["object_type"].lower() != "linear peptide":
        return "not_linear_peptide"
    peptide = record["peptide"].upper()
    if not peptide or not CANONICAL_PEPTIDE.fullmatch(peptide):
        return "noncanonical_or_empty_peptide"
    if "mass spectrometry" not in record["assay_method"].lower():
        return "not_mass_spectrometry"
    if not record["qualitative"].lower().startswith("positive"):
        return "not_positive"
    if record["mhc_class"].upper() not in {"I", "II"}:
        return "missing_or_invalid_mhc_class"
    if not record["mhc_restriction_name"]:
        return "missing_mhc_restriction"
    return None


def detail_tuple(record: dict[str, str]) -> tuple[str, ...]:
    return (
        record["peptide"].upper(),
        "Homo sapiens (human)",
        record["disease"],
        record["assay_method"],
        record["culture_condition"],
        record["mhc_class"].upper(),
        record["mhc_restriction_name"],
    )


def is_strict_normal(record: dict[str, str]) -> bool:
    return (
        record["disease"].strip().lower() in NORMAL_DISEASE_VALUES
        and record["culture_condition"].strip().lower() == DIRECT_EX_VIVO.lower()
    )


def open_database(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(path)
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.execute("PRAGMA temp_store=FILE")
    connection.execute(
        "CREATE TABLE evidence (peptide TEXT, host TEXT, disease TEXT, assay_method TEXT, "
        "culture_condition TEXT, mhc_class TEXT, mhc_restriction_name TEXT, "
        "PRIMARY KEY (peptide, host, disease, assay_method, culture_condition, mhc_class, mhc_restriction_name)) WITHOUT ROWID"
    )
    connection.execute("CREATE TABLE peptides (peptide TEXT PRIMARY KEY) WITHOUT ROWID")
    connection.execute("CREATE TABLE strict_normal (peptide TEXT PRIMARY KEY) WITHOUT ROWID")
    return connection


def write_query(connection: sqlite3.Connection, query: str, path: Path, fields: tuple[str, ...]) -> dict[str, object]:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = 0
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        for row in connection.execute(query):
            writer.writerow(row)
            rows += 1
    return {"path": path.name, "rows": rows, "size": path.stat().st_size, "sha256": sha256(path)}


def build(args: argparse.Namespace) -> dict[str, object]:
    archive_path = args.archive.resolve()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    build_dir = outdir / "build"
    temp_dir = args.temp_dir.resolve() if args.temp_dir else outdir
    temp_dir.mkdir(parents=True, exist_ok=True)
    db_path = temp_dir / f".iedb_builder.{os.getpid()}.sqlite3"
    if db_path.exists():
        db_path.unlink()
    connection = open_database(db_path)
    counts: Counter[str] = Counter()
    rejection_counts: Counter[str] = Counter()
    methods: Counter[str] = Counter()
    classes: Counter[str] = Counter()

    try:
        with zipfile.ZipFile(archive_path) as archive:
            member = choose_member(archive, args.member)
            with archive.open(member) as binary:
                text = io.TextIOWrapper(binary, encoding="utf-8-sig", newline="")
                reader = csv.reader(text, delimiter="\t")
                groups = next(reader)
                fields = next(reader)
                headers = qualified_headers(groups, fields)
                index = resolve_columns(headers)
                pending_evidence: list[tuple[str, ...]] = []
                pending_peptides: list[tuple[str]] = []
                pending_normal: list[tuple[str]] = []
                for row in reader:
                    counts["raw_data_rows"] += 1
                    record = normalize_row(row, index)
                    reason = rejection_reason(record)
                    if reason:
                        rejection_counts[reason] += 1
                        continue
                    counts["accepted_assay_rows"] += 1
                    methods[record["assay_method"]] += 1
                    classes[record["mhc_class"].upper()] += 1
                    peptide = record["peptide"].upper()
                    pending_evidence.append(detail_tuple(record))
                    pending_peptides.append((peptide,))
                    if is_strict_normal(record):
                        pending_normal.append((peptide,))
                        counts["strict_normal_assay_rows"] += 1
                    if len(pending_evidence) >= args.batch_size:
                        connection.executemany("INSERT OR IGNORE INTO evidence VALUES (?, ?, ?, ?, ?, ?, ?)", pending_evidence)
                        connection.executemany("INSERT OR IGNORE INTO peptides VALUES (?)", pending_peptides)
                        connection.executemany("INSERT OR IGNORE INTO strict_normal VALUES (?)", pending_normal)
                        connection.commit()
                        pending_evidence.clear()
                        pending_peptides.clear()
                        pending_normal.clear()
                if pending_evidence:
                    connection.executemany("INSERT OR IGNORE INTO evidence VALUES (?, ?, ?, ?, ?, ?, ?)", pending_evidence)
                    connection.executemany("INSERT OR IGNORE INTO peptides VALUES (?)", pending_peptides)
                    connection.executemany("INSERT OR IGNORE INTO strict_normal VALUES (?)", pending_normal)
                    connection.commit()

            outputs = {
                "detail": write_query(
                    connection,
                    "SELECT peptide, host, disease, assay_method, culture_condition, mhc_class, mhc_restriction_name "
                    "FROM evidence ORDER BY peptide, mhc_class, mhc_restriction_name, disease, assay_method, culture_condition",
                    build_dir / "iedb_human_ms_ligands_detail.tsv",
                    DETAIL_FIELDS,
                ),
                "all_human_ms_peptides": write_query(
                    connection,
                    "SELECT peptide FROM peptides ORDER BY peptide",
                    build_dir / "iedb_human_ms_ligands.tsv",
                    SIMPLE_FIELDS,
                ),
                "strict_normal_direct_ex_vivo": write_query(
                    connection,
                    "SELECT peptide FROM strict_normal ORDER BY peptide",
                    build_dir / "iedb_human_normal_direct_ex_vivo_ligands.tsv",
                    SIMPLE_FIELDS,
                ),
            }
            raw_record = {
                "path": archive_path.name,
                "manifest_path": f"raw/{archive_path.name}",
                "size": archive_path.stat().st_size,
                "sha256": sha256(archive_path),
                "zip_member": member.filename,
                "zip_member_uncompressed_size": member.file_size,
                "zip_member_compressed_size": member.compress_size,
                "zip_member_crc32": f"{member.CRC:08x}",
            }
    finally:
        connection.close()
        if db_path.exists() and not args.keep_database:
            db_path.unlink()

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "reference_name": "IEDB human MHC ligand mass-spectrometry reference",
        "source_database": "Immune Epitope Database and Analysis Resource (IEDB)",
        "source_export": "MHC ligand full TSV",
        "source_release": args.source_release,
        "source_url": "https://www.iedb.org/database_export_v3.php",
        "archive_local_timestamp": datetime.fromtimestamp(archive_path.stat().st_mtime, timezone.utc).isoformat(),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "builder": {"name": Path(__file__).name, "version": BUILDER_VERSION},
        "raw_archive": raw_record,
        "filters": {
            "host": "Host::Name contains Homo sapiens (case-insensitive)",
            "object_type": "Epitope::Object Type equals Linear peptide",
            "peptide": "uppercase; non-empty; only 20 canonical amino-acid letters ACDEFGHIKLMNPQRSTVWY",
            "assay_method": "Assay::Method contains mass spectrometry (case-insensitive)",
            "qualitative_measurement": "Assay::Qualitative Measurement starts with Positive (case-insensitive)",
            "mhc_class": "MHC Restriction::Class is I or II",
            "mhc_restriction": "MHC Restriction::Name is non-empty",
            "strict_normal_subset": "all preceding filters AND disease is blank or healthy AND culture condition is Direct Ex Vivo",
        },
        "deduplication": {
            "detail": "unique seven-column detail record; deterministic lexical ordering",
            "all_human_ms_peptides": "unique peptide sequence; deterministic lexical ordering",
            "strict_normal_direct_ex_vivo": "unique peptide sequence; deterministic lexical ordering",
        },
        "counts": {
            **dict(counts),
            "rejected_rows": sum(rejection_counts.values()),
            "rejection_reasons": dict(sorted(rejection_counts.items())),
            "accepted_assay_methods": dict(sorted(methods.items())),
            "accepted_mhc_classes": dict(sorted(classes.items())),
        },
        "outputs": outputs,
        "interpretation": {
            "all_human_ms": "supplemental observed human immunopeptidome evidence; not a pure normal-tissue safety reference",
            "strict_normal_direct_ex_vivo": "preferred IEDB-derived safety subset; exact allele and tissue context should still be retained during candidate review",
            "disease_or_cell_line_records": "review/caution evidence only; must not be promoted to normal-tissue exclusion evidence",
        },
    }
    manifest_path = outdir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
    checksums_path = outdir / "SHA256SUMS"
    with checksums_path.open("w", encoding="utf-8") as handle:
        handle.write(f"{raw_record['sha256']}  raw/{archive_path.name}\n")
        for output in outputs.values():
            handle.write(f"{output['sha256']}  build/{output['path']}\n")
    result = {
        "status": "PASS",
        "manifest": str(manifest_path),
        "checksums": str(checksums_path),
        "raw_sha256": raw_record["sha256"],
        "outputs": outputs,
    }
    print(json.dumps(result, indent=2))
    return result


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True, help="IEDB mhc_ligand_full_tsv.zip")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--member", help="TSV member name when the ZIP contains multiple TSV files")
    parser.add_argument("--source-release", default="export-date-unknown")
    parser.add_argument("--batch-size", type=int, default=10000)
    parser.add_argument("--temp-dir", type=Path, help="Local directory for the temporary SQLite database")
    parser.add_argument("--keep-database", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if not args.archive.is_file():
        raise FileNotFoundError(args.archive)
    build(args)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, zipfile.BadZipFile, csv.Error) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
