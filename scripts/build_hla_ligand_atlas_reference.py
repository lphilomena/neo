#!/usr/bin/env python3
"""Build the NeoAg benign-tissue ligand reference from HLA Ligand Atlas tables."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import shutil
import time
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


BUILDER_VERSION = "1.0.0"
REQUIRED_TABLES = ("peptides", "sample_hits")
OUTPUT_FIELDS = ("peptide", "source_tissue", "hla_allele", "hla_class", "source", "note")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, destination: Path, retries: int = 8) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".part")
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        offset = temporary.stat().st_size if temporary.exists() else 0
        headers = {"User-Agent": "neoag-reference-builder/1.0"}
        if offset:
            headers["Range"] = f"bytes={offset}-"
        request = urllib.request.Request(url, headers=headers)
        try:
            with urllib.request.urlopen(request, timeout=600) as response:
                resumed = offset > 0 and getattr(response, "status", None) == 206
                mode = "ab" if resumed else "wb"
                with temporary.open(mode) as out:
                    shutil.copyfileobj(response, out, length=1024 * 1024)
            temporary.replace(destination)
            return
        except Exception as exc:  # network failures are retried without hiding the final error
            last_error = exc
            if attempt < retries:
                time.sleep(min(5 * attempt, 30))
    raise RuntimeError(f"Failed to download {url} after {retries} attempts") from last_error


def ensure_raw_tables(raw_dir: Path, release: str, base_url: str, allow_download: bool) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for table in REQUIRED_TABLES:
        filename = f"{table}.tsv.gz"
        path = raw_dir / filename
        url = f"{base_url.rstrip('/')}/{release}/{filename}"
        if not path.is_file():
            if not allow_download:
                raise FileNotFoundError(f"Missing {path}; rerun with --download")
            download(url, path)
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            header = next(csv.reader(handle, delimiter="\t"), [])
        records.append({
            "table": table,
            "filename": filename,
            "url": url,
            "path": f"raw/{filename}",
            "size": path.stat().st_size,
            "sha256": sha256(path),
            "header": header,
        })
    return records


def build(peptides_path: Path, hits_path: Path, output: Path, release: str) -> dict[str, object]:
    peptide_order: list[tuple[str, str]] = []
    with gzip.open(peptides_path, "rt", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            peptide_id = str(row.get("peptide_sequence_id") or "").strip()
            sequence = str(row.get("peptide_sequence") or "").strip().upper()
            if peptide_id and sequence:
                peptide_order.append((peptide_id, sequence))

    tissues: dict[str, set[str]] = defaultdict(set)
    classes: dict[str, set[str]] = defaultdict(set)
    with gzip.open(hits_path, "rt", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            peptide_id = str(row.get("peptide_sequence_id") or "").strip()
            tissue = str(row.get("tissue") or "").strip()
            hla_class = str(row.get("hla_class") or "").strip()
            if peptide_id and tissue:
                tissues[peptide_id].add(tissue)
            if peptide_id and hla_class:
                classes[peptide_id].add(hla_class)

    output.parent.mkdir(parents=True, exist_ok=True)
    class_counts: Counter[str] = Counter()
    length_counts: Counter[int] = Counter()
    rows_written = 0
    with output.open("w", encoding="utf-8", newline="") as handle:
        # Keep CRLF for byte-level compatibility with the canonical 2020.12 build.
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\r\n")
        writer.writeheader()
        for peptide_id, sequence in peptide_order:
            observed_classes = classes.get(peptide_id, set())
            if observed_classes == {"HLA-I", "HLA-II"}:
                hla_class = "HLA-I+II"
            elif len(observed_classes) == 1:
                hla_class = next(iter(observed_classes))
            else:
                hla_class = "+".join(sorted(observed_classes))
            writer.writerow({
                "peptide": sequence,
                "source_tissue": ",".join(sorted(tissues.get(peptide_id, set()))),
                "hla_allele": "",
                "hla_class": hla_class,
                "source": f"HLA_Ligand_Atlas_{release}",
                "note": "peptide-level normal ligand reference; no direct HLA restriction claim",
            })
            rows_written += 1
            class_counts[hla_class] += 1
            length_counts[len(sequence)] += 1
    return {
        "rows": rows_written,
        "class_counts": dict(sorted(class_counts.items())),
        "peptide_length_counts": {str(key): value for key, value in sorted(length_counts.items())},
        "size": output.stat().st_size,
        "sha256": sha256(output),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--release", default="2020.12")
    parser.add_argument("--base-url", default="http://hla-ligand-atlas.org/rel")
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--checksums", type=Path, required=True)
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--compare-to", type=Path)
    args = parser.parse_args()

    raw_records = ensure_raw_tables(args.raw_dir, args.release, args.base_url, args.download)
    raw_by_name = {record["table"]: args.raw_dir / str(record["filename"]) for record in raw_records}
    output_stats = build(raw_by_name["peptides"], raw_by_name["sample_hits"], args.output, args.release)

    comparison: dict[str, object] = {"requested": False}
    if args.compare_to:
        comparison = {
            "requested": True,
            "path": str(args.compare_to.resolve()),
            "exists": args.compare_to.is_file(),
        }
        if args.compare_to.is_file():
            expected_hash = sha256(args.compare_to)
            comparison.update({
                "size": args.compare_to.stat().st_size,
                "sha256": expected_hash,
                "byte_identical": expected_hash == output_stats["sha256"],
            })

    args.checksums.parent.mkdir(parents=True, exist_ok=True)
    with args.checksums.open("w", encoding="utf-8") as handle:
        for record in raw_records:
            handle.write(f"{record['sha256']}  raw/{record['filename']}\n")
        handle.write(f"{output_stats['sha256']}  build/{args.output.name}\n")

    manifest = {
        "schema_version": "1.0",
        "reference_name": "NeoAg normal HLA ligandome",
        "source_database": "HLA Ligand Atlas",
        "source_release": args.release,
        "source_url": f"https://hla-ligand-atlas.org/rel/{args.release}/",
        "license": "CC-BY-4.0",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "builder": {
            "name": Path(__file__).name,
            "version": BUILDER_VERSION,
            "required_tables": list(REQUIRED_TABLES),
        },
        "method": {
            "join_key": "peptide_sequence_id",
            "one_row_per": "unique peptide sequence from peptides.tsv.gz",
            "tissue_aggregation": "sorted unique sample_hits.tissue values",
            "class_aggregation": "HLA-I, HLA-II, or HLA-I+II",
            "hla_allele": "left blank because sample-level donor alleles do not establish peptide-specific restriction",
        },
        "raw_files": raw_records,
        "output": {"path": f"build/{args.output.name}", **output_stats},
        "canonical_comparison": comparison,
        "checksums_file": args.checksums.name,
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(manifest, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "output": output_stats, "comparison": comparison}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
