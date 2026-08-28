#!/usr/bin/env python3
"""Build a SQLite membership index for normal splice junction panels."""

from __future__ import annotations

import argparse
import csv
import gzip
import sqlite3
from pathlib import Path
from typing import TextIO


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Normal junction TSV/TSV.GZ with a junction_id column")
    parser.add_argument("--output", required=True, help="SQLite index path")
    parser.add_argument("--force", action="store_true", help="Overwrite an existing index")
    parser.add_argument("--check", action="store_true", help="Exit 0 only when an existing index uses the current schema")
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)
    schema_version = "2"
    if args.check:
        if not output_path.is_file():
            return 1
        try:
            with sqlite3.connect(output_path) as conn:
                value = conn.execute("select value from meta where key='schema_version'").fetchone()
                columns = {row[1] for row in conn.execute("pragma table_info(junction_ids)")}
            return 0 if value and value[0] == schema_version and "normal_samples" in columns else 1
        except sqlite3.Error:
            return 1
    if not input_path.is_file():
        raise SystemExit(f"input missing: {input_path}")
    if output_path.exists() and not args.force:
        raise SystemExit(f"output exists, use --force to overwrite: {output_path}")

    tmp_path = output_path.with_suffix(output_path.suffix + ".tmp")
    if tmp_path.exists():
        tmp_path.unlink()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows = 0
    unique_rows = 0
    with sqlite3.connect(tmp_path) as conn:
        conn.execute("pragma journal_mode=off")
        conn.execute("pragma synchronous=off")
        conn.execute("pragma temp_store=memory")
        conn.execute(
            "create table junction_ids ("
            "junction_id text primary key, normal_samples integer, normal_reads integer, "
            "normal_total_reads integer, normal_tissues integer, tissue text, source text, dataset text"
            ") without rowid"
        )
        conn.execute("create table meta (key text primary key, value text)")

        batch: list[tuple[str, int, int, int, int, str, str, str]] = []
        with open_text(input_path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if "junction_id" not in (reader.fieldnames or []):
                raise SystemExit(f"input lacks junction_id column: {input_path}")
            for row in reader:
                junction_id = (row.get("junction_id") or "").strip()
                if not junction_id:
                    continue
                rows += 1
                def integer(name: str) -> int:
                    try:
                        return int(float(str(row.get(name) or 0)))
                    except ValueError:
                        return 0
                batch.append((
                    junction_id, integer("normal_samples"), integer("normal_reads"),
                    integer("normal_total_reads"), integer("normal_tissues"),
                    str(row.get("tissue") or ""), str(row.get("source") or ""),
                    str(row.get("dataset") or ""),
                ))
                if len(batch) >= 100_000:
                    conn.executemany(
                        "insert into junction_ids values (?,?,?,?,?,?,?,?) on conflict(junction_id) do update set "
                        "normal_samples=max(normal_samples,excluded.normal_samples), "
                        "normal_reads=max(normal_reads,excluded.normal_reads), "
                        "normal_total_reads=max(normal_total_reads,excluded.normal_total_reads), "
                        "normal_tissues=max(normal_tissues,excluded.normal_tissues)",
                        batch,
                    )
                    batch.clear()
        if batch:
            conn.executemany(
                "insert into junction_ids values (?,?,?,?,?,?,?,?) on conflict(junction_id) do update set "
                "normal_samples=max(normal_samples,excluded.normal_samples), "
                "normal_reads=max(normal_reads,excluded.normal_reads), "
                "normal_total_reads=max(normal_total_reads,excluded.normal_total_reads), "
                "normal_tissues=max(normal_tissues,excluded.normal_tissues)",
                batch,
            )
        unique_rows = int(conn.execute("select count(*) from junction_ids").fetchone()[0])
        conn.executemany(
            "insert into meta values (?, ?)",
            [
                ("source", str(input_path)),
                ("resolvable_rows", str(rows)),
                ("unique_junction_ids", str(unique_rows)),
                ("schema_version", schema_version),
            ],
        )
        conn.execute("pragma optimize")

    tmp_path.replace(output_path)
    print(f"[OK] indexed {unique_rows} unique junctions from {rows} rows -> {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
