#!/usr/bin/env python3
"""Validate a tumor/normal RNA cohort and prepare LeafCutter junction inputs."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


REQUIRED = {"sample_id", "role", "star_sj"}


def read_samples(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = REQUIRED - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"sample sheet missing columns: {', '.join(sorted(missing))}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
    if not rows:
        raise ValueError("sample sheet is empty")
    return rows


def convert_star(path: Path, output: Path, min_unique_reads: int) -> tuple[int, int]:
    retained = total = 0
    output.parent.mkdir(parents=True, exist_ok=True)
    with path.open(encoding="utf-8") as source, output.open("w", encoding="utf-8") as target:
        for line in source:
            if not line.strip():
                continue
            total += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                raise ValueError(f"invalid STAR SJ.out.tab row in {path}: {line[:80]!r}")
            unique_reads = int(fields[6])
            if unique_reads < min_unique_reads or int(fields[4]) == 0:
                continue
            strand = {0: "*", 1: "+", 2: "-"}.get(int(fields[3]), "*")
            # Match the official SpliceMutr STAR_to_leaf.R coordinate convention.
            values = [fields[0], str(int(fields[1]) - 1), str(int(fields[2]) + 1),
                      f"JUNC{retained + 1}", str(unique_reads), strand]
            target.write("\t".join(values) + "\n")
            retained += 1
    return total, retained


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--min-normal-samples", type=int, default=2)
    parser.add_argument("--min-tumor-samples", type=int, default=1)
    parser.add_argument("--min-unique-reads", type=int, default=10)
    parser.add_argument("--allow-low-power", action="store_true")
    args = parser.parse_args()

    rows = read_samples(args.samples)
    sample_ids = [row["sample_id"] for row in rows]
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("sample_id values must be unique")
    for row in rows:
        role = row["role"].lower()
        if role not in {"normal", "tumor"}:
            raise ValueError(f"invalid role for {row['sample_id']}: {row['role']}")
        row["role"] = role
        source = Path(row["star_sj"]).expanduser().resolve()
        if not source.is_file():
            raise FileNotFoundError(f"STAR junction file not found: {source}")
        row["star_sj"] = str(source)

    normals = [row for row in rows if row["role"] == "normal"]
    tumors = [row for row in rows if row["role"] == "tumor"]
    low_power = len(normals) < args.min_normal_samples or len(tumors) < args.min_tumor_samples
    if low_power and not args.allow_low_power:
        raise ValueError(
            f"cohort requires >= {args.min_normal_samples} normal and >= "
            f"{args.min_tumor_samples} tumor samples; found {len(normals)} and {len(tumors)}"
        )

    outdir = args.outdir.resolve()
    junction_dir = outdir / "junctions"
    junction_dir.mkdir(parents=True, exist_ok=True)
    normalized: list[dict[str, str]] = []
    for row in normals + tumors:  # LeafCutter reference group must be first.
        junc = junction_dir / f"{row['sample_id']}.junc"
        total, retained = convert_star(Path(row["star_sj"]), junc, args.min_unique_reads)
        normalized.append({**row, "junc_file": str(junc), "junction_rows": str(total),
                           "retained_junction_rows": str(retained)})

    fields = ["sample_id", "role", "star_sj", "junc_file", "junction_rows",
              "retained_junction_rows"]
    with (outdir / "cohort_samples.normalized.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader(); writer.writerows(normalized)
    with (outdir / "groups_file.txt").open("w", encoding="utf-8") as handle:
        for row in normalized:
            handle.write(f"{Path(row['junc_file']).name}\t{1 if row['role'] == 'normal' else 2}\n")
    with (outdir / "juncfiles.txt").open("w", encoding="utf-8") as handle:
        handle.writelines(f"{row['junc_file']}\n" for row in normalized)
    status = "LOW_POWER" if low_power else "READY"
    summary = {"status": status, "normal_samples": len(normals), "tumor_samples": len(tumors),
               "normal_is_reference_group": True, "missing_normal_is_not_negative": True,
               "min_unique_reads": args.min_unique_reads, "sample_sheet": str(args.samples.resolve())}
    (outdir / "cohort_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
