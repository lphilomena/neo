#!/usr/bin/env python3
"""Run NetChop 3.1d once on merged candidate peptides."""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


def read_peptides(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    by_id, sequence_ids = {}, {}
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            peptide, peptide_id = (row.get("peptide") or "").strip(), (row.get("peptide_id") or "").strip()
            if peptide and peptide_id:
                by_id[peptide_id] = peptide
                sequence_ids.setdefault(peptide, f"SEQ{len(sequence_ids):07d}")
    return by_id, sequence_ids


def parse_outputs(paths: list[Path]) -> dict[str, tuple[str, str, str, str]]:
    scores: dict[str, list[float]] = defaultdict(list)
    sites: dict[str, int] = defaultdict(int)
    cterm: dict[str, tuple[int, float]] = {}
    pattern = re.compile(r"^\s*(\d+)\s+[A-Z]\s+([CS.])\s+([0-9.]+)\s+(\S+)")
    for path in paths:
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            match = pattern.match(line)
            if not match:
                continue
            position = int(match.group(1))
            score = float(match.group(3)); ident = match.group(4)
            scores[ident].append(score)
            if position >= cterm.get(ident, (0, 0.0))[0]:
                cterm[ident] = (position, score)
            if match.group(2) == "S":
                sites[ident] += 1
    return {
        key: (
            f"{max(values):.6g}",
            f"{sum(values)/len(values):.6g}",
            f"{cterm[key][1]:.6g}",
            str(sites.get(key, 0)),
        )
        for key, values in scores.items()
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-peptides", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--binary", default=os.environ.get("NEOAG_NETCHOP_BIN", "netChop"))
    parser.add_argument("--home", default=os.environ.get("NETCHOP_HOME", ""))
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--chunk-size", type=int, default=100)
    parser.add_argument(
        "--reuse-existing-raw",
        action="store_true",
        help="reparse complete netchop.*.out batches without rerunning NetChop",
    )
    args = parser.parse_args()
    by_id, sequence_ids = read_peptides(args.raw_peptides)
    if not by_id:
        raise SystemExit("No peptide_id/peptide records found")
    work = args.output.parent / "netchop_raw"
    work.mkdir(parents=True, exist_ok=True)
    items = [(sid, sequence) for sequence, sid in sequence_ids.items()]
    chunks = [items[i:i + args.chunk_size] for i in range(0, len(items), args.chunk_size)]

    def run_chunk(index: int, chunk: list[tuple[str, str]]) -> Path:
        fasta, output = work / f"candidate.{index:05d}.fa", work / f"netchop.{index:05d}.out"
        with fasta.open("w", encoding="ascii") as handle:
            for ident, sequence in chunk:
                handle.write(f">{ident}\n{sequence}\n")
        env = os.environ.copy()
        if args.home:
            env["NETCHOP"] = str(Path(args.home) / "Linux_x86_64")
        with output.open("w", encoding="utf-8") as handle:
            subprocess.run([args.binary, "-tdir", "/tmp/netChopXXXXXX", str(fasta)], stdout=handle, stderr=subprocess.STDOUT, env=env, check=True)
        return output

    outputs = [work / f"netchop.{index:05d}.out" for index in range(len(chunks))]
    if args.reuse_existing_raw:
        missing = [path for path in outputs if not path.is_file() or path.stat().st_size == 0]
        if missing:
            raise SystemExit(f"Cannot reuse incomplete NetChop output: {missing[0]}")
    else:
        with ThreadPoolExecutor(max_workers=max(1, min(args.threads, len(chunks)))) as pool:
            outputs = list(pool.map(lambda item: run_chunk(*item), enumerate(chunks)))
    parsed = parse_outputs(outputs)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = ["peptide_id", "peptide", "netchop_31d_max_score", "netchop_31d_mean_score", "netchop_31d_cterm_score", "netchop_31d_cleavage_sites", "netchop_processing_status", "netchop_model"]
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for peptide_id, sequence in by_id.items():
            values = parsed.get(sequence_ids[sequence])
            writer.writerow({"peptide_id": peptide_id, "peptide": sequence, "netchop_31d_max_score": values[0] if values else "", "netchop_31d_mean_score": values[1] if values else "", "netchop_31d_cterm_score": values[2] if values else "", "netchop_31d_cleavage_sites": values[3] if values else "", "netchop_processing_status": "ASSESSED" if values else "UNASSESSED", "netchop_model": "NetChop 3.1d C-term"})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
