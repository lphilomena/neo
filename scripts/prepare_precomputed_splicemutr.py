#!/usr/bin/env python3
"""Convert precomputed SpliceMutr transcript output to splice evidence rows.

SpliceMutr's ``peptide`` column contains a translated protein, not a presented
short peptide. This converter deliberately imports transcript reconstruction
as event-level evidence and leaves peptide/HLA fields empty unless a separate,
auditable junction-crossing peptide mapping is available.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


FIELDS = [
    "sample_id",
    "event_id",
    "source_junction_id",
    "gene",
    "chrom",
    "start",
    "end",
    "strand",
    "transcript_id",
    "peptide",
    "hla_allele",
    "binding_rank",
    "crosses_junction",
    "source_tool",
    "evidence_status",
    "transcript_reconstruction_status",
    "presentation_status",
    "translated_protein_length",
    "translated_protein_sha256",
    "source_record_id",
]


def _value(row: dict[str, str], *names: str) -> str:
    for name in names:
        value = str(row.get(name) or "").strip()
        if value:
            return value
    return ""


def convert(source: Path, output: Path, sample_id: str, presentation: Path | None) -> int:
    presentation_status = (
        "ASSESSED_GLOBAL_NOT_JUNCTION_PEPTIDE_MAPPED"
        if presentation and presentation.is_file() and presentation.stat().st_size > 0
        else "UNASSESSED"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with source.open(encoding="utf-8", errors="replace", newline="") as handle, output.open(
        "w", encoding="utf-8", newline=""
    ) as target:
        reader = csv.DictReader(handle, delimiter="\t")
        writer = csv.DictWriter(target, delimiter="\t", fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        for index, row in enumerate(reader, 1):
            chrom = _value(row, "chr", "chrom", "chromosome")
            start = _value(row, "start", "junction_start")
            end = _value(row, "end", "junction_end")
            strand = _value(row, "strand") or "."
            junction_id = _value(row, "juncs", "junction_id", "event_id")
            if not junction_id and chrom and start and end:
                junction_id = f"{chrom}:{start}:{end}:{strand}"
            if not junction_id:
                continue
            protein = _value(row, "peptide", "protein", "translated_protein").rstrip("*")
            digest = hashlib.sha256(protein.encode("utf-8")).hexdigest() if protein else ""
            writer.writerow(
                {
                    "sample_id": sample_id,
                    "event_id": junction_id,
                    "source_junction_id": junction_id,
                    "gene": _value(row, "gene", "gene_name"),
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "transcript_id": _value(row, "tx_id", "transcript_id"),
                    "peptide": "",
                    "hla_allele": "",
                    "binding_rank": "",
                    "crosses_junction": "UNASSESSED",
                    "source_tool": "SpliceMutr",
                    "evidence_status": "TRANSCRIPT_RECONSTRUCTION_ONLY",
                    "transcript_reconstruction_status": "ASSESSED",
                    "presentation_status": presentation_status,
                    "translated_protein_length": str(len(protein)),
                    "translated_protein_sha256": digest,
                    "source_record_id": f"SpliceMutr:{index}",
                }
            )
            count += 1
    if count == 0:
        raise ValueError(f"No resolvable SpliceMutr rows found in {source}")
    return count


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--presentation", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    count = convert(args.input, args.output, args.sample_id, args.presentation)
    print(f"wrote {count} SpliceMutr event-evidence rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
