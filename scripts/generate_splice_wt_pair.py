#!/usr/bin/env python3
"""Reconstruct a reference transcript protein and pair a splice-junction peptide with WT."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

import pysam
from Bio.Seq import Seq


ATTR_RE = re.compile(r'(\w+) "([^"]+)"')


def attributes(text: str) -> dict[str, str]:
    return dict(ATTR_RE.findall(text))


def transcript_cds(gtf: Path, transcript_id: str) -> tuple[list[dict[str, object]], dict[str, str]]:
    target = transcript_id.split(".")[0]
    segments: list[dict[str, object]] = []
    metadata: dict[str, str] = {}
    with gtf.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or f'transcript_id "{target}"' not in line:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            attrs = attributes(fields[8])
            if fields[2] == "transcript":
                metadata = attrs
            if fields[2] == "CDS":
                segments.append(
                    {
                        "chrom": fields[0],
                        "start": int(fields[3]),
                        "end": int(fields[4]),
                        "strand": fields[6],
                        "phase": fields[7],
                        "exon_number": attrs.get("exon_number", ""),
                    }
                )
    if not segments:
        raise ValueError(f"No CDS records found for {target}")
    strand = str(segments[0]["strand"])
    segments.sort(key=lambda row: int(row["start"]), reverse=strand == "-")
    return segments, metadata


def reconstruct_cds(fasta: Path, segments: list[dict[str, object]]) -> str:
    reference = pysam.FastaFile(str(fasta))
    parts: list[str] = []
    for segment in segments:
        sequence = reference.fetch(
            str(segment["chrom"]), int(segment["start"]) - 1, int(segment["end"])
        ).upper()
        if segment["strand"] == "-":
            sequence = str(Seq(sequence).reverse_complement())
        parts.append(sequence)
    return "".join(parts)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--transcript-id", required=True)
    parser.add_argument("--left-nt", required=True)
    parser.add_argument("--right-nt", required=True)
    parser.add_argument("--mutant-peptide", required=True)
    parser.add_argument("--hla", required=True)
    parser.add_argument("--event-id", required=True)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    segments, metadata = transcript_cds(args.gtf, args.transcript_id)
    cds = reconstruct_cds(args.reference_fasta, segments)
    normal_protein = str(Seq(cds).translate(to_stop=True))
    abnormal_nt = (args.left_nt + args.right_nt).upper()
    abnormal_protein = str(Seq(abnormal_nt[: len(abnormal_nt) // 3 * 3]).translate())
    mt = args.mutant_peptide.upper()
    start = abnormal_protein.find(mt)
    if start < 0:
        raise ValueError(f"Mutant peptide {mt} not found in reconstructed abnormal protein")
    junction_aa = len(args.left_nt) // 3
    if not (start <= junction_aa < start + len(mt)):
        raise ValueError("Mutant peptide does not span the supplied junction")
    wt = normal_protein[start : start + len(mt)]
    if len(wt) != len(mt):
        raise ValueError("Reference protein is too short for an equal-length WT peptide")

    args.outdir.mkdir(parents=True, exist_ok=True)
    pair_path = args.outdir / "splice_mt_wt_pair.tsv"
    fields = [
        "event_id", "gene", "transcript_id", "transcript_name", "hla_allele", "peptide_type",
        "peptide", "peptide_length", "protein_start_1based", "junction_aa_1based",
        "junction_position_in_peptide_1based", "protein_context", "derivation",
    ]
    rows = []
    for peptide_type, peptide, protein, derivation in (
        ("MT", mt, abnormal_protein, "SNAF abnormal junction sequence"),
        ("WT", wt, normal_protein, "GENCODE reference transcript CDS"),
    ):
        rows.append(
            {
                "event_id": args.event_id,
                "gene": args.gene,
                "transcript_id": args.transcript_id.split(".")[0],
                "transcript_name": metadata.get("transcript_name", ""),
                "hla_allele": args.hla,
                "peptide_type": peptide_type,
                "peptide": peptide,
                "peptide_length": len(peptide),
                "protein_start_1based": start + 1,
                "junction_aa_1based": junction_aa + 1,
                "junction_position_in_peptide_1based": junction_aa - start + 1,
                "protein_context": protein[max(0, start - 8) : start + len(peptide) + 8],
                "derivation": derivation,
            }
        )
    with pair_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    predictor_path = args.outdir / "splice_mt_wt_predictor_input.tsv"
    with predictor_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["peptide", "hla_allele", "gene", "event_id", "peptide_type"],
            delimiter="\t",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "peptide": row["peptide"],
                    "hla_allele": args.hla,
                    "gene": args.gene,
                    "event_id": args.event_id,
                    "peptide_type": row["peptide_type"],
                }
            )

    manifest = {
        "event_id": args.event_id,
        "gene": args.gene,
        "transcript_id": args.transcript_id.split(".")[0],
        "transcript_version": metadata.get("transcript_version", ""),
        "transcript_name": metadata.get("transcript_name", ""),
        "protein_id": next(
            (
                attributes(line.split("\t", 8)[8]).get("protein_id", "")
                for line in args.gtf.open(encoding="utf-8")
                if not line.startswith("#")
                and "\tCDS\t" in line
                and f'transcript_id "{args.transcript_id.split(".")[0]}"' in line
            ),
            "",
        ),
        "cds_length": len(cds),
        "normal_protein_length": len(normal_protein),
        "abnormal_sequence_length": len(abnormal_nt),
        "abnormal_protein": abnormal_protein,
        "mutant_peptide": mt,
        "wildtype_peptide": wt,
        "junction_position_in_peptide_1based": junction_aa - start + 1,
        "outputs": {"pair_table": str(pair_path), "predictor_input": str(predictor_path)},
    }
    (args.outdir / "splice_mt_wt_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
