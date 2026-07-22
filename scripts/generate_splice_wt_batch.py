#!/usr/bin/env python3
"""Batch reconstruct WT controls for selected SNAF splice peptide-HLA rows."""

from __future__ import annotations

import argparse
import csv
import json
import pickle
import re
from collections import defaultdict
from pathlib import Path

import pysam
from Bio.Seq import Seq


ATTR_RE = re.compile(r'(\w+) "([^"]+)"')
TRANSCRIPT_RE = re.compile(r"ENST\d+(?:\.\d+)?")
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


class SnafPlaceholder:
    def __new__(cls, *args: object, **kwargs: object) -> "SnafPlaceholder":
        return object.__new__(cls)


class SnafUnpickler(pickle.Unpickler):
    def find_class(self, module: str, name: str) -> object:
        if module == "snaf.snaf":
            return SnafPlaceholder
        return super().find_class(module, name)


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def attrs(text: str) -> dict[str, str]:
    return dict(ATTR_RE.findall(text))


def transcript_candidates(row: dict[str, str]) -> list[str]:
    candidates = []
    for value in [row.get("expressed_transcript_id", ""), row.get("supporting_transcript_ids", "")]:
        candidates.extend(match.split(".")[0] for match in TRANSCRIPT_RE.findall(value))
    return list(dict.fromkeys(candidates))


def load_cds_segments(gtf: Path, targets: set[str]) -> tuple[dict[str, list[dict[str, object]]], dict[str, dict[str, str]]]:
    segments: dict[str, list[dict[str, object]]] = defaultdict(list)
    metadata: dict[str, dict[str, str]] = {}
    with gtf.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"transcript", "CDS"}:
                continue
            attributes = attrs(fields[8])
            transcript = attributes.get("transcript_id", "").split(".")[0]
            if transcript not in targets:
                continue
            if fields[2] == "transcript":
                metadata[transcript] = attributes
            else:
                segments[transcript].append(
                    {
                        "chrom": fields[0],
                        "start": int(fields[3]),
                        "end": int(fields[4]),
                        "strand": fields[6],
                        "protein_id": attributes.get("protein_id", ""),
                    }
                )
    for transcript, rows in segments.items():
        rows.sort(key=lambda item: int(item["start"]), reverse=rows[0]["strand"] == "-")
    return dict(segments), metadata


def reconstruct_cds(reference: pysam.FastaFile, segments: list[dict[str, object]]) -> str:
    parts: list[str] = []
    for segment in segments:
        sequence = reference.fetch(
            str(segment["chrom"]), int(segment["start"]) - 1, int(segment["end"])
        ).upper()
        if segment["strand"] == "-":
            sequence = str(Seq(sequence).reverse_complement())
        parts.append(sequence)
    return "".join(parts)


def anchor_local_start(left: str, right: str, cds: str) -> tuple[int | None, str, int]:
    for length in range(min(len(left), 180), 17, -1):
        query = left[-length:]
        hits = [match.start() for match in re.finditer(f"(?={re.escape(query)})", cds)]
        if len(hits) == 1:
            return hits[0] - (len(left) - length), "LEFT_SUFFIX", length
    for length in range(min(len(right), 180), 17, -1):
        query = right[:length]
        hits = [match.start() for match in re.finditer(f"(?={re.escape(query)})", cds)]
        if len(hits) == 1:
            return hits[0] - len(left), "RIGHT_PREFIX", length
    return None, "NO_UNIQUE_CDS_ANCHOR", 0


def reconstruct_pair(
    row: dict[str, str],
    junction: str,
    cds_by_transcript: dict[str, str],
    metadata: dict[str, dict[str, str]],
    protein_by_transcript: dict[str, str],
) -> dict[str, object]:
    result: dict[str, object] = {
        "event_id": f"SNAF_SPLICE|{row['uid']}",
        "uid": row["uid"],
        "gene": row.get("gene", ""),
        "hla_allele": row["hla"],
        "mutant_peptide": row["peptide"].upper(),
        "wildtype_peptide": "",
        "transcript_id": "",
        "transcript_name": "",
        "protein_id": "",
        "status": "WT_NOT_RECONSTRUCTABLE",
        "reason": "",
        "anchor_type": "",
        "anchor_length_nt": "",
        "local_translation_frame": "",
        "protein_start_1based": "",
        "junction_position_in_peptide_1based": "",
        "mutant_context": "",
        "wildtype_context": "",
    }
    if "," not in junction:
        result["reason"] = "INVALID_SNAF_JUNCTION_SEQUENCE"
        return result
    left, right = (part.upper() for part in junction.split(",", 1))
    abnormal = left + right
    mt = str(result["mutant_peptide"])
    reasons: list[str] = []
    for transcript in transcript_candidates(row):
        cds = cds_by_transcript.get(transcript)
        protein = protein_by_transcript.get(transcript)
        if not cds or not protein:
            reasons.append(f"{transcript}:NO_CDS")
            continue
        local_start, anchor_type, anchor_length = anchor_local_start(left, right, cds)
        if local_start is None:
            reasons.append(f"{transcript}:NO_UNIQUE_CDS_ANCHOR")
            continue
        frame = (-local_start) % 3
        translated = str(Seq(abnormal[frame : len(abnormal) - ((len(abnormal) - frame) % 3)]).translate())
        occurrences = [match.start() for match in re.finditer(f"(?={re.escape(mt)})", translated)]
        spanning: list[tuple[int, int]] = []
        for aa_start in occurrences:
            nt_start = frame + aa_start * 3
            nt_end = nt_start + len(mt) * 3
            if nt_start < len(left) < nt_end:
                spanning.append((aa_start, nt_start))
        if len(spanning) != 1:
            reasons.append(f"{transcript}:MT_SPANNING_OCCURRENCES_{len(spanning)}")
            continue
        aa_start, nt_start = spanning[0]
        global_nt_start = local_start + nt_start
        if global_nt_start < 0 or global_nt_start % 3:
            reasons.append(f"{transcript}:CDS_FRAME_MISMATCH")
            continue
        global_aa_start = global_nt_start // 3
        wt = protein[global_aa_start : global_aa_start + len(mt)]
        if len(wt) != len(mt) or not set(wt) <= VALID_AA:
            reasons.append(f"{transcript}:INVALID_WT_WINDOW")
            continue
        junction_aa_position = (len(left) - nt_start + 2) // 3
        result.update(
            {
                "wildtype_peptide": wt,
                "transcript_id": transcript,
                "transcript_name": metadata.get(transcript, {}).get("transcript_name", ""),
                "protein_id": metadata.get(transcript, {}).get("protein_id", ""),
                "status": "PASS" if wt != mt else "NON_MUTANT_SEQUENCE",
                "reason": "WT_RECONSTRUCTED_FROM_TRANSCRIPT_CDS" if wt != mt else "MT_EQUALS_WT",
                "anchor_type": anchor_type,
                "anchor_length_nt": anchor_length,
                "local_translation_frame": frame,
                "protein_start_1based": global_aa_start + 1,
                "junction_position_in_peptide_1based": junction_aa_position,
                "mutant_context": translated[max(0, aa_start - 8) : aa_start + len(mt) + 8],
                "wildtype_context": protein[max(0, global_aa_start - 8) : global_aa_start + len(mt) + 8],
            }
        )
        return result
    result["reason"] = ";".join(reasons[:8]) or "NO_SUPPORTING_TRANSCRIPT"
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selected", required=True, type=Path)
    parser.add_argument("--snaf-pickle", required=True, type=Path)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    _, selected = read_tsv(args.selected)
    targets = {transcript for row in selected for transcript in transcript_candidates(row)}
    segments, metadata = load_cds_segments(args.gtf, targets)
    reference = pysam.FastaFile(str(args.reference_fasta))
    cds_by_transcript = {
        transcript: reconstruct_cds(reference, rows) for transcript, rows in segments.items()
    }
    protein_by_transcript = {
        transcript: str(Seq(cds).translate(to_stop=True))
        for transcript, cds in cds_by_transcript.items()
    }

    with args.snaf_pickle.open("rb") as handle:
        query = SnafUnpickler(handle).load()
    junction_index = {
        str(obj.__dict__.get("uid", "")): str(obj.__dict__.get("junction", ""))
        for obj in query.__dict__.get("translated", [])
    }

    pairs = [
        reconstruct_pair(
            row,
            junction_index.get(row["uid"], ""),
            cds_by_transcript,
            metadata,
            protein_by_transcript,
        )
        for row in selected
    ]
    args.outdir.mkdir(parents=True, exist_ok=True)
    fields = list(pairs[0]) if pairs else []
    pair_path = args.outdir / "splice_mt_wt_pairs.tsv"
    write_tsv(pair_path, fields, pairs)

    predictor_rows: list[dict[str, object]] = []
    for pair in pairs:
        if pair["status"] != "PASS":
            continue
        for peptide_type, peptide in (
            ("MT", pair["mutant_peptide"]),
            ("WT", pair["wildtype_peptide"]),
        ):
            predictor_rows.append(
                {
                    "peptide": peptide,
                    "hla_allele": pair["hla_allele"],
                    "gene": pair["gene"],
                    "event_id": pair["event_id"],
                    "peptide_type": peptide_type,
                }
            )
    predictor_path = args.outdir / "splice_mt_wt_predictor_input.tsv"
    write_tsv(
        predictor_path,
        ["peptide", "hla_allele", "gene", "event_id", "peptide_type"],
        predictor_rows,
    )

    counts: dict[str, int] = defaultdict(int)
    for pair in pairs:
        counts[str(pair["status"])] += 1
    summary = {
        "status": "PASS",
        "selected_rows": len(selected),
        "transcript_targets": len(targets),
        "transcripts_with_cds": len(cds_by_transcript),
        "reconstruction_status_counts": dict(sorted(counts.items())),
        "prediction_rows": len(predictor_rows),
        "pair_table": str(pair_path),
        "predictor_input": str(predictor_path),
    }
    (args.outdir / "splice_mt_wt_batch_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
