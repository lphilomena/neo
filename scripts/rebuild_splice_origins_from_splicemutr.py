#!/usr/bin/env python3
"""Build formal transcript/ORF/peptide-origin chains from SpliceMutr output.

Matching is deliberately strict: exact canonical junction, exact peptide
sequence in the translated protein, and an auditable junction boundary.  No
gene-only, nearest-coordinate, or caller-ID fallback is used.
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
from pathlib import Path

from neoag.splice.coordinates import CanonicalJunction, normalize_chromosome
from neoag.splice.identifiers import (
    link_id, orf_id, peptide_id, peptide_origin_id, sequence_sha256,
    splice_event_id, transcript_hypothesis_id,
)
from neoag.splice.schemas import (
    ORF_FIELDS, PEPTIDE_ORIGIN_FIELDS, PEPTIDE_ORIGIN_LINK_FIELDS,
    TRANSCRIPT_HYPOTHESIS_FIELDS,
)
from neoag.schemas import PEPTIDE_FIELDS
from neoag.utils import write_tsv


def _int(value: str) -> int | None:
    try:
        return int(float(str(value).strip()))
    except (TypeError, ValueError):
        return None


def _boundaries(value: str) -> list[int]:
    return [x for token in str(value or "").replace(";", ",").split(",") if (x := _int(token)) is not None]


def _protein(value: str) -> str:
    return str(value or "").strip().upper().rstrip("*")


def _event_map(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            event = row.get("splice_event_id", "")
            for junction in row.get("junction_ids", "").split(";"):
                if junction and junction not in result:
                    result[junction] = event
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--genome-build", default="GRCh38")
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--formal-events", type=Path)
    parser.add_argument("--splicemutr-glob", required=True)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--tool-version", default="UNASSESSED")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    event_by_junction = _event_map(args.formal_events) if args.formal_events else {}
    candidates: dict[str, list[dict[str, str]]] = {}
    with args.candidates.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            junction = row.get("canonical_junction_id") or row.get("event_id", "")
            peptide = row.get("peptide", "").upper()
            if junction.startswith("SJ|") and peptide:
                candidates.setdefault(junction, []).append(row)

    transcripts: dict[str, dict[str, str]] = {}
    orfs: dict[str, dict[str, str]] = {}
    origins: dict[str, dict[str, str]] = {}
    links: dict[str, dict[str, str]] = {}
    matched_candidates: set[tuple[str, str]] = set()
    best_origin: dict[tuple[str, str], dict[str, str]] = {}
    translated_rows = 0
    files = sorted(glob.glob(args.splicemutr_glob, recursive=True))

    for source in files:
        with open(source, newline="") as handle:
            for row_no, row in enumerate(csv.DictReader(handle, delimiter="\t"), start=2):
                start, end = _int(row.get("start", "")), _int(row.get("end", ""))
                strand = row.get("strand", "")
                if start is None or end is None or strand not in {"+", "-"}:
                    continue
                # SpliceMutr/SNAF corrected tables store the flanking exon
                # boundary coordinates. Convert them to the first/last
                # intronic base used by the canonical 1-based closed model.
                junction = CanonicalJunction(
                    args.genome_build, normalize_chromosome(row.get("chr", "")), start + 1, end - 1, strand
                ).junction_id
                relevant = candidates.get(junction, [])
                protein = _protein(row.get("peptide", ""))
                if not relevant or not protein:
                    continue
                event_id = event_by_junction.get(junction) or splice_event_id(
                    genome_build=args.genome_build, event_type="NOVEL_JUNCTION", strand=strand,
                    junction_ids=[junction], gene=next((x.get("gene", "") for x in relevant if x.get("gene")), row.get("gene", "")),
                )
                translated_rows += 1
                transcript_id = row.get("tx_id", "")
                gene_id = row.get("gene", "")
                gene = next((x.get("gene", "") for x in relevant if x.get("gene")), gene_id)
                modified = row.get("modified", "").upper() or "UNASSESSED"
                boundaries = _boundaries(row.get("pep_junc_loc", ""))
                source_record = f"SpliceMutr|{Path(source).name}|{row_no}"
                sth = transcript_hypothesis_id(
                    splice_event_id_value=event_id, reference_transcript_id=transcript_id,
                    junction_chain=[junction], path_role=f"SPLICEMUTR_{modified}",
                    sequence_sha256=sequence_sha256(protein), source_generator="SpliceMutr",
                )
                transcripts[sth] = {
                    "transcript_hypothesis_id": sth, "splice_event_id": event_id,
                    "sample_id": args.sample_id, "gene": gene, "gene_id": gene_id,
                    "reference_transcript_id": transcript_id, "mane_status": "UNASSESSED",
                    "path_id": row.get("cluster", ""), "path_role": f"SPLICEMUTR_{modified}",
                    "exon_chain": f"{row.get('start_exon','')};{row.get('end_exon','')}",
                    "junction_chain": junction, "cds_start": row.get("start_stop", "").split(":")[0],
                    "cds_stop": row.get("start_stop", "").split(":")[-1],
                    "cds_phase_before_event": "", "cds_phase_after_event": "",
                    "frame_status": "TRANSLATED", "translation_start_source": "SpliceMutr",
                    "transcript_expression_tpm": "", "full_length_status": "TRANSLATED_TRANSCRIPT_HYPOTHESIS",
                    "long_read_support": "UNASSESSED", "nucleotide_sequence_sha256": "",
                    "source_generator": "SpliceMutr", "source_generator_version": args.tool_version,
                    "source_file": source, "source_record_id": source_record,
                    "hypothesis_status": "RESOLVED_EXACT_JUNCTION_TRANSLATION",
                    "evidence_conflict_status": "NONE",
                }
                oid = orf_id(
                    transcript_hypothesis_id_value=sth,
                    protein_sequence_sha256=sequence_sha256(protein),
                    orf_start=1, orf_stop=len(protein), frame_status="TRANSLATED",
                )
                orfs[oid] = {
                    "orf_id": oid, "transcript_hypothesis_id": sth, "splice_event_id": event_id,
                    "sample_id": args.sample_id, "gene": gene, "protein_sequence": protein,
                    "protein_sequence_sha256": sequence_sha256(protein), "protein_length": str(len(protein)),
                    "orf_start": "1", "orf_stop": str(len(protein)), "frame_status": "TRANSLATED",
                    "frameshift_status": "UNASSESSED", "novel_aa_start": "", "novel_aa_end": "",
                    "premature_stop_status": "PRESENT" if str(row.get("peptide", "")).endswith("*") else "NOT_REPORTED",
                    "nmd_risk": "UNASSESSED", "nmd_reason": "", "orf_validity_status": "VALID_TRANSLATED_HYPOTHESIS",
                    "source_generator": "SpliceMutr", "source_generator_version": args.tool_version,
                    "source_file": source, "source_record_id": source_record,
                    "evidence_conflict_status": "NONE",
                }
                for candidate in relevant:
                    peptide = candidate.get("peptide", "").upper()
                    positions = []
                    offset = 0
                    while peptide and (pos := protein.find(peptide, offset)) >= 0:
                        positions.append(pos + 1)
                        offset = pos + 1
                    if len(positions) != 1:
                        continue
                    pstart, pend = positions[0], positions[0] + len(peptide) - 1
                    crosses = bool(boundaries) and pstart <= min(boundaries) and pend >= max(boundaries)
                    if not crosses:
                        continue
                    junction_offset = ";".join(str(x - pstart + 1) for x in boundaries)
                    pid = peptide_id(peptide)
                    por = peptide_origin_id(
                        orf_id_value=oid, splice_event_id_value=event_id, peptide_sequence=peptide,
                        protein_start=pstart, protein_end=pend, junction_offset=junction_offset,
                    )
                    normal_path = modified == "NORMAL"
                    origins[por] = {
                        "origin_peptide_id": por, "peptide_id": pid, "orf_id": oid,
                        "transcript_hypothesis_id": sth, "splice_event_id": event_id,
                        "sample_id": args.sample_id, "gene": gene, "peptide_sequence": peptide,
                        "peptide_length": str(len(peptide)), "protein_start": str(pstart),
                        "protein_end": str(pend), "crosses_junction": "true", "junction_ids": junction,
                        "junction_offset_in_peptide": junction_offset,
                        "contains_novel_aa": "false" if normal_path else "true",
                        "novel_aa_positions": "" if normal_path else f"{pstart}-{pend}",
                        "wildtype_counterpart_status": "SAME_SEQUENCE_IN_NORMAL_TRANSCRIPT" if normal_path else "UNRESOLVED",
                        "wildtype_peptide": peptide if normal_path else "", "reference_proteome_match": "UNASSESSED",
                        "generator_group": "RNA_DRIVEN", "source_generator": "SpliceMutr",
                        "source_generator_version": args.tool_version, "source_file": source,
                        "source_record_id": source_record,
                        "origin_status": "RESOLVED_EXACT_JUNCTION_TRANSLATED_ORIGIN",
                        "evidence_conflict_status": "NORMAL_TRANSCRIPT_SEQUENCE" if normal_path else "NONE",
                    }
                    key = (junction, peptide)
                    current = best_origin.get(key)
                    if current is None or (
                        origins[por]["contains_novel_aa"] == "true"
                        and current.get("contains_novel_aa") != "true"
                    ):
                        best_origin[key] = origins[por]
                    lid = link_id("POL", pid, por, oid, sth, event_id)
                    links[lid] = {
                        "peptide_origin_link_id": lid, "peptide_id": pid, "origin_peptide_id": por,
                        "orf_id": oid, "transcript_hypothesis_id": sth, "splice_event_id": event_id,
                        "sample_id": args.sample_id, "link_status": "RESOLVED_EXACT",
                    }
                    matched_candidates.add((junction, peptide))

    write_tsv(args.outdir / "splice_transcript_hypotheses.tsv", transcripts.values(), TRANSCRIPT_HYPOTHESIS_FIELDS)
    write_tsv(args.outdir / "splice_orfs.tsv", orfs.values(), ORF_FIELDS)
    write_tsv(args.outdir / "splice_peptide_origins.tsv", origins.values(), PEPTIDE_ORIGIN_FIELDS)
    write_tsv(args.outdir / "splice_peptide_origin_links.tsv", links.values(), PEPTIDE_ORIGIN_LINK_FIELDS)
    formal_candidates: list[dict[str, str]] = []
    for junction, rows in candidates.items():
        for candidate in rows:
            row = dict(candidate)
            origin = best_origin.get((junction, candidate.get("peptide", "").upper()))
            if origin:
                row.update({
                    "transcript_hypothesis_id": origin["transcript_hypothesis_id"],
                    "orf_id": origin["orf_id"], "origin_peptide_id": origin["origin_peptide_id"],
                    "junction_ids": junction, "crosses_junction": "yes",
                    "contains_novel_aa": "yes" if origin["contains_novel_aa"] == "true" else "no",
                    "wildtype_peptide": origin["wildtype_peptide"], "orf_evidence_grade": "O1",
                })
                if origin["contains_novel_aa"] == "true":
                    row["mutant_specificity_status"] = "MT_SPECIFIC"
                    row["mutant_specificity_reason"] = "Exact junction-spanning peptide from a changed SpliceMutr translation path."
                else:
                    row["mutant_specificity_status"] = "NON_MUTANT_SEQUENCE"
                    row["mutant_specificity_gate_status"] = "FAIL"
                    row["mutant_specificity_reason"] = "The identical junction-spanning peptide occurs in a SpliceMutr normal transcript path."
                    row["mutant_specificity_priority_cap"] = "R4"
            formal_candidates.append(row)
    write_tsv(args.outdir / "raw_peptides.formal_origins.tsv", formal_candidates, PEPTIDE_FIELDS)
    summary = {
        "sample_id": args.sample_id, "candidate_junctions": len(candidates),
        "candidate_peptide_junction_pairs": sum(len(x) for x in candidates.values()),
        "splicemutr_files": files, "translated_rows_linked_to_candidate_junctions": translated_rows,
        "transcript_hypotheses": len(transcripts), "orfs": len(orfs),
        "peptide_origins": len(origins), "matched_candidate_pairs": len(matched_candidates),
        "formal_candidate_rows": len(formal_candidates),
        "formal_candidate_rows_with_origin": sum(bool(x.get("origin_peptide_id")) for x in formal_candidates),
        "formal_candidate_rows_non_mutant_sequence": sum(x.get("mutant_specificity_status") == "NON_MUTANT_SEQUENCE" for x in formal_candidates),
        "matching_policy": "exact canonical junction + unique exact peptide occurrence + junction-spanning boundary",
    }
    (args.outdir / "rebuild_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
