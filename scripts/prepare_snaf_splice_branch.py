#!/usr/bin/env python3
"""Normalize, filter, and stage SNAF splice candidates for a phased NeoAg run."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import defaultdict
from pathlib import Path

from neoag.schemas import EVENT_FIELDS, PEPTIDE_FIELDS, RNA_JUNCTION_EVIDENCE_FIELDS
from neoag.utils import read_tsv, write_tsv


COORD_RE = re.compile(r"^(chr[^:]+):(\d+)-(\d+)\(([+-])\)$")
TRANSCRIPT_RE = re.compile(r"ENST\d+(?:\.\d+)?")
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def fnum(value: object, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def safe_id(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_coord(value: str) -> dict[str, object] | None:
    match = COORD_RE.match(value.strip())
    if not match:
        return None
    chrom, start, end, strand = match.groups()
    start_i, end_i = int(start), int(end)
    return {
        "chrom": chrom,
        "start": start_i,
        "end": end_i,
        "strand": strand,
        "raw": f"{chrom}:{start_i}-{end_i}:{strand}",
        "inner": f"{chrom}:{start_i + 1}-{end_i - 1}:{strand}",
    }


def load_normal_junctions(path: Path) -> set[str]:
    keys: set[str] = set()
    for row in read_tsv(path):
        if row.get("junction_class") in {"", "normal_splice_junction"}:
            key = str(row.get("junction_id") or "").strip()
            if key:
                keys.add(key)
    return keys


def load_gene_map(path: Path | None) -> dict[str, str]:
    mapping: dict[str, str] = {}
    if not path or not path.is_file():
        return mapping
    for row in read_tsv(path):
        gene = str(row.get("gene") or "").strip()
        for ensembl in str(row.get("ensembl_gene_id") or "").split(";"):
            ensembl = ensembl.strip().split(".")[0]
            if ensembl and gene:
                mapping[ensembl] = gene
    return mapping


def load_tpm_maps(
    gene_path: Path, transcript_path: Path
) -> tuple[dict[str, float], dict[str, float]]:
    gene_tpm: dict[str, float] = {}
    for row in read_tsv(gene_path):
        value = fnum(row.get("TPM"), -1.0)
        if value < 0:
            continue
        keys = {
            str(row.get("gene") or "").strip(),
            str(row.get("gene_symbol") or "").strip(),
            str(row.get("gene_id") or "").strip().split(".")[0],
            str(row.get("ensembl_gene_id") or "").strip().split(".")[0],
        }
        for key in keys:
            if key:
                gene_tpm[key] = max(value, gene_tpm.get(key, -1.0))

    transcript_tpm: dict[str, float] = {}
    for row in read_tsv(transcript_path):
        transcript = str(row.get("transcript_id") or "").strip().split(".")[0]
        value = fnum(row.get("TPM"), -1.0)
        if transcript and value >= 0:
            transcript_tpm[transcript] = value
    return gene_tpm, transcript_tpm


def supporting_transcripts(evidence: object) -> list[str]:
    return sorted({match.split(".")[0] for match in TRANSCRIPT_RE.findall(str(evidence or ""))})


def empty_row(fields: list[str]) -> dict[str, str]:
    return {field: "" for field in fields}


def score_candidate(row: dict[str, object]) -> float:
    reads = fnum(row["junction_count"])
    mean = fnum(row["tumor_specificity_mean"], 99.0)
    mle = fnum(row["tumor_specificity_mle"], 1.0)
    rank = fnum(row["binding_affinity"], 99.0)
    immuno = fnum(row["immunogenicity"])
    read_score = min(1.0, math.log1p(reads) / math.log1p(200.0))
    normal_score = max(0.0, 1.0 - min(mean, 1.0))
    mle_score = max(0.0, 1.0 - min(mle, 1.0))
    binding_score = max(0.0, 1.0 - min(rank / 2.0, 1.0))
    normal_ref_score = 1.0 if row["normal_junction_status"] == "ABSENT_GTEX_V11" else 0.0
    frame_score = 1.0 if row["frame_evidence_status"] == "TRANSCRIPT_EVIDENCE_PRESENT" else 0.0
    return (
        0.20 * read_score
        + 0.15 * normal_score
        + 0.15 * mle_score
        + 0.20 * binding_score
        + 0.15 * immuno
        + 0.10 * normal_ref_score
        + 0.05 * frame_score
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--snaf-candidates", required=True, type=Path)
    parser.add_argument("--base-run", required=True, type=Path)
    parser.add_argument("--normal-junctions", required=True, type=Path)
    parser.add_argument("--normal-expression", type=Path)
    parser.add_argument("--gene-expression", required=True, type=Path)
    parser.add_argument("--transcript-expression", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--profile", default="sarcoma_rna_supported_v2_provisional")
    parser.add_argument("--total-splice-reads", type=int, default=184007902)
    parser.add_argument("--min-junction-reads", type=int, default=20)
    parser.add_argument("--max-normal-mean", type=float, default=1.0)
    parser.add_argument("--max-mle", type=float, default=0.1)
    parser.add_argument("--max-binding-rank", type=float, default=0.5)
    parser.add_argument("--min-immunogenicity", type=float, default=0.5)
    parser.add_argument("--max-events", type=int, default=500)
    parser.add_argument("--max-pairs-per-event", type=int, default=2)
    args = parser.parse_args()

    splice_dir = args.outdir / "splice_snaf"
    input_dir = args.outdir / "inputs"
    qc_dir = args.outdir / "qc"
    for directory in (splice_dir, input_dir, qc_dir):
        directory.mkdir(parents=True, exist_ok=True)

    normal_keys = load_normal_junctions(args.normal_junctions)
    gene_map = load_gene_map(args.normal_expression)
    gene_tpm, transcript_tpm = load_tpm_maps(args.gene_expression, args.transcript_expression)
    raw_rows = list(csv.DictReader(args.snaf_candidates.open(encoding="utf-8"), delimiter="\t"))

    raw_hits = inner_hits = 0
    seen_coords: set[str] = set()
    for row in raw_rows:
        coord = str(row.get("coord") or "")
        if coord in seen_coords:
            continue
        seen_coords.add(coord)
        parsed = parse_coord(coord)
        if parsed:
            raw_hits += int(parsed["raw"] in normal_keys)
            inner_hits += int(parsed["inner"] in normal_keys)
    convention = "inner" if inner_hits >= raw_hits else "raw"

    annotated: list[dict[str, object]] = []
    for row in raw_rows:
        coord_info = parse_coord(str(row.get("coord") or ""))
        normal_key = str(coord_info[convention]) if coord_info else ""
        if not coord_info:
            normal_status = "UNMAPPED_COORDINATE"
        elif normal_key in normal_keys:
            normal_status = "SEEN_GTEX_V11"
        else:
            normal_status = "ABSENT_GTEX_V11"
        evidence_text = str(row.get("evidences") or "").strip()
        frame_status = (
            "TRANSCRIPT_EVIDENCE_PRESENT"
            if evidence_text not in {"", "()", "[]", "None"}
            else "FRAME_UNRESOLVED"
        )
        ensembl = str(row.get("symbol") or str(row.get("uid") or "").split(":", 1)[0]).split(".")[0]
        gene = gene_map.get(ensembl, ensembl)
        transcripts = supporting_transcripts(row.get("evidences"))
        expressed_transcripts = [
            (transcript, transcript_tpm[transcript])
            for transcript in transcripts
            if transcript in transcript_tpm
        ]
        best_transcript, best_transcript_tpm = (
            max(expressed_transcripts, key=lambda item: item[1])
            if expressed_transcripts
            else ("", None)
        )
        gene_tpm_value = gene_tpm.get(ensembl, gene_tpm.get(gene))
        if gene_tpm_value is not None and best_transcript_tpm is not None:
            expression_status = "GENE_AND_TRANSCRIPT_SUPPORTED"
        elif gene_tpm_value is not None:
            expression_status = "GENE_ONLY_PARTIAL"
        elif best_transcript_tpm is not None:
            expression_status = "TRANSCRIPT_ONLY_PARTIAL"
        else:
            expression_status = "UNASSESSED"
        reads = int(fnum(row.get("junction_count")))
        cpm = reads * 1_000_000.0 / max(args.total_splice_reads, 1)
        peptide = str(row.get("peptide") or "").upper()
        reasons: list[str] = []
        if reads < args.min_junction_reads:
            reasons.append("LOW_JUNCTION_READS")
        if fnum(row.get("tumor_specificity_mean"), 99.0) >= args.max_normal_mean:
            reasons.append("NORMAL_MEAN_HIGH")
        if fnum(row.get("tumor_specificity_mle"), 1.0) > args.max_mle:
            reasons.append("MLE_HIGH")
        if fnum(row.get("binding_affinity"), 99.0) > args.max_binding_rank:
            reasons.append("BINDING_RANK_WEAK")
        if fnum(row.get("immunogenicity")) < args.min_immunogenicity:
            reasons.append("IMMUNOGENICITY_LOW")
        if normal_status != "ABSENT_GTEX_V11":
            reasons.append(normal_status)
        if frame_status != "TRANSCRIPT_EVIDENCE_PRESENT":
            reasons.append("FRAME_UNRESOLVED")
        if len(peptide) not in {8, 9, 10, 11} or not peptide or not set(peptide) <= VALID_AA:
            reasons.append("INVALID_PEPTIDE")
        enriched: dict[str, object] = dict(row)
        enriched.update({
            "gene": gene,
            "ensembl_gene_id": ensembl,
            "junction_coord_normalized": normal_key,
            "junction_coord_convention": convention,
            "normal_junction_status": normal_status,
            "frame_evidence_status": frame_status,
            "junction_cpm": f"{cpm:.6f}",
            "gene_expression_tpm": "" if gene_tpm_value is None else f"{gene_tpm_value:.6f}",
            "supporting_transcript_ids": ";".join(transcripts),
            "expressed_transcript_id": best_transcript,
            "transcript_expression_tpm": (
                "" if best_transcript_tpm is None else f"{best_transcript_tpm:.6f}"
            ),
            "expression_evidence_status": expression_status,
            "filter_status": "PASS" if not reasons else "REVIEW",
            "filter_reasons": ";".join(reasons),
        })
        enriched["splice_composite_score"] = f"{score_candidate(enriched):.6f}"
        annotated.append(enriched)

    annotated_fields = list(annotated[0]) if annotated else []
    all_path = splice_dir / "snaf_candidates.annotated.tsv"
    write_tsv(all_path, annotated, annotated_fields)
    eligible = [row for row in annotated if row["filter_status"] == "PASS"]
    eligible.sort(key=lambda row: (-fnum(row["splice_composite_score"]), -fnum(row["junction_count"]), fnum(row["binding_affinity"])))
    eligible_path = splice_dir / "snaf_candidates.eligible.tsv"
    write_tsv(eligible_path, eligible, annotated_fields)

    by_uid: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in eligible:
        by_uid[str(row["uid"])].append(row)
    event_order = sorted(
        by_uid,
        key=lambda uid: (
            -max(fnum(row["splice_composite_score"]) for row in by_uid[uid]),
            -max(fnum(row["junction_count"]) for row in by_uid[uid]),
            uid,
        ),
    )[: args.max_events]
    selected: list[dict[str, object]] = []
    for uid in event_order:
        rows = sorted(
            by_uid[uid],
            key=lambda row: (-fnum(row["splice_composite_score"]), fnum(row["binding_affinity"]), -fnum(row["immunogenicity"])),
        )
        used: set[tuple[str, str]] = set()
        for row in rows:
            key = (str(row["peptide"]), str(row["hla"]))
            if key in used:
                continue
            used.add(key)
            selected.append(row)
            if len(used) >= args.max_pairs_per_event:
                break
    selected_path = splice_dir / "snaf_candidates.selected.tsv"
    write_tsv(selected_path, selected, annotated_fields)

    event_rows: list[dict[str, str]] = []
    peptide_rows: list[dict[str, str]] = []
    junction_rows: list[dict[str, str]] = []
    pair_rows: list[dict[str, str]] = []
    selected_by_uid: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in selected:
        selected_by_uid[str(row["uid"])].append(row)

    for uid in event_order:
        if uid not in selected_by_uid:
            continue
        representative = selected_by_uid[uid][0]
        coord = parse_coord(str(representative.get("coord") or "")) or {}
        event_id = f"SNAF_SPLICE|{uid}"
        tumor_specificity = 0.5 * (1.0 - min(fnum(representative["tumor_specificity_mle"], 1.0), 1.0)) + 0.5 * max(0.0, 1.0 - min(fnum(representative["tumor_specificity_mean"], 1.0), 1.0))
        event = empty_row(EVENT_FIELDS)
        event.update({
            "event_id": event_id,
            "sample_id": args.sample_id,
            "disease_profile": args.profile,
            "event_type": "Splice",
            "mutation_source": "SNAF",
            "peptide_consequence": "splice_junction",
            "evidence_scope": "RNA_SPLICE",
            "priority_cap": "C_CAUTION",
            "gene": str(representative["gene"]),
            "event_name": str(representative["junction_coord_normalized"]),
            "chrom": str(coord.get("chrom", "")),
            "pos": str(coord.get("start", "")),
            "consequence": "aberrant_splice_junction",
            "transcript_id": str(representative["expressed_transcript_id"]),
            "rna_junction_reads": str(representative["junction_count"]),
            "rna_junction_source": str(args.snaf_candidates),
            "rna_frame_status": "translated_frame_with_transcript_evidence",
            "event_confidence": str(representative["splice_composite_score"]),
            "event_expression": str(
                representative["gene_expression_tpm"]
                or representative["transcript_expression_tpm"]
            ),
            "gene_expression_tpm": str(representative["gene_expression_tpm"]),
            "transcript_expression_tpm": str(representative["transcript_expression_tpm"]),
            "expression_evidence_status": str(representative["expression_evidence_status"]),
            "rna_support_status": "RNA_JUNCTION_SUPPORTED",
            "rna_evidence_completeness": "COMPLETE",
            "rna_evidence_score": "1.0000",
            "tumor_specificity": f"{tumor_specificity:.4f}",
            "ccf_status": "RNA_ONLY_UNRESOLVED",
            "ccf_confidence": "unresolved",
            "ccf_warning": "RNA-only splice event; DNA clonality is not estimated",
            "ccf_method": "RNA_ONLY_UNRESOLVED",
            "clonality_multiplier": "0.8500",
            "safety_status": "SAFETY_PARTIAL",
            "safety_reason": "RNA-only splice candidate requires normal-junction and peptide safety review",
            "source": "SNAF",
        })
        event_rows.append(event)
        for row in selected_by_uid[uid]:
            peptide = str(row["peptide"])
            hla = str(row["hla"])
            peptide_id = safe_id(f"{event_id}|{hla}|{peptide}")
            rank = fnum(row["binding_affinity"], 99.0)
            peptide_row = empty_row(PEPTIDE_FIELDS)
            peptide_row.update({
                "peptide_id": peptide_id,
                "event_id": event_id,
                "sample_id": args.sample_id,
                "event_type": "Splice",
                "mutation_source": "SNAF",
                "peptide_consequence": "splice_junction",
                "evidence_scope": "RNA_SPLICE",
                "priority_cap": "C_CAUTION",
                "gene": str(row["gene"]),
                "peptide": peptide,
                "crosses_junction": "yes",
                "contains_novel_aa": "yes",
                "rna_junction_reads": str(row["junction_count"]),
                "rna_junction_source": str(args.snaf_candidates),
                "rna_frame_status": "translated_frame_with_transcript_evidence",
                "gene_expression_tpm": str(row["gene_expression_tpm"]),
                "transcript_expression_tpm": str(row["transcript_expression_tpm"]),
                "expression_evidence_status": str(row["expression_evidence_status"]),
                "rna_support_status": "RNA_JUNCTION_SUPPORTED",
                "rna_evidence_completeness": "COMPLETE",
                "rna_evidence_score": "1.0000",
                "hla_allele": hla,
                "mhc_class": "I",
                "source_tool": "SNAF",
                "binding_rank": f"{rank:.6g}",
                "el_rank": f"{rank:.6g}",
                "netmhcpan_mt_rank_el": f"{rank:.6g}",
                "immunogenicity_score": str(row["immunogenicity"]),
                "deepimmuno_score": str(row["immunogenicity"]),
                "ccf_status": "RNA_ONLY_UNRESOLVED",
                "ccf_confidence": "unresolved",
                "ccf_warning": "RNA-only splice event",
                "ccf_method": "RNA_ONLY_UNRESOLVED",
                "ccf_multiplier": "0.8500",
                "safety_status": "SAFETY_PARTIAL",
                "safety_reason": "normal-junction and peptide safety review required",
                "review_required": "yes",
                "normal_junction_assessment_status": str(row["normal_junction_status"]),
                "safety_priority_cap": "C_CAUTION",
            })
            peptide_rows.append(peptide_row)
            junction = empty_row(RNA_JUNCTION_EVIDENCE_FIELDS)
            junction.update({
                "evidence_id": safe_id(f"SNAF|{uid}|{peptide}|{hla}"),
                "event_id": event_id,
                "peptide_id": peptide_id,
                "sample_id": args.sample_id,
                "gene": str(row["gene"]),
                "junction_reads": str(row["junction_count"]),
                "junction_source": "SNAF",
                "mutation_source": "SNAF",
                "peptide_consequence": "splice_junction",
                "rna_frame_status": "translated_frame_with_transcript_evidence",
                "rna_support_status": "RNA_JUNCTION_SUPPORTED",
                "rna_evidence_completeness": "COMPLETE",
                "rna_evidence_score": "1.0000",
                "targeted_validation_status": "NOT_PERFORMED",
            })
            junction_rows.append(junction)
            pair_rows.append({
                "peptide": peptide,
                "hla_allele": hla,
                "gene": str(row["gene"]),
                "event_id": event_id,
            })

    snaf_events = splice_dir / "snaf_raw_events.tsv"
    snaf_peptides = splice_dir / "snaf_raw_peptides.tsv"
    rna_junctions = splice_dir / "snaf_rna_junction_evidence.tsv"
    pair_input = splice_dir / "snaf_selected_peptide_hla.tsv"
    write_tsv(snaf_events, event_rows, EVENT_FIELDS)
    write_tsv(snaf_peptides, peptide_rows, PEPTIDE_FIELDS)
    write_tsv(rna_junctions, junction_rows, RNA_JUNCTION_EVIDENCE_FIELDS)
    write_tsv(pair_input, pair_rows, ["peptide", "hla_allele", "gene", "event_id"])

    base_events = read_tsv(args.base_run / "parsed" / "raw_events.tsv")
    base_peptides = read_tsv(args.base_run / "parsed" / "raw_peptides.tsv")
    combined_events = input_dir / "combined_raw_events.tsv"
    combined_peptides = input_dir / "combined_raw_peptides.tsv"
    write_tsv(combined_events, base_events + event_rows, EVENT_FIELDS)
    write_tsv(combined_peptides, base_peptides + peptide_rows, PEPTIDE_FIELDS)

    summary = {
        "status": "PASS",
        "source": str(args.snaf_candidates),
        "base_run": str(args.base_run),
        "coordinate_convention": convention,
        "coordinate_match_counts": {"raw": raw_hits, "inner_start_plus1_end_minus1": inner_hits},
        "total_candidate_rows": len(annotated),
        "eligible_rows": len(eligible),
        "eligible_junctions": len(by_uid),
        "selected_junctions": len(event_rows),
        "selected_peptide_hla_rows": len(peptide_rows),
        "thresholds": {
            "min_junction_reads": args.min_junction_reads,
            "max_normal_mean": args.max_normal_mean,
            "max_mle": args.max_mle,
            "max_binding_rank": args.max_binding_rank,
            "min_immunogenicity": args.min_immunogenicity,
            "frame_transcript_evidence_required": True,
            "normal_gtex_v11_absence_required": True,
            "max_events": args.max_events,
            "max_pairs_per_event": args.max_pairs_per_event,
        },
        "expression_note": (
            "gene TPM is loaded from the matched WTS RSEM table; transcript TPM is the highest "
            "TPM among SNAF frame-supporting transcripts; junction CPM remains junction evidence only"
        ),
        "expression_sources": {
            "gene_tpm": {
                "path": str(args.gene_expression),
                "sha256": sha256(args.gene_expression),
                "evidence_tool": "RSEM",
                "annotation": "GENCODE_v49",
            },
            "transcript_tpm": {
                "path": str(args.transcript_expression),
                "sha256": sha256(args.transcript_expression),
                "evidence_tool": "RSEM",
                "annotation": "GENCODE_v49",
            },
        },
        "expression_coverage": {
            "selected_events_with_gene_tpm": sum(bool(row.get("gene_expression_tpm")) for row in event_rows),
            "selected_events_with_transcript_tpm": sum(
                bool(row.get("transcript_expression_tpm")) for row in event_rows
            ),
        },
        "ccf_note": "all SNAF splice events are RNA_ONLY_UNRESOLVED",
        "outputs": {},
    }
    for name, path in {
        "annotated": all_path,
        "eligible": eligible_path,
        "selected": selected_path,
        "snaf_raw_events": snaf_events,
        "snaf_raw_peptides": snaf_peptides,
        "rna_junction_evidence": rna_junctions,
        "peptide_hla_input": pair_input,
        "combined_raw_events": combined_events,
        "combined_raw_peptides": combined_peptides,
    }.items():
        summary["outputs"][name] = {"path": str(path), "sha256": sha256(path), "size": path.stat().st_size}
    summary_path = qc_dir / "snaf_splice_prepare_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    (qc_dir / "snaf_splice_prepare_summary.md").write_text(
        "# SNAF splice branch preparation\n\n"
        f"- Candidate rows: {len(annotated):,}\n"
        f"- Eligible rows: {len(eligible):,}\n"
        f"- Eligible junctions: {len(by_uid):,}\n"
        f"- Selected junctions: {len(event_rows):,}\n"
        f"- Selected peptide-HLA pairs: {len(peptide_rows):,}\n"
        f"- GTEx coordinate convention: `{convention}` (raw matches={raw_hits}, inner matches={inner_hits})\n"
        "- CCF: `RNA_ONLY_UNRESOLVED`\n"
        "- Priority cap before orthogonal validation: `C_CAUTION`\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
