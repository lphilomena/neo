"""moPepGen variant-peptide FASTA adapter.

moPepGen FASTA records are peptide products, not guaranteed full ORFs. This
adapter therefore distinguishes exact peptide-level consensus from full ORF
consensus. Exact splice provenance requires either a junction-bearing GVF INFO
field or the explicit provenance map emitted by the NeoAg preparation wrapper.
"""
from __future__ import annotations

import gzip
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping

from neoag.splice.identifiers import (
    link_id, orf_id, peptide_id, peptide_origin_id, sequence_sha256,
    stable_id, transcript_hypothesis_id, unresolved_splice_event_id,
)

from ..variants import ExactEventIndex, parse_junction_token
from .base import clean, get, join_tokens, parse_attributes, read_delimited, row_hash, source_record_id

_AA = set("ACDEFGHIKLMNPQRSTVWY")


def _open(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open("r", encoding="utf-8")


def _fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header = ""
    chunks: list[str] = []
    with _open(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    records.append((header, "".join(chunks).upper()))
                header, chunks = line[1:].strip(), []
            else:
                chunks.append(line)
    if header:
        records.append((header, "".join(chunks).upper()))
    return records


def _valid_peptide(value: str) -> str:
    seq = "".join(clean(value).split()).upper().replace("*", "")
    return seq if seq and set(seq) <= _AA else ""


def _parse_info(value: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in clean(value).split(";"):
        if not item:
            continue
        if "=" in item:
            key, val = item.split("=", 1)
            attrs[key.strip()] = val.strip().strip('"')
    return attrs


def parse_mopepgen_gvf(paths: Iterable[str | Path], *, genome_build: str = "GRCh38") -> dict[str, dict[str, str]]:
    """Return GVF variant metadata keyed by moPepGen variant ID."""
    result: dict[str, dict[str, str]] = {}
    for path in paths:
        p = Path(path)
        header: list[str] | None = None
        with _open(p) as handle:
            for raw in handle:
                line = raw.rstrip("\n")
                if not line or line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    header = line.lstrip("#").split("\t")
                    continue
                parts = line.split("\t")
                if header and len(parts) >= len(header):
                    row = dict(zip(header, parts))
                elif len(parts) >= 8:
                    row = dict(zip(["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"], parts[:8]))
                else:
                    continue
                variant_id = clean(row.get("ID"))
                if not variant_id:
                    continue
                info = _parse_info(row.get("INFO", ""))
                junction_token = (
                    info.get("JUNCTION_ID") or info.get("JUNC_ID") or info.get("GENOMIC_JUNCTION")
                    or info.get("SPLICE_JUNCTION") or ""
                )
                result[variant_id] = {
                    "variant_id": variant_id,
                    "gene_id": clean(row.get("CHROM")),
                    "transcript_id": info.get("TRANSCRIPT_ID", "").split(",")[0],
                    "gene": info.get("GENE_SYMBOL", ""),
                    "junction_id": junction_token,
                    "source_file": str(p),
                    "ref": clean(row.get("REF")), "alt": clean(row.get("ALT")),
                    "pos": clean(row.get("POS")),
                    "info": clean(row.get("INFO")),
                }
    return result


def _map_rows(paths: Iterable[str | Path]) -> tuple[dict[str, list[dict[str, str]]], dict[str, list[dict[str, str]]]]:
    by_header: dict[str, list[dict[str, str]]] = defaultdict(list)
    by_variant: dict[str, list[dict[str, str]]] = defaultdict(list)
    for path in paths:
        for row in read_delimited(path):
            header = get(row, "mopepgen_header", "fasta_header", "header", "record_id")
            variant = get(row, "variant_id", "mopepgen_variant_id")
            if header:
                by_header[header].append(row)
            if variant:
                by_variant[variant].append(row)
    return by_header, by_variant


def _header_origins(header: str) -> list[tuple[str, list[str], str]]:
    origins: list[tuple[str, list[str], str]] = []
    for token in header.split():
        parts = [part for part in token.split("|") if part]
        if not parts:
            continue
        transcript = parts[0]
        index = parts[-1] if parts[-1].isdigit() else ""
        variants = parts[1:-1] if index else parts[1:]
        origins.append((transcript, variants, index))
    return origins


def parse_mopepgen(
    fasta_path: str | Path,
    *,
    sample_id: str,
    genome_build: str = "GRCh38",
    gvf_paths: Iterable[str | Path] = (),
    provenance_maps: Iterable[str | Path] = (),
    source_tool_version: str = "UNASSESSED",
    entity_bundle: Mapping[str, list[Mapping[str, str]]] | None = None,
    strict: bool = False,
) -> dict[str, list[dict[str, str]]]:
    p = Path(fasta_path)
    gvf = parse_mopepgen_gvf(gvf_paths, genome_build=genome_build)
    by_header, by_variant = _map_rows(provenance_maps)
    event_index = ExactEventIndex.from_tables(entity_bundle or {})
    result: dict[str, list[dict[str, str]]] = {
        "transcripts": [], "orfs": [], "peptide_origins": [], "peptide_origin_links": [],
        "tool_evidence": [], "conflicts": [],
    }
    for row_no, (header, raw_sequence) in enumerate(_fasta(p), start=1):
        peptide = _valid_peptide(raw_sequence)
        record_id = source_record_id("moPepGen", p, row_no, {"header": header, "sequence": raw_sequence})
        if not peptide:
            result["conflicts"].append({
                "conflict_id": stable_id("CFL", "moPepGen", record_id, "sequence"),
                "entity_type": "MOPEPGEN_FASTA", "entity_id": record_id, "sample_id": sample_id,
                "conflict_type": "MOPEPGEN_INVALID_PEPTIDE_SEQUENCE", "field_name": "sequence",
                "observed_values": raw_sequence, "source_tools": "moPepGen", "source_record_ids": record_id,
                "severity": "ERROR", "resolution_status": "UNRESOLVED",
                "resolution_reason": "FASTA sequence is not a valid amino-acid peptide.",
            })
            continue
        for origin_index, (transcript, variant_ids, ordinal) in enumerate(_header_origins(header), start=1):
            explicit_rows = list(by_header.get(header, []))
            for variant_id_value in variant_ids:
                explicit_rows.extend(by_variant.get(variant_id_value, []))
            # Deterministic de-duplication of identical map rows.
            seen_payloads: set[str] = set()
            explicit_rows = [row for row in explicit_rows if not (row_hash(row) in seen_payloads or seen_payloads.add(row_hash(row)))]
            candidate_jids: set[str] = set()
            candidate_event_ids: set[str] = set()
            genes: set[str] = set()
            map_protein = ""
            map_frame = ""
            map_crosses = "UNASSESSED"
            map_novel = "true"
            map_offset = ""
            for row in explicit_rows:
                token = get(row, "junction_id", "canonical_junction_id", "junc_id")
                junction = parse_junction_token(token, genome_build=genome_build) if token else None
                if junction:
                    candidate_jids.add(junction.junction_id)
                event_value = get(row, "splice_event_id", "event_id")
                if event_value:
                    candidate_event_ids.add(event_value)
                gene = get(row, "gene", "gene_name", "gene_symbol")
                if gene:
                    genes.add(gene)
                map_protein = _valid_peptide(get(row, "protein_sequence", "translated_sequence", "orf_sequence")) or map_protein
                map_frame = get(row, "frame_status", "reading_frame", default=map_frame)
                map_crosses = get(row, "crosses_junction", default=map_crosses)
                map_novel = get(row, "contains_novel_aa", default=map_novel)
                map_offset = get(row, "junction_offset_in_peptide", default=map_offset)
            for variant_id_value in variant_ids:
                meta = gvf.get(variant_id_value, {})
                genes.update([meta.get("gene", "")])
                token = meta.get("junction_id", "")
                junction = parse_junction_token(token, genome_build=genome_build) if token else None
                if junction:
                    candidate_jids.add(junction.junction_id)
            genes.discard("")
            event_id = ""
            resolution = ""
            if len(candidate_event_ids) == 1:
                event_id = next(iter(candidate_event_ids))
                resolution = "RESOLVED_EXPLICIT_EVENT_MAP"
            elif len(candidate_event_ids) > 1:
                resolution = "AMBIGUOUS_EXPLICIT_EVENT_MAP"
            elif candidate_jids:
                event_id, resolution = event_index.resolve(sorted(candidate_jids), gene=next(iter(genes)) if len(genes) == 1 else "")
            if not event_id:
                event_id = unresolved_splice_event_id(
                    source_tool="moPepGen", source_record_id=f"{record_id}|{origin_index}",
                    event_type="MOPEPGEN_PEPTIDE_UNRESOLVED_EVENT",
                )
            gene = next(iter(genes)) if len(genes) == 1 else ""
            protein = map_protein or peptide
            full_length_status = "EXPLICIT_TRANSLATED_SEQUENCE" if map_protein else "PEPTIDE_SEQUENCE_ONLY"
            frame = map_frame or ("UNASSESSED" if not map_protein else "FRAME_REPORTED")
            sth = transcript_hypothesis_id(
                splice_event_id_value=event_id, reference_transcript_id=transcript,
                junction_chain=sorted(candidate_jids), path_role="MOPEPGEN_VARIANT_PEPTIDE",
                source_generator="moPepGen",
            )
            result["transcripts"].append({
                "transcript_hypothesis_id": sth, "splice_event_id": event_id, "sample_id": sample_id,
                "gene": gene, "gene_id": "", "reference_transcript_id": transcript,
                "mane_status": "UNASSESSED", "path_id": header, "path_role": "MOPEPGEN_VARIANT_PEPTIDE",
                "exon_chain": "", "junction_chain": join_tokens(candidate_jids), "cds_start": "", "cds_stop": "",
                "cds_phase_before_event": "", "cds_phase_after_event": "", "frame_status": frame,
                "translation_start_source": "MOPEPGEN_VARIANT_GRAPH", "transcript_expression_tpm": "",
                "full_length_status": full_length_status, "long_read_support": "UNASSESSED",
                "nucleotide_sequence_sha256": "", "source_generator": "moPepGen",
                "source_generator_version": source_tool_version, "source_file": str(p),
                "source_record_id": record_id, "hypothesis_status": resolution or "EVENT_UNRESOLVED",
                "evidence_conflict_status": "NONE" if resolution.startswith("RESOLVED") else "EVENT_UNRESOLVED",
            })
            oid = orf_id(
                transcript_hypothesis_id_value=sth, protein_sequence_sha256=sequence_sha256(protein),
                frame_status=frame,
            )
            result["orfs"].append({
                "orf_id": oid, "transcript_hypothesis_id": sth, "splice_event_id": event_id,
                "sample_id": sample_id, "gene": gene, "protein_sequence": protein,
                "protein_sequence_sha256": sequence_sha256(protein), "protein_length": str(len(protein)),
                "orf_start": "", "orf_stop": "", "frame_status": frame,
                "frameshift_status": "UNASSESSED", "novel_aa_start": "1" if not map_protein else "",
                "novel_aa_end": str(len(protein)) if not map_protein else "", "premature_stop_status": "UNASSESSED",
                "nmd_risk": "UNASSESSED", "nmd_reason": "",
                "orf_validity_status": "VALID_TRANSLATED_SEQUENCE" if map_protein else "VALID_PEPTIDE_PRODUCT_ONLY",
                "source_generator": "moPepGen", "source_generator_version": source_tool_version,
                "source_file": str(p), "source_record_id": record_id,
                "evidence_conflict_status": "NONE" if resolution.startswith("RESOLVED") else "EVENT_UNRESOLVED",
            })
            start = protein.find(peptide) + 1 if peptide in protein else 1
            pid = peptide_id(peptide)
            por = peptide_origin_id(
                orf_id_value=oid, splice_event_id_value=event_id, peptide_sequence=peptide,
                protein_start=start, protein_end=start + len(peptide) - 1, junction_offset=map_offset,
            )
            origin_status = "RESOLVED_EXACT_EVENT_AND_PEPTIDE" if resolution.startswith("RESOLVED") else "EVENT_UNRESOLVED"
            result["peptide_origins"].append({
                "origin_peptide_id": por, "peptide_id": pid, "orf_id": oid,
                "transcript_hypothesis_id": sth, "splice_event_id": event_id, "sample_id": sample_id,
                "gene": gene, "peptide_sequence": peptide, "peptide_length": str(len(peptide)),
                "protein_start": str(start), "protein_end": str(start + len(peptide) - 1),
                "crosses_junction": map_crosses, "junction_ids": join_tokens(candidate_jids),
                "required_junction_ids": join_tokens(candidate_jids),
                "junction_offset_in_peptide": map_offset, "contains_novel_aa": map_novel,
                "novel_aa_positions": "", "wildtype_counterpart_status": "UNRESOLVED",
                "wildtype_peptide": "", "reference_proteome_match": "UNASSESSED",
                "generator_group": "RNA_DRIVEN", "source_generator": "moPepGen",
                "source_generator_version": source_tool_version, "source_file": str(p),
                "source_record_id": record_id, "origin_status": origin_status,
                "evidence_conflict_status": "NONE" if resolution.startswith("RESOLVED") else "EVENT_UNRESOLVED",
            })
            result["peptide_origin_links"].append({
                "peptide_origin_link_id": link_id("POL", pid, por), "peptide_id": pid,
                "origin_peptide_id": por, "orf_id": oid, "transcript_hypothesis_id": sth,
                "splice_event_id": event_id, "sample_id": sample_id,
                "link_status": "RESOLVED_EXACT" if resolution.startswith("RESOLVED") else "EVENT_UNRESOLVED",
            })
            result["tool_evidence"].append({
                "evidence_id": stable_id("EVD", "moPepGen", event_id, pid, record_id, origin_index),
                "entity_type": "PEPTIDE_ORIGIN", "entity_id": por, "sample_id": sample_id,
                "evidence_group": "RNA_DRIVEN_TRANSLATION", "evidence_type": "MOPEPGEN_VARIANT_PEPTIDE",
                "source_tool": "moPepGen", "source_tool_version": source_tool_version,
                "source_file": str(p), "source_row_number": str(row_no), "source_record_id": record_id,
                "provided_value": header, "verified_value": peptide,
                "resolution_status": origin_status,
                "resolution_reason": (
                    "Exact peptide sequence and exact event/junction map were retained."
                    if resolution.startswith("RESOLVED") else
                    "Peptide is retained but does not contribute cross-generator event consensus until an exact map is supplied."
                ),
                "raw_payload_sha256": row_hash({"header": header, "sequence": peptide}),
            })
            if not resolution.startswith("RESOLVED"):
                result["conflicts"].append({
                    "conflict_id": stable_id("CFL", "moPepGen", record_id, origin_index, "event"),
                    "entity_type": "PEPTIDE_ORIGIN", "entity_id": por, "sample_id": sample_id,
                    "conflict_type": "MOPEPGEN_EVENT_PROVENANCE_UNRESOLVED", "field_name": "splice_event_id",
                    "observed_values": join_tokens(candidate_event_ids or candidate_jids),
                    "source_tools": "moPepGen", "source_record_ids": record_id,
                    "severity": "ERROR" if strict else "WARNING", "resolution_status": "UNRESOLVED",
                    "resolution_reason": "Provide an exact moPepGen provenance map with splice_event_id or canonical junction_id.",
                })
    return result
