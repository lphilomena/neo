"""ImmunoPepper adapter.

ImmunoPepper emits event-local translated sequences rather than guaranteed
full-length transcripts.  The adapter therefore labels its transcript and ORF
entities as hypotheses/partial translated segments and never upgrades them to
full-length evidence without an independent source.
"""
from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path
from typing import Any

from neoag.splice.coordinates import CanonicalJunction, normalize_chromosome, normalize_genome_build, normalize_strand
from neoag.splice.identifiers import (
    link_id,
    orf_id,
    peptide_id,
    peptide_origin_id,
    sequence_sha256,
    splice_event_id,
    transcript_hypothesis_id,
    unresolved_splice_event_id,
)

from .base import as_int, clean, get, join_tokens, read_delimited, row_hash, source_record_id, truth_text

_COORD_RE = re.compile(r"(?:(?P<chrom>chr(?:[0-9]+|X|Y|M|MT))[:_])?(?P<start>\d+)\s*[-:]\s*(?P<end>\d+)", re.I)


def _parse_exons(value: str, default_chrom: str) -> list[tuple[str, int, int]]:
    """Parse ImmunoPepper exon coordinates without guessing across records.

    Official outputs commonly encode ``modifiedExonsCoord`` as four or six
    semicolon-separated integers (two or three exon intervals).  Some site
    workflows instead emit explicit ``chr:start-end`` intervals.  Both forms
    are accepted; odd or incomplete coordinate lists are rejected.
    """
    raw = clean(value)
    exons: list[tuple[str, int, int]] = []
    explicit = list(_COORD_RE.finditer(raw))
    if explicit:
        for match in explicit:
            chrom = normalize_chromosome(match.group("chrom") or default_chrom)
            start, end = int(match.group("start")), int(match.group("end"))
            exons.append((chrom, min(start, end), max(start, end)))
    else:
        tokens = [token for token in re.split(r"[;,:|\s]+", raw.strip("[](){}")) if token]
        numeric = [int(token) for token in tokens if re.fullmatch(r"-?\d+", token)]
        # ImmunoPepper emits two or three exon intervals.  Never silently drop
        # an unpaired coordinate because it would alter the junction identity.
        if len(numeric) >= 4 and len(numeric) % 2 == 0:
            chrom = normalize_chromosome(default_chrom)
            for idx in range(0, len(numeric), 2):
                start, end = numeric[idx], numeric[idx + 1]
                exons.append((chrom, min(start, end), max(start, end)))
    # Preserve source path order while removing exact duplicates.
    seen: set[tuple[str, int, int]] = set()
    return [x for x in exons if not (x in seen or seen.add(x))]


def _junction_aa_boundaries(exons: list[tuple[str, int, int]], strand: str) -> list[str]:
    """Return event-local amino-acid boundary tokens for a translated segment.

    ``N-M`` means the junction lies between amino acids N and M. ``N`` means
    the codon for amino acid N itself spans the exon junction.  These offsets
    are derived from ImmunoPepper's CDS-adjusted exon intervals and are used
    only to decide whether a pVACbind sub-peptide overlaps the junction.
    """
    ordered = list(reversed(exons)) if strand == "-" else list(exons)
    cumulative_nt = 0
    boundaries: list[str] = []
    for _, start, end in ordered[:-1]:
        cumulative_nt += end - start + 1
        quotient, remainder = divmod(cumulative_nt, 3)
        if remainder == 0:
            boundaries.append(f"{quotient}-{quotient + 1}")
        else:
            boundaries.append(str(quotient + 1))
    return boundaries


def _junctions_from_exons(build: str, strand: str, exons: list[tuple[str, int, int]]) -> list[CanonicalJunction]:
    if len(exons) < 2 or len({x[0] for x in exons}) != 1:
        return []
    genomic = sorted(exons, key=lambda x: (x[1], x[2]))
    result: list[CanonicalJunction] = []
    for (_, _, left_end), (chrom, right_start, _) in zip(genomic, genomic[1:]):
        if left_end + 1 <= right_start - 1:
            result.append(CanonicalJunction(build, chrom, left_end + 1, right_start - 1, strand))
    return result


def _valid_aa(sequence: str) -> str:
    seq = "".join(clean(sequence).split()).upper().replace("*", "")
    return seq if seq and all(ch in "ACDEFGHIKLMNPQRSTVWY" for ch in seq) else ""


def _event_type_from_mutation_mode(value: str) -> str:
    key = re.sub(r"[^a-z0-9]+", "_", clean(value).casefold()).strip("_")
    if "exitron" in key or "exonic_intron" in key:
        return "EXITRON"
    if "cryptic" in key and "exon" in key:
        return "CRYPTIC_EXON"
    if "novel" in key and "junction" in key:
        return "NOVEL_JUNCTION"
    if "intron" in key and ("ret" in key or "retain" in key):
        return "RI"
    return "COMPLEX_SPLICE"


def parse_immunopepper_meta(
    path: str | Path,
    *,
    sample_id: str,
    genome_build: str = "GRCh38",
    source_tool_version: str = "UNASSESSED",
) -> dict[str, list[dict[str, str]]]:
    p = Path(path)
    build = normalize_genome_build(genome_build)
    result: dict[str, list[dict[str, str]]] = {
        "junctions": [], "events": [], "event_junction_links": [], "transcripts": [],
        "orfs": [], "peptide_origins": [], "peptide_origin_links": [],
        "tool_evidence": [], "conflicts": [],
    }
    for row_no, row in enumerate(read_delimited(p), start=2):
        record_id = source_record_id("ImmunoPepper", p, row_no, row)
        gene = get(row, "geneName", "gene", "gene_name")
        gene_id = get(row, "geneId", "gene_id")
        chrom = normalize_chromosome(get(row, "geneChr", "chrom", "chr"))
        strand = normalize_strand(get(row, "geneStrand", "strand", default="."))
        exons = _parse_exons(get(row, "modifiedExonsCoord", "exon_chain", "modified_exons_coord"), chrom)
        junctions = _junctions_from_exons(build, strand, exons)
        junction_ids = [j.junction_id for j in junctions]
        source_id = get(row, "id", "transcript_id", "index", default=f"row_{row_no}")
        mutation_mode = get(row, "mutationMode", "mutation_mode", "kmerType", default="SPLICE")
        event_type = _event_type_from_mutation_mode(mutation_mode)
        event_id = (
            splice_event_id(
                genome_build=build, event_type=event_type, strand=strand,
                junction_ids=junction_ids, gene=gene or gene_id,
                affected_exons=[f"{c}:{a}-{b}:{strand}" for c, a, b in exons],
            )
            if junction_ids
            else unresolved_splice_event_id(source_tool="ImmunoPepper", source_record_id=record_id, event_type=event_type)
        )
        stranded = strand in {"+", "-"}
        for junction in junctions:
            result["junctions"].append({
                "junction_id": junction.junction_id, "sample_id": sample_id, "genome_build": build,
                "chrom": junction.chrom, "intron_start_1based": str(junction.intron_start_1based),
                "intron_end_1based": str(junction.intron_end_1based), "strand": junction.strand,
                "donor_1based": str(junction.donor_1based), "acceptor_1based": str(junction.acceptor_1based),
                "splice_motif": "", "annotation_status": "IMMUNOPEPPER_TRANSLATED_PATH",
                "unique_split_reads": "0", "multi_split_reads": "0", "total_split_reads": "0", "max_overhang": "",
                "source_coordinate_systems": "modified_exons_1based_closed", "source_tools": "ImmunoPepper",
                "source_tool_versions": source_tool_version, "source_files": str(p), "source_record_ids": record_id,
                "provenance_record_count": "1",
                "junction_resolution_status": "RESOLVED_FROM_EXON_CHAIN" if stranded else "RESOLVED_UNSTRANDED",
                "evidence_conflict_status": "NONE" if stranded else "JUNCTION_STRAND_UNRESOLVED",
            })
        result["events"].append({
            "splice_event_id": event_id, "sample_id": sample_id, "genome_build": build,
            "event_type": event_type, "gene": gene, "gene_id": gene_id, "strand": strand,
            "junction_ids": join_tokens(junction_ids), "reference_junction_ids": "",
            "alternative_junction_ids": join_tokens(junction_ids),
            "affected_exons": join_tokens(f"{c}:{a}-{b}:{strand}" for c, a, b in exons),
            "annotation_status": "IMMUNOPEPPER_TRANSLATED",
            "cryptic_exon_status": "PRESENT" if event_type == "CRYPTIC_EXON" else "NOT_APPLICABLE",
            "psi": "", "delta_psi": "", "qvalue": "", "outlier_score": "",
            "event_expression": get(row, "variantSegExpr", "expression"),
            "event_confidence": "TRANSLATED_PATH", "reference_path_status": "UNRESOLVED",
            "cohort_analysis_status": "NOT_APPLICABLE", "source_tools": "ImmunoPepper",
            "source_tool_versions": source_tool_version, "source_files": str(p),
            "source_record_ids": record_id, "provenance_record_count": "1",
            "event_resolution_status": ("RESOLVED" if stranded else "RESOLVED_UNSTRANDED") if junction_ids else "JUNCTION_UNRESOLVED",
            "evidence_conflict_status": ("NONE" if stranded else "JUNCTION_STRAND_UNRESOLVED") if junction_ids else "JUNCTION_UNRESOLVED",
        })
        if junction_ids and not stranded:
            result["conflicts"].append({
                "entity_type": "SPLICE_EVENT", "entity_id": event_id, "sample_id": sample_id,
                "conflict_type": "JUNCTION_STRAND_UNRESOLVED", "field_name": "geneStrand",
                "observed_values": strand, "source_tools": "ImmunoPepper", "source_record_ids": record_id,
                "severity": "WARNING", "resolution_status": "UNRESOLVED",
                "resolution_reason": "The translated path is retained for review but cannot contribute exact stranded junction support.",
            })
        exon_chain = [f"{c}:{a}-{b}:{strand}" for c, a, b in (reversed(exons) if strand == "-" else exons)]
        jchain = list(reversed(junction_ids)) if strand == "-" else junction_ids
        frame_raw = get(row, "readFrame", "reading_frame", "frame")
        frame_annotated = truth_text(get(row, "readFrameAnnotated", "frame_annotated"))
        frame_status = f"FRAME_{frame_raw}" if frame_raw else "UNASSESSED"
        if frame_annotated == "true":
            frame_status += "_ANNOTATED"
        sth = transcript_hypothesis_id(
            splice_event_id_value=event_id, reference_transcript_id=source_id,
            exon_chain=exon_chain, junction_chain=jchain, path_role="TRANSLATED_PATH",
            source_generator="ImmunoPepper",
        )
        result["transcripts"].append({
            "transcript_hypothesis_id": sth, "splice_event_id": event_id, "sample_id": sample_id,
            "gene": gene, "gene_id": gene_id, "reference_transcript_id": source_id,
            "mane_status": "UNASSESSED", "path_id": source_id, "path_role": "TRANSLATED_PATH",
            "exon_chain": ";".join(exon_chain), "junction_chain": ";".join(jchain),
            "cds_start": "", "cds_stop": "", "cds_phase_before_event": frame_raw,
            "cds_phase_after_event": "", "frame_status": frame_status,
            "translation_start_source": "IMMUNOPEPPER_FRAME", "transcript_expression_tpm": get(row, "variantSegExpr"),
            "full_length_status": "PARTIAL_JUNCTION_TRANSLATION", "long_read_support": "UNASSESSED",
            "nucleotide_sequence_sha256": "", "source_generator": "ImmunoPepper",
            "source_generator_version": source_tool_version, "source_file": str(p), "source_record_id": record_id,
            "hypothesis_status": "EVENT_LOCAL_TRANSLATED_HYPOTHESIS", "evidence_conflict_status": "NONE",
        })
        for edge_index, jid in enumerate(jchain, start=1):
            result["event_junction_links"].append({
                "event_junction_link_id": link_id("EJL", event_id, jid, source_id, edge_index),
                "splice_event_id": event_id, "junction_id": jid, "sample_id": sample_id,
                "path_id": source_id, "path_role": "TRANSLATED_PATH", "edge_index": str(edge_index),
                "junction_role": "TRANSLATED_EDGE", "source_tool": "ImmunoPepper",
                "source_record_id": record_id, "link_status": "RESOLVED",
            })
        protein = _valid_aa(get(row, "peptide", "protein_sequence", "translated_sequence"))
        is_isolated = truth_text(get(row, "isIsolated", "is_isolated"))
        crosses = "false" if is_isolated == "true" else ("true" if len(exons) >= 2 else "UNASSESSED")
        boundaries = _junction_aa_boundaries(exons, strand)
        if protein:
            psha = sequence_sha256(protein)
            oid = orf_id(
                transcript_hypothesis_id_value=sth, protein_sequence_sha256=psha,
                orf_start=1, orf_stop=len(protein), frame_status=frame_status,
            )
            has_stop = truth_text(get(row, "hasStopCodon", "has_stop_codon"))
            result["orfs"].append({
                "orf_id": oid, "transcript_hypothesis_id": sth, "splice_event_id": event_id,
                "sample_id": sample_id, "gene": gene, "protein_sequence": protein,
                "protein_sequence_sha256": psha, "protein_length": str(len(protein)), "orf_start": "1",
                "orf_stop": str(len(protein)), "frame_status": frame_status,
                "frameshift_status": "UNASSESSED", "novel_aa_start": "1", "novel_aa_end": str(len(protein)),
                "premature_stop_status": "STOP_PRESENT" if has_stop == "true" else "UNASSESSED",
                "nmd_risk": "UNASSESSED", "nmd_reason": "Full transcript context is unavailable.",
                "orf_validity_status": "PARTIAL_TRANSLATED_SEGMENT", "source_generator": "ImmunoPepper",
                "source_generator_version": source_tool_version, "source_file": str(p), "source_record_id": record_id,
                "evidence_conflict_status": "NONE",
            })
            pid = peptide_id(protein)
            por = peptide_origin_id(
                orf_id_value=oid, splice_event_id_value=event_id, peptide_sequence=protein,
                protein_start=1, protein_end=len(protein),
            )
            result["peptide_origins"].append({
                "origin_peptide_id": por, "peptide_id": pid, "orf_id": oid,
                "transcript_hypothesis_id": sth, "splice_event_id": event_id, "sample_id": sample_id,
                "gene": gene, "peptide_sequence": protein, "peptide_length": str(len(protein)),
                "protein_start": "1", "protein_end": str(len(protein)), "crosses_junction": crosses,
                "junction_ids": join_tokens(jchain), "junction_offset_in_peptide": ";".join(boundaries),
                "contains_novel_aa": "UNASSESSED", "novel_aa_positions": "",
                "wildtype_counterpart_status": "WT_UNRESOLVED", "wildtype_peptide": "",
                "reference_proteome_match": "UNASSESSED", "generator_group": "TRANSLATION",
                "source_generator": "ImmunoPepper", "source_generator_version": source_tool_version,
                "source_file": str(p), "source_record_id": record_id,
                "origin_status": "PARTIAL_TRANSLATED_CANDIDATE", "evidence_conflict_status": "NONE",
            })
            result["peptide_origin_links"].append({
                "peptide_origin_link_id": link_id("POL", pid, por, oid, sth, event_id),
                "peptide_id": pid, "origin_peptide_id": por, "orf_id": oid,
                "transcript_hypothesis_id": sth, "splice_event_id": event_id,
                "sample_id": sample_id, "link_status": "RESOLVED",
            })
        else:
            result["conflicts"].append({
                "entity_type": "TRANSCRIPT_HYPOTHESIS", "entity_id": sth, "sample_id": sample_id,
                "conflict_type": "TRANSLATED_SEQUENCE_MISSING_OR_INVALID", "field_name": "peptide",
                "observed_values": get(row, "peptide"), "source_tools": "ImmunoPepper",
                "source_record_ids": record_id, "severity": "WARNING", "resolution_status": "UNRESOLVED",
                "resolution_reason": "No valid amino-acid sequence was available for ORF/peptide registration.",
            })
        result["tool_evidence"].append({
            "entity_type": "TRANSCRIPT_HYPOTHESIS", "entity_id": sth, "sample_id": sample_id,
            "evidence_group": "TRANSLATION", "evidence_type": "IMMUNOPEPPER_TRANSLATED_PATH",
            "source_tool": "ImmunoPepper", "source_tool_version": source_tool_version,
            "source_file": str(p), "source_row_number": str(row_no), "source_record_id": record_id,
            "provided_value": protein, "verified_value": protein,
            "resolution_status": "PARTIAL_TRANSLATED_SEGMENT" if protein else "SEQUENCE_UNRESOLVED",
            "resolution_reason": (
                f"mutationMode={mutation_mode}; frame={frame_status}; source_id={source_id}; "
                f"isJunctionList={get(row, 'isJunctionList') or 'UNASSESSED'}; "
                f"isIsolated={is_isolated}; aa_boundaries={','.join(boundaries) or 'UNRESOLVED'}"
            ),
            "raw_payload_sha256": row_hash(row),
        })
    return result


def parse_immunopepper_kmers(
    path: str | Path,
    *,
    sample_id: str,
    source_tool_version: str = "UNASSESSED",
    meta_bundle: dict[str, list[dict[str, str]]] | None = None,
) -> dict[str, list[dict[str, str]]]:
    """Register ImmunoPepper k-mers and map them only when the ORF link is unique."""
    p = Path(path)
    result = {"peptide_origins": [], "peptide_origin_links": [], "tool_evidence": [], "conflicts": []}
    orfs = (meta_bundle or {}).get("orfs", [])
    transcripts = {r["transcript_hypothesis_id"]: r for r in (meta_bundle or {}).get("transcripts", [])}
    by_sequence: dict[str, list[dict[str, str]]] = defaultdict(list)
    for orf in orfs:
        by_sequence[orf.get("protein_sequence", "")].append(orf)
    for row_no, row in enumerate(read_delimited(p), start=2):
        record_id = source_record_id("ImmunoPepper", p, row_no, row)
        kmer = _valid_aa(get(row, "kmer", "peptide", "sequence"))
        if not kmer:
            continue
        candidates = [orf for protein, rows in by_sequence.items() if kmer in protein for orf in rows]
        if len(candidates) != 1:
            result["conflicts"].append({
                "entity_type": "PEPTIDE", "entity_id": peptide_id(kmer), "sample_id": sample_id,
                "conflict_type": "IMMUNOPEPPER_KMER_ORF_MAPPING_AMBIGUOUS" if candidates else "IMMUNOPEPPER_KMER_ORF_UNRESOLVED",
                "field_name": "kmer", "observed_values": kmer, "source_tools": "ImmunoPepper",
                "source_record_ids": record_id, "severity": "WARNING", "resolution_status": "UNRESOLVED",
                "resolution_reason": f"Candidate ORF count={len(candidates)}; no approximate mapping was used.",
            })
            continue
        orf = candidates[0]
        protein = orf["protein_sequence"]
        start = protein.index(kmer) + 1
        end = start + len(kmer) - 1
        sth = transcripts.get(orf["transcript_hypothesis_id"], {})
        event_id = orf["splice_event_id"]
        pid = peptide_id(kmer)
        por = peptide_origin_id(
            orf_id_value=orf["orf_id"], splice_event_id_value=event_id,
            peptide_sequence=kmer, protein_start=start, protein_end=end,
        )
        crosses = truth_text(get(row, "isCrossJunction", "is_cross_junction"))
        junction_annotated = truth_text(get(row, "junctionAnnotated", "junction_annotated"))
        result["peptide_origins"].append({
            "origin_peptide_id": por, "peptide_id": pid, "orf_id": orf["orf_id"],
            "transcript_hypothesis_id": orf["transcript_hypothesis_id"], "splice_event_id": event_id,
            "sample_id": sample_id, "gene": orf.get("gene", ""), "peptide_sequence": kmer,
            "peptide_length": str(len(kmer)), "protein_start": str(start), "protein_end": str(end),
            "crosses_junction": crosses, "junction_ids": sth.get("junction_chain", ""),
            "junction_offset_in_peptide": get(row, "coord", "junction_offset"),
            "contains_novel_aa": "UNASSESSED", "novel_aa_positions": "",
            "wildtype_counterpart_status": "WT_UNRESOLVED", "wildtype_peptide": "",
            "reference_proteome_match": "UNASSESSED", "generator_group": "TRANSLATION",
            "source_generator": "ImmunoPepper", "source_generator_version": source_tool_version,
            "source_file": str(p), "source_record_id": record_id,
            "origin_status": "KMER_MAPPED_TO_UNIQUE_ORF", "evidence_conflict_status": "NONE",
        })
        result["peptide_origin_links"].append({
            "peptide_origin_link_id": link_id("POL", pid, por, orf["orf_id"], event_id),
            "peptide_id": pid, "origin_peptide_id": por, "orf_id": orf["orf_id"],
            "transcript_hypothesis_id": orf["transcript_hypothesis_id"], "splice_event_id": event_id,
            "sample_id": sample_id, "link_status": "RESOLVED",
        })
        result["tool_evidence"].append({
            "entity_type": "PEPTIDE_ORIGIN", "entity_id": por, "sample_id": sample_id,
            "evidence_group": "TRANSLATION", "evidence_type": "IMMUNOPEPPER_KMER",
            "source_tool": "ImmunoPepper", "source_tool_version": source_tool_version,
            "source_file": str(p), "source_row_number": str(row_no), "source_record_id": record_id,
            "provided_value": kmer, "verified_value": kmer, "resolution_status": "MAPPED_UNIQUE_ORF",
            "resolution_reason": f"cross_junction={crosses}; junction_annotated={junction_annotated}", "raw_payload_sha256": row_hash(row),
        })
    return result
