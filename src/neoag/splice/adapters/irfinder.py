"""IRFinder-S adapter for retained-intron events."""
from __future__ import annotations

from pathlib import Path

from neoag.splice.coordinates import CanonicalJunction, JunctionNormalizationError, convert_interval, normalize_chromosome, normalize_genome_build, normalize_strand
from neoag.splice.identifiers import link_id, splice_event_id, transcript_hypothesis_id

from .base import as_float_text, as_int, clean, get, join_tokens, read_delimited, row_hash, source_record_id


def parse_irfinder(
    path: str | Path,
    *,
    sample_id: str,
    genome_build: str = "GRCh38",
    coordinate_system: str = "intron_1based_closed",
    source_tool_version: str = "UNASSESSED",
    strict: bool = False,
) -> dict[str, list[dict[str, str]]]:
    p = Path(path)
    build = normalize_genome_build(genome_build)
    result: dict[str, list[dict[str, str]]] = {
        "junctions": [], "events": [], "event_junction_links": [], "transcripts": [],
        "tool_evidence": [], "conflicts": [],
    }
    for row_no, row in enumerate(read_delimited(p), start=2):
        record_id = source_record_id("IRFinder-S", p, row_no, row)
        chrom = normalize_chromosome(get(row, "Chr", "chrom", "chromosome"))
        strand = normalize_strand(get(row, "Strand", "strand", default="."))
        start = as_int(get(row, "Start", "start", "intron_start"), 0)
        end = as_int(get(row, "End", "end", "intron_end"), 0)
        gene = get(row, "Name", "GeneName", "gene", "gene_name")
        gene_id = get(row, "GeneID", "gene_id")
        try:
            start1, end1 = convert_interval(start, end, coordinate_system)
            junction = CanonicalJunction(build, chrom, start1, end1, strand)
            if strict and strand not in {"+", "-"}:
                raise JunctionNormalizationError("strand is required in strict mode")
        except (ValueError, JunctionNormalizationError) as exc:
            result["conflicts"].append({
                "entity_type": "SPLICE_EVENT", "entity_id": record_id, "sample_id": sample_id,
                "conflict_type": "IR_COORDINATE_UNRESOLVED", "field_name": "Start/End",
                "observed_values": f"{chrom}:{start}-{end}:{strand}", "source_tools": "IRFinder-S",
                "source_record_ids": record_id, "severity": "ERROR" if strict else "WARNING",
                "resolution_status": "UNRESOLVED", "resolution_reason": str(exc),
            })
            if strict:
                raise
            continue
        stranded = strand in {"+", "-"}
        event_id = splice_event_id(
            genome_build=build, event_type="RI", strand=strand,
            junction_ids=[junction.junction_id], gene=gene or gene_id,
            affected_exons=[f"{chrom}:{start1}-{end1}:{strand}"],
        )
        ir_ratio = as_float_text(get(row, "IRratio", "IRRatio", "ir_ratio"))
        warnings = get(row, "Warnings", "Warning", "warnings")
        cnn = get(row, "CNN", "CNNscore", "Validation", "validated", "Pass")
        depth = get(row, "IntronDepth", "Depth", "Coverage", "intron_depth")
        splice_left = get(row, "SpliceLeft", "splice_left")
        splice_right = get(row, "SpliceRight", "splice_right")
        splice_exact = get(row, "SpliceExact", "MaxSplice", "splice_exact")
        annotation = "IRFINDER_VALIDATED" if clean(cnn).casefold() in {"pass", "true", "1", "validated"} else "IRFINDER_QUANTIFIED"
        result["junctions"].append({
            "junction_id": junction.junction_id, "sample_id": sample_id, "genome_build": build,
            "chrom": chrom, "intron_start_1based": str(start1), "intron_end_1based": str(end1),
            "strand": strand, "donor_1based": str(junction.donor_1based), "acceptor_1based": str(junction.acceptor_1based),
            "splice_motif": "", "annotation_status": annotation,
            "unique_split_reads": splice_exact, "multi_split_reads": "", "total_split_reads": splice_exact,
            "max_overhang": "", "source_coordinate_systems": coordinate_system,
            "source_tools": "IRFinder-S", "source_tool_versions": source_tool_version,
            "source_files": str(p), "source_record_ids": record_id, "provenance_record_count": "1",
            "junction_resolution_status": "RESOLVED_EXACT" if stranded else "RESOLVED_UNSTRANDED",
            "evidence_conflict_status": "NONE" if stranded else "JUNCTION_STRAND_UNRESOLVED",
        })
        result["events"].append({
            "splice_event_id": event_id, "sample_id": sample_id, "genome_build": build,
            "event_type": "RI", "gene": gene, "gene_id": gene_id, "strand": strand,
            "junction_ids": junction.junction_id, "reference_junction_ids": junction.junction_id,
            "alternative_junction_ids": "", "affected_exons": f"{chrom}:{start1}-{end1}:{strand}",
            "annotation_status": annotation, "cryptic_exon_status": "NOT_APPLICABLE",
            "psi": ir_ratio, "delta_psi": "", "qvalue": "", "outlier_score": "",
            "event_expression": depth, "event_confidence": "PASS" if not warnings else "REVIEW",
            "reference_path_status": "SPLICED_PATH_ANNOTATED", "cohort_analysis_status": "NOT_APPLICABLE",
            "source_tools": "IRFinder-S", "source_tool_versions": source_tool_version,
            "source_files": str(p), "source_record_ids": record_id, "provenance_record_count": "1",
            "event_resolution_status": "RESOLVED" if stranded else "RESOLVED_UNSTRANDED",
            "evidence_conflict_status": ("NONE" if not warnings else "IRFINDER_WARNING") if stranded else "JUNCTION_STRAND_UNRESOLVED",
        })
        paths = [("SPLICED", [junction.junction_id]), ("RETAINED", [])]
        for idx, (role, chain) in enumerate(paths, start=1):
            path_id = f"{event_id}:{role}"
            sth = transcript_hypothesis_id(
                splice_event_id_value=event_id, junction_chain=chain,
                path_role=role, source_generator="IRFinder-S",
            )
            result["transcripts"].append({
                "transcript_hypothesis_id": sth, "splice_event_id": event_id, "sample_id": sample_id,
                "gene": gene, "gene_id": gene_id, "reference_transcript_id": "", "mane_status": "UNASSESSED",
                "path_id": path_id, "path_role": role, "exon_chain": "", "junction_chain": ";".join(chain),
                "cds_start": "", "cds_stop": "", "cds_phase_before_event": "", "cds_phase_after_event": "",
                "frame_status": "UNASSESSED", "translation_start_source": "UNASSESSED",
                "transcript_expression_tpm": "", "full_length_status": "INTRON_LOCAL_EVENT",
                "long_read_support": "UNASSESSED", "nucleotide_sequence_sha256": "",
                "source_generator": "IRFinder-S", "source_generator_version": source_tool_version,
                "source_file": str(p), "source_record_id": record_id,
                "hypothesis_status": "IR_EVENT_PATH", "evidence_conflict_status": "NONE",
            })
            for edge_idx, jid in enumerate(chain, start=1):
                result["event_junction_links"].append({
                    "event_junction_link_id": link_id("EJL", event_id, jid, path_id),
                    "splice_event_id": event_id, "junction_id": jid, "sample_id": sample_id,
                    "path_id": path_id, "path_role": role, "edge_index": str(edge_idx),
                    "junction_role": "EXCISED_INTRON", "source_tool": "IRFinder-S",
                    "source_record_id": record_id, "link_status": "RESOLVED",
                })
        result["tool_evidence"].append({
            "entity_type": "SPLICE_EVENT", "entity_id": event_id, "sample_id": sample_id,
            "evidence_group": "INTRON_RETENTION", "evidence_type": "IRFINDER_IR_RATIO",
            "source_tool": "IRFinder-S", "source_tool_version": source_tool_version,
            "source_file": str(p), "source_row_number": str(row_no), "source_record_id": record_id,
            "provided_value": ir_ratio, "verified_value": ir_ratio if stranded else "",
            "resolution_status": "RESOLVED_EXACT" if stranded else "RESOLVED_UNSTRANDED",
            "resolution_reason": f"IRratio={ir_ratio}; depth={depth}; splice_left={splice_left}; splice_right={splice_right}; splice_exact={splice_exact}; warnings={warnings or 'NONE'}; validation={cnn or 'UNASSESSED'}; strand={'resolved' if stranded else 'unresolved'}",
            "raw_payload_sha256": row_hash(row),
        })
        if not stranded:
            result["conflicts"].append({
                "entity_type": "SPLICE_EVENT", "entity_id": event_id, "sample_id": sample_id,
                "conflict_type": "JUNCTION_STRAND_UNRESOLVED", "field_name": "Strand",
                "observed_values": strand, "source_tools": "IRFinder-S", "source_record_ids": record_id,
                "severity": "WARNING", "resolution_status": "UNRESOLVED",
                "resolution_reason": "The intron-retention record is retained for review but cannot contribute exact stranded support.",
            })
    return result
