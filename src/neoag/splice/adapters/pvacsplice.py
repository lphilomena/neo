"""pVACsplice 6.x/7.x adapter with exact variant+junction provenance.

pVACsplice coordinates are interpreted only when the report declares the
zero-based half-open fields documented by pVACtools. A strand must be present
in the report or supplied through an exact Junction→canonical junction map.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Mapping

from neoag.splice.coordinates import CanonicalJunction, normalize_chromosome, normalize_genome_build, normalize_strand
from neoag.splice.identifiers import (
    link_id, orf_id, peptide_id, peptide_origin_id, sequence_sha256,
    stable_id, transcript_hypothesis_id,
)

from ..variants import ExactEventIndex, canonical_variant_id, parse_junction_token, variant_from_row
from .base import as_float_text, as_int, clean, get, infer_mhc_class, read_delimited, row_hash, source_record_id

_AA = set("ACDEFGHIKLMNPQRSTVWY")


def _peptide(value: str) -> str:
    seq = "".join(clean(value).split()).upper().replace("*", "")
    return seq if seq and set(seq) <= _AA else ""


def _residue_positions(row: Mapping[str, str], peptide_length: int) -> list[int]:
    """Parse only explicit mutation/novel-residue positions local to epitope.

    pVACsplice versions differ in naming the mutation-position column.  Values
    outside the reported epitope cannot be interpreted as epitope-local and are
    deliberately discarded rather than guessed from the transcript position.
    """
    raw = get(
        row,
        "Novel AA Positions", "novel_aa_positions", "Mutation Position(s)",
        "Mutation Position", "mutation_position", "Pos",
    )
    positions = sorted({int(x) for x in re.findall(r"\d+", raw) if 1 <= int(x) <= peptide_length})
    return positions


def _junction_boundary(row: Mapping[str, str], peptide_length: int) -> str:
    raw = get(
        row,
        "Junction Offset in Peptide", "junction_offset_in_peptide",
        "Junction AA Position", "junction_aa_position", "Junction Position",
    )
    values = [int(x) for x in re.findall(r"\d+", raw)]
    if len(values) >= 2:
        left, right = values[0], values[1]
    elif len(values) == 1:
        left, right = values[0], values[0] + 1
    else:
        return ""
    if 1 <= left < right <= peptide_length:
        return f"{left}-{right}"
    return ""


def _load_junction_map(path: str | Path | None) -> dict[str, str]:
    if not path:
        return {}
    mapping: dict[str, str] = {}
    for row in read_delimited(path):
        source = get(row, "junction", "source_junction_id", "name", "Junction")
        canonical = get(row, "junction_id", "canonical_junction_id")
        if source and canonical:
            mapping[source] = canonical
    return mapping


def _resolve_junction(row: Mapping[str, str], *, genome_build: str, exact_map: Mapping[str, str]) -> tuple[str, str]:
    token = get(row, "Junction", "junction", "junction_id")
    if token in exact_map:
        return exact_map[token], "MAPPED_EXACT_JUNCTION_MAP"
    parsed = parse_junction_token(token, genome_build=genome_build)
    if parsed:
        return parsed.junction_id, "MAPPED_EXACT_JUNCTION_TOKEN"
    chrom = get(row, "Junction Chromosome", "Chromosome", "chromosome", "chrom")
    start_text = get(row, "Junction Start", "junction_start")
    stop_text = get(row, "Junction Stop", "junction_stop")
    strand = normalize_strand(get(row, "Junction Strand", "junction_strand", "Strand", "strand", default="."))
    if chrom and start_text and stop_text and strand in {"+", "-"}:
        try:
            # pVACsplice report coordinates are zero-based, half-open.
            junction = CanonicalJunction(
                normalize_genome_build(genome_build), normalize_chromosome(chrom),
                int(float(start_text)) + 1, int(float(stop_text)), strand,
            )
            return junction.junction_id, "MAPPED_ZERO_BASED_HALF_OPEN_WITH_EXPLICIT_STRAND"
        except Exception:
            pass
    return "", "JUNCTION_UNRESOLVED"


def parse_pvacsplice(
    path: str | Path,
    *,
    sample_id: str,
    genome_build: str = "GRCh38",
    source_tool_version: str = "UNASSESSED",
    junction_map: str | Path | None = None,
    entity_bundle: Mapping[str, list[Mapping[str, str]]] | None = None,
    strict: bool = False,
) -> dict[str, list[dict[str, str]]]:
    p = Path(path)
    build = normalize_genome_build(genome_build)
    exact_map = _load_junction_map(junction_map)
    event_index = ExactEventIndex.from_tables(entity_bundle or {})
    causal_index = {
        (row.get("variant_id", ""), row.get("junction_id", "")): row
        for row in (entity_bundle or {}).get("causal_links", [])
    }
    result: dict[str, list[dict[str, str]]] = {
        "variants": [], "causal_links": [], "transcripts": [], "orfs": [],
        "peptide_origins": [], "peptide_origin_links": [], "pvacsplice_predictions": [],
        "presentation": [], "tool_evidence": [], "conflicts": [],
    }
    for row_no, row in enumerate(read_delimited(p), start=2):
        rid = source_record_id("pVACsplice", p, row_no, row)
        parsed_variant = variant_from_row(row, genome_build=build, zero_based_half_open=True)
        jid, junction_resolution = _resolve_junction(row, genome_build=build, exact_map=exact_map)
        epitope = _peptide(get(row, "Epitope Seq", "epitope_sequence", "MT Epitope Seq"))
        if not parsed_variant or not jid or not epitope:
            missing = []
            if not parsed_variant:
                missing.append("variant")
            if not jid:
                missing.append("junction_with_strand")
            if not epitope:
                missing.append("epitope")
            result["conflicts"].append({
                "conflict_id": stable_id("CFL", "pVACsplice", rid, missing),
                "entity_type": "PVACSPLICE_ROW", "entity_id": rid, "sample_id": sample_id,
                "conflict_type": "PVACSPLICE_EXACT_PROVENANCE_UNRESOLVED", "field_name": ";".join(missing),
                "observed_values": get(row, "Junction", "junction"), "source_tools": "pVACsplice",
                "source_record_ids": rid, "severity": "ERROR" if strict else "WARNING",
                "resolution_status": "UNRESOLVED",
                "resolution_reason": "pVACsplice rows require exact variant, strand-aware junction, and valid epitope sequence.",
            })
            continue
        vb, vc, vp, vr, va = parsed_variant
        variant_id = canonical_variant_id(vb, vc, vp, vr, va)
        causal = causal_index.get((variant_id, jid), {})
        event_id = causal.get("splice_event_id", "")
        if not event_id:
            event_id, event_resolution = event_index.resolve([jid], gene=get(row, "Gene Name", "gene"))
        else:
            event_resolution = "RESOLVED_EXISTING_CAUSAL_LINK"
        if not event_id:
            result["conflicts"].append({
                "conflict_id": stable_id("CFL", "pVACsplice", rid, "event"),
                "entity_type": "PVACSPLICE_ROW", "entity_id": rid, "sample_id": sample_id,
                "conflict_type": "PVACSPLICE_EVENT_UNRESOLVED", "field_name": "splice_event_id",
                "observed_values": f"{variant_id};{jid}", "source_tools": "pVACsplice",
                "source_record_ids": rid, "severity": "ERROR" if strict else "WARNING",
                "resolution_status": "UNRESOLVED",
                "resolution_reason": "No exact event exists for the canonical junction. Build the formal event layer first.",
            })
            continue
        gene = get(row, "Gene Name", "gene", "gene_name")
        gene_id = get(row, "Ensembl Gene ID", "gene_id")
        transcript = get(row, "Transcript", "transcript_id")
        novel_positions = _residue_positions(row, len(epitope))
        junction_boundary = _junction_boundary(row, len(epitope))
        novelty_status = "true" if novel_positions else "UNASSESSED"
        crossing_status = "true" if junction_boundary else "UNASSESSED"
        result["variants"].append({
            "variant_id": variant_id, "sample_id": sample_id, "genome_build": vb, "chrom": vc,
            "pos_1based": str(vp), "ref": vr, "alt": va,
            "variant_type": get(row, "Variant Type", "variant_type"), "gene": gene, "gene_id": gene_id,
            "transcript_ids": transcript, "hgvsc": get(row, "HGVSc", "hgvsc"),
            "hgvsp": get(row, "HGVSp", "hgvsp"), "spliceai_score": "", "pangolin_score": "",
            "mmsplice_score": "", "ci_spliceai_score": "", "source_tools": "pVACsplice",
            "source_tool_versions": source_tool_version, "source_files": str(p), "source_record_ids": rid,
            "variant_resolution_status": "RESOLVED_EXACT", "evidence_conflict_status": "NONE",
        })
        causal_link_id = stable_id("DCL", variant_id, jid, event_id)
        jreads = as_int(get(row, "Junction Score", "junction_score"), 0)
        result["causal_links"].append({
            "causal_link_id": causal_link_id, "variant_id": variant_id, "junction_id": jid,
            "splice_event_id": event_id, "sample_id": sample_id, "gene": gene, "gene_id": gene_id,
            "transcript_id": transcript, "causal_status": "PVACSPLICE_SUPPORTED",
            "prediction_status": "PVACSPLICE_TRANSLATED_AND_PREDICTED",
            "rna_junction_status": "EXACT_RNA_SUPPORTED" if jreads > 0 else "UNASSESSED",
            "targeted_requant_status": causal.get("targeted_requant_status", "UNASSESSED"),
            "pvacsplice_status": "PVACSPLICE_SUPPORTED", "junction_reads": str(jreads),
            "easyquant_junction_reads": causal.get("easyquant_junction_reads", ""),
            "easyquant_spanning_pairs": causal.get("easyquant_spanning_pairs", ""),
            "spliceai_score": causal.get("spliceai_score", ""), "pangolin_score": causal.get("pangolin_score", ""),
            "mmsplice_score": causal.get("mmsplice_score", ""), "ci_spliceai_score": causal.get("ci_spliceai_score", ""),
            "source_tools": "pVACsplice", "source_tool_versions": source_tool_version,
            "source_files": str(p), "source_record_ids": rid,
            "link_resolution_status": "RESOLVED_EXACT_VARIANT_AND_JUNCTION",
            "resolution_reason": f"{junction_resolution}; {event_resolution}", "evidence_conflict_status": "NONE",
        })

        sth = transcript_hypothesis_id(
            splice_event_id_value=event_id, reference_transcript_id=transcript,
            junction_chain=[jid], path_role="PVACSPLICE_TRANSLATED_EPITOPE",
            source_generator="pVACsplice",
        )
        result["transcripts"].append({
            "transcript_hypothesis_id": sth, "splice_event_id": event_id, "sample_id": sample_id,
            "gene": gene, "gene_id": gene_id, "reference_transcript_id": transcript,
            "mane_status": get(row, "MANE Select", default="UNASSESSED"), "path_id": get(row, "Index", "Fasta Key"),
            "path_role": "PVACSPLICE_TRANSLATED_EPITOPE", "exon_chain": "", "junction_chain": jid,
            "cds_start": "", "cds_stop": "", "cds_phase_before_event": "", "cds_phase_after_event": "",
            "frame_status": "FRAMESHIFT" if get(row, "Frameshift Event").lower() == "yes" else "PVACSPLICE_REPORTED",
            "translation_start_source": "PVACSPLICE", "transcript_expression_tpm": get(row, "Gene Expression", "Transcript Expression"),
            "full_length_status": "EPITOPE_LOCAL_PRODUCT", "long_read_support": "UNASSESSED",
            "nucleotide_sequence_sha256": "", "source_generator": "pVACsplice",
            "source_generator_version": source_tool_version, "source_file": str(p), "source_record_id": rid,
            "hypothesis_status": "RESOLVED_EXACT_CAUSAL_EPITOPE", "evidence_conflict_status": "NONE",
        })
        frame = "FRAMESHIFT" if get(row, "Frameshift Event").lower() == "yes" else "PVACSPLICE_REPORTED"
        oid = orf_id(transcript_hypothesis_id_value=sth, protein_sequence_sha256=sequence_sha256(epitope), frame_status=frame)
        result["orfs"].append({
            "orf_id": oid, "transcript_hypothesis_id": sth, "splice_event_id": event_id,
            "sample_id": sample_id, "gene": gene, "protein_sequence": epitope,
            "protein_sequence_sha256": sequence_sha256(epitope), "protein_length": str(len(epitope)),
            "orf_start": "", "orf_stop": "", "frame_status": frame, "frameshift_status": get(row, "Frameshift Event"),
            "novel_aa_start": str(min(novel_positions)) if novel_positions else "",
            "novel_aa_end": str(max(novel_positions)) if novel_positions else "",
            "premature_stop_status": "UNASSESSED",
            "nmd_risk": "UNASSESSED", "nmd_reason": "", "orf_validity_status": "VALID_EPITOPE_PRODUCT_ONLY",
            "source_generator": "pVACsplice", "source_generator_version": source_tool_version,
            "source_file": str(p), "source_record_id": rid, "evidence_conflict_status": "NONE",
        })
        pid = peptide_id(epitope)
        por = peptide_origin_id(
            orf_id_value=oid, splice_event_id_value=event_id, peptide_sequence=epitope,
            protein_start=1, protein_end=len(epitope),
        )
        result["peptide_origins"].append({
            "origin_peptide_id": por, "peptide_id": pid, "orf_id": oid,
            "transcript_hypothesis_id": sth, "splice_event_id": event_id, "sample_id": sample_id,
            "gene": gene, "peptide_sequence": epitope, "peptide_length": str(len(epitope)),
            "protein_start": "1", "protein_end": str(len(epitope)), "crosses_junction": crossing_status,
            "junction_ids": jid, "required_junction_ids": jid,
            "junction_offset_in_peptide": junction_boundary or "UNASSESSED",
            "contains_novel_aa": novelty_status,
            "novel_aa_positions": ";".join(str(x) for x in novel_positions),
            "wildtype_counterpart_status": "UNRESOLVED",
            "wildtype_peptide": "", "reference_proteome_match": get(row, "Reference Match", default="UNASSESSED"),
            "generator_group": "DNA_CAUSAL", "source_generator": "pVACsplice",
            "source_generator_version": source_tool_version, "source_file": str(p), "source_record_id": rid,
            "origin_status": (
                "RESOLVED_EXACT_VARIANT_JUNCTION_EPITOPE_AND_RESIDUES"
                if novel_positions and junction_boundary
                else "RESOLVED_EXACT_VARIANT_JUNCTION_EPITOPE_NOVELTY_UNASSESSED"
            ),
            "evidence_conflict_status": "NONE",
        })
        result["peptide_origin_links"].append({
            "peptide_origin_link_id": link_id("POL", pid, por), "peptide_id": pid,
            "origin_peptide_id": por, "orf_id": oid, "transcript_hypothesis_id": sth,
            "splice_event_id": event_id, "sample_id": sample_id, "link_status": "RESOLVED_EXACT",
        })
        hla = get(row, "HLA Allele", "hla_allele")
        prediction_id = stable_id("PVS", causal_link_id, por, hla, epitope, rid)
        prediction = {
            "pvacsplice_prediction_id": prediction_id, "causal_link_id": causal_link_id,
            "variant_id": variant_id, "junction_id": jid, "splice_event_id": event_id,
            "origin_peptide_id": por, "peptide_id": pid, "sample_id": sample_id, "chrom": vc,
            "variant_start_1based": str(vp), "variant_stop_1based": str(vp + max(len(vr.replace('-', '')), 1) - 1),
            "ref": vr, "alt": va, "junction_score": str(jreads), "junction_anchor": get(row, "Junction Anchor"),
            "transcript_id": transcript, "gene": gene, "gene_id": gene_id, "hla_allele": hla,
            "mhc_class": infer_mhc_class(hla), "epitope_sequence": epitope,
            "epitope_length": str(len(epitope)), "protein_position": get(row, "Protein Position"),
            "best_ic50": as_float_text(get(row, "Best IC50 Score", "best_ic50")),
            "best_percentile": as_float_text(get(row, "Best Percentile", "best_percentile")),
            "best_binding_percentile": as_float_text(get(row, "Best IC50 Percentile", "Best Binding Percentile")),
            "best_presentation_percentile": as_float_text(get(row, "Best Presentation Percentile")),
            "presentation_score": as_float_text(get(row, "MHCflurryEL Presentation Score", "Best MHCflurryEL Presentation Score")),
            "immunogenicity_score": as_float_text(get(row, "Immunogenicity Score", "BigMHC_IM Score")),
            "prediction_methods": get(row, "Best IC50 Score Method", "Best Percentile Method"),
            "aggregate_tier": get(row, "Tier", "Evaluation"), "source_tool": "pVACsplice",
            "source_tool_version": source_tool_version, "source_file": str(p), "source_record_id": rid,
            "mapping_status": "MAPPED_EXACT_VARIANT_JUNCTION_EVENT",
            "evidence_conflict_status": "NONE",
        }
        result["pvacsplice_predictions"].append(prediction)
        presentation_id = stable_id("PRE", "pVACsplice", prediction_id)
        result["presentation"].append({
            "presentation_id": presentation_id, "origin_peptide_id": por, "peptide_id": pid,
            "orf_id": oid, "transcript_hypothesis_id": sth, "splice_event_id": event_id,
            "sample_id": sample_id, "index": get(row, "Index", "Fasta Key"), "hla_allele": hla,
            "mhc_class": infer_mhc_class(hla), "epitope_sequence": epitope,
            "epitope_length": str(len(epitope)), "sub_peptide_position": get(row, "Protein Position"),
            "median_ic50": as_float_text(get(row, "Median IC50 Score")), "best_ic50": prediction["best_ic50"],
            "median_percentile": as_float_text(get(row, "Median Percentile")), "best_percentile": prediction["best_percentile"],
            "median_binding_percentile": as_float_text(get(row, "Median IC50 Percentile", "Median Binding Percentile")),
            "best_binding_percentile": prediction["best_binding_percentile"],
            "median_presentation_percentile": as_float_text(get(row, "Median Presentation Percentile")),
            "best_presentation_percentile": prediction["best_presentation_percentile"],
            "median_immunogenicity_percentile": as_float_text(get(row, "Median Immunogenicity Percentile")),
            "best_immunogenicity_percentile": as_float_text(get(row, "Best Immunogenicity Percentile")),
            "presentation_score": prediction["presentation_score"], "immunogenicity_score": prediction["immunogenicity_score"],
            "prediction_methods": prediction["prediction_methods"], "aggregate_tier": prediction["aggregate_tier"],
            "result_scope": "DNA_CAUSAL_SPLICE", "source_tool": "pVACsplice",
            "source_tool_version": source_tool_version, "source_file": str(p), "source_record_id": rid,
            "mapping_status": "MAPPED_EXACT_VARIANT_JUNCTION_EVENT", "evidence_conflict_status": "NONE",
        })
        result["tool_evidence"].append({
            "evidence_id": stable_id("EVD", "pVACsplice", causal_link_id, prediction_id),
            "entity_type": "CAUSAL_LINK", "entity_id": causal_link_id, "sample_id": sample_id,
            "evidence_group": "DNA_CAUSAL", "evidence_type": "PVACSPLICE_SUPPORTED",
            "source_tool": "pVACsplice", "source_tool_version": source_tool_version,
            "source_file": str(p), "source_row_number": str(row_no), "source_record_id": rid,
            "provided_value": f"{variant_id}|{jid}|{epitope}|{hla}", "verified_value": prediction_id,
            "resolution_status": "RESOLVED_EXACT_VARIANT_JUNCTION_EVENT",
            "resolution_reason": f"{junction_resolution}; {event_resolution}", "raw_payload_sha256": row_hash(row),
        })
    return result
