"""SplAdder event adapters.

The parser accepts confirmed-event GFF3 and the documented event TXT family.
It preserves both isoform paths and derives exact intron edges from adjacent
1-based closed exons.  It deliberately does not claim which path is canonical
unless the source explicitly labels it.
"""
from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

from neoag.splice.coordinates import CanonicalJunction, normalize_chromosome, normalize_genome_build, normalize_strand
from neoag.splice.identifiers import link_id, splice_event_id, transcript_hypothesis_id
from neoag.utils import open_text_maybe_gz

from .base import as_float_text, as_int, clean, get, join_tokens, parse_attributes, read_delimited, row_hash, source_record_id

_EVENT_ALIASES = {
    "exon_skip": "SE", "exon_skipping": "SE", "skip": "SE", "se": "SE",
    "alt_3prime": "A3SS", "alt3": "A3SS", "a3ss": "A3SS",
    "alt_5prime": "A5SS", "alt5": "A5SS", "a5ss": "A5SS",
    "intron_retention": "RI", "retained_intron": "RI", "ri": "RI",
    "mutex_exons": "MXE", "mutually_exclusive": "MXE", "mxe": "MXE",
    "mult_exon_skip": "MULTI_SE", "multi_exon_skip": "MULTI_SE",
    "cryptic_exon": "CRYPTIC_EXON", "cryptic": "CRYPTIC_EXON",
    "exitron": "EXITRON", "exonic_intron": "EXITRON",
    "novel_junction": "NOVEL_JUNCTION", "novel_splice_junction": "NOVEL_JUNCTION",
    "complex_splice": "COMPLEX_SPLICE", "complex": "COMPLEX_SPLICE",
}


def infer_event_type(path: str | Path, explicit: str = "") -> str:
    if clean(explicit):
        key = re.sub(r"[^a-z0-9]+", "_", clean(explicit).casefold()).strip("_")
        return _EVENT_ALIASES.get(key, clean(explicit).upper())
    name = re.sub(r"[^a-z0-9]+", "_", Path(path).name.casefold()).strip("_")
    for token, mapped in sorted(_EVENT_ALIASES.items(), key=lambda kv: -len(kv[0])):
        if token in name:
            return mapped
    return "COMPLEX_SPLICE"


def _exon_token(chrom: str, start: int, end: int, strand: str) -> str:
    return f"{chrom}:{start}-{end}:{strand}"


def _path_junctions(genome_build: str, chrom: str, strand: str, exons: list[tuple[int, int]]) -> list[CanonicalJunction]:
    if len(exons) < 2:
        return []
    genomic = sorted(exons)
    result: list[CanonicalJunction] = []
    for (_, left_end), (right_start, _) in zip(genomic, genomic[1:]):
        intron_start = left_end + 1
        intron_end = right_start - 1
        if intron_start <= intron_end:
            result.append(CanonicalJunction(genome_build, chrom, intron_start, intron_end, strand))
    return result


def _classify_path_role(value: str) -> str:
    """Classify a path only when the source role is explicit and unambiguous."""
    token = re.sub(r"[^a-z0-9]+", "_", clean(value).casefold()).strip("_")
    if token in {"reference", "ref", "canonical", "normal", "reference_path", "canonical_path"}:
        return "REFERENCE"
    if token in {"alternative", "alt", "variant", "event", "alternative_path", "variant_path"}:
        return "ALTERNATIVE"
    return "UNRESOLVED"


def _junction_row(j: CanonicalJunction, *, sample_id: str, path: Path, record_id: str, tool_version: str) -> dict[str, str]:
    stranded = j.strand in {"+", "-"}
    return {
        "junction_id": j.junction_id, "sample_id": sample_id, "genome_build": j.genome_build,
        "chrom": j.chrom, "intron_start_1based": str(j.intron_start_1based),
        "intron_end_1based": str(j.intron_end_1based), "strand": j.strand,
        "donor_1based": str(j.donor_1based), "acceptor_1based": str(j.acceptor_1based),
        "splice_motif": "", "annotation_status": "SPLADDER_GRAPH_EDGE",
        "unique_split_reads": "0", "multi_split_reads": "0", "total_split_reads": "0", "max_overhang": "",
        "source_coordinate_systems": "gff3_exon_1based_closed", "source_tools": "SplAdder",
        "source_tool_versions": tool_version, "source_files": str(path), "source_record_ids": record_id,
        "provenance_record_count": "1",
        "junction_resolution_status": "RESOLVED_FROM_EXON_CHAIN" if stranded else "RESOLVED_UNSTRANDED",
        "evidence_conflict_status": "NONE" if stranded else "JUNCTION_STRAND_UNRESOLVED",
    }


def _build_event_bundle(
    *,
    sample_id: str,
    genome_build: str,
    chrom: str,
    strand: str,
    gene: str,
    gene_id: str,
    event_type: str,
    source_event_id: str,
    paths: list[dict[str, Any]],
    source_path: Path,
    source_record: str,
    tool_version: str,
    annotation_status: str = "SPLADDER_CONFIRMED",
    psi: str = "",
) -> dict[str, list[dict[str, str]]]:
    all_junctions: dict[str, CanonicalJunction] = {}
    path_junctions: dict[str, list[str]] = {}
    affected_exons: list[str] = []
    for path in paths:
        exons = [(int(a), int(b)) for a, b in path.get("exons", [])]
        affected_exons.extend(_exon_token(chrom, a, b, strand) for a, b in exons)
        js = _path_junctions(genome_build, chrom, strand, exons)
        path_junctions[path["path_id"]] = [j.junction_id for j in js]
        for j in js:
            all_junctions[j.junction_id] = j
    event_id = splice_event_id(
        genome_build=genome_build,
        event_type=event_type,
        strand=strand,
        junction_ids=all_junctions,
        gene=gene or gene_id,
        affected_exons=affected_exons,
    ) if all_junctions else link_id("SEV0", source_event_id, event_type, gene, affected_exons)
    stranded = strand in {"+", "-"}

    reference_junctions: set[str] = set()
    alternative_junctions: set[str] = set()
    classified_roles: list[str] = []
    for path in paths:
        classification = _classify_path_role(clean(path.get("path_role")))
        classified_roles.append(classification)
        if classification == "REFERENCE":
            reference_junctions.update(path_junctions.get(path["path_id"], []))
        elif classification == "ALTERNATIVE":
            alternative_junctions.update(path_junctions.get(path["path_id"], []))
    if reference_junctions and alternative_junctions:
        reference_path_status = "RESOLVED_EXPLICIT_SOURCE_ROLE"
    elif alternative_junctions:
        reference_path_status = "ALTERNATIVE_ONLY_REFERENCE_UNRESOLVED"
    elif reference_junctions:
        reference_path_status = "REFERENCE_ONLY_ALTERNATIVE_UNRESOLVED"
    else:
        reference_path_status = "UNRESOLVED"

    junction_rows = [_junction_row(j, sample_id=sample_id, path=source_path, record_id=source_record, tool_version=tool_version) for j in all_junctions.values()]
    event_row = {
        "splice_event_id": event_id, "sample_id": sample_id, "genome_build": genome_build,
        "event_type": event_type, "gene": gene, "gene_id": gene_id, "strand": strand,
        "junction_ids": join_tokens(all_junctions), "reference_junction_ids": join_tokens(reference_junctions),
        "alternative_junction_ids": join_tokens(alternative_junctions), "affected_exons": join_tokens(affected_exons),
        "annotation_status": annotation_status, "cryptic_exon_status": "UNASSESSED",
        "psi": psi, "delta_psi": "", "qvalue": "", "outlier_score": "", "event_expression": "",
        "event_confidence": "SPLADDER_GRAPH_CONFIRMED", "reference_path_status": reference_path_status,
        "cohort_analysis_status": "NOT_APPLICABLE", "source_tools": "SplAdder",
        "source_tool_versions": tool_version, "source_files": str(source_path),
        "source_record_ids": source_record, "provenance_record_count": "1",
        "event_resolution_status": ("RESOLVED" if stranded else "RESOLVED_UNSTRANDED") if all_junctions else "NO_JUNCTION_EDGE",
        "evidence_conflict_status": "NONE" if stranded else "JUNCTION_STRAND_UNRESOLVED",
    }
    links: list[dict[str, str]] = []
    transcripts: list[dict[str, str]] = []
    for idx, path in enumerate(paths, start=1):
        path_id = clean(path.get("path_id")) or f"path_{idx}"
        role = clean(path.get("path_role")) or f"ISOFORM_{idx}"
        classified_role = _classify_path_role(role)
        exons = [(int(a), int(b)) for a, b in path.get("exons", [])]
        exon_chain = [_exon_token(chrom, a, b, strand) for a, b in (sorted(exons, reverse=(strand == "-")))]
        junction_chain = path_junctions.get(path_id, [])
        if strand == "-":
            junction_chain = list(reversed(junction_chain))
        sth = transcript_hypothesis_id(
            splice_event_id_value=event_id,
            reference_transcript_id=clean(path.get("reference_transcript_id")),
            exon_chain=exon_chain,
            junction_chain=junction_chain,
            path_role=role,
            source_generator="SplAdder",
        )
        transcripts.append({
            "transcript_hypothesis_id": sth, "splice_event_id": event_id, "sample_id": sample_id,
            "gene": gene, "gene_id": gene_id, "reference_transcript_id": clean(path.get("reference_transcript_id")),
            "mane_status": "UNASSESSED", "path_id": path_id, "path_role": role,
            "exon_chain": ";".join(exon_chain), "junction_chain": ";".join(junction_chain),
            "cds_start": "", "cds_stop": "", "cds_phase_before_event": "", "cds_phase_after_event": "",
            "frame_status": "UNASSESSED", "translation_start_source": "UNASSESSED",
            "transcript_expression_tpm": "", "full_length_status": "EVENT_LOCAL_PATH",
            "long_read_support": "UNASSESSED", "nucleotide_sequence_sha256": "",
            "source_generator": "SplAdder", "source_generator_version": tool_version,
            "source_file": str(source_path), "source_record_id": source_record,
            "hypothesis_status": "GRAPH_PATH_HYPOTHESIS", "evidence_conflict_status": "NONE",
        })
        for edge_index, jid in enumerate(junction_chain, start=1):
            links.append({
                "event_junction_link_id": link_id("EJL", event_id, jid, path_id, edge_index),
                "splice_event_id": event_id, "junction_id": jid, "sample_id": sample_id,
                "path_id": path_id, "path_role": role, "edge_index": str(edge_index),
                "junction_role": f"{classified_role}_PATH_EDGE" if classified_role != "UNRESOLVED" else "PATH_EDGE_UNRESOLVED",
                "source_record_id": source_record, "link_status": "RESOLVED",
            })
    evidence = [{
        "entity_type": "SPLICE_EVENT", "entity_id": event_id, "sample_id": sample_id,
        "evidence_group": "SPLICE_GRAPH", "evidence_type": "SPLADDER_EVENT_GRAPH",
        "source_tool": "SplAdder", "source_tool_version": tool_version, "source_file": str(source_path),
        "source_row_number": "", "source_record_id": source_record, "provided_value": source_event_id,
        "verified_value": join_tokens(all_junctions), "resolution_status": event_row["event_resolution_status"],
        "resolution_reason": "Exact intron edges derived from adjacent GFF3/TXT exon coordinates.",
        "raw_payload_sha256": row_hash({"paths": paths, "source_event_id": source_event_id}),
    }]
    evidence.extend({
        "entity_type": "JUNCTION", "entity_id": jid, "sample_id": sample_id,
        "evidence_group": "SPLICE_GRAPH", "evidence_type": "SPLADDER_GRAPH_EDGE",
        "source_tool": "SplAdder", "source_tool_version": tool_version, "source_file": str(source_path),
        "source_row_number": "", "source_record_id": source_record, "provided_value": source_event_id,
        "verified_value": jid if stranded else "",
        "resolution_status": "RESOLVED_EXACT" if stranded else "RESOLVED_UNSTRANDED",
        "resolution_reason": "Junction is an exact stranded edge in a SplAdder event path." if stranded else "SplAdder edge has no resolved strand; exact support was withheld.",
        "raw_payload_sha256": row_hash({"junction_id": jid, "source_event_id": source_event_id}),
    } for jid in all_junctions)
    conflicts = [] if stranded else [{
        "entity_type": "SPLICE_EVENT", "entity_id": event_id, "sample_id": sample_id,
        "conflict_type": "JUNCTION_STRAND_UNRESOLVED", "field_name": "strand",
        "observed_values": strand, "source_tools": "SplAdder", "source_record_ids": source_record,
        "severity": "WARNING", "resolution_status": "UNRESOLVED",
        "resolution_reason": "The event graph is retained for review but cannot contribute exact stranded junction support.",
    }]
    return {
        "junctions": junction_rows, "events": [event_row], "event_junction_links": links,
        "transcripts": transcripts, "tool_evidence": evidence, "conflicts": conflicts,
    }


def _merge_bundles(bundles: Iterable[dict[str, list[dict[str, str]]]]) -> dict[str, list[dict[str, str]]]:
    result: dict[str, list[dict[str, str]]] = defaultdict(list)
    for bundle in bundles:
        for key, rows in bundle.items():
            result[key].extend(rows)
    return dict(result)


def parse_spladder_gff3(
    path: str | Path,
    *,
    sample_id: str,
    genome_build: str = "GRCh38",
    event_type: str = "",
    source_tool_version: str = "UNASSESSED",
) -> dict[str, list[dict[str, str]]]:
    p = Path(path)
    build = normalize_genome_build(genome_build)
    features: dict[str, dict[str, Any]] = {}
    children: dict[str, list[str]] = defaultdict(list)
    with open_text_maybe_gz(p) as fh:
        for line_no, line in enumerate(fh, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            cells = line.rstrip("\n").split("\t")
            if len(cells) < 9:
                continue
            chrom, source, feature, start, end, score, strand, phase, attr_text = cells[:9]
            attrs = parse_attributes(attr_text)
            fid = attrs.get("ID") or f"line_{line_no}"
            parents = [x for x in attrs.get("Parent", "").split(",") if x]
            features[fid] = {
                "id": fid, "parents": parents, "chrom": normalize_chromosome(chrom), "source": source,
                "feature": feature.casefold(), "start": as_int(start), "end": as_int(end),
                "strand": normalize_strand(strand), "attrs": attrs, "line_no": line_no,
            }
            for parent in parents:
                children[parent].append(fid)

    roots = [fid for fid, item in features.items() if item["feature"] in {"gene", "event"}]
    if not roots:
        roots = [fid for fid, item in features.items() if not item["parents"] and item["feature"] not in {"exon", "mrna", "transcript"}]
    bundles: list[dict[str, list[dict[str, str]]]] = []
    inferred_type = infer_event_type(p, event_type)
    for root_id in roots:
        root = features[root_id]
        transcript_ids = [cid for cid in children.get(root_id, []) if features.get(cid, {}).get("feature") in {"mrna", "transcript"}]
        if not transcript_ids:
            transcript_ids = [fid for fid, item in features.items() if item["feature"] in {"mrna", "transcript"} and root_id in item["parents"]]
        paths: list[dict[str, Any]] = []
        for idx, txid in enumerate(transcript_ids, start=1):
            tx = features[txid]
            exon_ids = [cid for cid in children.get(txid, []) if features.get(cid, {}).get("feature") == "exon"]
            exons = [(features[eid]["start"], features[eid]["end"]) for eid in exon_ids if features[eid]["start"] > 0]
            paths.append({
                "path_id": txid,
                "path_role": tx["attrs"].get("role") or tx["attrs"].get("type") or f"ISOFORM_{idx}",
                "reference_transcript_id": tx["attrs"].get("transcript_id", ""),
                "exons": exons,
            })
        if not paths:
            continue
        attrs = root["attrs"]
        gene = attrs.get("gene_name") or attrs.get("Name") or attrs.get("gene") or ""
        gene_id = attrs.get("gene_id") or root_id
        record_id = source_record_id("SplAdder", p, root["line_no"], {"root": root_id, "paths": paths})
        bundles.append(_build_event_bundle(
            sample_id=sample_id, genome_build=build, chrom=root["chrom"], strand=root["strand"],
            gene=gene, gene_id=gene_id, event_type=inferred_type, source_event_id=root_id,
            paths=paths, source_path=p, source_record=record_id, tool_version=source_tool_version,
        ))
    return _merge_bundles(bundles)


def _parse_exon_pairs(row: dict[str, str]) -> list[tuple[int, int]]:
    pairs: list[tuple[int, int]] = []
    # Headered forms: exon1_start/exon1_end, e1_start/e1_end, etc.
    lower = {k.casefold(): k for k in row}
    for idx in range(1, 20):
        starts = [f"exon{idx}_start", f"exon_{idx}_start", f"e{idx}_start"]
        ends = [f"exon{idx}_end", f"exon_{idx}_end", f"e{idx}_end"]
        s = next((as_int(row[lower[x]]) for x in starts if x in lower), 0)
        e = next((as_int(row[lower[x]]) for x in ends if x in lower), 0)
        if s and e:
            pairs.append((min(s, e), max(s, e)))
    if pairs:
        return pairs
    raw = get(row, "exons", "exon_chain", "event_coordinates")
    for a, b in re.findall(r"(\d+)\D+(\d+)", raw):
        pairs.append((min(int(a), int(b)), max(int(a), int(b))))
    return pairs


def _read_spladder_txt_rows(path: Path) -> list[dict[str, str]]:
    """Read both headered test outputs and headerless build-mode event TXT."""
    with open_text_maybe_gz(path) as fh:
        lines = [line.rstrip("\n") for line in fh if line.strip() and not line.startswith("#")]
    if not lines:
        return []
    first_cells = lines[0].split("\t")
    known_headers = {"event_id", "chrm", "chrom", "chr", "strand", "gene_name", "exon_pos"}
    if any(cell.strip().casefold() in known_headers for cell in first_cells):
        return read_delimited(path)
    rows: list[dict[str, str]] = []
    for line_no, line in enumerate(lines, start=1):
        cells = line.split("\t")
        if len(cells) < 6:
            continue
        chrom, strand, event_id, gene = cells[:4]
        coordinates: list[str] = []
        features: list[str] = []
        feature_started = False
        for cell in cells[4:]:
            token = cell.strip()
            if not feature_started and re.fullmatch(r"-?\d+", token):
                coordinates.append(token)
            else:
                feature_started = True
                features.append(token)
        exon_pairs = [f"{coordinates[i]}-{coordinates[i+1]}" for i in range(0, len(coordinates) - 1, 2)]
        rows.append({
            "chrom": chrom, "strand": strand, "event_id": event_id, "gene": gene,
            "exons": ";".join(exon_pairs), "sample_features": ";".join(features),
            "_source_line_number": str(line_no),
        })
    return rows


def parse_spladder_txt(
    path: str | Path,
    *,
    sample_id: str,
    genome_build: str = "GRCh38",
    event_type: str = "",
    source_tool_version: str = "UNASSESSED",
) -> dict[str, list[dict[str, str]]]:
    p = Path(path)
    build = normalize_genome_build(genome_build)
    rows = _read_spladder_txt_rows(p)
    bundles: list[dict[str, list[dict[str, str]]]] = []
    inferred_type = infer_event_type(p, event_type)
    for row_no, row in enumerate(rows, start=2):
        chrom = normalize_chromosome(get(row, "chrom", "chr", "gene_chr"))
        strand = normalize_strand(get(row, "strand", "gene_strand", default="."))
        source_event = get(row, "event_id", "id", "event", default=f"row_{row_no}")
        gene = get(row, "gene", "gene_name")
        gene_id = get(row, "gene_id")
        exons = _parse_exon_pairs(row)
        if not chrom or not exons:
            continue
        # When only one path is available, retain it as an event-local hypothesis.
        paths = [{"path_id": source_event, "path_role": "OBSERVED_PATH", "exons": exons}]
        record_id = source_record_id("SplAdder", p, row_no, row)
        bundle = _build_event_bundle(
            sample_id=sample_id, genome_build=build, chrom=chrom, strand=strand,
            gene=gene, gene_id=gene_id, event_type=inferred_type, source_event_id=source_event,
            paths=paths, source_path=p, source_record=record_id, tool_version=source_tool_version,
            psi=as_float_text(get(row, "psi", "event_psi")),
        )
        bundles.append(bundle)
    return _merge_bundles(bundles)
