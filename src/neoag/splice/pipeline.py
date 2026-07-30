"""Formal v0.5.0 Splice Provenance Layer orchestration."""
from __future__ import annotations

import json
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

from neoag.schemas import EVENT_FIELDS as LEGACY_EVENT_FIELDS, PEPTIDE_FIELDS as LEGACY_PEPTIDE_FIELDS, RNA_JUNCTION_EVIDENCE_FIELDS
from neoag.splice.coordinates import file_sha256
from neoag.splice.identifiers import link_id, splice_event_id, stable_digest, stable_id
from neoag.utils import write_json, write_tsv

from .adapters.immunopepper import parse_immunopepper_kmers, parse_immunopepper_meta
from .adapters.high_order import parse_high_order_evidence
from .adapters.irfinder import parse_irfinder
from .adapters.pvacbind import parse_pvacbind
from .adapters.regtools import parse_junction_source
from .adapters.spladder import parse_spladder_gff3, parse_spladder_txt
from .consensus import build_consensus, consensus_reason_conflicts
from .normal_background import parse_normal_coverage, parse_normal_junctions
from .projection import project_legacy
from .schemas import (
    CONFLICT_FIELDS, OUTPUT_FILENAMES, PVACBIND_FASTA_MAP_FIELDS,
    SPLICE_PROVENANCE_SCHEMA_VERSION, TABLE_FIELDS, TOOL_EVIDENCE_FIELDS,
)

_ID_FIELDS = {
    "junctions": "junction_id", "events": "splice_event_id",
    "event_junction_links": "event_junction_link_id", "transcripts": "transcript_hypothesis_id",
    "orfs": "orf_id", "peptide_origins": "origin_peptide_id",
    "peptide_origin_links": "peptide_origin_link_id", "presentation": "presentation_id",
    "normal_background": "normal_background_id", "tool_evidence": "evidence_id",
    "consensus": "consensus_id", "conflicts": "conflict_id",
}
_SET_FIELDS = {
    "source_tools", "source_tool_versions", "source_files", "source_record_ids",
    "junction_ids", "reference_junction_ids", "alternative_junction_ids", "affected_exons",
}
_MAX_FIELDS = {"unique_split_reads", "multi_split_reads", "total_split_reads", "max_overhang"}
_IDENTITY_FIELDS = {
    "genome_build", "chrom", "intron_start_1based", "intron_end_1based", "strand",
    "donor_1based", "acceptor_1based", "event_type", "junction_ids", "gene", "gene_id",
    "splice_event_id", "transcript_hypothesis_id", "orf_id", "protein_sequence",
    "protein_sequence_sha256", "frame_status", "peptide_sequence",
}


def _tokens(value: Any) -> set[str]:
    return {x.strip() for x in str(value or "").replace(",", ";").split(";") if x.strip()}


def _best(values: set[str]) -> str:
    non_missing = [x for x in values if x and x.upper() not in {"UNASSESSED", "UNRESOLVED", "NA", "NONE"}]
    return sorted(non_missing or [x for x in values if x])[0] if (non_missing or values) else ""


def _numeric_max(values: set[str]) -> str:
    parsed: list[float] = []
    for value in values:
        try:
            parsed.append(float(value))
        except Exception:
            pass
    if not parsed:
        return ""
    maximum = max(parsed)
    return str(int(maximum)) if maximum.is_integer() else f"{maximum:.12g}"


@dataclass
class SpliceLayer:
    sample_id: str
    genome_build: str = "GRCh38"
    disease_profile: str = "default"
    tables: dict[str, list[dict[str, str]]] = field(default_factory=lambda: defaultdict(list))
    input_files: list[dict[str, str]] = field(default_factory=list)

    def extend(self, bundle: dict[str, list[dict[str, str]]] | None) -> None:
        if not bundle:
            return
        for key, rows in bundle.items():
            if key in {"pvacbind_fasta", "manifest"}:
                continue
            self.tables.setdefault(key, []).extend(dict(row) for row in rows)

    def register_input(self, path: str | Path, *, role: str, tool: str, version: str = "UNASSESSED") -> None:
        p = Path(path)
        self.input_files.append({
            "path": str(p), "role": role, "tool": tool, "version": version,
            "sha256": file_sha256(p) if p.is_file() else "MISSING",
        })

    def _ensure_ids(self) -> None:
        for table, rows in self.tables.items():
            id_field = _ID_FIELDS.get(table)
            if not id_field:
                continue
            for idx, row in enumerate(rows, start=1):
                if row.get(id_field):
                    continue
                row[id_field] = stable_id(table[:3].upper(), self.sample_id, idx, row)

    def _merge_table(self, table: str, rows: list[dict[str, str]]) -> list[dict[str, str]]:
        id_field = _ID_FIELDS.get(table)
        if not id_field:
            return rows
        grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in rows:
            grouped[row[id_field]].append(row)
        merged: list[dict[str, str]] = []
        for entity_id, group in sorted(grouped.items()):
            fields = TABLE_FIELDS.get(table, [])
            output: dict[str, str] = {id_field: entity_id}
            for field_name in fields:
                if field_name == id_field:
                    continue
                values = {str(row.get(field_name, "")).strip() for row in group if str(row.get(field_name, "")).strip()}
                if field_name in _SET_FIELDS:
                    output[field_name] = ";".join(sorted({token for value in values for token in _tokens(value)}))
                elif field_name in _MAX_FIELDS:
                    output[field_name] = _numeric_max(values)
                elif field_name == "provenance_record_count":
                    records = {token for row in group for token in _tokens(row.get("source_record_ids", ""))}
                    output[field_name] = str(len(records) or len(group))
                else:
                    output[field_name] = _best(values)
                if field_name in _IDENTITY_FIELDS and len(values) > 1:
                    conflict = {
                        "entity_type": table.upper(), "entity_id": entity_id, "sample_id": self.sample_id,
                        "conflict_type": "IDENTITY_FIELD_CONFLICT", "field_name": field_name,
                        "observed_values": ";".join(sorted(values)),
                        "source_tools": ";".join(sorted({x for row in group for x in _tokens(row.get("source_tools") or row.get("source_generator") or row.get("source_tool"))})),
                        "source_record_ids": ";".join(sorted({x for row in group for x in _tokens(row.get("source_record_ids") or row.get("source_record_id"))})),
                        "severity": "ERROR", "resolution_status": "DETERMINISTIC_VALUE_SELECTED",
                        "resolution_reason": f"Exact {id_field} matched, but {field_name} differed; selected deterministic non-missing value.",
                    }
                    conflict["conflict_id"] = stable_id("CFL", conflict)
                    self.tables.setdefault("conflicts", []).append(conflict)
                    output["evidence_conflict_status"] = "IDENTITY_FIELD_CONFLICT"
            merged.append(output)
        return merged

    def consolidate(self) -> None:
        self._ensure_ids()
        # Merge entity tables only by their exact stable identifiers.
        for table in [
            "junctions", "events", "event_junction_links", "transcripts", "orfs",
            "peptide_origins", "peptide_origin_links", "presentation", "normal_background",
            "tool_evidence", "conflicts",
        ]:
            self.tables[table] = self._merge_table(table, self.tables.get(table, []))
        self._add_fallback_events()
        self.tables["events"] = self._merge_table("events", self.tables.get("events", []))
        self.tables["event_junction_links"] = self._merge_table("event_junction_links", self.tables.get("event_junction_links", []))

    def _add_fallback_events(self) -> None:
        linked = {row.get("junction_id", "") for row in self.tables.get("event_junction_links", [])}
        existing_events = {row.get("splice_event_id", "") for row in self.tables.get("events", [])}
        for junction in self.tables.get("junctions", []):
            jid = junction.get("junction_id", "")
            if not jid or jid in linked:
                continue
            event_id = splice_event_id(
                genome_build=junction.get("genome_build", self.genome_build), event_type="NOVEL_JUNCTION",
                strand=junction.get("strand", "."), junction_ids=[jid], gene="",
            )
            if event_id not in existing_events:
                self.tables["events"].append({
                    "splice_event_id": event_id, "sample_id": self.sample_id,
                    "genome_build": junction.get("genome_build", self.genome_build),
                    "event_type": "NOVEL_JUNCTION", "gene": "", "gene_id": "",
                    "strand": junction.get("strand", "."), "junction_ids": jid,
                    "reference_junction_ids": "", "alternative_junction_ids": jid, "affected_exons": "",
                    "annotation_status": "JUNCTION_ONLY", "cryptic_exon_status": "UNASSESSED",
                    "psi": "", "delta_psi": "", "qvalue": "", "outlier_score": "",
                    "event_expression": "", "event_confidence": "JUNCTION_ONLY",
                    "reference_path_status": "UNRESOLVED", "cohort_analysis_status": "NOT_APPLICABLE",
                    "source_tools": junction.get("source_tools", ""),
                    "source_tool_versions": junction.get("source_tool_versions", ""),
                    "source_files": junction.get("source_files", ""),
                    "source_record_ids": junction.get("source_record_ids", ""),
                    "provenance_record_count": junction.get("provenance_record_count", "1"),
                    "event_resolution_status": "RESOLVED_JUNCTION_ONLY", "evidence_conflict_status": "NONE",
                })
                existing_events.add(event_id)
            self.tables["event_junction_links"].append({
                "event_junction_link_id": link_id("EJL", event_id, jid, "junction_only"),
                "splice_event_id": event_id, "junction_id": jid, "sample_id": self.sample_id,
                "path_id": "junction_only", "path_role": "OBSERVED_JUNCTION", "edge_index": "1",
                "junction_role": "ALTERNATIVE_EDGE", "source_tool": junction.get("source_tools", ""),
                "source_record_id": junction.get("source_record_ids", ""), "link_status": "RESOLVED",
            })
            linked.add(jid)

    def write_pvacbind_fasta(self, outdir: str | Path) -> tuple[Path, Path]:
        out = Path(outdir)
        out.mkdir(parents=True, exist_ok=True)
        fasta_path = out / OUTPUT_FILENAMES["pvacbind_fasta"]
        map_path = out / OUTPUT_FILENAMES["pvacbind_fasta_map"]
        lines: list[str] = []
        map_rows: list[dict[str, str]] = []
        for orf in sorted(self.tables.get("orfs", []), key=lambda x: x.get("orf_id", "")):
            sequence = "".join(orf.get("protein_sequence", "").split()).upper()
            if not sequence or orf.get("orf_validity_status") in {"INVALID", "UNRESOLVED"}:
                continue
            index = f"SPORF_{stable_digest(orf.get('orf_id'), sequence, length=20)}"
            lines.append(f">{index}")
            lines.extend(sequence[i:i+60] for i in range(0, len(sequence), 60))
            map_rows.append({
                "index": index, "orf_id": orf.get("orf_id", ""),
                "transcript_hypothesis_id": orf.get("transcript_hypothesis_id", ""),
                "splice_event_id": orf.get("splice_event_id", ""), "sample_id": self.sample_id,
                "gene": orf.get("gene", ""), "sequence_type": "ALTERED_PROTEIN_OR_PARTIAL_ORF",
                "sequence_sha256": orf.get("protein_sequence_sha256", ""),
                "source_generator": orf.get("source_generator", ""),
                "source_record_id": orf.get("source_record_id", ""),
            })
        fasta_path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
        self.tables["pvacbind_fasta_map"] = map_rows
        write_tsv(map_path, map_rows, PVACBIND_FASTA_MAP_FIELDS)
        return fasta_path, map_path

    def _referential_qc(self) -> list[dict[str, str]]:
        events = {x.get("splice_event_id", "") for x in self.tables.get("events", [])}
        junctions = {x.get("junction_id", "") for x in self.tables.get("junctions", [])}
        transcripts = {x.get("transcript_hypothesis_id", "") for x in self.tables.get("transcripts", [])}
        orfs = {x.get("orf_id", "") for x in self.tables.get("orfs", [])}
        origins = {x.get("origin_peptide_id", "") for x in self.tables.get("peptide_origins", [])}
        unstranded_junctions = sum(
            1 for row in self.tables.get("junctions", []) if row.get("strand") not in {"+", "-"}
        )
        checks = {
            "event_junction_links_missing_event": sum(1 for r in self.tables.get("event_junction_links", []) if r.get("splice_event_id") not in events),
            "event_junction_links_missing_junction": sum(1 for r in self.tables.get("event_junction_links", []) if r.get("junction_id") not in junctions),
            "transcripts_missing_event": sum(1 for r in self.tables.get("transcripts", []) if r.get("splice_event_id") not in events),
            "orfs_missing_transcript": sum(1 for r in self.tables.get("orfs", []) if r.get("transcript_hypothesis_id") not in transcripts),
            "orfs_missing_event": sum(1 for r in self.tables.get("orfs", []) if r.get("splice_event_id") not in events),
            "origins_missing_orf": sum(1 for r in self.tables.get("peptide_origins", []) if r.get("orf_id") not in orfs),
            "origins_missing_event": sum(1 for r in self.tables.get("peptide_origins", []) if r.get("splice_event_id") not in events),
            "presentation_missing_origin": sum(1 for r in self.tables.get("presentation", []) if r.get("origin_peptide_id") not in origins),
        }
        rows = [{"metric": key, "value": str(value), "status": "PASS" if value == 0 else "FAIL", "detail": "Exact-ID referential integrity check."} for key, value in checks.items()]
        rows.extend([
            {
                "metric": "unstranded_canonical_junctions", "value": str(unstranded_junctions),
                "status": "PASS" if unstranded_junctions == 0 else "FAIL",
                "detail": "Canonical junctions require an explicit + or - strand for production use.",
            },
            {"metric": "junction_count", "value": str(len(junctions)), "status": "INFO", "detail": "Canonical junction entities."},
            {"metric": "splice_event_count", "value": str(len(events)), "status": "INFO", "detail": "Biological splice event entities."},
            {"metric": "transcript_hypothesis_count", "value": str(len(transcripts)), "status": "INFO", "detail": "Transcript hypotheses."},
            {"metric": "orf_count", "value": str(len(orfs)), "status": "INFO", "detail": "ORF/translated segment entities."},
            {"metric": "peptide_origin_count", "value": str(len(origins)), "status": "INFO", "detail": "Peptide origins."},
            {"metric": "gene_or_nearest_locus_fallbacks_used", "value": "0", "status": "PASS", "detail": "Approximate provenance linking is forbidden."},
        ])
        return rows

    def finalise(self) -> None:
        self.consolidate()
        self.tables["consensus"] = build_consensus(self.tables, sample_id=self.sample_id)
        self.tables.setdefault("conflicts", []).extend(
            consensus_reason_conflicts(self.tables["consensus"], sample_id=self.sample_id)
        )
        self.tables["conflicts"] = self._merge_table("conflicts", self.tables.get("conflicts", []))
        projection = project_legacy(self.tables, sample_id=self.sample_id, disease_profile=self.disease_profile)
        self.tables.update(projection)
        self.tables["qc"] = self._referential_qc()

    def write(self, outdir: str | Path, *, strict: bool = False) -> dict[str, Path]:
        out = Path(outdir)
        out.mkdir(parents=True, exist_ok=True)
        outputs: dict[str, Path] = {}
        for table, filename in OUTPUT_FILENAMES.items():
            if table in {"manifest", "pvacbind_fasta"}:
                continue
            path = out / filename
            if table == "raw_events":
                fields = LEGACY_EVENT_FIELDS
            elif table == "raw_peptides":
                fields = LEGACY_PEPTIDE_FIELDS
            elif table == "rna_junction_evidence":
                fields = RNA_JUNCTION_EVIDENCE_FIELDS
            else:
                fields = TABLE_FIELDS.get(table)
            write_tsv(path, self.tables.get(table, []), fields)
            outputs[table] = path
        fasta, fasta_map = self.write_pvacbind_fasta(out)
        outputs["pvacbind_fasta"] = fasta
        outputs["pvacbind_fasta_map"] = fasta_map
        manifest = {
            "schema_version": SPLICE_PROVENANCE_SCHEMA_VERSION,
            "software_version": "0.5.0", "sample_id": self.sample_id,
            "genome_build": self.genome_build, "disease_profile": self.disease_profile,
            "matching_policy": {
                "junction_identity": "genome_build+chrom+1based_closed_intron_start+intron_end+strand",
                "approximate_gene_or_nearest_locus_fallback": False,
                "pvacbind_mapping": "exact_fasta_index_only",
                "negative_normal_evidence_requires_coverage": True,
            },
            "inputs": self.input_files,
            "outputs": {key: {"path": path.name, "sha256": file_sha256(path)} for key, path in outputs.items() if path.is_file()},
            "counts": {key: len(value) for key, value in self.tables.items() if isinstance(value, list)},
        }
        manifest_path = out / OUTPUT_FILENAMES["manifest"]
        write_json(manifest_path, manifest)
        outputs["manifest"] = manifest_path
        failures = [row for row in self.tables.get("qc", []) if row.get("status") == "FAIL"]
        error_conflicts = [row for row in self.tables.get("conflicts", []) if row.get("severity") == "ERROR" and row.get("resolution_status") == "UNRESOLVED"]
        if strict and (failures or error_conflicts):
            raise ValueError(f"Splice Provenance Layer strict validation failed: qc_failures={len(failures)}, unresolved_errors={len(error_conflicts)}")
        return outputs


def _as_paths(values: Iterable[str | Path] | None) -> list[Path]:
    return [Path(x) for x in (values or [])]


def _missing_locked_versions(versions: dict[str, str], required_tools: Iterable[str]) -> list[str]:
    missing_values = {"", "UNASSESSED", "UNKNOWN", "NA", "N/A", "NONE"}
    return sorted({
        tool for tool in required_tools
        if str(versions.get(tool, "")).strip().upper() in missing_values
    })


def build_splice_provenance_layer(
    *,
    sample_id: str,
    outdir: str | Path,
    genome_build: str = "GRCh38",
    disease_profile: str = "default",
    junctions: str | Path | None = None,
    junction_coordinate_system: str = "auto",
    junction_source_assay_id: str = "",
    star_junctions: str | Path | None = None,
    star_junction_source_assay_id: str = "",
    spladder_gff3: Iterable[str | Path] | None = None,
    spladder_txt: Iterable[str | Path] | None = None,
    irfinder: Iterable[str | Path] | None = None,
    irfinder_coordinate_system: str = "UNSPECIFIED",
    immunopepper_meta: Iterable[str | Path] | None = None,
    immunopepper_kmers: Iterable[str | Path] | None = None,
    pvacbind: Iterable[str | Path] | None = None,
    pvacbind_fasta_map: str | Path | None = None,
    normal_junctions: Iterable[str | Path] | None = None,
    normal_coordinate_system: str = "auto",
    normal_coverage: Iterable[str | Path] | None = None,
    high_order_evidence: Iterable[str | Path] | None = None,
    tool_versions: dict[str, str] | None = None,
    strict: bool = False,
) -> dict[str, Path]:
    versions = tool_versions or {}
    if irfinder and str(irfinder_coordinate_system).strip().upper() in {"", "UNSPECIFIED", "AUTO"}:
        raise ValueError("IRFinder-S inputs require an explicit --irfinder-coordinate-system declaration")
    required_tools: list[str] = []
    if junctions:
        required_tools.append("RegTools")
    if star_junctions:
        required_tools.append("STAR")
    if spladder_gff3 or spladder_txt:
        required_tools.append("SplAdder")
    if irfinder:
        required_tools.append("IRFinder-S")
    if immunopepper_meta or immunopepper_kmers:
        required_tools.append("ImmunoPepper")
    if pvacbind:
        required_tools.append("pVACbind")
    if strict:
        missing_versions = _missing_locked_versions(versions, required_tools)
        if missing_versions:
            raise ValueError(
                "strict mode requires explicit TOOL=VERSION locks for all executed tools: "
                + ", ".join(missing_versions)
            )
    layer = SpliceLayer(sample_id=sample_id, genome_build=genome_build, disease_profile=disease_profile)
    if junctions:
        layer.register_input(junctions, role="primary_rna_junctions", tool="RegTools", version=versions.get("RegTools", "UNASSESSED"))
        layer.extend(parse_junction_source(
            junctions, sample_id=sample_id, source_tool="RegTools", source_tool_version=versions.get("RegTools", "UNASSESSED"),
            source_assay_id=junction_source_assay_id,
            genome_build=genome_build, coordinate_system=junction_coordinate_system, strict=strict,
        ))
    if star_junctions:
        layer.register_input(star_junctions, role="primary_rna_junctions", tool="STAR-SJ", version=versions.get("STAR", "UNASSESSED"))
        layer.extend(parse_junction_source(
            star_junctions, sample_id=sample_id, source_tool="STAR-SJ", source_tool_version=versions.get("STAR", "UNASSESSED"),
            source_assay_id=star_junction_source_assay_id,
            genome_build=genome_build, coordinate_system="star_sj", strict=strict,
        ))
    for path in _as_paths(spladder_gff3):
        layer.register_input(path, role="splice_event_graph", tool="SplAdder", version=versions.get("SplAdder", "UNASSESSED"))
        layer.extend(parse_spladder_gff3(path, sample_id=sample_id, genome_build=genome_build, source_tool_version=versions.get("SplAdder", "UNASSESSED")))
    for path in _as_paths(spladder_txt):
        layer.register_input(path, role="splice_event_table", tool="SplAdder", version=versions.get("SplAdder", "UNASSESSED"))
        layer.extend(parse_spladder_txt(path, sample_id=sample_id, genome_build=genome_build, source_tool_version=versions.get("SplAdder", "UNASSESSED")))
    for path in _as_paths(irfinder):
        layer.register_input(path, role="intron_retention", tool="IRFinder-S", version=versions.get("IRFinder-S", "UNASSESSED"))
        layer.extend(parse_irfinder(
            path, sample_id=sample_id, genome_build=genome_build,
            coordinate_system=irfinder_coordinate_system, source_tool_version=versions.get("IRFinder-S", "UNASSESSED"), strict=strict,
        ))
    meta_bundle: dict[str, list[dict[str, str]]] = defaultdict(list)
    for path in _as_paths(immunopepper_meta):
        layer.register_input(path, role="translated_splice_paths", tool="ImmunoPepper", version=versions.get("ImmunoPepper", "UNASSESSED"))
        bundle = parse_immunopepper_meta(path, sample_id=sample_id, genome_build=genome_build, source_tool_version=versions.get("ImmunoPepper", "UNASSESSED"))
        layer.extend(bundle)
        for key, rows in bundle.items():
            meta_bundle[key].extend(rows)
    for path in _as_paths(immunopepper_kmers):
        layer.register_input(path, role="translated_splice_kmers", tool="ImmunoPepper", version=versions.get("ImmunoPepper", "UNASSESSED"))
        layer.extend(parse_immunopepper_kmers(path, sample_id=sample_id, source_tool_version=versions.get("ImmunoPepper", "UNASSESSED"), meta_bundle=dict(meta_bundle)))
    for path in _as_paths(normal_junctions):
        layer.register_input(path, role="normal_junction_background", tool="NormalPanel")
        layer.extend(parse_normal_junctions(
            path, sample_id=sample_id, genome_build=genome_build,
            coordinate_system=normal_coordinate_system, strict=strict,
        ))
    for path in _as_paths(normal_coverage):
        layer.register_input(path, role="normal_coverage_background", tool="NormalCoverage")
        layer.extend(parse_normal_coverage(path, sample_id=sample_id))

    # Preliminary consolidation is required to create an exact ORF→FASTA index.
    layer.consolidate()
    for path in _as_paths(high_order_evidence):
        layer.register_input(path, role="high_order_validation", tool="ValidatedEvidence")
        layer.extend(parse_high_order_evidence(
            path, sample_id=sample_id, entity_bundle=layer.tables, strict=strict,
        ))
    _, generated_map = layer.write_pvacbind_fasta(outdir)
    map_for_pvac = Path(pvacbind_fasta_map) if pvacbind_fasta_map else generated_map
    for path in _as_paths(pvacbind):
        layer.register_input(path, role="presentation_prediction", tool="pVACbind", version=versions.get("pVACbind", "UNASSESSED"))
        layer.extend(parse_pvacbind(
            path, sample_id=sample_id, fasta_map=map_for_pvac,
            source_tool_version=versions.get("pVACbind", "UNASSESSED"), entity_bundle=layer.tables,
        ))
    layer.finalise()
    return layer.write(outdir, strict=strict)
