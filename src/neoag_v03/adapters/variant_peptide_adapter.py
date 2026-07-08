"""Convert extract-variant-peptides TSV rows to standard raw_events/raw_peptides."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

from ..config import load_profile
from ..input_router import resolve_entry_mode
from ..model_layers import enrich_event_layers, enrich_peptide_layers, infer_mutation_source, infer_peptide_consequence
from ..schemas import EVENT_FIELDS, PEPTIDE_FIELDS
from ..utils import read_tsv, safe_id, write_tsv


def _event_type_from_row(row: dict[str, str]) -> str:
    cons = str(row.get("consequence") or "").lower()
    flag = str(row.get("multi_aa_flag") or "").lower()
    if flag == "fusion_neo" or "fusion" in cons:
        return "Fusion"
    if flag == "frameshift" or "frameshift" in cons:
        return "InDel"
    if "inframe" in cons or flag == "inframe_multi":
        return "InDel"
    if "missense" in cons or flag in {"single_aa", "multi_aa_substitution"}:
        return "SNV"
    return "SNV"


def _mhc_class_for_hla(hla: str) -> str:
    h = str(hla or "").upper()
    return "II" if any(x in h for x in ("DR", "DQ", "DP")) else "I"


def event_from_variant_row(row: dict[str, str], sample_id: str, profile_name: str) -> dict[str, str]:
    gene = row.get("gene") or "UNKNOWN"
    event_type = _event_type_from_row(row)
    consequence = row.get("consequence") or ""
    event_name = row.get("hgvsp") or row.get("variant_key") or consequence
    event_id = row.get("variant_key") or safe_id(
        f"{sample_id}_{gene}_{row.get('chrom')}_{row.get('pos')}_{row.get('ref')}_{row.get('alt')}"
    )
    tool = "extract-variant-peptides"
    base = {
        "event_id": event_id,
        "sample_id": sample_id,
        "disease_profile": profile_name,
        "event_type": event_type,
        "gene": gene,
        "event_name": event_name,
        "chrom": row.get("chrom", ""),
        "pos": row.get("pos", ""),
        "ref": row.get("ref", ""),
        "alt": row.get("alt", ""),
        "transcript_id": row.get("transcript_id", ""),
        "consequence": consequence,
        "rna_junction_reads": "",
        "event_confidence": "0.7",
        "event_expression": "0.0",
        "driver_relevance": "0.0",
        "tumor_vaf": row.get("vaf") or "0.0",
        "tumor_depth": row.get("tumor_depth") or "",
        "tumor_alt_count": row.get("tumor_alt_count") or "",
        "rna_vaf": row.get("rna_vaf") or "",
        "rna_alt_reads": row.get("rna_alt_reads") or "",
        "rna_depth": row.get("rna_depth") or "",
        "clonality": "0.5",
        "persistence": "0.5",
        "tumor_specificity": "0.7",
        "source": tool,
        "mutation_source": infer_mutation_source(event_type=event_type, tool=tool, consequence=consequence),
        "peptide_consequence": infer_peptide_consequence(event_type=event_type, consequence=consequence, tool=tool),
    }
    return enrich_event_layers(base)


def peptide_from_variant_row(
    row: dict[str, str],
    sample_id: str,
    event: dict[str, str],
    hla: str,
) -> dict[str, str]:
    peptide = row.get("mutant_peptide") or ""
    wt = row.get("wildtype_peptide") or ""
    pid = safe_id(f"{event['event_id']}_{hla}_{peptide}")
    base = {
        "peptide_id": pid,
        "event_id": event["event_id"],
        "sample_id": sample_id,
        "event_type": event["event_type"],
        "mutation_source": event.get("mutation_source", ""),
        "peptide_consequence": event.get("peptide_consequence", ""),
        "gene": event["gene"],
        "peptide": peptide,
        "wildtype_peptide": wt,
        "crosses_junction": "",
        "contains_novel_aa": "",
        "rna_junction_reads": "",
        "hla_allele": hla,
        "mhc_class": _mhc_class_for_hla(hla),
        "source_tool": "extract-variant-peptides",
        "binding_rank": "99",
        "el_rank": "99",
        "presentation_score": "0.0",
        "immunogenicity_score": "0.5",
        "wildtype_binding_rank": "99",
        "self_similarity_score": "0.0",
        "normal_hla_ligand_overlap": "no",
    }
    return enrich_peptide_layers(base, event)


def variant_peptide_rows_to_raw_tables(
    rows: list[dict[str, str]],
    *,
    sample_id: str,
    profile_name: str,
    hla_alleles: list[str],
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    if not hla_alleles:
        raise ValueError("variant_peptide_rows_to_raw_tables requires at least one HLA allele")

    events: dict[str, dict[str, str]] = {}
    peptides: list[dict[str, str]] = []
    for row in rows:
        ev = event_from_variant_row(row, sample_id, profile_name)
        if ev["event_id"] not in events:
            events[ev["event_id"]] = ev
        event = events[ev["event_id"]]
        for hla in hla_alleles:
            peptides.append(peptide_from_variant_row(row, sample_id, event, hla))
    return list(events.values()), peptides


def resolve_variant_peptide_options(cfg: dict[str, Any], *, root: Path | None = None) -> dict[str, Any]:
    from ..vep.extract_peptides import parse_peptide_lengths

    inputs = cfg.get("inputs") or {}
    lengths = parse_peptide_lengths(
        str(inputs.get("variant_peptide_lengths") or inputs.get("peptide_lengths") or ""),
        length_min=inputs.get("variant_peptide_length_min"),
        length_max=inputs.get("variant_peptide_length_max"),
    )
    normal_fasta = (
        inputs.get("normal_proteome_fasta")
        or os.environ.get("NEOAG_NORMAL_PROTEOME_FASTA")
        or ""
    )
    annotate_only = bool(inputs.get("variant_peptide_annotate_normal_only"))
    filter_normal = inputs.get("variant_peptide_filter_normal_proteome")
    if filter_normal is None:
        filter_normal = bool(normal_fasta) and not annotate_only
    return {
        "lengths": lengths,
        "mini_len": int(inputs.get("variant_peptide_mini_len") or inputs.get("mini_len") or 10),
        "minigene_total_len": (
            int(inputs.get("variant_peptide_minigene_total_len") or inputs.get("minigene_total_len"))
            if (inputs.get("variant_peptide_minigene_total_len") or inputs.get("minigene_total_len"))
            else None
        ),
        "exclude_multi_aa": bool(inputs.get("variant_peptide_exclude_multi_aa")),
        "single_aa_only": bool(inputs.get("variant_peptide_single_aa_only")),
        "normal_proteome_fasta": str(normal_fasta) if normal_fasta else None,
        "filter_normal_proteome": bool(filter_normal),
        "annotate_normal_proteome_only": annotate_only,
    }


def variant_peptide_extraction_enabled(cfg: dict[str, Any], variants_vcf: Path | None) -> bool:
    inputs = cfg.get("inputs") or {}
    explicit = inputs.get("variant_peptide_extraction")
    if explicit is not None:
        return bool(explicit)
    if inputs.get("pvac_files"):
        return False
    if inputs.get("raw_peptides") or inputs.get("raw_events"):
        return False
    mode = resolve_entry_mode(cfg)
    if mode not in {"snv_indel", "e2e"}:
        return False
    return variants_vcf is not None and variants_vcf.is_file()


def _catalog_rows_to_raw_peptides(
    rows: list[dict[str, str]],
    events_by_id: dict[str, dict[str, str]],
    *,
    sample_id: str,
    hla_alleles: list[str],
    source_tool: str = "extract-variant-peptides",
) -> list[dict[str, str]]:
    peptides: list[dict[str, str]] = []
    for row in rows:
        event = events_by_id.get(row.get("variant_key", ""))
        if event is None:
            continue
        for hla in hla_alleles:
            pep = peptide_from_variant_row(row, sample_id, event, hla)
            pep["source_tool"] = source_tool
            if event.get("event_type") == "Fusion":
                pep["crosses_junction"] = "yes"
            peptides.append(pep)
    return peptides


def _merge_event_catalog(
    primary: list[dict[str, str]],
    extra: list[dict[str, str]],
) -> list[dict[str, str]]:
    merged: dict[str, dict[str, str]] = {row["event_id"]: row for row in primary if row.get("event_id")}
    for row in extra:
        eid = row.get("event_id")
        if eid and eid not in merged:
            merged[eid] = row
    return list(merged.values())


def _coerce_path_list(value: Any) -> list[Path]:
    if not value:
        return []
    if isinstance(value, (str, Path)):
        text = str(value)
        parts = [p.strip() for p in text.replace(";", ",").split(",") if p.strip()]
        return [Path(p) for p in parts]
    if isinstance(value, (list, tuple)):
        return [Path(str(p)) for p in value if str(p).strip()]
    return []


def _read_extra_variant_peptide_rows(inputs: dict[str, Any], tools_dir: Path) -> tuple[list[dict[str, str]], list[str]]:
    paths: list[Path] = []
    for key in (
        "extra_variant_peptide_tsvs",
        "extra_variant_peptides_tsv",
        "arriba_fusion_peptides_tsv",
        "diagnostic_fusion_rescue_peptides_tsv",
    ):
        paths.extend(_coerce_path_list(inputs.get(key)))

    rows: list[dict[str, str]] = []
    used: list[str] = []
    for path in paths:
        if not path.is_file():
            continue
        chunk = read_tsv(path)
        if not chunk:
            continue
        rows.extend(chunk)
        used.append(str(path))

    if rows:
        sidecar = tools_dir / "extra_variant_peptides.tsv"
        fields: list[str] = []
        for row in rows:
            for key in row:
                if key not in fields:
                    fields.append(key)
        write_tsv(sidecar, rows, fields)
    return rows, used


def run_variant_peptide_upstream(
    cfg: dict[str, Any],
    *,
    variants_vcf: Path,
    parsed_dir: Path,
    tools_dir: Path,
    sample_id: str,
    profile_name: str,
    hla_alleles: list[str],
) -> dict[str, str]:
    from ..vep.extract_peptides import OUTPUT_FIELDS, extract_variant_peptides_from_vcf
    from .easyfuse_variant_peptide import (
        build_easyfuse_catalog,
        resolve_easyfuse_filter_config,
        resolve_easyfuse_peptide_config,
        write_easyfuse_qc_tables,
    )

    opts = resolve_variant_peptide_options(cfg)
    inputs = cfg.get("inputs") or {}
    tumor_sample_name = inputs.get("tumor_sample_name") or sample_id
    rna_sample_name = inputs.get("rna_sample_name") or ""
    variant_tsv = tools_dir / "variant_peptides.tsv"
    easyfuse_path = inputs.get("easyfuse_tsv") or inputs.get("easyfuse_pass_csv")
    easyfuse_file = Path(str(easyfuse_path)) if easyfuse_path else None
    if easyfuse_file and not easyfuse_file.is_file():
        easyfuse_file = None

    summary: dict[str, Any] = {}
    vcf_rows: list[dict[str, str]] = []
    if variants_vcf.is_file():
        summary = extract_variant_peptides_from_vcf(
            variants_vcf,
            variant_tsv,
            lengths=opts["lengths"],
            sample_id=sample_id,
            exclude_multi_aa=opts["exclude_multi_aa"],
            single_aa_only=opts["single_aa_only"],
            mini_len=opts["mini_len"],
            minigene_total_len=opts["minigene_total_len"],
            normal_proteome_fasta=opts["normal_proteome_fasta"],
            filter_normal_proteome=opts["filter_normal_proteome"],
            hla_alleles=hla_alleles,
            tumor_sample_name=str(tumor_sample_name) if tumor_sample_name else None,
            rna_sample_name=str(rna_sample_name) if rna_sample_name else None,
        )
        vcf_rows = read_tsv(variant_tsv)

    ef_rows: list[dict[str, str]] = []
    ef_events: list[dict[str, str]] = []
    ef_qc_paths: dict[str, str] = {}
    if easyfuse_file:
        ef_result = build_easyfuse_catalog(
            easyfuse_file,
            sample_id,
            profile_name,
            lengths=opts["lengths"],
            filter_cfg=resolve_easyfuse_filter_config(cfg),
            peptide_cfg=resolve_easyfuse_peptide_config(cfg),
        )
        ef_rows = ef_result.catalog_rows
        ef_events = ef_result.events
        ef_qc_paths = write_easyfuse_qc_tables(ef_result, tools_dir)
        ef_sidecar = tools_dir / "easyfuse_variant_peptides.tsv"
        write_tsv(ef_sidecar, ef_rows, list(ef_rows[0].keys()) if ef_rows else [])

    extra_rows, extra_sources = _read_extra_variant_peptide_rows(inputs, tools_dir)

    catalog_rows = vcf_rows + ef_rows + extra_rows
    if not catalog_rows:
        raise ValueError(
            f"variant_peptide_extraction produced no peptides from {variants_vcf}"
            f"{f' and {easyfuse_file}' if easyfuse_file else ''}"
        )

    write_tsv(variant_tsv, catalog_rows, OUTPUT_FIELDS)

    events, peptides = variant_peptide_rows_to_raw_tables(
        vcf_rows,
        sample_id=sample_id,
        profile_name=profile_name,
        hla_alleles=hla_alleles,
    )
    events = _merge_event_catalog(events, ef_events)
    ef_event_map = {e["event_id"]: e for e in ef_events}
    peptides = peptides + _catalog_rows_to_raw_peptides(
        ef_rows,
        ef_event_map,
        sample_id=sample_id,
        hla_alleles=hla_alleles,
        source_tool="EasyFuse",
    )
    extra_events, _ = variant_peptide_rows_to_raw_tables(
        extra_rows,
        sample_id=sample_id,
        profile_name=profile_name,
        hla_alleles=hla_alleles,
    ) if extra_rows else ([], [])
    events = _merge_event_catalog(events, extra_events)
    extra_event_map = {e["event_id"]: e for e in extra_events}
    peptides = peptides + _catalog_rows_to_raw_peptides(
        extra_rows,
        extra_event_map,
        sample_id=sample_id,
        hla_alleles=hla_alleles,
        source_tool="extra-fusion-peptides",
    )

    if not peptides:
        raise ValueError(
            f"variant_peptide_extraction produced no peptide/HLA pairs from {variants_vcf}"
            f"{f' and {easyfuse_file}' if easyfuse_file else ''} "
            f"(catalog_rows={len(catalog_rows)})"
        )

    raw_events = parsed_dir / "raw_events.tsv"
    raw_peptides = parsed_dir / "raw_peptides.tsv"
    write_tsv(raw_events, events, EVENT_FIELDS)
    write_tsv(raw_peptides, peptides, PEPTIDE_FIELDS)
    result = {
        "variant_peptides": str(variant_tsv),
        "easyfuse_variant_peptides": str(tools_dir / "easyfuse_variant_peptides.tsv") if ef_rows else "",
        "raw_events": str(raw_events),
        "raw_peptides": str(raw_peptides),
        "variant_peptide_rows": str(len(catalog_rows)),
        "variant_peptide_rows_vcf": str(len(vcf_rows)),
        "variant_peptide_rows_easyfuse": str(len(ef_rows)),
        "variant_peptide_rows_extra": str(len(extra_rows)),
        "extra_variant_peptides": str(tools_dir / "extra_variant_peptides.tsv") if extra_rows else "",
        "extra_variant_peptide_sources": ",".join(extra_sources),
        "variant_peptides_filtered_normal": str(summary.get("peptides_filtered_normal_proteome", 0)),
        "sample_hla_alleles": summary.get("sample_hla_alleles", ",".join(hla_alleles)),
    }
    result.update(ef_qc_paths)
    return result
