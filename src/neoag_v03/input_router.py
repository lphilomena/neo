"""Multi-entry input router: map entry modes A–F to standard raw_events/raw_peptides."""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any

from .adapters.easyfuse_adapter import is_easyfuse_table, write_fusion_evidence
from .adapters.easyfuse_variant_peptide import (
    build_easyfuse_catalog,
    resolve_easyfuse_filter_config,
    resolve_easyfuse_peptide_config,
    write_easyfuse_qc_tables,
)
from .vep.extract_peptides import parse_peptide_lengths
from .adapters.event_catalog import parse_fusion_catalog, parse_splice_catalog, write_event_catalog
from .adapters.peptide_input import build_raw_events_from_peptides, convert_peptide_input
from .adapters.pvactools_parser import parse_pvactools_outputs
from .config import load_profile
from .schemas import INPUT_MODES, EVENT_FIELDS, PEPTIDE_FIELDS
from .utils import read_tsv, write_tsv


def resolve_entry_mode(cfg: dict[str, Any]) -> str:
    mode = (cfg.get("inputs") or {}).get("entry_mode") or cfg.get("entry_mode") or "pvac"
    mode = str(mode).strip().lower()
    aliases = {
        "a": "snv_indel",
        "b": "fusion",
        "c": "splice_junction",
        "d": "sv",
        "e": "peptide_only",
        "f": "e2e",
    }
    return aliases.get(mode, mode)


def _path_or_none(val, root: Path | None = None) -> Path | None:
    if not val:
        return None
    p = Path(val)
    if not p.is_absolute() and root:
        p = (root / p).resolve()
    return p if p.is_file() else None


def _copy_if_exists(src: Path | None, dst: Path) -> bool:
    if src and src.is_file():
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        return True
    return False


def _merge_peptides(existing: list[dict[str, str]], extra: list[dict[str, str]]) -> list[dict[str, str]]:
    by_id = {p["peptide_id"]: p for p in existing if p.get("peptide_id")}
    for pep in extra:
        pid = pep.get("peptide_id")
        if pid and pid not in by_id:
            by_id[pid] = pep
    return list(by_id.values())


def _ingest_fusion_table(
    fusion_path: Path,
    sample_id: str,
    profile_name: str,
    events: list[dict[str, str]],
    peptides: list[dict[str, str]],
    fusion_evidence: list[dict[str, str]],
    *,
    cfg: dict[str, Any] | None = None,
    tools_dir: Path | None = None,
    hla_alleles: list[str] | None = None,
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    if is_easyfuse_table(fusion_path):
        cfg = cfg or {}
        inputs = cfg.get("inputs") or {}
        lengths = parse_peptide_lengths(str(inputs.get("variant_peptide_lengths", "8,9,10,11")))
        result = build_easyfuse_catalog(
            fusion_path,
            sample_id,
            profile_name,
            lengths=lengths,
            filter_cfg=resolve_easyfuse_filter_config(cfg),
            peptide_cfg=resolve_easyfuse_peptide_config(cfg),
        )
        if tools_dir is not None:
            write_easyfuse_qc_tables(result, tools_dir)
        events = _merge_events(events, result.events)
        fusion_evidence = fusion_evidence + result.fusion_evidence
        if hla_alleles:
            from .adapters.variant_peptide_adapter import _catalog_rows_to_raw_peptides

            event_map = {e["event_id"]: e for e in result.events}
            peptides = _merge_peptides(
                peptides,
                _catalog_rows_to_raw_peptides(
                    result.catalog_rows,
                    event_map,
                    sample_id=sample_id,
                    hla_alleles=hla_alleles,
                    source_tool="EasyFuse",
                ),
            )
    else:
        events = _merge_events(events, parse_fusion_catalog(fusion_path, sample_id, profile_name))
    return events, peptides, fusion_evidence


def _merge_events(existing: list[dict[str, str]], extra: list[dict[str, str]]) -> list[dict[str, str]]:
    by_id = {e["event_id"]: e for e in existing if e.get("event_id")}
    for ev in extra:
        eid = ev.get("event_id")
        if eid and eid not in by_id:
            by_id[eid] = ev
    return list(by_id.values())


def build_raw_intermediates(
    cfg: dict[str, Any],
    outdir: str | Path,
    *,
    root: Path | None = None,
) -> dict[str, str]:
    """Build or passthrough parsed/raw_events.tsv + raw_peptides.tsv for any entry mode."""
    outdir = Path(outdir)
    parsed = outdir / "parsed"
    parsed.mkdir(parents=True, exist_ok=True)

    sample = cfg.get("sample") or {}
    inputs = cfg.get("inputs") or {}
    sample_id = str(sample.get("id") or "SAMPLE001")
    profile_name = load_profile(sample.get("profile") or "default")["_profile_name"]
    mode = resolve_entry_mode(cfg)

    raw_events = parsed / "raw_events.tsv"
    raw_peptides = parsed / "raw_peptides.tsv"

    pre_events = _path_or_none(inputs.get("raw_events"), root)
    pre_peptides = _path_or_none(inputs.get("raw_peptides"), root)
    if pre_events and pre_peptides:
        shutil.copy2(pre_events, raw_events)
        shutil.copy2(pre_peptides, raw_peptides)
        return {
            "entry_mode": "intermediates",
            "raw_events": str(raw_events),
            "raw_peptides": str(raw_peptides),
        }

    if mode == "intermediates":
        raise ValueError("entry_mode=intermediates requires inputs.raw_events and inputs.raw_peptides")

    if mode == "peptide_only":
        peptide_table = _path_or_none(inputs.get("peptide_table") or inputs.get("raw_peptide_table"), root)
        if not peptide_table:
            raise ValueError("peptide_only mode requires inputs.peptide_table")
        summary = convert_peptide_input(
            peptide_table,
            outdir,
            sample_id=sample_id,
            require_hla=bool(inputs.get("require_hla", True)),
        )
        build_raw_events_from_peptides(raw_peptides, raw_events, sample_id, profile_name)
        return {
            "entry_mode": mode,
            "raw_events": str(raw_events),
            "raw_peptides": str(raw_peptides),
            "peptide_input_summary": summary.raw_peptides_tsv,
        }

    if mode == "sv":
        sv_events = _path_or_none(inputs.get("sv_raw_events"), root)
        sv_peptides = _path_or_none(inputs.get("sv_raw_peptides"), root)
        if not (sv_events and sv_peptides):
            raise ValueError(
                "sv mode requires inputs.sv_raw_events + inputs.sv_raw_peptides "
                "(from neoag-v03 sv-build-raw / sv-run-full)"
            )
        shutil.copy2(sv_events, raw_events)
        shutil.copy2(sv_peptides, raw_peptides)
        return {"entry_mode": mode, "raw_events": str(raw_events), "raw_peptides": str(raw_peptides)}

    # Modes A/B/C/F and legacy pvac: pVAC parse + optional fusion/splice catalogs
    pvac_paths: list[str] = []
    for key in ("pvac_files",):
        for p in inputs.get(key) or []:
            path = _path_or_none(p, root)
            if path:
                pvac_paths.append(str(path))

    if mode in {"e2e", "snv_indel", "fusion", "splice_junction", "pvac"}:
        pass  # pvac_paths may be filled by upstream caller

    if pvac_paths:
        parse_pvactools_outputs(pvac_paths, sample_id, profile_name, raw_events, raw_peptides)
        events = read_tsv(raw_events)
        peptides = read_tsv(raw_peptides)
    else:
        events, peptides = [], []
        write_tsv(raw_events, events, EVENT_FIELDS)
        write_tsv(raw_peptides, peptides, PEPTIDE_FIELDS)

    fusion_evidence_rows: list[dict[str, str]] = []
    tools_dir = outdir / "tools"
    tools_dir.mkdir(parents=True, exist_ok=True)
    hla_alleles = [str(a) for a in (inputs.get("hla_alleles") or []) if str(a).strip()]

    fusion_path = _path_or_none(
        inputs.get("easyfuse_tsv")
        or inputs.get("easyfuse_pass_csv")
        or inputs.get("fusion_tsv"),
        root,
    )
    if fusion_path and mode in {"fusion", "e2e", "pvac", "snv_indel"}:
        events, peptides, fusion_evidence_rows = _ingest_fusion_table(
            fusion_path,
            sample_id,
            profile_name,
            events,
            peptides,
            fusion_evidence_rows,
            cfg=cfg,
            tools_dir=tools_dir,
            hla_alleles=hla_alleles,
        )

    splice_path = _path_or_none(
        inputs.get("splice_junction_tsv") or inputs.get("regtools_tsv"),
        root,
    )
    if splice_path and mode in {"splice_junction", "e2e", "pvac", "snv_indel"}:
        from .adapters.splice_junction_adapter import merge_splice_into_catalog

        variants_vcf = _path_or_none(
            inputs.get("variants_vcf") or inputs.get("tumor_vcf"),
            root,
        )
        events, peptides = merge_splice_into_catalog(
            splice_path,
            sample_id,
            profile_name,
            events,
            peptides,
            variants_vcf=variants_vcf,
            hla_alleles=inputs.get("hla_alleles") or [],
            cfg=cfg,
            tools_dir=outdir / "tools",
        )

    if events:
        write_event_catalog(events, raw_events)
    if not raw_peptides.is_file() or peptides:
        write_tsv(raw_peptides, peptides if peptides else read_tsv(raw_peptides) if raw_peptides.is_file() else [], PEPTIDE_FIELDS)

    fusion_evidence_out = parsed / "fusion_evidence.tsv"
    if fusion_evidence_rows:
        write_fusion_evidence(fusion_evidence_rows, fusion_evidence_out)

    if not events and not raw_peptides.is_file():
        raise ValueError(
            f"entry_mode={mode}: no raw tables produced. "
            f"Provide pvac_files, fusion_tsv, splice_junction_tsv, peptide_table, or pre-built raw TSVs. "
            f"Supported modes: {', '.join(INPUT_MODES)}"
        )

    result = {
        "entry_mode": mode,
        "raw_events": str(raw_events),
        "raw_peptides": str(raw_peptides),
    }
    if fusion_evidence_out.is_file():
        result["fusion_evidence"] = str(fusion_evidence_out)
    return result
