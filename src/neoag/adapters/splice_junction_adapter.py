"""Mode C splice junction adapter — RegTools RNA support + pVACsplice / VCF peptide calling."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ..model_layers import enrich_peptide_layers
from ..utils import first, read_tsv, to_float
from .event_catalog import parse_splice_catalog

SPLICE_CONSEQUENCES = frozenset({
    "splice_donor_variant",
    "splice_acceptor_variant",
    "splice_donor_region_variant",
    "splice_acceptor_region_variant",
    "splice_region_variant",
    "splice_polypyrimidine_tract_variant",
})

_VARIANT_INFO_RE = re.compile(
    r"^(?P<chrom>[^:]+):(?P<start>\d+)(?:-(?P<end>\d+))?$"
)


@dataclass
class SpliceJunctionFilterConfig:
  min_junction_reads: int = 3


def _consequences(cons_field: str) -> set[str]:
    return {c.strip() for c in str(cons_field or "").split("&") if c.strip()}


def is_splice_consequence(cons_field: str) -> bool:
    return bool(_consequences(cons_field) & SPLICE_CONSEQUENCES)


def parse_regtools_variant_info(value: str) -> tuple[str, str, str] | None:
    match = _VARIANT_INFO_RE.match(str(value or "").strip())
    if not match:
        return None
    chrom = match.group("chrom")
    start = match.group("start")
    end = match.group("end") or start
    return chrom, start, end


def build_junction_support_index(splice_path: str | Path) -> dict[str, Any]:
    """Index RegTools junction support by gene and variant locus."""
    by_gene: dict[str, int] = {}
    by_locus: dict[str, int] = {}
    by_junction: dict[str, int] = {}

    for row in read_tsv(splice_path):
        gene = first(row, ["gene", "Gene", "gene_name", "gene_names", "symbol"], "")
        if gene.upper() in {"", "NA", "N/A", "."}:
            gene = ""
        reads = int(
            to_float(
                first(row, ["counts", "junction_reads", "reads", "split_reads", "score"], "0"),
                0.0,
            )
        )
        junction_id = first(row, ["event_id", "name", "junction_id"], "")
        if junction_id:
            by_junction[junction_id] = max(by_junction.get(junction_id, 0), reads)
        if gene:
            key = gene.upper()
            by_gene[key] = max(by_gene.get(key, 0), reads)
        variant_info = parse_regtools_variant_info(first(row, ["variant_info"], ""))
        if variant_info:
            chrom, start, end = variant_info
            by_locus[f"{chrom}:{start}"] = max(by_locus.get(f"{chrom}:{start}", 0), reads)
            by_locus[f"{chrom}:{end}"] = max(by_locus.get(f"{chrom}:{end}", 0), reads)
    return {"by_gene": by_gene, "by_locus": by_locus, "by_junction": by_junction}


def junction_reads_for_peptide(
    peptide: dict[str, str],
    index: dict[str, Any],
    *,
    event: dict[str, str] | None = None,
) -> int:
    reads = int(to_float(peptide.get("rna_junction_reads"), 0.0))
    gene = str(peptide.get("gene") or (event or {}).get("gene") or "").upper()
    if gene:
        reads = max(reads, index["by_gene"].get(gene, 0))
    chrom = str((event or {}).get("chrom") or peptide.get("chrom") or "")
    pos = str((event or {}).get("pos") or peptide.get("pos") or "")
    if chrom and pos:
        reads = max(reads, index["by_locus"].get(f"{chrom}:{pos}", 0))
    event_name = str((event or {}).get("event_name") or "")
    if event_name:
        reads = max(reads, index["by_junction"].get(event_name, 0))
    return reads


def enrich_splice_peptide_layers(
    peptide: dict[str, str],
    *,
    index: dict[str, Any] | None = None,
    event: dict[str, str] | None = None,
) -> dict[str, str]:
    out = dict(peptide)
    out["peptide_consequence"] = "splice_junction"
    if not out.get("crosses_junction"):
        out["crosses_junction"] = "yes"
    if index is not None:
        out["rna_junction_reads"] = str(junction_reads_for_peptide(out, index, event=event))
    if not out.get("source_tool"):
        out["source_tool"] = "splice-junction-adapter"
    return enrich_peptide_layers(out, event)


def filter_splice_variant_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [row for row in rows if is_splice_consequence(row.get("consequence", ""))]


def build_splice_peptides_from_vcf(
    variants_vcf: str | Path,
    splice_path: str | Path,
    *,
    sample_id: str,
    profile_name: str,
    hla_alleles: list[str],
    cfg: dict[str, Any] | None = None,
    tools_dir: Path | None = None,
    min_junction_reads: int = 0,
) -> dict[str, Any]:
    """Extract splice-affecting variant peptides and attach RegTools junction support."""
    from ..vep.extract_peptides import extract_variant_peptides_from_vcf
    from .variant_peptide_adapter import (
        resolve_variant_peptide_options,
        variant_peptide_rows_to_raw_tables,
    )

    if not hla_alleles:
        raise ValueError("build_splice_peptides_from_vcf requires inputs.hla_alleles")

    cfg = cfg or {}
    inputs = cfg.get("inputs") or {}
    opts = resolve_variant_peptide_options(cfg)
    out_tsv = (tools_dir or Path(".")) / "splice_variant_peptides.tsv"
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    summary = extract_variant_peptides_from_vcf(
        variants_vcf,
        out_tsv,
        lengths=opts["lengths"],
        sample_id=sample_id,
        exclude_multi_aa=opts["exclude_multi_aa"],
        single_aa_only=opts["single_aa_only"],
        mini_len=opts["mini_len"],
        normal_proteome_fasta=opts["normal_proteome_fasta"],
        filter_normal_proteome=opts["filter_normal_proteome"],
        hla_alleles=hla_alleles,
        tumor_sample_name=str(inputs.get("tumor_sample_name") or sample_id) or None,
        rna_sample_name=str(inputs.get("rna_sample_name") or "") or None,
        consequence_filter="splice",
    )
    splice_rows = filter_splice_variant_rows(read_tsv(out_tsv))
    events, peptides = variant_peptide_rows_to_raw_tables(
        splice_rows,
        sample_id=sample_id,
        profile_name=profile_name,
        hla_alleles=hla_alleles,
    )
    index = build_junction_support_index(splice_path)
    events_by_id = {e["event_id"]: e for e in events}
    enriched: list[dict[str, str]] = []
    for pep in peptides:
        event = events_by_id.get(pep["event_id"])
        pep = enrich_splice_peptide_layers(pep, index=index, event=event)
        if min_junction_reads and int(pep.get("rna_junction_reads") or 0) < min_junction_reads:
            continue
        enriched.append(pep)

    return {
        "events": events,
        "peptides": enriched,
        "splice_variant_peptides": str(out_tsv),
        "summary": summary,
        "splice_variant_rows": len(splice_rows),
    }


def merge_splice_into_catalog(
    splice_path: str | Path,
    sample_id: str,
    profile_name: str,
    events: list[dict[str, str]],
    peptides: list[dict[str, str]],
    *,
    variants_vcf: str | Path | None = None,
    hla_alleles: list[str] | None = None,
    cfg: dict[str, Any] | None = None,
    tools_dir: Path | None = None,
    filter_cfg: SpliceJunctionFilterConfig | None = None,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Merge RegTools junction events and optional splice peptides into raw tables."""
    filter_cfg = filter_cfg or SpliceJunctionFilterConfig()
    splice_events = parse_splice_catalog(splice_path, sample_id, profile_name)
    event_by_id = {e["event_id"]: e for e in events if e.get("event_id")}
    for ev in splice_events:
        event_by_id.setdefault(ev["event_id"], ev)
    merged_events = list(event_by_id.values())

    index = build_junction_support_index(splice_path)
    if peptides:
        merged_peptides = [
            enrich_splice_peptide_layers(
                pep,
                index=index,
                event=event_by_id.get(pep.get("event_id", "")),
            )
            for pep in peptides
        ]
    else:
        merged_peptides = []

    if not merged_peptides and variants_vcf and hla_alleles:
        built = build_splice_peptides_from_vcf(
            variants_vcf,
            splice_path,
            sample_id=sample_id,
            profile_name=profile_name,
            hla_alleles=list(hla_alleles),
            cfg=cfg,
            tools_dir=tools_dir,
            min_junction_reads=filter_cfg.min_junction_reads,
        )
        for ev in built["events"]:
            event_by_id.setdefault(ev["event_id"], ev)
        merged_events = list(event_by_id.values())
        merged_peptides = built["peptides"]

    if filter_cfg.min_junction_reads:
        merged_peptides = [
            pep
            for pep in merged_peptides
            if int(pep.get("rna_junction_reads") or 0) >= filter_cfg.min_junction_reads
        ]

    return merged_events, merged_peptides


def run_splice_junction_upstream(
    cfg: dict[str, Any],
    *,
    splice_path: Path,
    variants_vcf: Path | None,
    parsed_dir: Path,
    tools_dir: Path,
    sample_id: str,
    profile_name: str,
    hla_alleles: list[str],
    pvacsplice_tsv: Path | None = None,
) -> dict[str, str]:
    """Build parsed/raw tables for splice_junction entry mode."""
    from ..schemas import EVENT_FIELDS, PEPTIDE_FIELDS
    from ..utils import write_tsv
    from .pvactools_parser import parse_pvactools_outputs

    events: list[dict[str, str]] = []
    peptides: list[dict[str, str]] = []
    outputs: dict[str, str] = {"splice_junction_tsv": str(splice_path)}

    if pvacsplice_tsv and pvacsplice_tsv.is_file():
        events, peptides = parse_pvactools_outputs(
            [pvacsplice_tsv],
            sample_id,
            profile_name,
        )
        outputs["pvacsplice"] = str(pvacsplice_tsv)
        outputs["peptide_source"] = "pvacsplice"

    events, peptides = merge_splice_into_catalog(
        splice_path,
        sample_id,
        profile_name,
        events,
        peptides,
        variants_vcf=variants_vcf,
        hla_alleles=hla_alleles,
        cfg=cfg,
        tools_dir=tools_dir,
    )

    if not peptides:
        raise ValueError(
            "splice_junction mode produced no peptides. Provide pvacsplice aggregated output "
            "or variants_vcf + hla_alleles with splice-affecting PASS variants."
        )

    raw_events = parsed_dir / "raw_events.tsv"
    raw_peptides = parsed_dir / "raw_peptides.tsv"
    write_tsv(raw_events, events, EVENT_FIELDS)
    write_tsv(raw_peptides, peptides, PEPTIDE_FIELDS)
    outputs["raw_events"] = str(raw_events)
    outputs["raw_peptides"] = str(raw_peptides)
    if outputs.get("peptide_source") is None:
        outputs["peptide_source"] = "splice-variant-peptides"
    return outputs
