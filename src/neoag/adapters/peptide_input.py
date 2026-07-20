"""Adapt arbitrary peptide tables to strict peptide–HLA pair inputs.

Migrated from neoantigen2 ``run_peptide_predict.sh`` input conversion, with
pair-level deduplication (never collapse rows that share a peptide but differ
in HLA allele).
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ..schemas import EVENT_FIELDS, PEPTIDE_FIELDS
from ..model_layers import enrich_event_layers, infer_mutation_source, infer_peptide_consequence
from ..utils import MISSING, safe_id, write_tsv

PEPTIDE_COLUMN_ALIASES = frozenset({
    "peptide",
    "peptide_seq",
    "seq",
    "sequence",
    "epitope",
    "mer",
    "mt epitope seq",
    "mt epitope",
})
HLA_COLUMN_ALIASES = frozenset({
    "hla",
    "hla_allele",
    "allele",
    "mhc",
    "hla_type",
    "hla allele",
})
SAMPLE_COLUMN_ALIASES = frozenset({"sample_id", "sample", "sample_name"})
GENE_COLUMN_ALIASES = frozenset({"gene_name", "gene", "genename"})
VARIANT_COLUMN_ALIASES = frozenset({"variant_id", "variant", "var_id", "event_id"})
MHC_CLASS_ALIASES = frozenset({"mhc_class", "class", "mhc"})
STANDARD_AAS = frozenset("ACDEFGHIKLMNPQRSTVWY")


@dataclass
class PeptideHlaRecord:
    peptide: str
    hla_allele: str
    sample_id: str = "SAMPLE001"
    peptide_id: str = ""
    event_id: str = ""
    gene: str = ""
    variant_id: str = ""
    mhc_class: str = "I"
    source_tool: str = "peptide_input"
    extra: dict[str, str] = field(default_factory=dict)


@dataclass
class PeptideInputSummary:
    input_path: str
    delimiter: str
    peptide_column: str
    hla_column: str
    input_rows: int
    pair_rows: int
    skipped_rows: int
    unique_peptides: int
    unique_hla: int
    raw_peptides_tsv: str
    pairs_tsv: str
    hla_alleles_txt: str


def sniff_delimiter(path: str | Path, sample_bytes: int = 8192) -> str:
    text = Path(path).read_text(encoding="utf-8")[:sample_bytes]
    tab_count = text.count("\t")
    comma_count = text.count(",")
    return "\t" if tab_count > comma_count else ","


def read_peptide_table(path: str | Path) -> tuple[str, list[str], list[dict[str, str]]]:
    path = Path(path)
    delimiter = sniff_delimiter(path)
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        headers = [h.strip() for h in (reader.fieldnames or []) if h]
        rows = [dict(row) for row in reader]
    return delimiter, headers, rows


def detect_column(headers: list[str], aliases: frozenset[str]) -> str | None:
    for header in headers:
        if header.strip().lower() in aliases:
            return header
    return None


def normalize_peptide(raw: str, *, require_standard_aa: bool = True) -> str | None:
    pep = raw.strip().upper()
    if not pep or pep in MISSING:
        return None
    if require_standard_aa and not set(pep).issubset(STANDARD_AAS):
        return None
    return pep


def normalize_hla_allele(raw: str) -> str:
    s = raw.strip().upper()
    if not s or s in MISSING:
        return ""
    compact = s.replace("HLA-", "").replace("*", "").replace(":", "")
    m = re.match(r"^([ABC])(\d{2})(\d{2})$", compact)
    if m:
        return f"HLA-{m.group(1)}*{m.group(2)}:{m.group(3)}"
    m = re.match(r"^(?:HLA-)?([ABC])\*(\d{2}):(\d{2})$", s)
    if m:
        return f"HLA-{m.group(1)}*{m.group(2)}:{m.group(3)}"
    if s.startswith("HLA-"):
        return s
    if re.match(r"^[ABC]\*", s):
        return f"HLA-{s}"
    return s


def pair_key(peptide: str, hla_allele: str) -> tuple[str, str]:
    return peptide.strip().upper(), normalize_hla_allele(hla_allele)


def extract_peptide_hla_records(
    rows: list[dict[str, str]],
    headers: list[str],
    *,
    sample_id: str = "SAMPLE001",
    require_hla: bool = True,
    require_standard_aa: bool = True,
) -> tuple[list[PeptideHlaRecord], int]:
    peptide_col = detect_column(headers, PEPTIDE_COLUMN_ALIASES)
    if peptide_col is None:
        raise ValueError(
            "Cannot detect peptide column. Supported names: "
            + ", ".join(sorted(PEPTIDE_COLUMN_ALIASES))
            + f". Found headers: {headers}"
        )

    hla_col = detect_column(headers, HLA_COLUMN_ALIASES)
    if require_hla and hla_col is None:
        raise ValueError(
            "Cannot detect HLA column. Supported names: "
            + ", ".join(sorted(HLA_COLUMN_ALIASES))
            + f". Found headers: {headers}"
        )

    sample_col = detect_column(headers, SAMPLE_COLUMN_ALIASES)
    gene_col = detect_column(headers, GENE_COLUMN_ALIASES)
    variant_col = detect_column(headers, VARIANT_COLUMN_ALIASES)
    mhc_col = detect_column(headers, MHC_CLASS_ALIASES)

    records: list[PeptideHlaRecord] = []
    skipped = 0
    for i, row in enumerate(rows):
        peptide = normalize_peptide(row.get(peptide_col, ""), require_standard_aa=require_standard_aa)
        if not peptide:
            skipped += 1
            continue

        hla_raw = row.get(hla_col, "").strip() if hla_col else ""
        hla = normalize_hla_allele(hla_raw)
        if require_hla and not hla:
            skipped += 1
            continue

        sample = row.get(sample_col, "").strip() if sample_col else sample_id
        gene = row.get(gene_col, "").strip() if gene_col else ""
        variant = row.get(variant_col, "").strip() if variant_col else ""
        mhc_class = row.get(mhc_col, "").strip() if mhc_col else "I"
        if not mhc_class or mhc_class in MISSING:
            mhc_class = "I"

        event_id = variant or (f"NA_{i + 1}" if not gene else safe_id(f"{gene}_{variant or i + 1}"))
        peptide_id = safe_id(f"pep_{peptide}_{hla or 'NA'}_{i + 1}")

        records.append(
            PeptideHlaRecord(
                peptide=peptide,
                hla_allele=hla,
                sample_id=sample or sample_id,
                peptide_id=peptide_id,
                event_id=event_id,
                gene=gene or "NA",
                variant_id=variant,
                mhc_class=mhc_class,
            )
        )
    return records, skipped


def unique_peptide_hla_records(records: list[PeptideHlaRecord]) -> list[PeptideHlaRecord]:
    """Keep first occurrence of each (peptide, HLA) pair."""
    seen: set[tuple[str, str]] = set()
    out: list[PeptideHlaRecord] = []
    for rec in records:
        key = pair_key(rec.peptide, rec.hla_allele)
        if key in seen:
            continue
        seen.add(key)
        out.append(rec)
    return out


def records_to_raw_peptide_rows(records: list[PeptideHlaRecord]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for rec in records:
        row = {field: "" for field in PEPTIDE_FIELDS}
        row.update(
            {
                "peptide_id": rec.peptide_id,
                "event_id": rec.event_id,
                "sample_id": rec.sample_id,
                "event_type": "SNV",
                "mutation_source": "SNV",
                "peptide_consequence": "missense",
                "gene": rec.gene,
                "peptide": rec.peptide,
                "hla_allele": rec.hla_allele,
                "mhc_class": rec.mhc_class,
                "source_tool": rec.source_tool,
            }
        )
        rows.append(row)
    return rows


def write_hla_allele_list(records: list[PeptideHlaRecord], path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    alleles: list[str] = []
    seen: set[str] = set()
    for rec in records:
        hla = normalize_hla_allele(rec.hla_allele)
        if hla and hla not in seen:
            seen.add(hla)
            alleles.append(hla)
    path.write_text("\n".join(alleles) + ("\n" if alleles else ""), encoding="utf-8")
    return path


def write_pairs_tsv(records: list[PeptideHlaRecord], path: str | Path) -> Path:
    path = Path(path)
    rows = [
        {
            "sample_id": rec.sample_id,
            "peptide": rec.peptide,
            "hla_allele": rec.hla_allele,
            "peptide_hla_key": safe_id(f"{rec.peptide}_{rec.hla_allele}"),
            "gene": rec.gene,
            "event_id": rec.event_id,
            "peptide_id": rec.peptide_id,
        }
        for rec in records
    ]
    write_tsv(
        path,
        rows,
        [
            "sample_id",
            "peptide",
            "hla_allele",
            "peptide_hla_key",
            "gene",
            "event_id",
            "peptide_id",
        ],
    )
    return path


def build_raw_events_from_peptides(
    raw_peptides_path: str | Path,
    raw_events_path: str | Path,
    sample_id: str,
    profile_name: str,
) -> list[dict[str, str]]:
    """Synthesize one Layer-1 event per unique event_id in a peptide-only table."""
    peptides = read_tsv(raw_peptides_path) if Path(raw_peptides_path).is_file() else []
    events: dict[str, dict[str, str]] = {}
    for pep in peptides:
        eid = pep.get("event_id") or safe_id(f"{sample_id}_PEP_{pep.get('peptide_id', '')}")
        if eid in events:
            continue
        base = {
            "event_id": eid,
            "sample_id": sample_id,
            "disease_profile": profile_name,
            "event_type": pep.get("event_type") or "SNV",
            "mutation_source": pep.get("mutation_source") or infer_mutation_source(event_type=pep.get("event_type", "SNV")),
            "peptide_consequence": pep.get("peptide_consequence") or infer_peptide_consequence(event_type=pep.get("event_type", "SNV")),
            "gene": pep.get("gene") or "UNKNOWN",
            "event_name": pep.get("gene") or eid,
            "chrom": "",
            "pos": "",
            "ref": "",
            "alt": "",
            "transcript_id": "",
            "consequence": "",
            "rna_junction_reads": pep.get("rna_junction_reads", ""),
            "event_confidence": "0.5",
            "event_expression": "0.0",
            "driver_relevance": "0.0",
            "tumor_vaf": "0.0",
            "tumor_depth": "",
            "tumor_alt_count": "",
            "rna_vaf": "",
            "rna_alt_reads": "",
            "rna_depth": "",
            "clonality": "0.5",
            "persistence": "0.5",
            "tumor_specificity": "0.7",
            "source": pep.get("source_tool") or "peptide_input",
        }
        events[eid] = enrich_event_layers(base)
    rows = list(events.values())
    write_tsv(raw_events_path, rows, EVENT_FIELDS)
    return rows


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    from ..utils import read_tsv as _read

    return _read(path)


def convert_peptide_input(
    input_path: str | Path,
    outdir: str | Path,
    *,
    sample_id: str = "SAMPLE001",
    require_hla: bool = True,
) -> PeptideInputSummary:
    """Convert CSV/TSV peptide table into pair-level pipeline inputs."""
    input_path = Path(input_path)
    outdir = Path(outdir)
    parsed = outdir / "parsed"
    parsed.mkdir(parents=True, exist_ok=True)
    prep = parsed  # standard intermediate layout (legacy 00_input removed)

    delimiter, headers, rows = read_peptide_table(input_path)
    if not rows:
        raise ValueError(f"Input file is empty: {input_path}")

    peptide_col = detect_column(headers, PEPTIDE_COLUMN_ALIASES)
    hla_col = detect_column(headers, HLA_COLUMN_ALIASES)
    if peptide_col is None:
        raise ValueError(f"Cannot detect peptide column in {input_path}")

    records, skipped = extract_peptide_hla_records(
        rows,
        headers,
        sample_id=sample_id,
        require_hla=require_hla,
    )
    pair_records = unique_peptide_hla_records(records)
    if not pair_records:
        raise ValueError(f"No valid peptide–HLA pairs found in {input_path}")

    raw_peptides = prep / "raw_peptides.tsv"
    pairs_tsv = prep / "peptide_hla_pairs.tsv"
    hla_txt = prep / "hla_alleles.txt"

    write_tsv(raw_peptides, records_to_raw_peptide_rows(pair_records), PEPTIDE_FIELDS)
    write_pairs_tsv(pair_records, pairs_tsv)
    write_hla_allele_list(pair_records, hla_txt)

    unique_peptides = len({rec.peptide for rec in pair_records})
    unique_hla = len({rec.hla_allele for rec in pair_records})

    return PeptideInputSummary(
        input_path=str(input_path.resolve()),
        delimiter=delimiter,
        peptide_column=peptide_col,
        hla_column=hla_col or "",
        input_rows=len(rows),
        pair_rows=len(pair_records),
        skipped_rows=skipped,
        unique_peptides=unique_peptides,
        unique_hla=unique_hla,
        raw_peptides_tsv=str(raw_peptides.resolve()),
        pairs_tsv=str(pairs_tsv.resolve()),
        hla_alleles_txt=str(hla_txt.resolve()),
    )


def unique_peptide_hla_pairs_from_table(path: str | Path) -> list[tuple[str, str]]:
    """Load any supported peptide table and return unique (peptide, HLA) pairs."""
    _, headers, rows = read_peptide_table(path)
    records, _ = extract_peptide_hla_records(rows, headers)
    return [(rec.peptide, rec.hla_allele) for rec in unique_peptide_hla_records(records)]
