"""Canonical splice-junction coordinates and source-row parsing.

v0.4.4 establishes a single internal representation for junction evidence:

* genome build is explicit;
* chromosome names are normalized;
* intron coordinates are 1-based, closed;
* strand is part of the identity;
* source coordinates and the conversion rule are retained;
* unresolved records never receive another junction's evidence.

The module intentionally refuses gene-level or variant-locus fallback matching.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping

from neoag.utils import MISSING, first, to_float

CANONICAL_JUNCTION_SCHEMA_VERSION = "neoag.splice.junction.v0.4.4"

CANONICAL_JUNCTION_FIELDS = [
    "sample_id",
    "junction_id",
    "genome_build",
    "chrom",
    "intron_start_1based",
    "intron_end_1based",
    "strand",
    "donor_1based",
    "acceptor_1based",
    "source_coordinate_system",
    "source_chrom",
    "source_start",
    "source_end",
    "source_tool",
    "source_tool_version",
    "source_file",
    "source_row_number",
    "source_record_id",
    "source_junction_id",
    "gene",
    "gene_id",
    "transcript_ids",
    "unique_split_reads",
    "multi_split_reads",
    "total_split_reads",
    "splice_motif",
    "known_donor",
    "known_acceptor",
    "known_junction",
    "resolution_status",
    "resolution_method",
    "coordinate_warning",
    "record_sha256",
]

PROVENANCE_EVIDENCE_FIELDS = [
    *CANONICAL_JUNCTION_FIELDS,
    "evidence_role",
    "peptide",
    "hla_allele",
    "crosses_junction",
    "generation_status",
    "raw_event_id",
]

_COORD_PATTERNS = (
    re.compile(
        r"^(?P<chrom>[^:|]+):(?P<start>\d+)-(?P<end>\d+)(?::|\()(?P<strand>[+\-.])\)?$"
    ),
    re.compile(r"^(?P<chrom>[^:|]+):(?P<start>\d+)-(?P<end>\d+)$"),
)
_CANONICAL_ID_RE = re.compile(
    r"^SJ\|(?P<build>[^|]+)\|(?P<chrom>[^|]+)\|(?P<start>\d+)\|(?P<end>\d+)\|(?P<strand>[+\-.])$"
)


class JunctionNormalizationError(ValueError):
    """Raised when strict coordinate normalization cannot be completed."""


@dataclass(frozen=True)
class CanonicalJunction:
    genome_build: str
    chrom: str
    intron_start_1based: int
    intron_end_1based: int
    strand: str

    def __post_init__(self) -> None:
        if self.intron_start_1based <= 0:
            raise JunctionNormalizationError("intron_start_1based must be > 0")
        if self.intron_end_1based < self.intron_start_1based:
            raise JunctionNormalizationError("intron_end_1based precedes intron_start_1based")
        if self.strand not in {"+", "-", "."}:
            raise JunctionNormalizationError(f"unsupported strand: {self.strand!r}")

    @property
    def junction_id(self) -> str:
        return canonical_junction_id(
            self.genome_build,
            self.chrom,
            self.intron_start_1based,
            self.intron_end_1based,
            self.strand,
        )

    @property
    def donor_1based(self) -> int:
        # Coordinates denote the first/last intronic base. On the minus strand,
        # transcription proceeds from high to low genomic coordinates.
        return self.intron_end_1based if self.strand == "-" else self.intron_start_1based

    @property
    def acceptor_1based(self) -> int:
        return self.intron_start_1based if self.strand == "-" else self.intron_end_1based

    @property
    def stranded_key(self) -> tuple[str, str, int, int, str]:
        return (
            self.genome_build,
            self.chrom,
            self.intron_start_1based,
            self.intron_end_1based,
            self.strand,
        )

    @property
    def unstranded_key(self) -> tuple[str, str, int, int]:
        return (
            self.genome_build,
            self.chrom,
            self.intron_start_1based,
            self.intron_end_1based,
        )


@dataclass
class JunctionSourceRecord:
    row: dict[str, str]
    junction: CanonicalJunction | None
    sample_id: str
    source_tool: str
    source_tool_version: str
    source_file: str
    source_row_number: int
    source_record_id: str
    source_junction_id: str
    source_coordinate_system: str
    source_chrom: str
    source_start: str
    source_end: str
    gene: str
    gene_id: str
    transcript_ids: str
    unique_split_reads: int
    multi_split_reads: int
    total_split_reads: int
    splice_motif: str
    known_donor: str
    known_acceptor: str
    known_junction: str
    resolution_status: str
    resolution_method: str
    coordinate_warning: str
    record_sha256: str

    @property
    def junction_id(self) -> str:
        return self.junction.junction_id if self.junction else ""

    def as_dict(self) -> dict[str, str]:
        junction = self.junction
        result = {
            "sample_id": self.sample_id,
            "junction_id": junction.junction_id if junction else "",
            "genome_build": junction.genome_build if junction else "",
            "chrom": junction.chrom if junction else "",
            "intron_start_1based": str(junction.intron_start_1based) if junction else "",
            "intron_end_1based": str(junction.intron_end_1based) if junction else "",
            "strand": junction.strand if junction else "",
            "donor_1based": str(junction.donor_1based) if junction else "",
            "acceptor_1based": str(junction.acceptor_1based) if junction else "",
            "source_coordinate_system": self.source_coordinate_system,
            "source_chrom": self.source_chrom,
            "source_start": self.source_start,
            "source_end": self.source_end,
            "source_tool": self.source_tool,
            "source_tool_version": self.source_tool_version,
            "source_file": self.source_file,
            "source_row_number": str(self.source_row_number),
            "source_record_id": self.source_record_id,
            "source_junction_id": self.source_junction_id,
            "gene": self.gene,
            "gene_id": self.gene_id,
            "transcript_ids": self.transcript_ids,
            "unique_split_reads": str(self.unique_split_reads),
            "multi_split_reads": str(self.multi_split_reads),
            "total_split_reads": str(self.total_split_reads),
            "splice_motif": self.splice_motif,
            "known_donor": self.known_donor,
            "known_acceptor": self.known_acceptor,
            "known_junction": self.known_junction,
            "resolution_status": self.resolution_status,
            "resolution_method": self.resolution_method,
            "coordinate_warning": self.coordinate_warning,
            "record_sha256": self.record_sha256,
        }
        return result


def normalize_genome_build(value: Any, default: str = "GRCh38") -> str:
    raw = str(value or default).strip()
    folded = raw.lower().replace("_", "").replace("-", "")
    mapping = {
        "grch38": "GRCh38",
        "hg38": "GRCh38",
        "b38": "GRCh38",
        "grch37": "GRCh37",
        "hg19": "GRCh37",
        "b37": "GRCh37",
        "t2tchm13": "T2T-CHM13",
        "chm13": "T2T-CHM13",
    }
    return mapping.get(folded, raw or default)


def normalize_chromosome(value: Any) -> str:
    raw = str(value or "").strip()
    if not raw or raw in MISSING:
        return ""
    if raw.lower().startswith("chr"):
        suffix = raw[3:]
    else:
        suffix = raw
    if suffix.upper() in {"M", "MT", "MITO"}:
        suffix = "M"
    elif suffix.upper() in {"X", "Y"}:
        suffix = suffix.upper()
    return f"chr{suffix}"


def normalize_strand(value: Any) -> str:
    raw = str(value or "").strip()
    if raw in {"+", "1", "F", "f", "forward", "Forward"}:
        return "+"
    if raw in {"-", "-1", "R", "r", "reverse", "Reverse"}:
        return "-"
    return "."


def canonical_junction_id(
    genome_build: str,
    chrom: str,
    intron_start_1based: int,
    intron_end_1based: int,
    strand: str,
) -> str:
    return (
        f"SJ|{normalize_genome_build(genome_build)}|{normalize_chromosome(chrom)}|"
        f"{int(intron_start_1based)}|{int(intron_end_1based)}|{normalize_strand(strand)}"
    )


def parse_canonical_junction_id(value: Any) -> CanonicalJunction | None:
    match = _CANONICAL_ID_RE.match(str(value or "").strip())
    if not match:
        return None
    return CanonicalJunction(
        normalize_genome_build(match.group("build")),
        normalize_chromosome(match.group("chrom")),
        int(match.group("start")),
        int(match.group("end")),
        normalize_strand(match.group("strand")),
    )


def _integer(value: Any) -> int | None:
    raw = str(value or "").strip().replace(",", "")
    if not raw or raw in MISSING:
        return None
    try:
        return int(float(raw))
    except (TypeError, ValueError):
        return None


def _row_hash(row: Mapping[str, Any]) -> str:
    payload = json.dumps(
        {str(k): "" if v is None else str(v) for k, v in sorted(row.items(), key=lambda item: str(item[0]))},
        sort_keys=True,
        ensure_ascii=False,
        separators=(",", ":"),
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_coordinate_token(
    value: Any,
    *,
    genome_build: str = "GRCh38",
    coordinate_system: str = "intron_1based_closed",
) -> CanonicalJunction | None:
    canonical = parse_canonical_junction_id(value)
    if canonical:
        return canonical
    raw = str(value or "").strip()
    for pattern in _COORD_PATTERNS:
        match = pattern.match(raw)
        if not match:
            continue
        chrom = normalize_chromosome(match.group("chrom"))
        start = int(match.group("start"))
        end = int(match.group("end"))
        strand = normalize_strand(match.groupdict().get("strand", "."))
        start1, end1 = convert_interval(start, end, coordinate_system)
        return CanonicalJunction(normalize_genome_build(genome_build), chrom, start1, end1, strand)
    return None


def convert_interval(start: int, end: int, coordinate_system: str) -> tuple[int, int]:
    """Convert a source interval to 1-based closed intron coordinates.

    Supported systems:
      * intron_1based_closed / star_sj
      * bed0_half_open / regtools_annotated / splice_boundaries_0based
      * splice_boundaries_1based / snaf_uid (outer exon-boundary coordinates)
    """

    system = str(coordinate_system or "").strip().lower()
    if system in {"intron_1based_closed", "1based", "one_based", "star_sj", "canonical"}:
        start1, end1 = start, end
    elif system in {
        "bed0_half_open",
        "bed0",
        "regtools_annotated",
        "splice_boundaries_0based",
    }:
        start1, end1 = start + 1, end
    elif system in {"splice_boundaries_1based", "snaf_uid"}:
        start1, end1 = start + 1, end - 1
    else:
        raise JunctionNormalizationError(f"unsupported coordinate system: {coordinate_system!r}")
    if start1 <= 0 or end1 < start1:
        raise JunctionNormalizationError(
            f"invalid converted intron interval {start1}-{end1} from {start}-{end} ({coordinate_system})"
        )
    return start1, end1


def _bed12_intron(row: Mapping[str, Any]) -> tuple[int, int]:
    chrom_start0 = _integer(first(row, ["start", "chromStart", "chrom_start"], ""))
    block_count = _integer(first(row, ["blockCount", "block_count"], ""))
    sizes_raw = first(row, ["blockSizes", "block_sizes"], "")
    starts_raw = first(row, ["blockStarts", "block_starts"], "")
    sizes = [_integer(value) for value in str(sizes_raw).strip(",").split(",") if value != ""]
    starts = [_integer(value) for value in str(starts_raw).strip(",").split(",") if value != ""]
    if chrom_start0 is None or not block_count or block_count < 2:
        raise JunctionNormalizationError("BED12 row lacks a valid two-block junction")
    if len(sizes) < 2 or len(starts) < 2 or any(value is None for value in sizes + starts):
        raise JunctionNormalizationError("BED12 blockSizes/blockStarts are incomplete")
    left_exon_end0 = chrom_start0 + int(starts[0]) + int(sizes[0])
    right_exon_start0 = chrom_start0 + int(starts[-1])
    return left_exon_end0 + 1, right_exon_start0


def detect_coordinate_system(
    row: Mapping[str, Any],
    *,
    source_tool: str,
    requested: str = "auto",
) -> str:
    requested_folded = str(requested or "auto").strip().lower()
    if requested_folded not in {"", "auto"}:
        return requested_folded
    keys = {str(key).lower() for key in row}
    tool = str(source_tool or "").strip().lower()
    explicit = first(row, ["source_coordinate_system", "coordinate_system"], "")
    if explicit:
        return explicit.strip().lower()
    if {"intron_start_1based", "intron_end_1based"} <= keys:
        return "intron_1based_closed"
    if {"blockcount", "blocksizes", "blockstarts"} <= keys or {
        "block_count",
        "block_sizes",
        "block_starts",
    } <= keys:
        return "bed12"
    if tool == "regtools" or {"splice_site", "known_junction"} <= keys:
        return "regtools_annotated"
    # SNAF commonly emits chr:start-end(strand) UIDs as outer splice boundaries.
    uid = first(row, ["uid", "junction", "junction_id"], "")
    if tool == "snaf" and re.search(r":\d+-\d+\([+-]\)$", uid):
        return "snaf_uid"
    # Generic normalized tables are treated as direct intron coordinates. The
    # source system remains explicit in the emitted provenance table.
    return "intron_1based_closed"


def _source_junction_alias(row: Mapping[str, Any]) -> str:
    return first(
        row,
        [
            "source_junction_id",
            "junction_id",
            "Splice Junction",
            "splice_junction",
            "event_id",
            "UID",
            "uid",
            "name",
            "junction",
        ],
        "",
    ).strip()


def _peptide_from_row(row: Mapping[str, Any]) -> str:
    return first(
        row,
        ["peptide", "junction_peptide", "mutant_peptide", "neoepitope", "MT Epitope Seq", "epitope"],
        "",
    ).strip()


def _hla_from_row(row: Mapping[str, Any]) -> str:
    return first(row, ["hla_allele", "HLA Allele", "hla", "allele"], "").strip()


def junction_record_from_row(
    row: Mapping[str, Any],
    *,
    sample_id: str,
    source_tool: str,
    source_file: str | Path,
    source_row_number: int,
    genome_build: str = "GRCh38",
    coordinate_system: str = "auto",
    source_tool_version: str = "UNASSESSED",
    strict: bool = False,
) -> JunctionSourceRecord:
    raw = {str(key): "" if value is None else str(value) for key, value in row.items()}
    build = normalize_genome_build(
        first(raw, ["genome_build", "reference_build", "assembly", "build"], genome_build),
        genome_build,
    )
    alias = _source_junction_alias(raw)
    record_id = first(raw, ["source_record_id", "row_id", "record_id"], "") or alias
    if not record_id:
        record_id = f"{source_tool}:{source_row_number}"
    system = detect_coordinate_system(raw, source_tool=source_tool, requested=coordinate_system)
    source_chrom = first(raw, ["source_chrom", "chrom1", "chrom", "Chromosome", "chr", "seqnames"], "")
    source_start = first(
        raw,
        ["source_start", "intron_start_1based", "start1", "junction_start", "intron_start", "start", "Start", "donor"],
        "",
    )
    source_end = first(
        raw,
        ["source_end", "intron_end_1based", "end1", "junction_end", "intron_end", "end", "End", "acceptor"],
        "",
    )
    strand = normalize_strand(first(raw, ["strand", "Strand"], "."))
    warning = ""
    method = ""
    junction: CanonicalJunction | None = None

    # A canonical junction ID is authoritative and already carries build/strand.
    for token_key in ("junction_id", "canonical_junction_id", "event_id"):
        parsed = parse_canonical_junction_id(raw.get(token_key, ""))
        if parsed:
            junction = parsed
            system = "canonical"
            method = f"canonical_id:{token_key}"
            source_chrom = source_chrom or parsed.chrom
            source_start = source_start or str(parsed.intron_start_1based)
            source_end = source_end or str(parsed.intron_end_1based)
            break

    if junction is None and system == "bed12":
        try:
            start1, end1 = _bed12_intron(raw)
            chrom = normalize_chromosome(source_chrom)
            if not chrom:
                raise JunctionNormalizationError("BED12 row lacks chromosome")
            junction = CanonicalJunction(build, chrom, start1, end1, strand)
            method = "bed12_blocks"
        except JunctionNormalizationError as exc:
            warning = str(exc)

    if junction is None:
        chrom = normalize_chromosome(source_chrom)
        start = _integer(source_start)
        end = _integer(source_end)
        if chrom and start is not None and end is not None:
            try:
                start1, end1 = convert_interval(start, end, system)
                junction = CanonicalJunction(build, chrom, start1, end1, strand)
                method = f"explicit_coordinates:{system}"
            except JunctionNormalizationError as exc:
                warning = str(exc)

    if junction is None:
        # Coordinate-bearing UID is permitted, but only under an explicit or
        # tool-specific conversion convention. It is never matched by gene.
        token = alias or first(raw, ["junction", "uid", "event_name"], "")
        token_system = system if system not in {"auto", "bed12"} else "intron_1based_closed"
        try:
            parsed_token = parse_coordinate_token(token, genome_build=build, coordinate_system=token_system)
        except JunctionNormalizationError as exc:
            parsed_token = None
            warning = str(exc)
        if parsed_token:
            junction = parsed_token
            method = f"coordinate_token:{token_system}"
            source_chrom = source_chrom or parsed_token.chrom
            source_start = source_start or str(parsed_token.intron_start_1based)
            source_end = source_end or str(parsed_token.intron_end_1based)

    if strict and junction is None:
        raise JunctionNormalizationError(
            f"unable to normalize {source_tool} row {source_row_number} from {source_file}: {warning or 'coordinates missing'}"
        )
    if strict and junction is not None and junction.strand == ".":
        raise JunctionNormalizationError(
            f"unable to normalize {source_tool} row {source_row_number} from {source_file}: strand is required in strict mode"
        )

    unique_reads = int(
        to_float(
            first(
                raw,
                [
                    "unique_split_reads",
                    "junction_reads",
                    "rna_junction_reads",
                    "RNA Junction Reads",
                    "read_count",
                    "counts",
                    "reads",
                    "split_reads",
                    "score",
                    "JunctionReadCount",
                ],
                "0",
            ),
            0.0,
        )
    )
    multi_reads = int(to_float(first(raw, ["multi_split_reads", "multi_reads"], "0"), 0.0))
    total_reads = int(
        to_float(
            first(raw, ["total_split_reads", "total_junction_reads"], str(unique_reads + multi_reads)),
            float(unique_reads + multi_reads),
        )
    )
    if total_reads < unique_reads:
        total_reads = unique_reads
    status = "RESOLVED" if junction else "UNRESOLVED"
    if junction and junction.strand == ".":
        status = "RESOLVED_UNSTRANDED"
        warning = "; ".join(value for value in (warning, "strand unavailable") if value)

    transcripts = first(raw, ["transcript_ids", "transcripts", "Transcript", "transcript_id"], "")
    record = JunctionSourceRecord(
        row=raw,
        junction=junction,
        sample_id=sample_id,
        source_tool=source_tool,
        source_tool_version=source_tool_version,
        source_file=str(source_file),
        source_row_number=source_row_number,
        source_record_id=record_id,
        source_junction_id=alias,
        source_coordinate_system=system,
        source_chrom=source_chrom,
        source_start=source_start,
        source_end=source_end,
        gene=first(raw, ["gene", "Gene", "gene_name", "gene_names", "symbol"], ""),
        gene_id=first(raw, ["gene_id", "gene_ids", "ensembl_gene_id"], ""),
        transcript_ids=transcripts,
        unique_split_reads=unique_reads,
        multi_split_reads=multi_reads,
        total_split_reads=total_reads,
        splice_motif=first(raw, ["splice_motif", "splice_site", "motif"], ""),
        known_donor=first(raw, ["known_donor"], ""),
        known_acceptor=first(raw, ["known_acceptor"], ""),
        known_junction=first(raw, ["known_junction"], ""),
        resolution_status=status,
        resolution_method=method or "unresolved",
        coordinate_warning=warning,
        record_sha256=_row_hash(raw),
    )
    return record


def _looks_like_header(cells: list[str]) -> bool:
    folded = {cell.strip().lower() for cell in cells}
    known = {
        "chrom",
        "chromosome",
        "chr",
        "start",
        "end",
        "name",
        "junction_id",
        "event_id",
        "gene",
        "uid",
        "splice junction",
    }
    return bool(folded & known)


def read_source_rows(path: str | Path) -> list[dict[str, str]]:
    """Read TSV/CSV or headerless RegTools BED/BED12 rows."""

    source = Path(path)
    if not source.is_file() or source.stat().st_size == 0:
        return []
    lines = [
        line
        for line in source.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.strip() and not line.startswith("#")
    ]
    if not lines:
        return []
    delimiter = "\t" if "\t" in lines[0] else ","
    first_cells = next(csv.reader([lines[0]], delimiter=delimiter))
    if _looks_like_header(first_cells):
        return [
            {str(key): "" if value is None else str(value) for key, value in row.items()}
            for row in csv.DictReader(lines, delimiter=delimiter)
        ]

    parsed: list[dict[str, str]] = []
    for line in lines:
        cells = next(csv.reader([line], delimiter=delimiter))
        if len(cells) >= 12:
            keys = [
                "chrom",
                "start",
                "end",
                "name",
                "score",
                "strand",
                "thickStart",
                "thickEnd",
                "itemRgb",
                "blockCount",
                "blockSizes",
                "blockStarts",
            ]
        elif len(cells) >= 6:
            keys = ["chrom", "start", "end", "name", "score", "strand"]
        elif len(cells) >= 5:
            keys = ["chrom", "start", "end", "name", "score"]
        else:
            continue
        parsed.append({key: cells[index] if index < len(cells) else "" for index, key in enumerate(keys)})
    return parsed


def iter_junction_records(
    path: str | Path,
    *,
    sample_id: str,
    source_tool: str,
    genome_build: str = "GRCh38",
    coordinate_system: str = "auto",
    source_tool_version: str = "UNASSESSED",
    strict: bool = False,
) -> Iterable[JunctionSourceRecord]:
    for index, row in enumerate(read_source_rows(path), 1):
        yield junction_record_from_row(
            row,
            sample_id=sample_id,
            source_tool=source_tool,
            source_file=path,
            source_row_number=index,
            genome_build=genome_build,
            coordinate_system=coordinate_system,
            source_tool_version=source_tool_version,
            strict=strict,
        )


def peptide_metadata(record: JunctionSourceRecord) -> dict[str, str]:
    row = record.row
    return {
        "peptide": _peptide_from_row(row),
        "hla_allele": _hla_from_row(row),
        "crosses_junction": first(row, ["crosses_junction", "Crosses Junction"], "true"),
        "generation_status": first(
            row,
            ["generation_status", "status"],
            "provided_by_splice_caller",
        ),
        "raw_event_id": first(row, ["event_id", "uid", "junction_id", "Splice Junction", "name"], ""),
    }
