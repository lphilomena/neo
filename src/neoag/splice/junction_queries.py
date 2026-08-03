"""Build exact EasyQuant query contexts from canonical splice junctions."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from neoag.splice.identifiers import stable_id

from .sequence_queries import make_sequence_query


_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _tokens(value: str) -> set[str]:
    return {token.strip() for token in str(value or "").split(";") if token.strip()}


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1].upper()


@dataclass(frozen=True)
class FaiEntry:
    length: int
    offset: int
    line_bases: int
    line_width: int


class IndexedFasta:
    """Minimal read-only FASTA index reader with no optional dependencies."""

    def __init__(self, fasta: str | Path):
        self.path = Path(fasta)
        fai = Path(f"{self.path}.fai")
        if not self.path.is_file() or not fai.is_file():
            raise ValueError(f"reference FASTA and .fai are required: {self.path}")
        self.index: dict[str, FaiEntry] = {}
        for line in fai.read_text(encoding="utf-8").splitlines():
            cells = line.split("\t")
            if len(cells) >= 5:
                self.index[cells[0]] = FaiEntry(*(int(value) for value in cells[1:5]))

    def resolve_contig(self, chrom: str) -> str:
        candidates = [chrom]
        if chrom.startswith("chr"):
            candidates.append(chrom[3:])
        else:
            candidates.append(f"chr{chrom}")
        if chrom in {"chrM", "M"}:
            candidates.extend(["MT", "chrMT"])
        for candidate in candidates:
            if candidate in self.index:
                return candidate
        return ""

    def fetch(self, chrom: str, start_1based: int, end_1based: int) -> str:
        contig = self.resolve_contig(chrom)
        if not contig:
            return ""
        entry = self.index[contig]
        start = max(1, start_1based)
        end = min(entry.length, end_1based)
        if end < start:
            return ""
        zero = start - 1
        byte_offset = entry.offset + (zero // entry.line_bases) * entry.line_width + (zero % entry.line_bases)
        wanted = end - start + 1
        chunks = bytearray()
        with self.path.open("rb") as handle:
            handle.seek(byte_offset)
            while len(chunks) < wanted:
                block = handle.read(max(4096, wanted - len(chunks) + 128))
                if not block:
                    break
                chunks.extend(byte for byte in block if byte not in b"\r\n\t ")
        return bytes(chunks[:wanted]).decode("ascii").upper()


def build_canonical_junction_queries(
    tables: dict[str, list[dict[str, str]]],
    *,
    sample_id: str,
    reference_fasta: str | Path,
    flank_bases: int = 31,
) -> dict[str, list[dict[str, str]]]:
    """Generate one strand-aware query per peptide-bearing event/junction pair."""

    if flank_bases < 8:
        raise ValueError("junction query flank must be at least 8 bases")
    fasta = IndexedFasta(reference_fasta)
    junction_by_id = {row.get("junction_id", ""): row for row in tables.get("junctions", [])}
    targets: set[tuple[str, str]] = set()
    for origin in tables.get("peptide_origins", []):
        if str(origin.get("crosses_junction", "")).lower() not in {"1", "true", "yes"}:
            continue
        event_id = origin.get("splice_event_id", "")
        for junction_id in _tokens(origin.get("junction_ids", "")):
            if event_id and junction_id:
                targets.add((event_id, junction_id))

    queries: list[dict[str, str]] = []
    conflicts: list[dict[str, str]] = []
    for event_id, junction_id in sorted(targets):
        junction = junction_by_id.get(junction_id)
        if not junction:
            reason = "Candidate peptide references a junction absent from the canonical junction table."
        else:
            try:
                intron_start = int(junction.get("intron_start_1based", ""))
                intron_end = int(junction.get("intron_end_1based", ""))
            except ValueError:
                intron_start = intron_end = 0
            chrom = junction.get("chrom", "")
            strand = junction.get("strand", "")
            left = fasta.fetch(chrom, intron_start - flank_bases, intron_start - 1) if intron_start else ""
            right = fasta.fetch(chrom, intron_end + 1, intron_end + flank_bases) if intron_end else ""
            if strand == "+":
                sequence = left + right
            elif strand == "-":
                sequence = _reverse_complement(right) + _reverse_complement(left)
            else:
                sequence = ""
            reason = "" if len(left) == flank_bases and len(right) == flank_bases and sequence else "Reference flank or strand could not be resolved."
        if reason:
            conflict = {
                "entity_type": "JUNCTION", "entity_id": junction_id, "sample_id": sample_id,
                "conflict_type": "JUNCTION_QUERY_UNRESOLVED", "field_name": "nucleotide_sequence",
                "observed_values": junction_id, "source_tools": "NeoAgCanonicalJunctionContext",
                "source_record_ids": junction_id, "severity": "WARNING", "resolution_status": "UNRESOLVED",
                "resolution_reason": reason,
            }
            conflict["conflict_id"] = stable_id("CFL", conflict)
            conflicts.append(conflict)
            continue
        queries.append(make_sequence_query(
            sample_id=sample_id,
            query_type="TARGETED_CANONICAL_JUNCTION",
            nucleotide_sequence=sequence,
            position_1based=flank_bases,
            query_length=flank_bases * 2,
            splice_event_id=event_id,
            junction_id=junction_id,
            sequence_scope="CANONICAL_JUNCTION_CONTEXT",
            source_generator="NeoAgCanonicalJunctionContext",
            source_file=str(Path(reference_fasta)),
            source_record_id=junction_id,
        ))
    return {
        "sequence_queries": queries,
        "conflicts": conflicts,
        "manifest": [{
            "adapter": "canonical_junction_queries",
            "input_path": str(Path(reference_fasta)),
            "target_event_junction_pairs": str(len(targets)),
            "queries_generated": str(len(queries)),
            "queries_unresolved": str(len(conflicts)),
            "flank_bases": str(flank_bases),
            "selection_policy": "PEPTIDE_ORIGIN_EXACT_CANONICAL_JUNCTION",
        }],
    }
