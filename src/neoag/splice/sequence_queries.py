"""Exact nucleotide-context query registry for EasyQuant and k4neo."""
from __future__ import annotations

from pathlib import Path
from typing import Any

from neoag.splice.identifiers import stable_id, sequence_sha256
from neoag.utils import write_tsv

from .schemas import SEQUENCE_QUERY_FIELDS

_DNA = set("ACGTN")


def clean_dna(value: Any) -> str:
    seq = "".join(str(value or "").split()).upper().replace("U", "T")
    return seq if seq and set(seq) <= _DNA else ""


def make_sequence_query(
    *,
    sample_id: str,
    query_type: str,
    nucleotide_sequence: str,
    position_1based: int | str = "",
    query_length: int | str = "",
    splice_event_id: str = "",
    junction_id: str = "",
    variant_id: str = "",
    transcript_hypothesis_id: str = "",
    orf_id: str = "",
    origin_peptide_id: str = "",
    sequence_scope: str = "JUNCTION_CONTEXT",
    source_generator: str = "",
    source_file: str = "",
    source_record_id: str = "",
    query_name: str = "",
) -> dict[str, str]:
    sequence = clean_dna(nucleotide_sequence)
    try:
        pos = int(float(position_1based)) if str(position_1based).strip() else 0
    except Exception:
        pos = 0
    if sequence and pos and not 1 <= pos <= len(sequence):
        raise ValueError(f"Query position {pos} lies outside sequence of length {len(sequence)}")
    qid = stable_id(
        "SEQ",
        sample_id,
        query_type,
        sequence_sha256(sequence),
        pos,
        query_length,
        splice_event_id,
        junction_id,
        variant_id,
        transcript_hypothesis_id,
        origin_peptide_id,
    )
    return {
        "query_id": qid,
        "sample_id": sample_id,
        "query_name": query_name or qid,
        "query_type": query_type,
        "splice_event_id": splice_event_id,
        "junction_id": junction_id,
        "variant_id": variant_id,
        "transcript_hypothesis_id": transcript_hypothesis_id,
        "orf_id": orf_id,
        "origin_peptide_id": origin_peptide_id,
        "nucleotide_sequence": sequence,
        "sequence_sha256": sequence_sha256(sequence),
        "position_1based": str(pos) if pos else "",
        "query_length": str(query_length or ""),
        "sequence_scope": sequence_scope,
        "source_generator": source_generator,
        "source_file": source_file,
        "source_record_id": source_record_id,
        "query_status": "READY" if sequence and (pos or query_type.upper() == "K4NEO_FULL_SEQUENCE") else "INCOMPLETE",
        "evidence_conflict_status": "NONE" if sequence else "INVALID_NUCLEOTIDE_SEQUENCE",
    }


def write_external_query_files(outdir: str | Path, queries: list[dict[str, str]]) -> dict[str, Path]:
    """Write exact query inputs and maps for both external tools.

    Query names/cts_ids are always the content-derived ``query_id``. This is the
    only key accepted when importing EasyQuant or k4neo output.
    """
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    ready = [q for q in queries if q.get("query_status") == "READY" and q.get("nucleotide_sequence")]

    easy = [q for q in ready if q.get("position_1based")]
    easy_input = out / "splice_easyquant_input.tsv"
    write_tsv(easy_input, [
        {"name": q["query_id"], "sequence": q["nucleotide_sequence"], "position": q["position_1based"]}
        for q in easy
    ], ["name", "sequence", "position"])
    easy_map = out / "splice_easyquant_query_map.tsv"
    write_tsv(easy_map, easy, SEQUENCE_QUERY_FIELDS)

    k4 = ready
    k4_input = out / "splice_k4neo_input.tsv"
    write_tsv(k4_input, [
        {
            "cts_id": q["query_id"],
            "cts_seq": q["nucleotide_sequence"],
            "pos": q.get("position_1based", ""),
            "query_length": q.get("query_length", ""),
        }
        for q in k4
    ], ["cts_id", "cts_seq", "pos", "query_length"])
    k4_map = out / "splice_k4neo_query_map.tsv"
    write_tsv(k4_map, k4, SEQUENCE_QUERY_FIELDS)
    return {
        "easyquant_input": easy_input,
        "easyquant_query_map": easy_map,
        "k4neo_input": k4_input,
        "k4neo_query_map": k4_map,
    }
