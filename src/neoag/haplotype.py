from __future__ import annotations

import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from .utils import safe_id


def parse_variant(token: str) -> dict[str, Any]:
    match = re.fullmatch(r"(\d+):([A-Za-z]+)>([A-Za-z]+)", token.strip())
    if not match:
        raise ValueError(f"Invalid variant {token!r}; expected POS:REF>ALT")
    return {"pos": int(match.group(1)), "ref": match.group(2).upper(), "alt": match.group(3).upper()}


def read_bases_at_positions(
    fields: Sequence[str], positions: Sequence[int], *, min_baseq: int = 20
) -> dict[int, str]:
    start = int(fields[3])
    cigar, seq, qual = fields[5], fields[9], fields[10]
    ref_pos = start
    query_pos = 0
    bases: dict[int, str] = {}
    for length_text, operation in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        length = int(length_text)
        if operation in "M=X":
            for target in positions:
                if ref_pos <= target < ref_pos + length:
                    index = query_pos + target - ref_pos
                    if index < len(seq) and ord(qual[index]) - 33 >= min_baseq:
                        bases[target] = seq[index].upper()
            ref_pos += length
            query_pos += length
        elif operation in "DN":
            ref_pos += length
        elif operation in "IS":
            query_pos += length
    return bases


def phase_bam_region(
    bam: str | Path,
    chrom: str,
    variants: Sequence[Mapping[str, Any]],
    *,
    samtools: str = "samtools",
    min_mapq: int = 20,
    min_baseq: int = 20,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    positions = sorted(int(v["pos"]) for v in variants)
    region = f"{chrom}:{positions[0]}-{positions[-1]}"
    process = subprocess.Popen(
        [samtools, "view", str(bam), region], stdout=subprocess.PIPE, text=True
    )
    fragment_calls: dict[str, list[dict[str, Any]]] = defaultdict(list)
    assert process.stdout is not None
    for line in process.stdout:
        fields = line.rstrip("\n").split("\t")
        flag, mapq = int(fields[1]), int(fields[4])
        if flag & (0x100 | 0x200 | 0x400 | 0x800) or mapq < min_mapq:
            continue
        bases = read_bases_at_positions(fields, positions, min_baseq=min_baseq)
        if len(bases) != len(positions):
            continue
        fragment_calls[fields[0]].append({
            "haplotype": "".join(bases[pos] for pos in positions),
            "flag": flag,
            "mapq": mapq,
            "cigar": fields[5],
        })
    if process.wait() != 0:
        raise RuntimeError(f"samtools view failed for {bam} {region}")

    rows: list[dict[str, Any]] = []
    counts: Counter[str] = Counter()
    conflicts = 0
    for read_name, calls in sorted(fragment_calls.items()):
        haplotypes = {call["haplotype"] for call in calls}
        if len(haplotypes) != 1:
            conflicts += 1
            continue
        haplotype = next(iter(haplotypes))
        counts[haplotype] += 1
        rows.append({
            "read_name": read_name,
            "haplotype": haplotype,
            "n_alignments": len(calls),
            "max_mapq": max(call["mapq"] for call in calls),
            "flags": ",".join(str(call["flag"]) for call in calls),
        })

    ref_haplotype = "".join(str(v["ref"]) for v in variants)
    alt_haplotype = "".join(str(v["alt"]) for v in variants)
    alt_support = counts[alt_haplotype]
    nonref_support = sum(n for hap, n in counts.items() if hap != ref_haplotype)
    single_variant_support = nonref_support - alt_support
    if alt_support >= 3 and (nonref_support == 0 or alt_support / nonref_support >= 0.8):
        status = "PHASED_CIS"
    elif alt_support == 0 and single_variant_support >= 3:
        status = "PHASED_NOT_CIS"
    else:
        status = "PHASING_UNRESOLVED"
    summary = {
        "bam": str(bam), "chrom": chrom, "positions": positions,
        "reference_haplotype": ref_haplotype, "alternate_haplotype": alt_haplotype,
        "haplotype_counts": dict(sorted(counts.items())),
        "informative_fragments": sum(counts.values()), "conflicting_fragments": conflicts,
        "cis_alt_fragments": alt_support, "single_variant_fragments": single_variant_support,
        "haplotype_status": status,
        "phase_confidence": "high" if status == "PHASED_CIS" and alt_support >= 4 else "moderate" if status == "PHASED_CIS" else "low",
    }
    return rows, summary


def annotate_nearby_variant_groups(
    events: list[dict[str, Any]], *, max_distance_bp: int = 3
) -> list[dict[str, Any]]:
    candidates: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for event in events:
        try:
            pos = int(event.get("pos", ""))
        except (TypeError, ValueError):
            continue
        if len(str(event.get("ref", ""))) != 1 or len(str(event.get("alt", ""))) != 1:
            continue
        key = (str(event.get("chrom", "")), str(event.get("transcript_id", "")))
        candidates[key].append(event)
    for (chrom, transcript), group in candidates.items():
        group.sort(key=lambda row: int(row["pos"]))
        for left, right in zip(group, group[1:]):
            if int(right["pos"]) - int(left["pos"]) > max_distance_bp:
                continue
            members = [left, right]
            group_id = safe_id(f"phase_{chrom}_{left['pos']}_{right['pos']}_{transcript}")
            component_ids = ";".join(str(row.get("event_id", "")) for row in members)
            for event in members:
                event.setdefault("phase_group_id", group_id)
                event.setdefault("component_event_ids", component_ids)
                if not event.get("haplotype_status"):
                    event["haplotype_status"] = "PHASING_REQUIRED"
                    event["phase_confidence"] = "unresolved"
    return events


def combined_haplotype_peptides(
    component_rows: Sequence[Mapping[str, Any]],
    *, peptide_lengths: Iterable[int] = (8, 9, 10, 11),
) -> tuple[str, str, list[dict[str, Any]]]:
    representatives: dict[str, Mapping[str, Any]] = {}
    for row in component_rows:
        key = str(row.get("variant_key") or row.get("event_id") or "")
        if key and key not in representatives and row.get("minigene"):
            representatives[key] = row
    if len(representatives) < 2:
        raise ValueError("At least two component variants with minigene context are required")

    reference: dict[int, str] = {}
    alternate: dict[int, str] = {}
    event_ids: list[str] = []
    for key, row in representatives.items():
        before, alt_segment, after = str(row["minigene"]).split("|")
        protein_pos = int(str(row["protein_position"]).split("-")[0])
        ref_aa, alt_aa = str(row["amino_acids"]).split("/")
        if len(ref_aa) != 1 or len(alt_aa) != 1 or len(alt_segment) != 1:
            raise ValueError("Combined reconstruction currently supports adjacent single-AA substitutions")
        start = protein_pos - len(before)
        wt_context = before + ref_aa + after
        for offset, aa in enumerate(wt_context):
            coordinate = start + offset
            previous = reference.get(coordinate)
            if previous and previous != aa:
                raise ValueError(f"Conflicting reference context at protein position {coordinate}")
            reference[coordinate] = aa
        alternate[protein_pos] = alt_aa
        event_ids.append(key)

    start, end = min(reference), max(reference)
    wt_protein = "".join(reference[pos] for pos in range(start, end + 1))
    mt_protein = "".join(alternate.get(pos, reference[pos]) for pos in range(start, end + 1))
    changed_indices = [pos - start for pos in sorted(alternate)]
    windows: list[dict[str, Any]] = []
    for length in peptide_lengths:
        for offset in range(0, len(mt_protein) - length + 1):
            if not all(offset <= index < offset + length for index in changed_indices):
                continue
            windows.append({
                "peptide": mt_protein[offset:offset + length],
                "wildtype_peptide": wt_protein[offset:offset + length],
                "peptide_length": length,
                "peptide_start_aa": start + offset,
                "peptide_end_aa": start + offset + length - 1,
                "mutation_positions_in_peptide": ",".join(str(index - offset + 1) for index in changed_indices),
            })
    return ";".join(event_ids), f"{wt_protein}>{mt_protein}", windows
