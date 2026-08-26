#!/usr/bin/env python3
"""Build a provenance-preserving fusion event/peptide union from completed callers."""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

from neoag.input_router import build_raw_intermediates
from neoag.adapters.diagnostic_fusion_rescue import (
    DEFAULT_DIAGNOSTIC_FUSION_WHITELIST,
    diagnostic_rescue_rows_from_easyfuse,
    infer_unfiltered_easyfuse_path,
    normalize_fusion_label,
    write_diagnostic_fusion_rescue,
)
from neoag.model_layers import enrich_event_layers, enrich_peptide_layers, infer_mutation_source, infer_peptide_consequence
from neoag.provenance import merge_rows_preserving_provenance
from neoag.schemas import EVENT_FIELDS, PEPTIDE_FIELDS
from neoag.utils import first, safe_id, to_float, write_tsv

HLA_RE = re.compile(r"(?:HLA-)?(?:A|B|C)\*[0-9]{2,3}(?::[0-9A-Z]{2,3}){1,4}", re.I)
AA_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")

STAR_FUSION_PATTERNS = (
    "**/star-fusion.fusion_predictions.abridged.tsv",
    "**/star-fusion.fusion_predictions.tsv",
    "**/*fusion_predictions.abridged.tsv",
    "**/*fusion_predictions.tsv",
)
ARRIBA_PATTERNS = (
    "**/*.fusions.tsv",
    "**/fusions.tsv",
)
FUSIONCATCHER_PATTERNS = (
    "**/fusioncatcher.final-list.txt",
    "**/final-list_candidate-fusion-genes*.txt",
    "**/final-list_candidate-fusion-genes*",
)
JAFFAL_PATTERNS = (
    "**/jaffa_results.csv",
    "**/jaffal_results.csv",
)
EASYFUSE_PATTERNS = (
    "**/fusions.pass.csv",
)
STAR_CHIMERIC_PATTERNS = (
    "**/Chimeric.out.junction",
)

TARGETED_FUSION_REGIONS = {
    "EWSR1_WT1": {
        "EWSR1": ("chr22", 29268009, 29300525),
        "WT1": ("chr11", 32387775, 32435564),
    },
}


def read_hla(path: Path) -> list[str]:
    values = HLA_RE.findall(path.read_text(encoding="utf-8", errors="replace"))
    return list(dict.fromkeys(value.upper() if value.upper().startswith("HLA-") else "HLA-" + value.upper() for value in values))


def read_table(path: Path) -> list[dict[str, str]]:
    if not path or not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        header = handle.readline()
        handle.seek(0)
        delimiter = "\t" if "\t" in header else ","
        return [{str(key): str(value or "") for key, value in row.items()} for row in csv.DictReader(handle, delimiter=delimiter)]


def clean_gene(value: str) -> str:
    return str(value or "").split("^", 1)[0].strip()


def first(row: dict[str, str], names: list[str], default: str = "") -> str:
    normalized = {str(key).lower().replace("#", ""): str(value or "") for key, value in row.items()}
    for name in names:
        value = normalized.get(name.lower().replace("#", ""), "").strip()
        if value:
            return value
    return default


def normalize_breakpoint(value: str) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    text = text.replace(";", ":").replace("|", ":")
    parts = [part for part in text.split(":") if part]
    if len(parts) >= 2 and parts[-1] in {"+", "-"}:
        parts = parts[:-1]
    return ":".join(parts[:2]) if len(parts) >= 2 else text


def parse_breakpoint(value: str) -> tuple[str, int | None, str]:
    text = str(value or "").strip()
    parts = text.split(":")
    if len(parts) < 2:
        return "", None, ""
    chrom = parts[0]
    if chrom and not chrom.startswith("chr"):
        chrom = "chr" + chrom
    try:
        pos = int(parts[1])
    except ValueError:
        pos = None
    strand = parts[2] if len(parts) > 2 else ""
    return chrom, pos, strand


def _same_breakpoint_pair(
    observed: tuple[str, int, str, int],
    expected: tuple[str, int, str, int],
    tolerance: int,
) -> bool:
    ochr1, opos1, ochr2, opos2 = observed
    echr1, epos1, echr2, epos2 = expected
    direct = ochr1 == echr1 and ochr2 == echr2 and abs(opos1 - epos1) <= tolerance and abs(opos2 - epos2) <= tolerance
    reverse = ochr1 == echr2 and ochr2 == echr1 and abs(opos1 - epos2) <= tolerance and abs(opos2 - epos1) <= tolerance
    return direct or reverse


def verify_event_junction_reads(
    audit: list[dict[str, str]],
    star_chimeric_files: list[Path],
    rna_bam: Path | None,
    *,
    samtools: str = "samtools",
    tolerance: int = 3,
) -> tuple[dict[str, dict[str, object]], list[dict[str, str]]]:
    """Match caller events to STAR junctions, then optionally confirm read names in the RNA BAM."""
    targets: dict[str, set[tuple[str, int, str, int]]] = defaultdict(set)
    for row in audit:
        event_id = str(row.get("event_id") or "").strip()
        chrom1, pos1, _ = parse_breakpoint(row.get("left_breakpoint", ""))
        chrom2, pos2, _ = parse_breakpoint(row.get("right_breakpoint", ""))
        if event_id and chrom1 and pos1 is not None and chrom2 and pos2 is not None:
            targets[event_id].add((chrom1, pos1, chrom2, pos2))
    target_index: dict[tuple[str, str], set[str]] = defaultdict(set)
    for event_id, expected_pairs in targets.items():
        for chrom1, _, chrom2, _ in expected_pairs:
            target_index[tuple(sorted((chrom1, chrom2)))].add(event_id)

    matched_names: dict[str, set[str]] = defaultdict(set)
    matched_sources: dict[str, set[str]] = defaultdict(set)
    for path in star_chimeric_files:
        if not path.is_file() or path.stat().st_size == 0:
            continue
        with path.open(encoding="utf-8", errors="replace") as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 10:
                    continue
                chrom1 = parts[0] if parts[0].startswith("chr") else "chr" + parts[0]
                chrom2 = parts[3] if parts[3].startswith("chr") else "chr" + parts[3]
                try:
                    observed = (chrom1, int(parts[1]), chrom2, int(parts[4]))
                except ValueError:
                    continue
                read_name = parts[9].split()[0]
                if not read_name:
                    continue
                candidate_ids = target_index.get(tuple(sorted((chrom1, chrom2))), set())
                for event_id in candidate_ids:
                    expected_pairs = targets[event_id]
                    if any(_same_breakpoint_pair(observed, expected, tolerance) for expected in expected_pairs):
                        matched_names[event_id].add(read_name)
                        matched_sources[event_id].add(str(path))

    all_names = sorted({name for names in matched_names.values() for name in names})
    bam_names: set[str] = set()
    bam_error = ""
    if rna_bam and rna_bam.is_file() and rna_bam.stat().st_size > 0 and all_names:
        name_file = None
        try:
            samtools_path = shutil.which(samtools)
            sibling_samtools = Path(sys.executable).with_name("samtools")
            if not samtools_path and sibling_samtools.is_file():
                samtools_path = str(sibling_samtools)
            if not samtools_path:
                raise FileNotFoundError(f"samtools executable not found: {samtools}")
            with tempfile.NamedTemporaryFile("w", encoding="utf-8", delete=False) as handle:
                handle.write("\n".join(all_names) + "\n")
                name_file = Path(handle.name)
            result = subprocess.run(
                [samtools_path, "view", "-N", str(name_file), str(rna_bam)],
                check=False,
                capture_output=True,
                text=True,
            )
            if result.returncode == 0:
                bam_names = {line.split("\t", 1)[0] for line in result.stdout.splitlines() if line}
            else:
                bam_error = (result.stderr or f"samtools exited {result.returncode}").strip()
        except (OSError, subprocess.SubprocessError) as exc:
            bam_error = str(exc)
        finally:
            if name_file:
                name_file.unlink(missing_ok=True)

    measurements: dict[str, dict[str, object]] = {}
    sidecar: list[dict[str, str]] = []
    for event_id in sorted(targets):
        star_names = matched_names.get(event_id, set())
        verified_names = star_names & bam_names if bam_names else set()
        if verified_names:
            status = "BAM_VERIFIED"
            method = "star_chimeric_breakpoint_plus_bam_qname"
            count = len(verified_names)
        elif star_names:
            status = "STAR_JUNCTION_VERIFIED"
            method = "star_chimeric_breakpoint"
            count = len(star_names)
        else:
            status = "NO_EXACT_JUNCTION_MATCH"
            method = "star_chimeric_breakpoint"
            count = 0
        source_parts = sorted(matched_sources.get(event_id, set()))
        if verified_names and rna_bam:
            source_parts.append(str(rna_bam))
        note = f"coordinate_tolerance_bp={tolerance}; unique_STAR_read_names={len(star_names)}"
        if bam_error:
            note += f"; BAM verification unavailable: {bam_error}"
        measurements[event_id] = {
            "verified_count": count,
            "status": status,
            "method": method,
            "source": ";".join(source_parts),
            "note": note,
        }
        sidecar.append({
            "event_id": event_id,
            "verified_rna_junction_reads": str(count) if count else "",
            "caller_rna_junction_reads": "",
            "junction_match_status": status,
            "junction_match_method": method,
            "junction_verification_source": ";".join(source_parts),
            "junction_verification_note": note,
        })
    return measurements, sidecar


def apply_junction_measurements(
    rows: list[dict[str, str]], measurements: dict[str, dict[str, object]]
) -> None:
    for row in rows:
        measurement = measurements.get(str(row.get("event_id") or ""))
        if not measurement or not int(measurement.get("verified_count") or 0):
            continue
        caller_count = str(row.get("provided_rna_junction_reads") or row.get("rna_junction_reads") or "").strip()
        if caller_count:
            row["provided_rna_junction_reads"] = caller_count
        row["rna_junction_reads"] = str(measurement["verified_count"])
        row["junction_match_status"] = str(measurement["status"])
        row["junction_match_method"] = str(measurement["method"])
        row["rna_junction_source"] = str(measurement["source"])


def infer_star_chimeric_from_junctions(path: Path | None) -> Path | None:
    if not path:
        return None
    p = Path(path)
    if p.name == "Chimeric.out.junction" and p.is_file():
        return p
    candidates = [p.with_name("Chimeric.out.junction")]
    for candidate in candidates:
        if candidate.is_file() and candidate.stat().st_size > 0:
            return candidate
    return None


def star_chimeric_support_count(path: Path | None, fusion_label: str, bp1: str, bp2: str) -> int:
    if not path or not path.is_file() or path.stat().st_size == 0:
        return 0
    norm = normalize_fusion_label(fusion_label)
    regions = TARGETED_FUSION_REGIONS.get(norm)
    if not regions:
        return 0
    genes = [gene for gene in norm.split("_") if gene in regions]
    if len(genes) != 2:
        return 0
    left_region = regions[genes[0]]
    right_region = regions[genes[1]]

    def in_region(chrom: str, pos: int, region: tuple[str, int, int]) -> bool:
        return chrom == region[0] and region[1] <= pos <= region[2]

    count = 0
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom1 = parts[0] if parts[0].startswith("chr") else "chr" + parts[0]
            chrom2 = parts[3] if parts[3].startswith("chr") else "chr" + parts[3]
            try:
                pos1 = int(parts[1])
                pos2 = int(parts[4])
            except ValueError:
                continue
            if (
                in_region(chrom1, pos1, left_region) and in_region(chrom2, pos2, right_region)
            ) or (
                in_region(chrom1, pos1, right_region) and in_region(chrom2, pos2, left_region)
            ):
                count += 1
    return count


def peptide_windows(sequence: str, lengths: tuple[int, ...] = (8, 9, 10, 11)) -> list[str]:
    seq = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", str(sequence or "").upper())
    if not AA_RE.fullmatch(seq or ""):
        return []
    if min(lengths) <= len(seq) <= max(lengths):
        return [seq]
    return list(dict.fromkeys(seq[start:start + length] for length in lengths for start in range(max(0, len(seq) - length + 1))))


def breakpoint_window_records(
    sequence: str,
    breakpoint: int,
    *,
    left_gene: str = "",
    right_gene: str = "",
    lengths: tuple[int, ...] = (8, 9, 10, 11),
) -> list[dict[str, str]]:
    """Return peptide windows with at least one residue from each fusion side."""
    seq = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", str(sequence or "").upper())
    if not seq or breakpoint <= 0 or breakpoint >= len(seq):
        return []
    records: list[dict[str, str]] = []
    seen: set[tuple[str, int]] = set()
    for length in lengths:
        first_start = max(0, breakpoint - length + 1)
        last_start = min(breakpoint - 1, len(seq) - length)
        for start in range(first_start, last_start + 1):
            peptide = seq[start:start + length]
            left_count = breakpoint - start
            if len(peptide) != length or not (1 <= left_count < length):
                continue
            key = (peptide, left_count)
            if key in seen:
                continue
            seen.add(key)
            left = peptide[:left_count]
            right = peptide[left_count:]
            records.append({
                "peptide": peptide,
                "crosses_junction": "yes",
                "contains_novel_aa": "yes",
                "junction_position_in_peptide_1based": str(left_count),
                "fusion_left_gene": left_gene,
                "fusion_right_gene": right_gene,
                "fusion_left_peptide": left,
                "fusion_right_peptide": right,
                "fusion_junction_display": f"{left}|{right}",
                "fusion_peptide_classification": "JUNCTION_SPANNING",
            })
    return records


def targeted_rescue_rows(
    easyfuse_files: list[Path],
    *,
    star_chimeric_files: list[Path],
    sample_id: str,
    profile: str,
    hla: list[str],
    whitelist: list[str],
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    events: list[dict[str, str]] = []
    peptides: list[dict[str, str]] = []
    audit: list[dict[str, str]] = []
    rescue_sidecar: list[dict[str, str]] = []
    seen: set[str] = set()
    norm_whitelist = [normalize_fusion_label(item) for item in whitelist]
    for source in easyfuse_files:
        rows = diagnostic_rescue_rows_from_easyfuse(
            source,
            sample_id=sample_id,
            whitelist=norm_whitelist,
            pass_keys=set(),
        )
        for row in rows:
            fusion_norm = row.get("fusion_gene_normalized", "")
            bp1 = row.get("breakpoint1", "")
            bp2 = row.get("breakpoint2", "")
            event_id = safe_id(f"TARGETED_RESCUE|{sample_id}|{fusion_norm}|{bp1}|{bp2}")
            if event_id in seen:
                continue
            seen.add(event_id)
            star_support = sum(star_chimeric_support_count(path, fusion_norm, bp1, bp2) for path in star_chimeric_files)
            rescue_status = "STAR-junction-supported" if star_support > 0 else "single-caller"
            reads = str(max(int(to_float(row.get("rna_junction_reads"), 0.0)), star_support))
            gene_pair_display = row.get("fusion_gene", "").replace("_", "::")
            chrom, pos, _strand = parse_breakpoint(bp1)
            confidence = "0.85" if star_support > 0 else "0.65"
            base = {field: "" for field in EVENT_FIELDS}
            base.update({
                "event_id": event_id,
                "sample_id": sample_id,
                "disease_profile": profile,
                "event_type": "Fusion",
                "mutation_source": infer_mutation_source(event_type="Fusion", tool="TARGETED_RESCUE"),
                "peptide_consequence": infer_peptide_consequence(event_type="Fusion", consequence="fusion", tool="TARGETED_RESCUE"),
                "gene": gene_pair_display,
                "event_name": gene_pair_display,
                "chrom": chrom,
                "pos": str(pos or ""),
                "transcript_id": row.get("ftid", ""),
                "consequence": row.get("frame_status", "") or "fusion",
                "rna_junction_reads": reads,
                "rna_junction_source": rescue_status,
                "event_confidence": confidence,
                "event_expression": "0.0",
                "driver_relevance": "1.0",
                "tumor_vaf": "0.0",
                "clonality": "0.5",
                "persistence": "0.5",
                "tumor_specificity": "0.8",
                "source": f"TARGETED_RESCUE:{source}",
                "source_file": str(source),
                "source_record_id": row.get("rescue_id", ""),
                "source_tools": "TARGETED_RESCUE,EasyFuseRaw" + (",STAR-Chimeric" if star_support > 0 else ""),
            })
            events.append(enrich_event_layers(base))
            breakpoint = int(to_float(row.get("neo_peptide_sequence_bp"), 0.0))
            side_genes = gene_pair_display.split("::", 1)
            windows = breakpoint_window_records(
                row.get("neo_peptide_sequence", ""), breakpoint,
                left_gene=side_genes[0] if side_genes else "",
                right_gene=side_genes[1] if len(side_genes) > 1 else "",
            )
            for window in windows:
                peptide = window["peptide"]
                for allele in hla:
                    pbase = {field: "" for field in PEPTIDE_FIELDS}
                    pbase.update({
                        "peptide_id": safe_id(f"{event_id}|{allele}|{peptide}"),
                        "event_id": event_id,
                        "sample_id": sample_id,
                        "event_type": "Fusion",
                        "mutation_source": base["mutation_source"],
                        "peptide_consequence": base["peptide_consequence"],
                        "gene": gene_pair_display,
                        "peptide": peptide,
                        "hla_allele": allele,
                        "mhc_class": "I",
                        "source_tool": "TARGETED_RESCUE",
                        "source_file": str(source),
                        "source_record_id": row.get("rescue_id", ""),
                        **window,
                        "transcript_id": row.get("ftid", ""),
                        "fusion_orf_comparison_status": "TRACEABLE_TO_CALLER_TRANSCRIPT" if row.get("ftid") else "ORF_TRANSCRIPT_UNASSESSED",
                        "rna_junction_reads": reads,
                        "rna_junction_source": rescue_status,
                        "binding_rank": "99",
                        "el_rank": "99",
                        "presentation_score": "0.0",
                        "immunogenicity_score": "0.5",
                        "wildtype_binding_rank": "99",
                        "self_similarity_score": "0.0",
                        "normal_hla_ligand_overlap": "no",
                    })
                    peptides.append(enrich_peptide_layers(pbase, base))
            peptide_status = "TARGETED_RESCUE:" + rescue_status + (":GENERATED_FROM_RESCUE_ORF" if windows else ":ORF_PEPTIDE_UNAVAILABLE_REVIEW_ONLY")
            audit.append({
                "event_id": event_id,
                "gene_pair": gene_pair_display,
                "left_breakpoint": bp1,
                "right_breakpoint": bp2,
                "direction": rescue_status,
                "source_tool": "TARGETED_RESCUE",
                "source_file": str(source),
                "source_row": row.get("rescue_id", ""),
                "peptide_status": peptide_status,
            })
            rescue_sidecar.append({
                **row,
                "rescue_reason": rescue_status,
                "peptide_generation_status": "generated_for_ranking" if windows else "not_generated_no_rescue_orf",
                "rna_junction_reads": reads,
                "notes": (row.get("notes", "") + f" TARGETED_RESCUE included in fusion peptide generation; STAR Chimeric support={star_support}.").strip(),
            })
    return events, peptides, audit, rescue_sidecar


def existing(paths: list[Path | None]) -> list[Path]:
    seen: set[Path] = set()
    result: list[Path] = []
    for path in paths:
        if not path:
            continue
        resolved = path.expanduser()
        if resolved.is_file() and resolved.stat().st_size > 0 and resolved not in seen:
            seen.add(resolved)
            result.append(resolved)
    return result


def discover_files(roots: list[Path], patterns: tuple[str, ...]) -> list[Path]:
    found: list[Path] = []
    seen: set[Path] = set()
    for root in roots:
        if not root or not root.exists():
            continue
        search_root = root if root.is_dir() else root.parent
        for pattern in patterns:
            for path in sorted(search_root.glob(pattern)):
                if path.is_file() and path.stat().st_size > 0 and path not in seen:
                    seen.add(path)
                    found.append(path)
    return found


def gene_pair(row: dict[str, str]) -> tuple[str, str]:
    combined = first(row, ["FusionName", "#FusionName", "fusion", "fusion_name", "Fusion_Gene", "fusion genes", "fusion_genes"], "")
    for sep in ("::", "--", "_"):
        if sep in combined:
            left, right = combined.split(sep, 1)
            return clean_gene(left), clean_gene(right)
    return (
        clean_gene(first(row, ["LeftGene", "left_gene", "gene1", "Gene1", "gene5"], "")),
        clean_gene(first(row, ["RightGene", "right_gene", "gene2", "Gene2", "gene3"], "")),
    )


def generic_caller_rows(path: Path, tool: str, sample_id: str, profile: str, hla: list[str]):
    events, peptides, audit = [], [], []
    for index, row in enumerate(read_table(path), 1):
        if tool == "JAFFAL":
            classification = first(row, ["classification", "confidence"], "")
            if classification and classification.lower() != "highconfidence":
                continue
        if tool == "Arriba" and first(row, ["confidence"], "").strip().lower() == "low":
            continue
        filter_status = first(row, ["filter", "filter_status", "status"], "").strip().lower()
        if filter_status and filter_status in {"fail", "failed", "reject", "rejected", "filtered"}:
            continue
        left_gene, right_gene = gene_pair(row)
        if not left_gene or not right_gene:
            continue
        left_bp = first(row, ["LeftBreakpoint", "breakpoint1", "left_breakpoint", "breakpoint_1"], "")
        right_bp = first(row, ["RightBreakpoint", "breakpoint2", "right_breakpoint", "breakpoint_2"], "")
        if not left_bp:
            chrom, base = first(row, ["chrom1", "chr1"], ""), first(row, ["base1", "position1", "pos1"], "")
            left_bp = f"{chrom}:{base}" if chrom and base else ""
        if not right_bp:
            chrom, base = first(row, ["chrom2", "chr2"], ""), first(row, ["base2", "position2", "pos2"], "")
            right_bp = f"{chrom}:{base}" if chrom and base else ""
        direction = first(row, ["direction", "strand", "Strand1", "strand1(gene/fusion)"], "") + "/" + first(row, ["Strand2", "strand2(gene/fusion)"], "")
        pair = f"{left_gene}::{right_gene}"
        norm_left_bp = normalize_breakpoint(left_bp)
        norm_right_bp = normalize_breakpoint(right_bp)
        event_id = safe_id(f"FUSION|{pair}|{norm_left_bp}|{norm_right_bp}")
        reads = first(row, ["JunctionReadCount", "junction_reads", "split_reads", "supporting_reads", "split_reads1", "spanning reads", "spanning_reads", "spanning pairs", "spanning_pairs"], "")
        frame = first(row, ["frame", "reading_frame", "in_frame", "reading_frame_status"], "")
        base = {field: "" for field in EVENT_FIELDS}
        base.update({
            "event_id": event_id, "sample_id": sample_id, "disease_profile": profile,
            "event_type": "Fusion", "gene": pair, "event_name": pair,
            "consequence": frame or "fusion_orf_unassessed", "rna_junction_reads": reads,
            "event_confidence": "0.7", "event_expression": "0.0", "driver_relevance": "0.0",
            "clonality": "0.5", "persistence": "0.5", "tumor_specificity": "0.7",
            "source": f"{tool}:{path}", "source_file": str(path), "source_row_number": str(index),
            "source_tools": tool, "mutation_source": infer_mutation_source(event_type="Fusion", tool=tool),
            "peptide_consequence": infer_peptide_consequence(event_type="Fusion", consequence="fusion", tool=tool),
        })
        events.append(enrich_event_layers(base))
        sequence = first(row, ["junction_peptide", "neo_peptide_sequence", "fusion_peptide", "mutant_peptide", "peptide", "peptide_sequence"], "").upper()
        sequence = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", sequence)
        explicit_hla = first(row, ["hla", "hla_allele", "allele"], "")
        row_hla = read_hla_text(explicit_hla) or hla
        window_records: list[dict[str, str]] = []
        breakpoint = int(to_float(first(row, ["neo_peptide_sequence_bp", "junction_position_in_peptide_1based", "junction_offset_in_peptide"], "0"), 0.0))
        if AA_RE.fullmatch(sequence or ""):
            if breakpoint:
                window_records = breakpoint_window_records(
                    sequence, breakpoint, left_gene=left_gene, right_gene=right_gene,
                )
            else:
                provided = [sequence] if 8 <= len(sequence) <= 12 else peptide_windows(sequence)
                window_records = [{
                    "peptide": peptide,
                    "crosses_junction": "UNASSESSED",
                    "contains_novel_aa": "UNASSESSED",
                    "fusion_left_gene": left_gene,
                    "fusion_right_gene": right_gene,
                    "fusion_peptide_classification": "BOUNDARY_UNASSESSED",
                } for peptide in provided]
        for window in window_records:
            peptide = window["peptide"]
            for allele in row_hla:
                pbase = {field: "" for field in PEPTIDE_FIELDS}
                pbase.update({
                    "peptide_id": safe_id(f"{event_id}|{allele}|{peptide}"), "event_id": event_id,
                    "sample_id": sample_id, "event_type": "Fusion", "gene": pair,
                    "peptide": peptide, "hla_allele": allele, "mhc_class": "I",
                    "source_tool": tool, "source_file": str(path), **window,
                    "rna_junction_reads": reads,
                    "transcript_id": first(row, ["transcript_id", "ftid", "transcript"], ""),
                    "orf_id": first(row, ["orf_id", "fusion_orf_id"], ""),
                    "fusion_orf_comparison_status": "TRACEABLE_TO_CALLER_ORF" if first(row, ["orf_id", "fusion_orf_id", "transcript_id", "ftid"], "") else "ORF_TRANSCRIPT_UNASSESSED",
                    "mutation_source": base["mutation_source"], "peptide_consequence": base["peptide_consequence"],
                    "binding_rank": "99", "el_rank": "99", "presentation_score": "0.0",
                    "immunogenicity_score": "0.5", "wildtype_binding_rank": "99", "self_similarity_score": "0.0",
                })
                peptides.append(enrich_peptide_layers(pbase))
        audit.append({"event_id": event_id, "gene_pair": pair, "left_breakpoint": left_bp, "right_breakpoint": right_bp, "direction": direction, "source_tool": tool, "source_file": str(path), "source_row": str(index), "peptide_status": ("JUNCTION_SPANNING_PEPTIDE" if breakpoint and window_records else "BOUNDARY_UNASSESSED_PEPTIDE" if window_records else "ORF_PEPTIDE_UNAVAILABLE_REVIEW_ONLY"), "admission_policy": "CALLER_PASS_OR_INDEPENDENT_CALLER", "rescue_reason": ""})
    return events, peptides, audit


def _breakpoint_windows(sequence: str, breakpoint: int, lengths: tuple[int, ...] = (8, 9, 10, 11)) -> list[str]:
    """Return only peptide windows that contain residues on both sides of a junction."""
    return [row["peptide"] for row in breakpoint_window_records(sequence, breakpoint, lengths=lengths)]


def diagnostic_rescue_entities(
    rows: list[dict[str, str]], sample_id: str, profile: str, hla: list[str]
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    """Promote exact-whitelist rescue rows into the formal union with an auditable cap."""
    events: list[dict[str, str]] = []
    peptides: list[dict[str, str]] = []
    audit: list[dict[str, str]] = []
    for row in rows:
        pair = str(row.get("fusion_gene") or "").strip()
        bp1, bp2 = str(row.get("breakpoint1") or ""), str(row.get("breakpoint2") or "")
        event_id = safe_id(f"FUSION|{pair}|{normalize_breakpoint(bp1)}|{normalize_breakpoint(bp2)}")
        chrom, pos = "", ""
        parts = normalize_breakpoint(bp1).split(":", 1)
        if len(parts) == 2:
            chrom, pos = parts
        event = {field: "" for field in EVENT_FIELDS}
        event.update({
            "event_id": event_id,
            "sample_id": sample_id,
            "disease_profile": profile,
            "event_type": "Fusion",
            "mutation_source": infer_mutation_source(event_type="Fusion", tool="EasyFuse"),
            "peptide_consequence": "fusion",
            "evidence_scope": "DIAGNOSTIC_WHITELIST_RESCUE",
            "priority_cap": "C_CAUTION",
            "gene": pair,
            "event_name": pair,
            "chrom": chrom,
            "pos": pos,
            "transcript_id": row.get("ftid", ""),
            "consequence": row.get("frame_status") or "fusion_orf_unassessed",
            "rna_junction_reads": row.get("rna_junction_reads", ""),
            "event_confidence": "0.600",
            "driver_relevance": "1.0",
            "clonality": "0.5",
            "persistence": "0.5",
            "tumor_specificity": "0.7",
            "source": f"EasyFuse:diagnostic_whitelist_rescue:{row.get('source_file', '')}",
        })
        events.append(enrich_event_layers(event))
        sequence = row.get("neo_peptide_sequence", "")
        breakpoint = int(to_float(row.get("neo_peptide_sequence_bp"), 0.0))
        side_genes = pair.split("::", 1)
        windows = breakpoint_window_records(
            sequence, breakpoint,
            left_gene=side_genes[0] if side_genes else "",
            right_gene=side_genes[1] if len(side_genes) > 1 else "",
        )
        for window in windows:
            peptide = window["peptide"]
            for allele in hla:
                pbase = {field: "" for field in PEPTIDE_FIELDS}
                pbase.update({
                    "peptide_id": safe_id(f"{event_id}|{allele}|{peptide}"),
                    "event_id": event_id,
                    "sample_id": sample_id,
                    "event_type": "Fusion",
                    "mutation_source": event["mutation_source"],
                    "peptide_consequence": "fusion",
                    "evidence_scope": "DIAGNOSTIC_WHITELIST_RESCUE",
                    "priority_cap": "C_CAUTION",
                    "gene": pair,
                    "peptide": peptide,
                    **window,
                    "transcript_id": row.get("ftid", ""),
                    "fusion_orf_comparison_status": "TRACEABLE_TO_CALLER_TRANSCRIPT" if row.get("ftid") else "ORF_TRANSCRIPT_UNASSESSED",
                    "rna_junction_reads": row.get("rna_junction_reads", ""),
                    "hla_allele": allele,
                    "mhc_class": "I",
                    "source_tool": "EasyFuse_diagnostic_whitelist_rescue",
                    "binding_rank": "99",
                    "el_rank": "99",
                    "presentation_score": "0.0",
                    "immunogenicity_score": "0.5",
                    "wildtype_binding_rank": "99",
                    "self_similarity_score": "0.0",
                })
                peptides.append(enrich_peptide_layers(pbase))
        audit.append({
            "event_id": event_id,
            "gene_pair": pair,
            "left_breakpoint": bp1,
            "right_breakpoint": bp2,
            "direction": "",
            "source_tool": "EasyFuse",
            "source_file": row.get("source_file", ""),
            "source_row": "",
            "peptide_status": "JUNCTION_WINDOWS_GENERATED" if windows else "ORF_PEPTIDE_UNAVAILABLE_REVIEW_ONLY",
            "admission_policy": "DIAGNOSTIC_WHITELIST_RESCUE",
            "rescue_reason": row.get("rescue_reason", ""),
        })
    return events, peptides, audit


def read_hla_text(value: str) -> list[str]:
    return list(dict.fromkeys(match.upper() if match.upper().startswith("HLA-") else "HLA-" + match.upper() for match in HLA_RE.findall(value or "")))


def write_consensus(path: Path, audit: list[dict[str, str]]) -> None:
    grouped: dict[str, dict[str, object]] = defaultdict(lambda: {
        "gene_pair": "",
        "tools": set(),
        "left_breakpoints": set(),
        "right_breakpoints": set(),
        "peptide_statuses": set(),
    })
    for row in audit:
        event_id = row.get("event_id", "")
        if not event_id:
            continue
        item = grouped[event_id]
        item["gene_pair"] = item["gene_pair"] or row.get("gene_pair", "")
        item["tools"].add(row.get("source_tool", ""))  # type: ignore[union-attr]
        if row.get("left_breakpoint"):
            item["left_breakpoints"].add(row["left_breakpoint"])  # type: ignore[union-attr]
        if row.get("right_breakpoint"):
            item["right_breakpoints"].add(row["right_breakpoint"])  # type: ignore[union-attr]
        if row.get("peptide_status"):
            item["peptide_statuses"].add(row["peptide_status"])  # type: ignore[union-attr]
    rows: list[dict[str, str]] = []
    for event_id, item in grouped.items():
        tools = sorted(tool for tool in item["tools"] if tool)  # type: ignore[union-attr]
        rows.append({
            "event_id": event_id,
            "fusion": str(item["gene_pair"]),
            "support_tools": ",".join(tools),
            "n_tools": str(len(tools)),
            "left_breakpoints": ";".join(sorted(item["left_breakpoints"])),  # type: ignore[arg-type]
            "right_breakpoints": ";".join(sorted(item["right_breakpoints"])),  # type: ignore[arg-type]
            "peptide_status": ";".join(sorted(item["peptide_statuses"])),  # type: ignore[arg-type]
            "status": "CROSS_VALIDATED" if len(tools) >= 2 else "SINGLE_TOOL",
        })
    rows.sort(key=lambda row: (-int(row["n_tools"]), row["fusion"]))
    write_tsv(path, rows, ["event_id", "fusion", "support_tools", "n_tools", "left_breakpoints", "right_breakpoints", "peptide_status", "status"])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True); parser.add_argument("--profile", default="default")
    parser.add_argument("--hla-file", required=True, type=Path); parser.add_argument("--easyfuse", type=Path)
    parser.add_argument("--easyfuse-unfiltered", type=Path)
    parser.add_argument("--diagnostic-fusion-whitelist", action="append", default=[])
    parser.add_argument("--disable-diagnostic-fusion-rescue", action="store_true")
    parser.add_argument("--star-fusion", type=Path); parser.add_argument("--arriba", type=Path)
    parser.add_argument("--fusioncatcher", type=Path)
    parser.add_argument("--easyfuse-all", action="append", type=Path, default=[])
    parser.add_argument("--star-chimeric", action="append", type=Path, default=[])
    parser.add_argument("--rna-bam", type=Path, help="RNA BAM used to confirm STAR junction read names")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--junction-coordinate-tolerance", type=int, default=3)
    parser.add_argument("--targeted-fusion-rescue", action="store_true", default=True)
    parser.add_argument("--no-targeted-fusion-rescue", action="store_false", dest="targeted_fusion_rescue")
    parser.add_argument("--targeted-fusion-whitelist", default=",".join(DEFAULT_DIAGNOSTIC_FUSION_WHITELIST))
    parser.add_argument("--jaffal", type=Path)
    parser.add_argument("--caller-root", action="append", type=Path, default=[])
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    hla = read_hla(args.hla_file)
    if not hla:
        raise SystemExit("HLA consensus has no class-I alleles")
    args.outdir.mkdir(parents=True, exist_ok=True)
    events, peptides, audit = [], [], []
    roots = [root for root in args.caller_root if root]
    easyfuse_files = existing([args.easyfuse]) or discover_files(roots, EASYFUSE_PATTERNS)
    easyfuse_all_files = existing(args.easyfuse_all)
    for easyfuse in easyfuse_files:
        inferred = infer_unfiltered_easyfuse_path(easyfuse)
        if inferred:
            easyfuse_all_files.extend(existing([inferred]))
    easyfuse_all_files.extend([
        path for path in discover_files(roots, ("**/fusions.csv",))
        if path.name != "fusions.pass.csv"
    ])
    star_fusion_files = existing([args.star_fusion]) + discover_files(roots, STAR_FUSION_PATTERNS)
    arriba_files = existing([args.arriba]) + [
        path for path in discover_files(roots, ARRIBA_PATTERNS)
        if path.name != "fusions.pass.csv"
    ]
    fusioncatcher_files = existing([args.fusioncatcher]) + discover_files(roots, FUSIONCATCHER_PATTERNS)
    star_chimeric_files = existing(args.star_chimeric) + discover_files(roots, STAR_CHIMERIC_PATTERNS)
    jaffal_files = existing([args.jaffal]) + discover_files(roots, JAFFAL_PATTERNS)

    for easyfuse in existing(easyfuse_files):
        cfg = {"sample": {"id": args.sample_id, "profile": args.profile}, "inputs": {"entry_mode": "fusion", "easyfuse_tsv": str(easyfuse.resolve()), "hla_alleles": hla}}
        easyfuse_out = args.outdir / "easyfuse"
        build_raw_intermediates(cfg, easyfuse_out, root=Path.cwd())
        easyfuse_events = read_table(easyfuse_out / "parsed/raw_events.tsv")
        easyfuse_source_rows = read_table(easyfuse)
        events.extend(easyfuse_events)
        peptides.extend(read_table(easyfuse_out / "parsed/raw_peptides.tsv"))
        for event_index, event in enumerate(easyfuse_events, 1):
            source_row_number = int(to_float(event.get("source_row_number"), event_index))
            source_row = easyfuse_source_rows[source_row_number - 1] if 0 < source_row_number <= len(easyfuse_source_rows) else {}
            audit.append({
                "event_id": event.get("event_id", ""),
                "gene_pair": event.get("gene", ""),
                "left_breakpoint": first(source_row, ["Breakpoint1", "breakpoint1", "LeftBreakpoint", "left_breakpoint"], ""),
                "right_breakpoint": first(source_row, ["Breakpoint2", "breakpoint2", "RightBreakpoint", "right_breakpoint"], ""),
                "direction": "",
                "source_tool": "EasyFuse",
                "source_file": str(easyfuse),
                "source_row": str(source_row_number),
                "peptide_status": "GENERATED_FROM_EASYFUSE_ORF",
                "admission_policy": "CALLER_PASS",
                "rescue_reason": "",
            })
    rescue_rows: list[dict[str, str]] = []
    rescue_sources: list[Path] = []
    if not args.disable_diagnostic_fusion_rescue:
        whitelist = args.diagnostic_fusion_whitelist or list(DEFAULT_DIAGNOSTIC_FUSION_WHITELIST)
        whitelist = sorted({normalize_fusion_label(value) for value in whitelist if normalize_fusion_label(value)})
        rescue_sources = existing([args.easyfuse_unfiltered])
        if not rescue_sources:
            rescue_sources = existing([infer_unfiltered_easyfuse_path(path) for path in easyfuse_files])
        pass_keys: set[tuple[str, str, str]] = set()
        for path in easyfuse_files:
            for row in read_table(path):
                left_gene, right_gene = gene_pair(row)
                label = normalize_fusion_label(f"{left_gene}::{right_gene}")
                pass_keys.add((label, first(row, ["Breakpoint1", "breakpoint1"], ""), first(row, ["Breakpoint2", "breakpoint2"], "")))
        for path in rescue_sources:
            rescue_rows.extend(diagnostic_rescue_rows_from_easyfuse(path, sample_id=args.sample_id, whitelist=whitelist, pass_keys=pass_keys))
        rescued_events, rescued_peptides, rescued_audit = diagnostic_rescue_entities(rescue_rows, args.sample_id, args.profile, hla)
        existing_event_ids = {row.get("event_id", "") for row in events}
        rescued_event_ids = {row.get("event_id", "") for row in rescued_events if row.get("event_id", "") not in existing_event_ids}
        events.extend(row for row in rescued_events if row.get("event_id", "") in rescued_event_ids)
        peptides.extend(row for row in rescued_peptides if row.get("event_id", "") in rescued_event_ids)
        audit.extend(row for row in rescued_audit if row.get("event_id", "") in rescued_event_ids)
    write_diagnostic_fusion_rescue(
        rescue_rows,
        args.outdir / "diagnostic_fusion_rescue.tsv",
        source_path=rescue_sources[0] if rescue_sources else None,
    )
    for paths, tool in ((star_fusion_files, "STAR-Fusion"), (arriba_files, "Arriba"), (fusioncatcher_files, "FusionCatcher"), (jaffal_files, "JAFFAL")):
        for path in existing(paths):
            e, p, a = generic_caller_rows(path, tool, args.sample_id, args.profile, hla)
            events.extend(e); peptides.extend(p); audit.extend(a)
    rescue_sidecar: list[dict[str, str]] = []
    if args.targeted_fusion_rescue:
        whitelist = [item.strip() for item in args.targeted_fusion_whitelist.replace(";", ",").split(",") if item.strip()]
        rescue_events, rescue_peptides, rescue_audit, rescue_sidecar = targeted_rescue_rows(
            existing(easyfuse_all_files),
            star_chimeric_files=existing(star_chimeric_files),
            sample_id=args.sample_id,
            profile=args.profile,
            hla=hla,
            whitelist=whitelist,
        )
        events.extend(rescue_events)
        peptides.extend(rescue_peptides)
        audit.extend(rescue_audit)
    if not events:
        raise SystemExit("No fusion events were parsed from supplied caller outputs")
    measurements, verification_rows = verify_event_junction_reads(
        audit,
        existing(star_chimeric_files),
        args.rna_bam,
        samtools=args.samtools,
        tolerance=max(0, args.junction_coordinate_tolerance),
    )
    caller_counts: dict[str, set[str]] = defaultdict(set)
    for row in events + peptides:
        event_id = str(row.get("event_id") or "")
        value = str(row.get("provided_rna_junction_reads") or row.get("rna_junction_reads") or "").strip()
        if event_id and value:
            caller_counts[event_id].add(value)
    for row in verification_rows:
        row["caller_rna_junction_reads"] = ";".join(sorted(caller_counts.get(row["event_id"], set())))
    apply_junction_measurements(events, measurements)
    apply_junction_measurements(peptides, measurements)
    merged_events, _, _ = merge_rows_preserving_provenance(events, EVENT_FIELDS, ("event_id",), entity_type="fusion_union_event")
    merged_peptides, _, _ = merge_rows_preserving_provenance(peptides, PEPTIDE_FIELDS, ("event_id", "peptide", "hla_allele"), entity_type="fusion_union_peptide")
    write_tsv(args.outdir / "raw_events.tsv", merged_events, EVENT_FIELDS)
    write_tsv(args.outdir / "raw_peptides.tsv", merged_peptides, PEPTIDE_FIELDS)
    write_tsv(args.outdir / "parsed/raw_events.tsv", merged_events, EVENT_FIELDS)
    write_tsv(args.outdir / "parsed/raw_peptides.tsv", merged_peptides, PEPTIDE_FIELDS)
    write_tsv(args.outdir / "fusion_caller_union.tsv", audit, ["event_id", "gene_pair", "left_breakpoint", "right_breakpoint", "direction", "source_tool", "source_file", "source_row", "peptide_status", "admission_policy", "rescue_reason"])
    write_tsv(
        args.outdir / "junction_read_verification.tsv",
        verification_rows,
        ["event_id", "verified_rna_junction_reads", "caller_rna_junction_reads", "junction_match_status", "junction_match_method", "junction_verification_source", "junction_verification_note"],
    )
    write_consensus(args.outdir / "fusion_consensus.tsv", audit)
    if rescue_sidecar:
        fields = [
            "rescue_id", "sample_id", "fusion_gene", "fusion_gene_raw", "fusion_gene_normalized",
            "gene5", "gene3", "breakpoint1", "breakpoint2", "ftid", "fusion_type",
            "frame_status", "neo_peptide_sequence", "neo_peptide_sequence_bp",
            "rna_junction_reads", "rna_spanning_reads", "anchor_size", "star_detected",
            "fusioncatcher_detected", "arriba_detected", "tools_detected", "tool_count",
            "prediction_class", "prediction_prob", "easyfuse_pass_status",
            "diagnostic_whitelist_status", "diagnostic_relevance", "rescue_reason",
            "peptide_generation_status", "source_file", "notes",
        ]
        write_tsv(args.outdir / "targeted_fusion_rescue.tsv", rescue_sidecar, fields)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
