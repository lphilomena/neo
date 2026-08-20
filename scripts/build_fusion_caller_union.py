#!/usr/bin/env python3
"""Build a provenance-preserving fusion event/peptide union from completed callers."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

from neoag.input_router import build_raw_intermediates
from neoag.adapters.diagnostic_fusion_rescue import (
    DEFAULT_DIAGNOSTIC_FUSION_WHITELIST,
    diagnostic_rescue_rows_from_easyfuse,
    infer_unfiltered_easyfuse_path,
    normalize_fusion_label,
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
            windows = peptide_windows(row.get("neo_peptide_sequence", ""))
            for peptide in windows:
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
                        "crosses_junction": "true",
                        "contains_novel_aa": "true",
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
    combined = first(row, ["FusionName", "#FusionName", "fusion", "fusion_name", "Fusion_Gene"], "")
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
        left_gene, right_gene = gene_pair(row)
        if not left_gene or not right_gene:
            continue
        left_bp = first(row, ["LeftBreakpoint", "breakpoint1", "left_breakpoint", "breakpoint_1"], "")
        right_bp = first(row, ["RightBreakpoint", "breakpoint2", "right_breakpoint", "breakpoint_2"], "")
        direction = first(row, ["direction", "strand", "Strand1", "strand1(gene/fusion)"], "") + "/" + first(row, ["Strand2", "strand2(gene/fusion)"], "")
        pair = f"{left_gene}::{right_gene}"
        norm_left_bp = normalize_breakpoint(left_bp)
        norm_right_bp = normalize_breakpoint(right_bp)
        event_id = safe_id(f"FUSION|{pair}|{norm_left_bp}|{norm_right_bp}")
        reads = first(row, ["JunctionReadCount", "junction_reads", "split_reads", "supporting_reads", "split_reads1"], "")
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
        windows: list[str] = []
        if AA_RE.fullmatch(sequence or ""):
            if 8 <= len(sequence) <= 12:
                windows = [sequence]
            else:
                windows = list(dict.fromkeys(sequence[start:start + length] for length in (8, 9, 10, 11) for start in range(max(0, len(sequence) - length + 1))))
        for peptide in windows:
            for allele in row_hla:
                pbase = {field: "" for field in PEPTIDE_FIELDS}
                pbase.update({
                    "peptide_id": safe_id(f"{event_id}|{allele}|{peptide}"), "event_id": event_id,
                    "sample_id": sample_id, "event_type": "Fusion", "gene": pair,
                    "peptide": peptide, "hla_allele": allele, "mhc_class": "I",
                    "source_tool": tool, "source_file": str(path), "crosses_junction": "true",
                    "contains_novel_aa": "true", "rna_junction_reads": reads,
                    "mutation_source": base["mutation_source"], "peptide_consequence": base["peptide_consequence"],
                    "binding_rank": "99", "el_rank": "99", "presentation_score": "0.0",
                    "immunogenicity_score": "0.5", "wildtype_binding_rank": "99", "self_similarity_score": "0.0",
                })
                peptides.append(enrich_peptide_layers(pbase))
        audit.append({"event_id": event_id, "gene_pair": pair, "left_breakpoint": left_bp, "right_breakpoint": right_bp, "direction": direction, "source_tool": tool, "source_file": str(path), "source_row": str(index), "peptide_status": "PROVIDED_ORF_PEPTIDE" if windows else "ORF_PEPTIDE_UNAVAILABLE_REVIEW_ONLY"})
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
    parser.add_argument("--star-fusion", type=Path); parser.add_argument("--arriba", type=Path)
    parser.add_argument("--fusioncatcher", type=Path)
    parser.add_argument("--easyfuse-all", action="append", type=Path, default=[])
    parser.add_argument("--star-chimeric", action="append", type=Path, default=[])
    parser.add_argument("--targeted-fusion-rescue", action="store_true", default=True)
    parser.add_argument("--no-targeted-fusion-rescue", action="store_false", dest="targeted_fusion_rescue")
    parser.add_argument("--targeted-fusion-whitelist", default=",".join(DEFAULT_DIAGNOSTIC_FUSION_WHITELIST))
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

    for easyfuse in existing(easyfuse_files):
        cfg = {"sample": {"id": args.sample_id, "profile": args.profile}, "inputs": {"entry_mode": "fusion", "easyfuse_tsv": str(easyfuse.resolve()), "hla_alleles": hla}}
        easyfuse_out = args.outdir / "easyfuse"
        build_raw_intermediates(cfg, easyfuse_out, root=Path.cwd())
        easyfuse_events = read_table(easyfuse_out / "parsed/raw_events.tsv")
        events.extend(easyfuse_events)
        peptides.extend(read_table(easyfuse_out / "parsed/raw_peptides.tsv"))
        for event in easyfuse_events:
            audit.append({"event_id": event.get("event_id", ""), "gene_pair": event.get("gene", ""), "left_breakpoint": "", "right_breakpoint": "", "direction": "", "source_tool": "EasyFuse", "source_file": str(easyfuse), "source_row": "", "peptide_status": "GENERATED_FROM_EASYFUSE_ORF"})
    for paths, tool in ((star_fusion_files, "STAR-Fusion"), (arriba_files, "Arriba"), (fusioncatcher_files, "FusionCatcher")):
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
    merged_events, _, _ = merge_rows_preserving_provenance(events, EVENT_FIELDS, ("event_id",), entity_type="fusion_union_event")
    merged_peptides, _, _ = merge_rows_preserving_provenance(peptides, PEPTIDE_FIELDS, ("event_id", "peptide", "hla_allele"), entity_type="fusion_union_peptide")
    write_tsv(args.outdir / "raw_events.tsv", merged_events, EVENT_FIELDS)
    write_tsv(args.outdir / "raw_peptides.tsv", merged_peptides, PEPTIDE_FIELDS)
    write_tsv(args.outdir / "parsed/raw_events.tsv", merged_events, EVENT_FIELDS)
    write_tsv(args.outdir / "parsed/raw_peptides.tsv", merged_peptides, PEPTIDE_FIELDS)
    write_tsv(args.outdir / "fusion_caller_union.tsv", audit, ["event_id", "gene_pair", "left_breakpoint", "right_breakpoint", "direction", "source_tool", "source_file", "source_row", "peptide_status"])
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
