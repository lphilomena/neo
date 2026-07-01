"""Extract variant short peptides from VEP-annotated VCF (sliding window, pVACseq-style)."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Iterator

from ..utils import open_text_maybe_gz, write_tsv

# Standard genetic code (Mutation2Peptide-compatible minigene_nt)
_CODON_TABLE = {
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "TGC": "C", "TGT": "C",
    "GAC": "D", "GAT": "D",
    "GAA": "E", "GAG": "E",
    "TTC": "F", "TTT": "F",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
    "CAC": "H", "CAT": "H",
    "ATA": "I", "ATC": "I", "ATT": "I",
    "AAA": "K", "AAG": "K",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L", "TTA": "L", "TTG": "L",
    "ATG": "M",
    "AAC": "N", "AAT": "N",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CAA": "Q", "CAG": "Q",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R", "AGA": "R", "AGG": "R",
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S", "AGC": "S", "AGT": "S",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "TGG": "W",
    "TAC": "Y", "TAT": "Y",
    "TAA": "X", "TAG": "X", "TGA": "X",
}
_REVERSE_CODON = {
    "A": "GCA", "C": "TGC", "D": "GAC", "E": "GAA", "F": "TTC", "G": "GGA",
    "H": "CAC", "I": "ATA", "K": "AAA", "L": "CTA", "M": "ATG", "N": "AAC",
    "P": "CCA", "Q": "CAA", "R": "CGA", "S": "TCA", "T": "ACA", "V": "GTA",
    "W": "TGG", "Y": "TAC", "X": "TAG",
}

MIN_PEPTIDE_LENGTH = 5
MAX_PEPTIDE_LENGTH = 50
DEFAULT_PEPTIDE_LENGTHS = (8, 9, 10, 11)
DEFAULT_MINIGENE_FLANK = 10
DEFAULT_MINIGENE_TOTAL_LEN: int | None = None

# pVACseq-compatible protein-altering consequences
PEPTIDE_CONSEQUENCES = frozenset({
    "missense_variant",
    "frameshift_variant",
    "inframe_insertion",
    "inframe_deletion",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "protein_altering_variant",
})

SPLICE_CONSEQUENCES = frozenset({
    "splice_donor_variant",
    "splice_acceptor_variant",
    "splice_donor_region_variant",
    "splice_acceptor_region_variant",
    "splice_region_variant",
    "splice_polypyrimidine_tract_variant",
})

# multi_aa_flag values written to output TSV
MULTI_AA_SINGLE = "single_aa"
MULTI_AA_SUBSTITUTION = "multi_aa_substitution"
MULTI_AA_INFRAME = "inframe_multi"
MULTI_AA_FRAMESHIFT = "frameshift"
MULTI_AA_COMPLEX = "complex"

# Categories excluded when exclude_multi_aa=True
MULTI_AA_EXCLUDED_FLAGS = frozenset({
    MULTI_AA_SUBSTITUTION,
    MULTI_AA_INFRAME,
    MULTI_AA_COMPLEX,
})

from ..adapters.peptide_netmhcpan import BINDING_ANNOTATION_FIELDS

OUTPUT_FIELDS = [
    "peptide_id",
    "gene",
    "ensembl_gene_id",
    "transcript_id",
    "hgvsc",
    "hgvsp",
    "chrom",
    "pos",
    "ref",
    "alt",
    "vaf",
    "tumor_depth",
    "tumor_alt_count",
    "rna_vaf",
    "rna_alt_reads",
    "rna_depth",
    "consequence",
    "protein_position",
    "amino_acids",
    "multi_aa_flag",
    "peptide_length",
    "peptide_start_aa",
    "peptide_end_aa",
    "mutation_position_in_peptide",
    "mutant_peptide",
    "wildtype_peptide",
    "minigene",
    "minigene_nt",
    "in_normal_proteome",
    *BINDING_ANNOTATION_FIELDS,
    "variant_key",
    "peptide_label",
    "peptide_source",
    "generation_method",
    "crosses_junction",
    "contains_novel_aa",
    "fusion_window_type",
    "fusion_breakpoint_position_raw",
    "fusion_generation_method",
]


def parse_peptide_lengths(
    lengths_str: str = "",
    *,
    length_min: int | None = None,
    length_max: int | None = None,
) -> tuple[int, ...]:
    """Parse custom peptide lengths from comma list or inclusive min–max range."""
    if length_min is not None or length_max is not None:
        if length_min is None or length_max is None:
            raise ValueError("length_min and length_max must be used together")
        if length_min > length_max:
            raise ValueError(f"length_min ({length_min}) > length_max ({length_max})")
        parsed = tuple(range(length_min, length_max + 1))
    elif lengths_str.strip():
        parsed = tuple(sorted({int(x.strip()) for x in lengths_str.split(",") if x.strip()}))
    else:
        parsed = DEFAULT_PEPTIDE_LENGTHS
    if not parsed:
        raise ValueError("No peptide lengths specified")
    for length in parsed:
        if length < MIN_PEPTIDE_LENGTH or length > MAX_PEPTIDE_LENGTH:
            raise ValueError(
                f"Peptide length {length} out of allowed range "
                f"[{MIN_PEPTIDE_LENGTH}, {MAX_PEPTIDE_LENGTH}]"
            )
    return parsed


def generation_method_for(
    lengths: tuple[int, ...],
    *,
    mini_len: int,
    filter_normal_proteome: bool,
    minigene_total_len: int | None = None,
) -> str:
    k = ",".join(str(x) for x in lengths)
    parts = [
        "sliding_window_on_mutant_protein",
        f"k={k}",
        "missense=windows_covering_changed_AA_region",
        "multi_aa=replace_ref_peptides_with_alt_by_length",
        "frameshift=all_windows_on_FrameshiftSequence",
        "exclude_identical_to_wt_window",
        f"minigene_flank={mini_len}",
        (f"minigene_total_len={minigene_total_len}" if minigene_total_len else "minigene_total_len=flank_based"),
        "snv_minigene=27mer_centered_when_total_len_set",
        "indel_minigene=peptide_centered",
    ]
    if filter_normal_proteome:
        parts.append("exclude_peptides_in_normal_proteome")
    return ";".join(parts) + ";"


GENERATION_METHOD = generation_method_for(
    DEFAULT_PEPTIDE_LENGTHS,
    mini_len=DEFAULT_MINIGENE_FLANK,
    filter_normal_proteome=False,
)


def _open_vcf(path: Path):
    return open_text_maybe_gz(path)


def _resolve_tumor_sample_index(header: list[str], tumor_sample_name: str | None) -> int | None:
    """Return tumor sample column index, or None when VCF has no FORMAT/samples."""
    if len(header) <= 9:
        return None
    if tumor_sample_name and tumor_sample_name in header:
        return header.index(tumor_sample_name)
    if len(header) > 10:
        return 10
    return 9


def _resolve_rna_sample_index(
    header: list[str],
    rna_sample_name: str | None,
    *,
    tumor_sample_name: str | None = None,
) -> int | None:
    """Return RNA-seq sample column index when present in a multi-sample VCF."""
    if len(header) <= 9:
        return None
    if rna_sample_name and rna_sample_name in header:
        return header.index(rna_sample_name)
    for col in header[9:]:
        upper = col.upper()
        if "RNA" in upper and (not tumor_sample_name or col != tumor_sample_name):
            return header.index(col)
    return None


def _parse_allele_metrics_from_format(
    fmt_str: str,
    sample_str: str,
    *,
    alt_index: int = 0,
) -> tuple[str, str, str]:
    """Extract (vaf, depth, alt_count) from VCF FORMAT/sample (AF preferred, AD fallback)."""
    fmt = fmt_str.split(":")
    vals = sample_str.split(":")
    if not fmt or not vals:
        return "", "", ""
    if len(vals) < len(fmt):
        vals = vals + [""] * (len(fmt) - len(vals))
    fmt_map = dict(zip(fmt, vals))

    vaf = ""
    depth = ""
    alt_count = ""

    af = fmt_map.get("AF", "")
    if af and af != ".":
        parts = af.split(",")
        if alt_index < len(parts):
            try:
                vaf = f"{float(parts[alt_index]):.4f}"
            except ValueError:
                pass

    ad = fmt_map.get("AD", "")
    if ad and ad != ".":
        depths = [x for x in ad.split(",") if x and x != "."]
        try:
            depth_vals = [int(float(x)) for x in depths]
        except ValueError:
            depth_vals = []
        if len(depth_vals) > alt_index + 1:
            total = sum(depth_vals)
            alt_n = depth_vals[alt_index + 1]
            alt_count = str(alt_n)
            if total > 0:
                depth = str(total)
                if not vaf:
                    vaf = f"{alt_n / total:.4f}"

    if not depth:
        dp = fmt_map.get("DP", "")
        if dp and dp != ".":
            try:
                depth = str(int(float(dp)))
            except ValueError:
                pass

    return vaf, depth, alt_count


def _parse_vaf_from_format(fmt_str: str, sample_str: str, *, alt_index: int = 0) -> str:
    """Extract tumor VAF from VCF FORMAT/sample (AF preferred, AD fallback)."""
    vaf, _, _ = _parse_allele_metrics_from_format(fmt_str, sample_str, alt_index=alt_index)
    return vaf


def parse_csq_header(header_line: str) -> dict[str, int]:
    marker = "Format: "
    if marker not in header_line:
        raise ValueError("VCF missing CSQ Format header (run vep-annotate first)")
    tail = header_line.split(marker, 1)[1].strip().rstrip(">").strip('"')
    fields = [f.strip().strip('"') for f in tail.split("|")]
    return {name: i for i, name in enumerate(fields)}


def _csq_entries(info_csq: str) -> list[list[str]]:
    return [e.split("|") for e in info_csq.split(",") if e]


def _field(parts: list[str], idx: int) -> str:
    return parts[idx] if idx < len(parts) else ""


def _consequences(cons_field: str) -> set[str]:
    return {c.strip() for c in cons_field.split("&") if c.strip()}


def _is_protein_altering(cons_field: str, *, consequence_filter: str | None = None) -> bool:
    cons = _consequences(cons_field)
    if consequence_filter == "splice":
        return bool(cons & SPLICE_CONSEQUENCES)
    return bool(cons & PEPTIDE_CONSEQUENCES)


def _csq_rank(parts: list[str], idx: dict[str, int]) -> tuple[int, int, int]:
    mane = _field(parts, idx.get("MANE_SELECT", -1))
    canonical = _field(parts, idx.get("CANONICAL", -1))
    return (
        0 if mane in {"MANE_Select", "MANE"} else 1,
        0 if canonical == "YES" else 1,
        0,
    )


def pick_csq_transcript(
    csq_entries: list[list[str]],
    idx: dict[str, int],
    *,
    consequence_filter: str | None = None,
) -> list[str] | None:
    candidates: list[tuple[tuple[int, int, int], list[str]]] = []
    for parts in csq_entries:
        cons = _field(parts, idx["Consequence"])
        if not _is_protein_altering(cons, consequence_filter=consequence_filter):
            continue
        wt = _field(parts, idx.get("WildtypeProtein", -1))
        fs = _field(parts, idx.get("FrameshiftSequence", -1))
        if not wt and not fs:
            continue
        candidates.append((_csq_rank(parts, idx), parts))
    if not candidates:
        return None
    candidates.sort(key=lambda x: x[0])
    return candidates[0][1]


def _clean_protein(seq: str) -> str:
    return seq.strip().upper().replace("*", "")


def _clean_aa_token(token: str) -> str:
    return token.strip().upper().replace("*", "")


def _cdna_translate(seq: str) -> str:
    """Encode an amino-acid string to concatenated DNA codons (minigene_nt)."""
    if not seq:
        return ""
    return "".join(_REVERSE_CODON.get(aa, _REVERSE_CODON["X"]) for aa in seq)


def _centered_minigene_segments(
    seq: str,
    start0: int,
    end0: int,
    total_len: int,
) -> tuple[str, str, str]:
    center = seq[start0:end0]
    if total_len <= 0 or len(center) >= total_len:
        return "", center[:total_len] if total_len > 0 else center, ""
    remaining = total_len - len(center)
    left_want = remaining // 2
    right_want = remaining - left_want
    left_start = max(0, start0 - left_want)
    right_end = min(len(seq), end0 + right_want)

    # Rebalance near protein edges to keep the long peptide as close to total_len as possible.
    left_short = left_want - (start0 - left_start)
    if left_short > 0:
        right_end = min(len(seq), right_end + left_short)
    right_short = right_want - (right_end - end0)
    if right_short > 0:
        left_start = max(0, left_start - right_short)

    return seq[left_start:start0], center, seq[end0:right_end]


def build_minigene(
    wt_protein: str,
    mut_protein: str,
    *,
    anchor_start: int | None,
    anchor_end: int | None,
    amino_acids: str,
    frameshift: bool = False,
    mini_len: int = DEFAULT_MINIGENE_FLANK,
    minigene_total_len: int | None = DEFAULT_MINIGENE_TOTAL_LEN,
    peptide_start_aa: int | None = None,
    peptide_end_aa: int | None = None,
    peptide_centered: bool = False,
) -> tuple[str, str]:
    """Build minigene and minigene_nt.

    SNV/missense keeps the historical Mutation2Peptide-style
    upstream|alt|downstream representation. For indel-derived peptide rows,
    callers can request peptide_centered so each row is represented as
    left_flank|mutant_peptide|right_flank around the actual short peptide.
    """
    wt = _clean_protein(wt_protein)
    mut = _clean_protein(mut_protein)
    if anchor_start is None or mini_len < 0:
        return "", ""

    if peptide_centered and peptide_start_aa is not None and peptide_end_aa is not None:
        start0 = max(0, int(peptide_start_aa) - 1)
        end0 = min(len(mut), int(peptide_end_aa))
        if start0 < end0:
            if minigene_total_len is not None:
                seq_f, center, seq_b = _centered_minigene_segments(
                    mut, start0, end0, int(minigene_total_len)
                )
            else:
                seq_f = mut[max(0, start0 - mini_len):start0]
                center = mut[start0:end0]
                seq_b = mut[end0:end0 + mini_len]
            minigene = f"{seq_f}|{center}|{seq_b}"
            minigene_nt = (
                f"{_cdna_translate(seq_f)}|{_cdna_translate(center)}|{_cdna_translate(seq_b)}"
            )
            return minigene, minigene_nt

    region_end = anchor_end if anchor_end is not None else anchor_start
    pos0 = anchor_start - 1

    if frameshift:
        seq_f = mut[max(0, pos0 - mini_len):pos0]
        novel = mut[pos0:]
        minigene = f"{seq_f}|{novel}"
        return minigene, f"{_cdna_translate(seq_f)}|{_cdna_translate(novel)}"

    alt_segment = ""
    if amino_acids and "/" in amino_acids:
        _ref_part, alt_part = amino_acids.split("/", 1)
        if alt_part not in {"", "-"}:
            alt_segment = _clean_aa_token(alt_part)

    if minigene_total_len is not None and alt_segment:
        center_end = min(len(mut), pos0 + len(alt_segment))
        if pos0 < center_end:
            seq_f, center, seq_b = _centered_minigene_segments(
                mut, pos0, center_end, int(minigene_total_len)
            )
            minigene = f"{seq_f}|{center}|{seq_b}"
            minigene_nt = (
                f"{_cdna_translate(seq_f)}|{_cdna_translate(center)}|{_cdna_translate(seq_b)}"
            )
            return minigene, minigene_nt

    seq_f = wt[max(0, pos0 - mini_len):pos0]
    tail_start = region_end
    seq_b = wt[tail_start:tail_start + mini_len]
    minigene = f"{seq_f}|{alt_segment}|{seq_b}"
    minigene_nt = (
        f"{_cdna_translate(seq_f)}|{_cdna_translate(alt_segment)}|{_cdna_translate(seq_b)}"
    )
    return minigene, minigene_nt


def _open_text_maybe_gz(path: Path):
    return open_text_maybe_gz(path)


def read_fasta_sequences(path: str | Path) -> list[str]:
    """Load protein sequences from FASTA (.fa/.fasta, optionally .gz)."""
    fasta = Path(path).resolve()
    sequences: list[str] = []
    current: list[str] = []
    with _open_text_maybe_gz(fasta) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(";"):
                continue
            if line.startswith(">"):
                if current:
                    sequences.append(_clean_protein("".join(current)))
                    current = []
                continue
            current.append(line)
    if current:
        sequences.append(_clean_protein("".join(current)))
    return [seq for seq in sequences if seq]


def build_proteome_kmer_index(
    sequences: list[str],
    lengths: tuple[int, ...],
) -> dict[int, set[str]]:
    """Index all k-mers of requested lengths from a proteome for fast lookup."""
    index: dict[int, set[str]] = {length: set() for length in lengths}
    for seq in sequences:
        for length in lengths:
            if len(seq) < length:
                continue
            for start in range(len(seq) - length + 1):
                index[length].add(seq[start:start + length])
    return index


def peptide_in_normal_proteome(peptide: str, proteome_index: dict[int, set[str]]) -> bool:
    length = len(peptide)
    bucket = proteome_index.get(length)
    return bool(bucket and peptide in bucket)


def parse_protein_position_range(raw: str) -> tuple[int, int] | None:
    """Parse VEP Protein_position (single value or inclusive range, 1-based)."""
    if not raw or raw == ".":
        return None
    if "-" in raw:
        start_s, end_s = raw.split("-", 1)
        try:
            return int(start_s), int(end_s)
        except ValueError:
            return None
    m = re.match(r"^(\d+)", raw)
    if not m:
        return None
    pos = int(m.group(1))
    return pos, pos


def classify_multi_aa_flag(
    consequence: str,
    amino_acids: str,
    protein_position_raw: str,
) -> str:
    """Classify whether one VEP record affects multiple amino acids."""
    cons = _consequences(consequence)
    if len(cons) > 1:
        return MULTI_AA_COMPLEX
    if "frameshift_variant" in cons:
        return MULTI_AA_FRAMESHIFT

    pos_range = parse_protein_position_range(protein_position_raw)
    if pos_range and pos_range[1] > pos_range[0]:
        if "missense_variant" in cons:
            return MULTI_AA_SUBSTITUTION
        if "inframe_deletion" in cons or "inframe_insertion" in cons:
            return MULTI_AA_INFRAME

    if amino_acids and "/" in amino_acids:
        ref_part, alt_part = amino_acids.split("/", 1)
        ref_clean = _clean_aa_token(ref_part)
        alt_clean = _clean_aa_token(alt_part)
        if ref_part == "-":
            ref_len = 0
        else:
            ref_len = len(ref_clean)
        alt_len = 0 if alt_part in {"", "-"} else len(alt_clean)
        if "missense_variant" in cons and (ref_len > 1 or alt_len > 1):
            return MULTI_AA_SUBSTITUTION
        if ("inframe_deletion" in cons or "inframe_insertion" in cons) and (ref_len > 1 or alt_len > 1):
            return MULTI_AA_INFRAME

    return MULTI_AA_SINGLE


def _truncate_alt_before_stop(alt_part: str) -> str:
    """Keep alt amino acids before the first stop codon token (*)."""
    raw = alt_part.strip().upper()
    if "*" in raw:
        raw = raw.split("*", 1)[0]
    return _clean_aa_token(raw)


def sliding_full_mutant_mode(consequence: str, multi_aa_flag: str) -> bool:
    """Use full-length mutant-protein sliding (no anchor overlap filter).

    Frameshift neo-proteins and stop-gained variants with in-frame insertions
    can carry novel epitopes entirely downstream of the VEP anchor interval.
    """
    if multi_aa_flag == MULTI_AA_FRAMESHIFT:
        return True
    if multi_aa_flag != MULTI_AA_COMPLEX:
        return False
    cons = _consequences(consequence)
    return "stop_gained" in cons


def _apply_aa_change(wt: str, pos_start: int, ref_part: str, alt_part: str) -> str:
    """Replace ref_part with alt_part at 1-based pos_start on wildtype protein."""
    pos0 = pos_start - 1
    if pos0 < 0:
        return ""

    ref_clean = _clean_aa_token(ref_part)
    alt_clean = (
        _truncate_alt_before_stop(alt_part)
        if "*" in alt_part
        else _clean_aa_token(alt_part)
    )
    ref_len = 0 if ref_part in {"", "-"} else len(ref_clean)
    alt_is_empty = alt_part in {"", "-"} or not alt_clean

    if ref_len == 0 and not alt_is_empty:
        return wt[:pos0] + alt_clean + wt[pos0:]
    if ref_len > 0 and alt_is_empty:
        return wt[:pos0] + wt[pos0 + ref_len :]
    if ref_len > 0:
        return wt[:pos0] + alt_clean + wt[pos0 + ref_len :]
    return wt


def build_mutant_protein(
    consequence: str,
    wt_protein: str,
    *,
    protein_position_raw: str,
    amino_acids: str,
    frameshift_sequence: str,
) -> tuple[str, int | None, int | None]:
    """Return (mutant_protein, anchor_start, anchor_end) for sliding-window overlap."""
    wt = _clean_protein(wt_protein)
    cons = _consequences(consequence)
    pos_range = parse_protein_position_range(protein_position_raw)
    if not pos_range:
        if frameshift_sequence:
            mut = _clean_protein(frameshift_sequence)
            return mut, 1, len(mut) if mut else None
        return "", None, None
    pos_start, pos_end = pos_range

    if "frameshift_variant" in cons and frameshift_sequence:
        mut = _clean_protein(frameshift_sequence)
        return mut, pos_start, pos_end

    if not wt:
        if frameshift_sequence:
            mut = _clean_protein(frameshift_sequence)
            return mut, pos_start, pos_end
        return "", None, None

    pos0 = pos_start - 1
    if pos0 < 0 or pos0 >= len(wt):
        if frameshift_sequence:
            mut = _clean_protein(frameshift_sequence)
            return mut, pos_start, pos_end
        return "", None, None

    if amino_acids and "/" in amino_acids:
        ref_part, alt_part = amino_acids.split("/", 1)
        if (
            "missense_variant" in cons
            or "inframe_deletion" in cons
            or "inframe_insertion" in cons
            or "protein_altering_variant" in cons
        ):
            mut = _apply_aa_change(wt, pos_start, ref_part, alt_part)
            if mut:
                anchor_end = pos_start + max(len(_clean_aa_token(ref_part)) if ref_part != "-" else 0, 1) - 1
                anchor_end = max(anchor_end, pos_end)
                return mut, pos_start, anchor_end

    if "stop_gained" in cons and len(cons) == 1:
        return wt[:pos0], pos_start, pos_start

    if frameshift_sequence:
        mut = _clean_protein(frameshift_sequence)
        return mut, pos_start, pos_end

    return wt, pos_start, pos_end


def sliding_variant_peptides(
    mutant_protein: str,
    wildtype_protein: str,
    *,
    anchor_start: int | None,
    anchor_end: int | None = None,
    lengths: tuple[int, ...] = (8, 9, 10, 11),
    frameshift_mode: bool = False,
) -> Iterator[dict[str, Any]]:
    """Enumerate k-mers on mutant protein; drop windows identical to wildtype."""
    mut = _clean_protein(mutant_protein)
    wt = _clean_protein(wildtype_protein)
    if not mut:
        return

    region_end = anchor_end if anchor_end is not None else anchor_start

    for length in lengths:
        if len(mut) < length:
            continue
        for start0 in range(0, len(mut) - length + 1):
            end0 = start0 + length
            mt_pep = mut[start0:end0]
            if not mt_pep or re.search(r"[^ACDEFGHIKLMNPQRSTVWY]", mt_pep):
                continue
            wt_pep = wt[start0:end0] if end0 <= len(wt) else ""
            if mt_pep == wt_pep and wt_pep:
                continue
            if not frameshift_mode and anchor_start is not None and region_end is not None:
                # Window must overlap the changed amino-acid region [anchor_start, region_end]
                if not (start0 + 1 <= region_end and end0 >= anchor_start):
                    continue
            mut_pos_in_pep = (anchor_start - start0) if anchor_start else ""
            yield {
                "peptide_length": length,
                "peptide_start_aa": start0 + 1,
                "peptide_end_aa": end0,
                "mutation_position_in_peptide": mut_pos_in_pep,
                "mutant_peptide": mt_pep,
                "wildtype_peptide": wt_pep,
            }


def _variant_key(gene: str, chrom: str, pos: str, ref: str, alt: str) -> str:
    return f"{gene}|{chrom}:{pos}{ref}>{alt}"


def _peptide_label(
    *,
    gene: str,
    transcript_id: str,
    variant_key: str,
    hgvsp: str,
    consequence: str,
    multi_aa_flag: str,
    mutant_peptide: str,
    wildtype_peptide: str,
) -> str:
    return (
        f"GENE={gene}|TX={transcript_id}|VAR={variant_key}|HGVSp={hgvsp or '.'}|"
        f"CONS={consequence}|MULTI_AA={multi_aa_flag}|"
        f"WT={wildtype_peptide or '.'}|MT={mutant_peptide}"
    )


def extract_variant_peptides_from_vcf(
    input_vcf: str | Path,
    output_tsv: str | Path,
    *,
    lengths: tuple[int, ...] = DEFAULT_PEPTIDE_LENGTHS,
    pass_only: bool = True,
    sample_id: str = "SAMPLE",
    exclude_multi_aa: bool = False,
    single_aa_only: bool = False,
    mini_len: int = DEFAULT_MINIGENE_FLANK,
    minigene_total_len: int | None = DEFAULT_MINIGENE_TOTAL_LEN,
    normal_proteome_fasta: str | Path | None = None,
    filter_normal_proteome: bool = False,
    hla_alleles: list[str] | None = None,
    netmhcpan_xls: str | Path | None = None,
    mhcflurry_csv: str | Path | None = None,
    annotate_netmhcpan: bool = False,
    tumor_sample_name: str | None = None,
    rna_sample_name: str | None = None,
    consequence_filter: str | None = None,
) -> dict[str, Any]:
    """Parse VEP-annotated VCF and write variant short peptides with metadata."""
    path = Path(input_vcf).resolve()
    out = Path(output_tsv).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)

    if mini_len < 0:
        raise ValueError(f"mini_len must be >= 0, got {mini_len}")
    if minigene_total_len is not None and int(minigene_total_len) <= 0:
        raise ValueError(f"minigene_total_len must be > 0, got {minigene_total_len}")

    proteome_index: dict[int, set[str]] | None = None
    if normal_proteome_fasta:
        proteome_seqs = read_fasta_sequences(normal_proteome_fasta)
        if not proteome_seqs:
            raise ValueError(f"No sequences loaded from normal proteome FASTA: {normal_proteome_fasta}")
        proteome_index = build_proteome_kmer_index(proteome_seqs, lengths)

    method = generation_method_for(
        lengths,
        mini_len=mini_len,
        filter_normal_proteome=filter_normal_proteome and proteome_index is not None,
        minigene_total_len=minigene_total_len,
    )

    csq_idx: dict[str, int] | None = None
    tumor_sample_idx: int | None = None
    rna_sample_idx: int | None = None
    header_cols: list[str] = []
    rows: list[dict[str, Any]] = []
    variants_seen = 0
    variants_with_peptides = 0
    variants_skipped_multi_aa = 0
    peptides_filtered_normal_proteome = 0
    peptide_counter = 0
    flag_counts: dict[str, int] = {}

    with _open_vcf(path) as fh:
        for line in fh:
            if line.startswith("##INFO=<ID=CSQ"):
                csq_idx = parse_csq_header(line)
            if line.startswith("#CHROM"):
                header_cols = line.rstrip("\n").split("\t")
                tumor_sample_idx = _resolve_tumor_sample_index(
                    header_cols,
                    tumor_sample_name or sample_id,
                )
                rna_sample_idx = _resolve_rna_sample_index(
                    header_cols,
                    rna_sample_name,
                    tumor_sample_name=tumor_sample_name,
                )
            if line.startswith("#"):
                continue
            if not csq_idx:
                raise ValueError(f"No CSQ annotations in {path}; run neoag-v03 vep-annotate first")

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, _vid, ref, alt, _qual, filt = parts[:7]
            if pass_only and filt != "PASS":
                continue

            vaf = ""
            tumor_depth = ""
            tumor_alt_count = ""
            rna_vaf = ""
            rna_alt_reads = ""
            rna_depth = ""
            if tumor_sample_idx is not None and len(parts) > tumor_sample_idx:
                vaf, tumor_depth, tumor_alt_count = _parse_allele_metrics_from_format(
                    parts[8], parts[tumor_sample_idx]
                )
            if rna_sample_idx is not None and len(parts) > rna_sample_idx:
                rna_vaf, rna_depth, rna_alt_reads = _parse_allele_metrics_from_format(
                    parts[8], parts[rna_sample_idx]
                )

            info = parts[7]
            csq_raw = ""
            for item in info.split(";"):
                if item.startswith("CSQ="):
                    csq_raw = item[4:]
                    break
            if not csq_raw:
                continue

            picked = pick_csq_transcript(
                _csq_entries(csq_raw),
                csq_idx,
                consequence_filter=consequence_filter,
            )
            if not picked:
                continue

            variants_seen += 1
            gene = _field(picked, csq_idx["SYMBOL"])
            ensembl_gene = _field(picked, csq_idx["Gene"])
            transcript = _field(picked, csq_idx["Feature"])
            hgvsc = _field(picked, csq_idx["HGVSc"])
            hgvsp = _field(picked, csq_idx["HGVSp"])
            consequence = _field(picked, csq_idx["Consequence"])
            protein_position_raw = _field(picked, csq_idx["Protein_position"])
            amino_acids = _field(picked, csq_idx["Amino_acids"])
            wt_protein = _field(picked, csq_idx.get("WildtypeProtein", -1))
            fs_protein = _field(picked, csq_idx.get("FrameshiftSequence", -1))

            multi_aa_flag = classify_multi_aa_flag(consequence, amino_acids, protein_position_raw)
            flag_counts[multi_aa_flag] = flag_counts.get(multi_aa_flag, 0) + 1

            if exclude_multi_aa and multi_aa_flag in MULTI_AA_EXCLUDED_FLAGS:
                variants_skipped_multi_aa += 1
                continue
            if single_aa_only and multi_aa_flag != MULTI_AA_SINGLE:
                variants_skipped_multi_aa += 1
                continue

            mut_protein, anchor_start, anchor_end = build_mutant_protein(
                consequence,
                wt_protein,
                protein_position_raw=protein_position_raw,
                amino_acids=amino_acids,
                frameshift_sequence=fs_protein,
            )
            if not mut_protein:
                continue

            fs_mode = sliding_full_mutant_mode(consequence, multi_aa_flag)
            variant_peps = list(
                sliding_variant_peptides(
                    mut_protein,
                    wt_protein,
                    anchor_start=anchor_start,
                    anchor_end=anchor_end,
                    lengths=lengths,
                    frameshift_mode=fs_mode,
                )
            )
            if not variant_peps:
                continue

            peptide_centered_minigene = multi_aa_flag in {
                MULTI_AA_INFRAME,
                MULTI_AA_FRAMESHIFT,
                MULTI_AA_COMPLEX,
            }

            variants_with_peptides += 1
            vkey = _variant_key(gene, chrom, pos, ref, alt)
            for pep in variant_peps:
                minigene, minigene_nt = build_minigene(
                    wt_protein,
                    mut_protein,
                    anchor_start=anchor_start,
                    anchor_end=anchor_end,
                    amino_acids=amino_acids,
                    frameshift=multi_aa_flag == MULTI_AA_FRAMESHIFT,
                    mini_len=mini_len,
                    minigene_total_len=minigene_total_len,
                    peptide_start_aa=pep.get("peptide_start_aa"),
                    peptide_end_aa=pep.get("peptide_end_aa"),
                    peptide_centered=peptide_centered_minigene,
                )
                in_normal = "no"
                if proteome_index is not None:
                    in_normal = (
                        "yes"
                        if peptide_in_normal_proteome(pep["mutant_peptide"], proteome_index)
                        else "no"
                    )
                    if filter_normal_proteome and in_normal == "yes":
                        peptides_filtered_normal_proteome += 1
                        continue

                peptide_counter += 1
                rows.append({
                    "peptide_id": f"{sample_id}.{peptide_counter}",
                    "gene": gene,
                    "ensembl_gene_id": ensembl_gene,
                    "transcript_id": transcript,
                    "hgvsc": hgvsc,
                    "hgvsp": hgvsp,
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "vaf": vaf,
                    "tumor_depth": tumor_depth,
                    "tumor_alt_count": tumor_alt_count,
                    "rna_vaf": rna_vaf,
                    "rna_alt_reads": rna_alt_reads,
                    "rna_depth": rna_depth,
                    "consequence": consequence,
                    "protein_position": protein_position_raw,
                    "amino_acids": amino_acids,
                    "multi_aa_flag": multi_aa_flag,
                    "minigene": minigene,
                    "minigene_nt": minigene_nt,
                    "in_normal_proteome": in_normal,
                    "variant_key": vkey,
                    "peptide_source": "snv",
                    "generation_method": method,
                    "peptide_label": _peptide_label(
                        gene=gene,
                        transcript_id=transcript,
                        variant_key=vkey,
                        hgvsp=hgvsp,
                        consequence=consequence,
                        multi_aa_flag=multi_aa_flag,
                        mutant_peptide=pep["mutant_peptide"],
                        wildtype_peptide=pep["wildtype_peptide"],
                    ),
                    **pep,
                })

    if hla_alleles and (annotate_netmhcpan or netmhcpan_xls or mhcflurry_csv):
        from ..adapters.peptide_netmhcpan import annotate_variant_peptide_rows

        rows = annotate_variant_peptide_rows(
            rows,
            list(hla_alleles),
            netmhcpan_xls=netmhcpan_xls,
            mhcflurry_csv=mhcflurry_csv,
            fetch_missing=annotate_netmhcpan and not netmhcpan_xls,
        )

    write_tsv(out, rows, OUTPUT_FIELDS)
    return {
        "input_vcf": str(path),
        "output_tsv": str(out),
        "variants_parsed": variants_seen,
        "variants_with_peptides": variants_with_peptides,
        "variants_skipped_multi_aa": variants_skipped_multi_aa,
        "peptide_rows": len(rows),
        "unique_mutant_peptides": len({r["mutant_peptide"] for r in rows}),
        "multi_aa_flag_counts": flag_counts,
        "lengths": ",".join(str(x) for x in lengths),
        "mini_len": mini_len,
        "minigene_total_len": minigene_total_len,
        "exclude_multi_aa": exclude_multi_aa,
        "single_aa_only": single_aa_only,
        "normal_proteome_fasta": str(normal_proteome_fasta) if normal_proteome_fasta else "",
        "filter_normal_proteome": filter_normal_proteome and proteome_index is not None,
        "peptides_filtered_normal_proteome": peptides_filtered_normal_proteome,
        "generation_method": method,
        "sample_hla_alleles": ",".join(hla_alleles) if hla_alleles else "",
        "netmhcpan_annotated": bool(hla_alleles and (annotate_netmhcpan or netmhcpan_xls)),
        "mhcflurry_annotated": bool(hla_alleles and mhcflurry_csv),
    }
