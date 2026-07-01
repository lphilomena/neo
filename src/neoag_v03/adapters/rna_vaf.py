"""Parse RNA allele count / RNA VAF support tables."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

from ..utils import first, read_tsv, to_float


@dataclass(frozen=True)
class RNAVafSupport:
    event_id: str = ""
    gene: str = ""
    chrom: str = ""
    pos: str = ""
    ref: str = ""
    alt: str = ""
    rna_alt_reads: str = ""
    rna_ref_reads: str = ""
    rna_depth: str = ""
    rna_vaf: str = ""
    source: str = ""

    def as_fields(self) -> dict[str, str]:
        return {
            "rna_alt_reads": self.rna_alt_reads,
            "rna_ref_reads": self.rna_ref_reads,
            "rna_depth": self.rna_depth,
            "rna_vaf": self.rna_vaf,
            "rna_vaf_source": self.source,
        }


def _fmt_int(value: str) -> str:
    if value == "":
        return ""
    return str(int(to_float(value, 0.0)))


def _fmt_vaf(value: str, alt: str = "", depth: str = "") -> str:
    if value != "":
        v = to_float(value, 0.0)
        if v > 1.0 and v <= 100.0:
            v = v / 100.0
        return f"{v:.4f}"
    a = to_float(alt, -1.0)
    d = to_float(depth, 0.0)
    if a >= 0 and d > 0:
        return f"{a / d:.4f}"
    return ""


def _variant_key(row: Mapping[str, str]) -> str:
    chrom = first(row, ["chrom", "chr", "chromosome", "contig"], "")
    pos = first(row, ["pos", "position", "start"], "")
    ref = first(row, ["ref", "reference", "reference_allele"], "")
    alt = first(row, ["alt", "alternate", "alt_allele", "variant_allele"], "")
    if chrom and pos and ref and alt:
        return f"{chrom}:{pos}:{ref}>{alt}"
    return ""


def parse_rna_vaf_table(path: str | Path) -> list[RNAVafSupport]:
    """Parse a generic RNA allele-count TSV.

    Accepted columns include event_id/gene/chrom/pos/ref/alt plus common aliases
    for RNA alt/ref/depth/VAF from bcftools mpileup, pysam counters, or custom
    variant support tables.
    """
    source = str(path)
    out: list[RNAVafSupport] = []
    for row in read_tsv(path):
        event_id = first(row, ["event_id", "variant_id", "id", "mutation_id"], "")
        gene = first(row, ["gene", "Gene", "symbol"], "")
        chrom = first(row, ["chrom", "chr", "chromosome", "contig"], "")
        pos = first(row, ["pos", "position", "start"], "")
        ref = first(row, ["ref", "reference", "reference_allele"], "")
        alt = first(row, ["alt", "alternate", "alt_allele", "variant_allele"], "")
        alt_reads = first(row, ["rna_alt_reads", "rna_alt_count", "RNA Alt Count", "alt_reads", "alt_count", "AD_ALT", "t_alt_count"], "")
        ref_reads = first(row, ["rna_ref_reads", "rna_ref_count", "RNA Ref Count", "ref_reads", "ref_count", "AD_REF", "t_ref_count"], "")
        depth = first(row, ["rna_depth", "RNA Depth", "depth", "dp", "DP", "coverage", "total_reads"], "")
        if depth == "" and (alt_reads or ref_reads):
            depth = str(int(to_float(alt_reads, 0.0) + to_float(ref_reads, 0.0)))
        vaf = _fmt_vaf(first(row, ["rna_vaf", "RNA VAF", "vaf", "VAF", "af", "AF", "allele_fraction"], ""), alt_reads, depth)
        out.append(RNAVafSupport(
            event_id=event_id,
            gene=gene,
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            rna_alt_reads=_fmt_int(alt_reads),
            rna_ref_reads=_fmt_int(ref_reads),
            rna_depth=_fmt_int(depth),
            rna_vaf=vaf,
            source=source,
        ))
    return out


def load_rna_vaf_support(path: str | Path | None) -> dict[str, RNAVafSupport]:
    """Return support indexed by event id, gene, and chrom:pos:ref>alt keys."""
    if not path or not Path(path).is_file():
        return {}
    index: dict[str, RNAVafSupport] = {}
    for support in parse_rna_vaf_table(path):
        keys = [support.event_id, support.gene]
        if support.chrom and support.pos and support.ref and support.alt:
            keys.append(f"{support.chrom}:{support.pos}:{support.ref}>{support.alt}")
        for key in keys:
            if key and key not in index:
                index[key] = support
    return index


def support_from_raw_row(row: Mapping[str, str], source: str) -> RNAVafSupport:
    alt_reads = first(row, ["rna_alt_reads", "rna_alt_count", "RNA Alt Count"], "")
    ref_reads = first(row, ["rna_ref_reads", "rna_ref_count", "RNA Ref Count"], "")
    depth = first(row, ["rna_depth", "RNA Depth"], "")
    vaf = _fmt_vaf(first(row, ["rna_vaf", "RNA VAF"], ""), alt_reads, depth)
    return RNAVafSupport(
        event_id=first(row, ["event_id"], ""),
        gene=first(row, ["gene"], ""),
        chrom=first(row, ["chrom"], ""),
        pos=first(row, ["pos"], ""),
        ref=first(row, ["ref"], ""),
        alt=first(row, ["alt"], ""),
        rna_alt_reads=_fmt_int(alt_reads),
        rna_ref_reads=_fmt_int(ref_reads),
        rna_depth=_fmt_int(depth),
        rna_vaf=vaf,
        source=source if any([alt_reads, ref_reads, depth, vaf]) else "",
    )


def choose_rna_vaf_support(row: Mapping[str, str], external: Mapping[str, RNAVafSupport], raw_source: str) -> RNAVafSupport:
    keys = [
        first(row, ["event_id"], ""),
        first(row, ["gene"], ""),
        _variant_key(row),
    ]
    for key in keys:
        if key and key in external:
            return external[key]
    return support_from_raw_row(row, raw_source)
