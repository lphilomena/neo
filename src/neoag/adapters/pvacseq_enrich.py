"""Post-process pVACseq aggregated output with minigene and normal-proteome annotation."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ..schemas import EVENT_FIELDS, PEPTIDE_FIELDS
from ..utils import first, read_tsv, write_tsv


from .peptide_netmhcpan import (
    BINDING_ANNOTATION_FIELDS,
    annotate_variant_peptide_row,
    build_mhcflurry_index,
    format_sample_hla_alleles,
)

ENRICHMENT_COLUMNS = (
    "minigene",
    "minigene_nt",
    "multi_aa_flag",
    "in_normal_proteome",
) + BINDING_ANNOTATION_FIELDS

# MHC I + II epitope lengths used by pVACseq runner (-e1 / -e2)
PROTEOME_INDEX_LENGTHS = tuple(range(8, 22))


@dataclass(frozen=True)
class VariantCsqRecord:
    gene: str
    transcript: str
    consequence: str
    protein_position: str
    amino_acids: str
    wt_protein: str
    fs_protein: str
    mut_protein: str
    anchor_start: int | None
    anchor_end: int | None
    multi_aa_flag: str
    minigene: str
    minigene_nt: str


def parse_pvac_variant_id(pvac_id: str) -> tuple[str, str, str, str] | None:
    """Parse pVACseq ``ID`` (chrom-start-stop-ref-alt) into chrom, pos, ref, alt."""
    raw = str(pvac_id or "").strip()
    if not raw:
        return None
    parts = raw.split("-")
    if len(parts) < 5:
        return None
    alt = parts[-1]
    ref = parts[-2]
    _stop = parts[-3]
    pos = parts[-4]
    chrom = "-".join(parts[:-4])
    if not chrom or not pos.isdigit():
        return None
    return chrom, pos, ref, alt


def normalize_chrom(chrom: str) -> str:
    c = str(chrom or "").strip()
    if c.lower().startswith("chr"):
        return c[3:]
    return c


def variant_lookup_keys(chrom: str, pos: str, ref: str, alt: str) -> list[str]:
    nc = normalize_chrom(chrom)
    prefixed = f"chr{nc}" if nc else chrom
    keys = [
        f"{nc}:{pos}:{ref}>{alt}",
        f"{prefixed}:{pos}:{ref}>{alt}",
        f"{chrom}:{pos}:{ref}>{alt}",
    ]
    seen: set[str] = set()
    ordered: list[str] = []
    for key in keys:
        if key not in seen:
            seen.add(key)
            ordered.append(key)
    return ordered


def resolve_pvacseq_enrich_options(cfg: dict[str, Any]) -> dict[str, Any]:
    inputs = cfg.get("inputs") or {}
    normal_fasta = (
        inputs.get("normal_proteome_fasta")
        or os.environ.get("NEOAG_NORMAL_PROTEOME_FASTA")
        or ""
    )
    annotate_only = bool(
        inputs.get("pvacseq_annotate_normal_proteome_only")
        or inputs.get("variant_peptide_annotate_normal_only")
    )
    enrich_minigene = inputs.get("pvacseq_enrich_minigene")
    if enrich_minigene is None:
        enrich_minigene = True

    filter_normal = inputs.get("pvacseq_filter_normal_proteome")
    if filter_normal is None:
        filter_normal = inputs.get("variant_peptide_filter_normal_proteome")
    if filter_normal is None:
        filter_normal = bool(normal_fasta) and not annotate_only

    return {
        "enrich_minigene": bool(enrich_minigene),
        "mini_len": int(
            inputs.get("pvacseq_mini_len")
            or inputs.get("variant_peptide_mini_len")
            or inputs.get("mini_len")
            or 10
        ),
        "normal_proteome_fasta": str(normal_fasta) if normal_fasta else None,
        "filter_normal_proteome": bool(filter_normal),
        "annotate_normal_proteome_only": annotate_only,
    }


def pvacseq_enrich_enabled(
    cfg: dict[str, Any],
    variants_vcf: Path | None,
    *,
    has_pvacseq_output: bool,
) -> bool:
    if not has_pvacseq_output:
        return False
    if variants_vcf is None or not variants_vcf.is_file():
        return False
    inputs = cfg.get("inputs") or {}
    if inputs.get("pvacseq_enrich") is False:
        return False
    opts = resolve_pvacseq_enrich_options(cfg)
    if inputs.get("pvacseq_enrich") is True:
        return True
    return opts["enrich_minigene"] or bool(opts["normal_proteome_fasta"])


def _required_keys_from_aggregated(rows: list[dict[str, str]]) -> set[str]:
    keys: set[str] = set()
    for row in rows:
        parsed = parse_pvac_variant_id(first(row, ["ID", "id"], ""))
        if not parsed:
            continue
        chrom, pos, ref, alt = parsed
        keys.update(variant_lookup_keys(chrom, pos, ref, alt))
    return keys


def _variant_matches_required(
    chrom: str,
    pos: str,
    ref: str,
    alt: str,
    required_keys: set[str] | None,
) -> bool:
    if not required_keys:
        return True
    return any(key in required_keys for key in variant_lookup_keys(chrom, pos, ref, alt))


def _clean_peptide_token(peptide: str) -> str:
    return str(peptide or "").strip().upper().replace("*", "")


def _indel_minigene_should_center(multi_aa_flag: str) -> bool:
    return multi_aa_flag in {"inframe_multi", "frameshift", "complex"}


def _peptide_window_on_mutant(meta: VariantCsqRecord, peptide: str) -> tuple[int | None, int | None]:
    pep = _clean_peptide_token(peptide)
    mut = _clean_peptide_token(meta.mut_protein)
    if not pep or not mut:
        return None, None
    starts: list[int] = []
    pos = mut.find(pep)
    while pos >= 0:
        starts.append(pos)
        pos = mut.find(pep, pos + 1)
    if not starts:
        return None, None

    anchor_start = meta.anchor_start
    anchor_end = meta.anchor_end if meta.anchor_end is not None else meta.anchor_start
    if anchor_start is not None and anchor_end is not None:
        for start0 in starts:
            end0 = start0 + len(pep)
            if start0 + 1 <= anchor_end and end0 >= anchor_start:
                return start0 + 1, end0
        anchor_center = (anchor_start + anchor_end) / 2
        best = min(starts, key=lambda x: abs((x + 1 + x + len(pep)) / 2 - anchor_center))
        return best + 1, best + len(pep)

    return starts[0] + 1, starts[0] + len(pep)


def _minigene_for_pvac_peptide(meta: VariantCsqRecord, peptide: str, mini_len: int) -> tuple[str, str]:
    if not _indel_minigene_should_center(meta.multi_aa_flag):
        return meta.minigene, meta.minigene_nt
    from ..vep.extract_peptides import build_minigene

    pep_start, pep_end = _peptide_window_on_mutant(meta, peptide)
    if pep_start is None or pep_end is None:
        return meta.minigene, meta.minigene_nt
    return build_minigene(
        meta.wt_protein,
        meta.mut_protein,
        anchor_start=meta.anchor_start,
        anchor_end=meta.anchor_end,
        amino_acids=meta.amino_acids,
        frameshift=meta.multi_aa_flag == "frameshift",
        mini_len=mini_len,
        peptide_start_aa=pep_start,
        peptide_end_aa=pep_end,
        peptide_centered=True,
    )


def _build_vcf_variant_index(
    vcf_path: Path,
    *,
    mini_len: int,
    pass_only: bool = True,
    required_keys: set[str] | None = None,
) -> dict[str, VariantCsqRecord]:
    from ..vep.extract_peptides import (
        _csq_entries,
        _field,
        _open_vcf,
        _variant_key,
        build_minigene,
        build_mutant_protein,
        classify_multi_aa_flag,
        parse_csq_header,
        pick_csq_transcript,
    )

    index: dict[str, VariantCsqRecord] = {}
    csq_idx: dict[str, int] | None = None

    with _open_vcf(vcf_path) as fh:
        for line in fh:
            if line.startswith("##INFO=<ID=CSQ"):
                csq_idx = parse_csq_header(line)
            if line.startswith("#"):
                continue
            if not csq_idx:
                raise ValueError(f"No CSQ annotations in {vcf_path}; run vep-annotate first")

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, _vid, ref, alt, _qual, filt = parts[:7]
            if pass_only and filt != "PASS":
                continue
            if not _variant_matches_required(chrom, pos, ref, alt, required_keys):
                continue

            info = parts[7]
            csq_raw = ""
            for item in info.split(";"):
                if item.startswith("CSQ="):
                    csq_raw = item[4:]
                    break
            if not csq_raw:
                continue

            picked = pick_csq_transcript(_csq_entries(csq_raw), csq_idx)
            if not picked:
                continue

            gene = _field(picked, csq_idx["SYMBOL"])
            transcript = _field(picked, csq_idx["Feature"])
            consequence = _field(picked, csq_idx["Consequence"])
            protein_position = _field(picked, csq_idx["Protein_position"])
            amino_acids = _field(picked, csq_idx["Amino_acids"])
            wt_protein = _field(picked, csq_idx.get("WildtypeProtein", -1))
            fs_protein = _field(picked, csq_idx.get("FrameshiftSequence", -1))
            multi_aa_flag = classify_multi_aa_flag(consequence, amino_acids, protein_position)

            mut_protein, anchor_start, anchor_end = build_mutant_protein(
                consequence,
                wt_protein,
                protein_position_raw=protein_position,
                amino_acids=amino_acids,
                frameshift_sequence=fs_protein,
            )
            minigene, minigene_nt = "", ""
            if mut_protein:
                minigene, minigene_nt = build_minigene(
                    wt_protein,
                    mut_protein,
                    anchor_start=anchor_start,
                    anchor_end=anchor_end,
                    amino_acids=amino_acids,
                    frameshift=multi_aa_flag == "frameshift",
                    mini_len=mini_len,
                )

            record = VariantCsqRecord(
                gene=gene,
                transcript=transcript,
                consequence=consequence,
                protein_position=protein_position,
                amino_acids=amino_acids,
                wt_protein=wt_protein,
                fs_protein=fs_protein,
                mut_protein=mut_protein,
                anchor_start=anchor_start,
                anchor_end=anchor_end,
                multi_aa_flag=multi_aa_flag,
                minigene=minigene,
                minigene_nt=minigene_nt,
            )
            for key in variant_lookup_keys(chrom, pos, ref, alt):
                index[key] = record
            gene_key = _variant_key(gene, chrom, pos, ref, alt)
            index[gene_key] = record
    return index


def _lookup_variant(
    row: dict[str, str],
    vcf_index: dict[str, VariantCsqRecord],
) -> VariantCsqRecord | None:
    pvac_id = first(row, ["ID", "id"], "")
    parsed = parse_pvac_variant_id(pvac_id)
    if parsed:
        chrom, pos, ref, alt = parsed
        for key in variant_lookup_keys(chrom, pos, ref, alt):
            hit = vcf_index.get(key)
            if hit:
                return hit
    return None


def _peptide_lengths_from_rows(rows: list[dict[str, str]]) -> tuple[int, ...]:
    lengths: set[int] = set()
    for row in rows:
        peptide = first(
            row,
            ["Best Peptide", "MT Epitope Seq", "MT Epitope", "Peptide", "peptide"],
            "",
        )
        if peptide:
            lengths.add(len(peptide))
    return tuple(sorted(lengths)) or PROTEOME_INDEX_LENGTHS


def _peptide_lengths_from_raw(rows: list[dict[str, str]]) -> tuple[int, ...]:
    lengths = {len(row.get("peptide") or "") for row in rows}
    lengths.discard(0)
    return tuple(sorted(lengths)) or PROTEOME_INDEX_LENGTHS


def _peptide_in_normal(
    peptide: str,
    proteome_index: dict[int, set[str]] | None,
) -> str:
    if not peptide or proteome_index is None:
        return "no"
    from ..vep.extract_peptides import peptide_in_normal_proteome

    return "yes" if peptide_in_normal_proteome(peptide, proteome_index) else "no"


def _load_pvac_all_epitopes_rows(aggregated_path: Path) -> list[dict[str, str]]:
    work = aggregated_path.parent
    if work.name == "tools":
        work = work / "pvacseq"
    elif (work / "pvacseq").is_dir():
        work = work / "pvacseq"
    candidates = sorted(work.rglob("*all_epitopes.tsv")) if work.is_dir() else []
    return read_tsv(candidates[0]) if candidates else []


def _wt_epitope_for_aggregated_row(
    row: dict[str, str],
    all_epitopes_rows: list[dict[str, str]],
) -> str:
    best = first(row, ["Best Peptide", "MT Epitope Seq"], "")
    allele = first(row, ["Allele", "HLA Allele"], "")
    for ep in all_epitopes_rows:
        if (
            first(ep, ["MT Epitope Seq", "MT Epitope"], "") == best
            and first(ep, ["HLA Allele", "Allele"], "") == allele
        ):
            return first(ep, ["WT Epitope Seq", "WT Epitope"], "")
    return ""


def enrich_pvacseq_aggregated(
    aggregated_tsv: str | Path,
    variants_vcf: str | Path,
    out_enriched_tsv: str | Path,
    *,
    enrich_minigene: bool = True,
    mini_len: int = 10,
    normal_proteome_fasta: str | Path | None = None,
    hla_alleles: list[str] | None = None,
    netmhcpan_xls: str | Path | None = None,
    mhcflurry_csv: str | Path | None = None,
) -> dict[str, Any]:
    """Annotate pVACseq aggregated rows with minigene / normal-proteome columns."""
    from ..vep.extract_peptides import build_proteome_kmer_index, read_fasta_sequences

    aggregated_path = Path(aggregated_tsv)
    vcf_path = Path(variants_vcf)
    out_path = Path(out_enriched_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows_in = read_tsv(aggregated_path)
    all_epitopes_rows = _load_pvac_all_epitopes_rows(aggregated_path)

    netmhcpan_index = None
    if netmhcpan_xls:
        from .peptide_netmhcpan import build_netmhcpan_index

        netmhcpan_index = build_netmhcpan_index(netmhcpan_xls)

    mhcflurry_index = build_mhcflurry_index(mhcflurry_csv) if mhcflurry_csv else None

    proteome_index: dict[int, set[str]] | None = None
    if normal_proteome_fasta:
        proteome_seqs = read_fasta_sequences(normal_proteome_fasta)
        if proteome_seqs:
            pep_lengths = _peptide_lengths_from_rows(rows_in)
            proteome_index = build_proteome_kmer_index(proteome_seqs, pep_lengths)
    required_keys = _required_keys_from_aggregated(rows_in) if enrich_minigene else None
    rows_out: list[dict[str, str]] = []
    matched_variants = 0
    in_normal_yes = 0

    vcf_index: dict[str, VariantCsqRecord] = {}
    if enrich_minigene:
        vcf_index = _build_vcf_variant_index(
            vcf_path,
            mini_len=mini_len,
            required_keys=required_keys,
        )

    for row in rows_in:
        out = dict(row)
        wt_ep = _wt_epitope_for_aggregated_row(row, all_epitopes_rows)
        if wt_ep:
            out["WT Epitope Seq"] = wt_ep
        peptide = first(
            row,
            ["Best Peptide", "MT Epitope Seq", "MT Epitope", "Peptide", "peptide"],
            "",
        )

        meta = _lookup_variant(row, vcf_index) if vcf_index else None
        if meta:
            matched_variants += 1
            out["multi_aa_flag"] = meta.multi_aa_flag
            if enrich_minigene:
                minigene, minigene_nt = _minigene_for_pvac_peptide(meta, peptide, mini_len)
                out["minigene"] = minigene
                out["minigene_nt"] = minigene_nt
        else:
            out.setdefault("multi_aa_flag", "")
            out.setdefault("minigene", "")
            out.setdefault("minigene_nt", "")

        in_normal = _peptide_in_normal(peptide, proteome_index)
        out["in_normal_proteome"] = in_normal
        if in_normal == "yes":
            in_normal_yes += 1

        if hla_alleles:
            out = annotate_variant_peptide_row(
                out,
                hla_alleles,
                netmhcpan_index,
                mhcflurry_index=mhcflurry_index,
                fetch_missing=False,
                prefer_local_netmhcpan=bool(netmhcpan_index),
                mt_key="Best Peptide",
                wt_key="WT Epitope Seq",
                peptide_key="Best Peptide",
            )
        rows_out.append(out)

    fieldnames = list(rows_in[0].keys()) if rows_in else []
    if all_epitopes_rows and "WT Epitope Seq" not in fieldnames:
        fieldnames.append("WT Epitope Seq")
    for col in ENRICHMENT_COLUMNS:
        if col not in fieldnames:
            fieldnames.append(col)
    write_tsv(out_path, rows_out, fieldnames)

    return {
        "input_aggregated": str(aggregated_path),
        "output_enriched": str(out_path),
        "rows_in": len(rows_in),
        "rows_out": len(rows_out),
        "variants_matched": matched_variants,
        "in_normal_proteome_yes": in_normal_yes,
        "enrich_minigene": enrich_minigene,
        "normal_proteome_fasta": str(normal_proteome_fasta) if normal_proteome_fasta else "",
    }


def filter_raw_peptides_normal_proteome(
    raw_peptides_tsv: str | Path,
    *,
    normal_proteome_fasta: str | Path,
    annotate_only: bool = False,
) -> dict[str, Any]:
    """Drop or annotate raw_peptides rows whose peptide appears in the reference proteome."""
    from ..vep.extract_peptides import (
        build_proteome_kmer_index,
        peptide_in_normal_proteome,
        read_fasta_sequences,
    )

    path = Path(raw_peptides_tsv)
    rows = read_tsv(path)
    proteome_seqs = read_fasta_sequences(normal_proteome_fasta)
    if not proteome_seqs:
        raise ValueError(f"No sequences loaded from normal proteome FASTA: {normal_proteome_fasta}")
    pep_lengths = _peptide_lengths_from_raw(rows)
    proteome_index = build_proteome_kmer_index(proteome_seqs, pep_lengths)
    kept: list[dict[str, str]] = []
    filtered = 0
    for row in rows:
        peptide = row.get("peptide") or ""
        if peptide and peptide_in_normal_proteome(peptide, proteome_index):
            filtered += 1
            if annotate_only:
                kept.append(row)
            continue
        kept.append(row)

    if not annotate_only:
        write_tsv(path, kept, PEPTIDE_FIELDS)
    return {
        "raw_peptides": str(path),
        "rows_in": len(rows),
        "rows_kept": len(kept) if not annotate_only else len(rows),
        "peptides_filtered_normal_proteome": filtered if not annotate_only else 0,
        "annotate_only": annotate_only,
    }


def refresh_raw_peptides_from_enriched(
    enriched_tsv: str | Path,
    *,
    sample_id: str,
    profile_name: str,
    raw_peptides_tsv: str | Path,
    raw_events_tsv: str | Path | None = None,
) -> int:
    """Re-parse pVACseq enriched rows so raw_peptides carries WT + NetMHCpan columns."""
    from .pvactools_parser import parse_pvactools_outputs

    _, peptides = parse_pvactools_outputs(
        [str(enriched_tsv)],
        sample_id,
        profile_name,
        out_events=raw_events_tsv,
        out_peptides=raw_peptides_tsv,
    )
    return len(peptides)


def enrich_pvacseq_outputs(
    cfg: dict[str, Any],
    *,
    aggregated_tsv: str | Path,
    variants_vcf: str | Path,
    raw_peptides_tsv: str | Path | None = None,
    raw_events_tsv: str | Path | None = None,
    out_enriched_tsv: str | Path,
    sample_id: str | None = None,
    profile_name: str | None = None,
) -> dict[str, str]:
    """Run pVACseq post-enrichment and optional raw_peptides filtering."""
    opts = resolve_pvacseq_enrich_options(cfg)
    inputs = cfg.get("inputs") or {}
    hla_alleles = list(inputs.get("hla_alleles") or [])
    netmhcpan_xls = inputs.get("netmhcpan")
    mhcflurry_csv = inputs.get("mhcflurry")
    summary = enrich_pvacseq_aggregated(
        aggregated_tsv,
        variants_vcf,
        out_enriched_tsv,
        enrich_minigene=opts["enrich_minigene"],
        mini_len=opts["mini_len"],
        normal_proteome_fasta=opts["normal_proteome_fasta"],
        hla_alleles=hla_alleles,
        netmhcpan_xls=netmhcpan_xls,
        mhcflurry_csv=mhcflurry_csv,
    )

    outputs: dict[str, str] = {
        "pvacseq_enriched": summary["output_enriched"],
        "pvacseq_enrich_rows": str(summary["rows_out"]),
        "pvacseq_enrich_variants_matched": str(summary["variants_matched"]),
        "pvacseq_enrich_in_normal_yes": str(summary["in_normal_proteome_yes"]),
    }

    sid = sample_id or (cfg.get("sample") or {}).get("id") or "SAMPLE001"
    prof = profile_name or (cfg.get("sample") or {}).get("profile") or "default"
    if raw_peptides_tsv:
        n = refresh_raw_peptides_from_enriched(
            summary["output_enriched"],
            sample_id=sid,
            profile_name=prof,
            raw_peptides_tsv=raw_peptides_tsv,
            raw_events_tsv=raw_events_tsv,
        )
        outputs["pvacseq_raw_peptides_refreshed"] = str(n)

    if (
        raw_peptides_tsv
        and opts["normal_proteome_fasta"]
        and opts["filter_normal_proteome"]
        and not opts["annotate_normal_proteome_only"]
    ):
        filt = filter_raw_peptides_normal_proteome(
            raw_peptides_tsv,
            normal_proteome_fasta=opts["normal_proteome_fasta"],
            annotate_only=False,
        )
        outputs["pvacseq_peptides_filtered_normal"] = str(filt["peptides_filtered_normal_proteome"])
        outputs["pvacseq_raw_peptides_rows"] = str(filt["rows_kept"])

    # raw_events unchanged; reserved for future event-level annotations
    _ = raw_events_tsv
    return outputs
