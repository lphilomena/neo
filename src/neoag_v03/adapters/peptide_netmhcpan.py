"""Annotate variant-peptide tables with HLA typing and NetMHCpan / MHCflurry MT/WT binding."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ..utils import first, read_tsv, to_float, write_tsv

NETMHCPAN_ANNOTATION_FIELDS = (
    "sample_hla_alleles",
    "hla_allele",
    "netmhcpan_mt_ic50",
    "netmhcpan_mt_rank_ba",
    "netmhcpan_mt_rank_el",
    "netmhcpan_wt_ic50",
    "netmhcpan_wt_rank_ba",
    "netmhcpan_wt_rank_el",
)

MHCFLURRY_ANNOTATION_FIELDS = (
    "mhcflurry_mt_affinity",
    "mhcflurry_mt_affinity_percentile",
    "mhcflurry_mt_processing_score",
    "mhcflurry_mt_presentation_score",
    "mhcflurry_wt_affinity",
    "mhcflurry_wt_affinity_percentile",
    "mhcflurry_wt_processing_score",
    "mhcflurry_wt_presentation_score",
)

NETMHCSTABPAN_ANNOTATION_FIELDS = (
    "netmhcstabpan_score",
    "netmhcstabpan_rank",
    "netmhcstabpan_wt_score",
    "netmhcstabpan_wt_rank",
)

PRIME_ANNOTATION_FIELDS = (
    "prime_score",
    "prime_rank",
    "prime_wt_score",
    "prime_wt_rank",
)

BIGMHC_IM_ANNOTATION_FIELDS = (
    "bigmhc_im_score",
    "bigmhc_im_wt_score",
)

IEDB_ANNOTATION_FIELDS = (
    "iedb_immunogenicity_score",
    "iedb_immunogenicity_wt_score",
)

BINDING_ANNOTATION_FIELDS = NETMHCPAN_ANNOTATION_FIELDS + MHCFLURRY_ANNOTATION_FIELDS

TOOL_ANNOTATION_FIELDS = (
    BINDING_ANNOTATION_FIELDS
    + NETMHCSTABPAN_ANNOTATION_FIELDS
    + PRIME_ANNOTATION_FIELDS
    + BIGMHC_IM_ANNOTATION_FIELDS
    + IEDB_ANNOTATION_FIELDS
)


def format_sample_hla_alleles(hla_alleles: list[str] | None) -> str:
    return ",".join(str(a).strip() for a in (hla_alleles or []) if str(a).strip())


def resolve_peptide_source(row: dict[str, str]) -> str:
    """Return catalog peptide origin: snv (VCF) or easyfuse (fusion neo)."""
    explicit = str(row.get("peptide_source") or "").strip().lower()
    if explicit in {"snv", "easyfuse"}:
        return explicit
    if str(row.get("multi_aa_flag") or "").lower() == "fusion_neo":
        return "easyfuse"
    if "easyfuse_neo_peptide_sliding_window" in str(row.get("generation_method") or ""):
        return "easyfuse"
    return "snv"


def _normalize_hla(allele: str) -> str:
    from .netmhcpan import _hla_with_star

    return _hla_with_star(allele)


def build_mhcflurry_index(csv_path: str | Path) -> dict[tuple[str, str], dict[str, str]]:
    """Index MHCflurry CSV rows by (peptide, normalized HLA)."""
    from .mhcflurry import parse_mhcflurry

    index: dict[tuple[str, str], dict[str, str]] = {}
    for row in parse_mhcflurry(csv_path):
        peptide = row.get("peptide") or ""
        allele = _normalize_hla(row.get("hla_allele") or "")
        if peptide and allele:
            index[(peptide, allele)] = row
    return index


def build_netmhcstabpan_index(tsv_path: str | Path) -> dict[tuple[str, str], dict[str, str]]:
    """Index NetMHCstabpan rows by (peptide, normalized HLA)."""
    from .netmhcstabpan import parse_netmhcstabpan

    index: dict[tuple[str, str], dict[str, str]] = {}
    for row in parse_netmhcstabpan(tsv_path):
        peptide = (row.get("peptide") or "").strip().upper()
        allele = _normalize_hla(row.get("hla_allele") or "")
        if peptide and allele:
            index[(peptide, allele)] = row
    return index


def build_prime_indexes(
    path: str | Path,
) -> tuple[dict[str, dict[str, str]], dict[tuple[str, str], dict[str, str]]]:
    """Return (wide-by-peptide, long-by-pair) PRIME indexes."""
    from .prime import parse_prime, read_prime_wide_rows

    wide: dict[str, dict[str, str]] = {}
    for row in read_prime_wide_rows(path):
        peptide = (row.get("Peptide") or row.get("peptide") or "").strip().upper()
        if peptide:
            wide[peptide] = row
    pair: dict[tuple[str, str], dict[str, str]] = {}
    for row in parse_prime(path):
        peptide = (row.get("peptide") or "").strip().upper()
        allele = _normalize_hla(row.get("hla_allele") or "")
        if peptide and allele:
            pair[(peptide, allele)] = row
    return wide, pair


def build_bigmhc_im_index(path: str | Path) -> dict[tuple[str, str], dict[str, str]]:
    """Index BigMHC_IM rows by (peptide, normalized HLA)."""
    from .bigmhc_im import bigmhc_by_pair, parse_bigmhc_im

    return bigmhc_by_pair(parse_bigmhc_im(path))


def build_iedb_immunogenicity_index(path: str | Path) -> dict[tuple[str, str], dict[str, str]]:
    """Index IEDB immunogenicity evidence by (peptide, normalized HLA)."""
    index: dict[tuple[str, str], dict[str, str]] = {}
    for row in read_tsv(path):
        peptide = (row.get("peptide") or row.get("mutant_peptide") or "").strip().upper()
        allele = _normalize_hla(row.get("hla_allele") or "")
        if peptide and allele:
            index[(peptide, allele)] = row
    return index


def _mhcflurry_columns_from_pred(
    pred: dict[str, str] | None,
    *,
    prefix: str,
) -> dict[str, str]:
    if not pred:
        return {
            f"mhcflurry_{prefix}_affinity": "",
            f"mhcflurry_{prefix}_affinity_percentile": "",
            f"mhcflurry_{prefix}_processing_score": "",
            f"mhcflurry_{prefix}_presentation_score": "",
        }
    return {
        f"mhcflurry_{prefix}_affinity": pred.get("mhcflurry_affinity", ""),
        f"mhcflurry_{prefix}_affinity_percentile": pred.get("mhcflurry_affinity_percentile", ""),
        f"mhcflurry_{prefix}_processing_score": pred.get("mhcflurry_processing_score", ""),
        f"mhcflurry_{prefix}_presentation_score": pred.get("mhcflurry_presentation_score", ""),
    }


def _lookup_mhcflurry(
    peptide: str,
    allele: str,
    index: dict[tuple[str, str], dict[str, str]] | None,
) -> dict[str, str] | None:
    if not peptide or not allele or not index:
        return None
    return index.get((peptide, _normalize_hla(allele)))


def pick_best_hla_mhcflurry_for_peptide(
    peptide: str,
    hla_alleles: list[str],
    index: dict[tuple[str, str], dict[str, str]] | None,
) -> tuple[str, dict[str, str]] | tuple[None, None]:
    best_hla = ""
    best_pred: dict[str, str] | None = None
    best_pct = 99.0
    for allele in hla_alleles:
        pred = _lookup_mhcflurry(peptide, allele, index)
        if not pred:
            continue
        pct = to_float(pred.get("mhcflurry_affinity_percentile"), 99.0)
        if pct < best_pct:
            best_pct = pct
            best_hla = _normalize_hla(allele)
            best_pred = pred
    if best_hla and best_pred:
        return best_hla, best_pred
    return None, None


def apply_mhcflurry_annotation(
    out: dict[str, str],
    mt_peptide: str,
    wt_peptide: str,
    hla_alleles: list[str],
    mhcflurry_index: dict[tuple[str, str], dict[str, str]] | None,
    *,
    bind_hla: str = "",
) -> dict[str, str]:
    if not mhcflurry_index:
        for col in MHCFLURRY_ANNOTATION_FIELDS:
            out.setdefault(col, "")
        return out

    allele = bind_hla or out.get("hla_allele") or ""
    if not allele and mt_peptide and hla_alleles:
        picked_hla, _ = pick_best_hla_mhcflurry_for_peptide(mt_peptide, hla_alleles, mhcflurry_index)
        allele = picked_hla or ""

    mt_pred = _lookup_mhcflurry(mt_peptide, allele, mhcflurry_index) if allele else None
    out.update(_mhcflurry_columns_from_pred(mt_pred, prefix="mt"))

    if wt_peptide and allele:
        wt_pred = _lookup_mhcflurry(wt_peptide, allele, mhcflurry_index)
        out.update(_mhcflurry_columns_from_pred(wt_pred, prefix="wt"))
    else:
        out.update(_mhcflurry_columns_from_pred(None, prefix="wt"))

    for col in MHCFLURRY_ANNOTATION_FIELDS:
        out.setdefault(col, "")
    return out


def build_netmhcpan_index(xls_path: str | Path) -> dict[tuple[str, str], dict[str, str]]:
    """Index NetMHCpan XLS/TSV rows by (peptide, normalized HLA)."""
    from .netmhcpan import _hla_with_star, parse_netmhcpan

    index: dict[tuple[str, str], dict[str, str]] = {}
    for row in parse_netmhcpan(xls_path):
        peptide = row.get("peptide") or ""
        allele = _hla_with_star(row.get("hla_allele") or "")
        if peptide and allele:
            index[(peptide, allele)] = row
    return index


def _prediction_from_index(
    peptide: str,
    allele: str,
    index: dict[tuple[str, str], dict[str, str]],
) -> dict[str, str] | None:
    if not peptide or not allele:
        return None
    return index.get((peptide, _normalize_hla(allele)))


def _fetch_iedb_prediction(peptide: str, allele: str) -> dict[str, str]:
    from ..tools.runner import _iedb_mhci_row

    ba = _iedb_mhci_row("netmhcpan_ba", peptide, _normalize_hla(allele))
    el = _iedb_mhci_row("netmhcpan_el", peptide, _normalize_hla(allele))
    return {
        "netmhcpan_ba_score": str(to_float(ba.get("ic50") or ba.get("Score_BA"), 0.0)),
        "netmhcpan_ba_rank": str(to_float(ba.get("percentile_rank") or ba.get("%Rank_BA"), 99.0)),
        "netmhcpan_el_score": str(to_float(el.get("score") or el.get("Score_EL"), 0.0)),
        "netmhcpan_el_rank": str(to_float(el.get("percentile_rank") or el.get("%Rank_EL"), 99.0)),
    }


def _lookup_prediction(
    peptide: str,
    allele: str,
    index: dict[tuple[str, str], dict[str, str]] | None,
    *,
    fetch_missing: bool,
    cache: dict[tuple[str, str], dict[str, str]],
) -> dict[str, str] | None:
    key = (peptide, _normalize_hla(allele))
    if key in cache:
        return cache[key]
    pred = _prediction_from_index(peptide, allele, index or {})
    if pred is None and fetch_missing:
        pred = _fetch_iedb_prediction(peptide, allele)
    if pred is not None:
        cache[key] = pred
    return pred


def pick_best_hla_for_peptide(
    peptide: str,
    hla_alleles: list[str],
    index: dict[tuple[str, str], dict[str, str]] | None,
    *,
    fetch_missing: bool = False,
    cache: dict[tuple[str, str], dict[str, str]] | None = None,
) -> tuple[str, dict[str, str]] | tuple[None, None]:
    best_hla = ""
    best_pred: dict[str, str] | None = None
    best_rank = 99.0
    local_cache = cache if cache is not None else {}
    for allele in hla_alleles:
        pred = _lookup_prediction(
            peptide,
            allele,
            index,
            fetch_missing=fetch_missing,
            cache=local_cache,
        )
        if not pred:
            continue
        rank = to_float(pred.get("netmhcpan_ba_rank"), 99.0)
        if rank < best_rank:
            best_rank = rank
            best_hla = _normalize_hla(allele)
            best_pred = pred
    if best_hla and best_pred:
        return best_hla, best_pred
    return None, None


def netmhcpan_columns_from_pvac_row(row: dict[str, str]) -> dict[str, str]:
    """Map pVACseq aggregated / all-epitopes NetMHCpan fields to standard names."""
    return {
        "hla_allele": first(row, ["Allele", "HLA Allele", "hla_allele"], ""),
        "netmhcpan_mt_ic50": first(
            row,
            ["IC50 MT", "NetMHCpan MT IC50 Score", "Best MT IC50 Score", "netmhcpan_mt_ic50"],
            "",
        ),
        "netmhcpan_wt_ic50": first(
            row,
            ["IC50 WT", "NetMHCpan WT IC50 Score", "Corresponding WT IC50 Score", "netmhcpan_wt_ic50"],
            "",
        ),
        "netmhcpan_mt_rank_ba": first(
            row,
            ["%ile MT", "NetMHCpan MT Percentile", "Best MT Percentile", "netmhcpan_mt_rank_ba"],
            "",
        ),
        "netmhcpan_wt_rank_ba": first(
            row,
            ["%ile WT", "NetMHCpan WT Percentile", "Corresponding WT Percentile", "netmhcpan_wt_rank_ba"],
            "",
        ),
        "netmhcpan_mt_rank_el": first(
            row,
            ["EL Rank", "Best MT EL Score", "netmhcpan_mt_rank_el"],
            "",
        ),
        "netmhcpan_wt_rank_el": first(row, ["netmhcpan_wt_rank_el"], ""),
    }


def _lookup_stabpan(
    peptide: str,
    allele: str,
    index: dict[tuple[str, str], dict[str, str]] | None,
) -> dict[str, str] | None:
    if not peptide or not allele or not index:
        return None
    return index.get((peptide.strip().upper(), _normalize_hla(allele)))


def _lookup_prime(
    peptide: str,
    allele: str,
    wide_index: dict[str, dict[str, str]] | None,
    pair_index: dict[tuple[str, str], dict[str, str]] | None,
) -> tuple[str, str]:
    from .prime import extract_prime_pair

    pep = peptide.strip().upper()
    hla = _normalize_hla(allele)
    if pep and hla and pair_index:
        rec = pair_index.get((pep, hla))
        if rec:
            return rec.get("prime_score", ""), rec.get("prime_rank", "")
    if pep and hla and wide_index and pep in wide_index:
        return extract_prime_pair(wide_index[pep], pep, hla)
    return "", ""


def _lookup_bigmhc_im(
    peptide: str,
    allele: str,
    index: dict[tuple[str, str], dict[str, str]] | None,
) -> str:
    if not peptide or not allele or not index:
        return ""
    rec = index.get((peptide.strip().upper(), _normalize_hla(allele)))
    if not rec:
        return ""
    return rec.get("bigmhc_im_score", "")


def _lookup_iedb_immunogenicity(
    peptide: str,
    allele: str,
    index: dict[tuple[str, str], dict[str, str]] | None,
) -> str:
    pep = peptide.strip().upper()
    hla = _normalize_hla(allele)
    if not pep or not hla:
        return ""
    if index:
        rec = index.get((pep, hla))
        if rec:
            score = str(rec.get("iedb_immunogenicity_score", "")).strip()
            if score:
                return score
    from .iedb_immunogenicity import predict_immunogenicity

    return f"{predict_immunogenicity(pep, hla):.5f}"


def apply_extra_tool_annotations(
    out: dict[str, str],
    mt_peptide: str,
    wt_peptide: str = "",
    *,
    stabpan_index: dict[tuple[str, str], dict[str, str]] | None = None,
    prime_wide_index: dict[str, dict[str, str]] | None = None,
    prime_pair_index: dict[tuple[str, str], dict[str, str]] | None = None,
    bigmhc_im_index: dict[tuple[str, str], dict[str, str]] | None = None,
    iedb_immunogenicity_index: dict[tuple[str, str], dict[str, str]] | None = None,
) -> dict[str, str]:
    allele = out.get("hla_allele") or ""
    wt = (wt_peptide or "").strip()

    stab = _lookup_stabpan(mt_peptide, allele, stabpan_index) if allele else None
    out["netmhcstabpan_score"] = stab.get("netmhcstabpan_score", "") if stab else ""
    out["netmhcstabpan_rank"] = stab.get("netmhcstabpan_rank", "") if stab else ""
    if wt and allele:
        wt_stab = _lookup_stabpan(wt, allele, stabpan_index)
        out["netmhcstabpan_wt_score"] = wt_stab.get("netmhcstabpan_score", "") if wt_stab else ""
        out["netmhcstabpan_wt_rank"] = wt_stab.get("netmhcstabpan_rank", "") if wt_stab else ""
    else:
        out["netmhcstabpan_wt_score"] = "NA"
        out["netmhcstabpan_wt_rank"] = "NA"

    prime_score, prime_rank = _lookup_prime(
        mt_peptide,
        allele,
        prime_wide_index,
        prime_pair_index,
    )
    out["prime_score"] = prime_score
    out["prime_rank"] = prime_rank
    if wt and allele:
        prime_wt_score, prime_wt_rank = _lookup_prime(
            wt,
            allele,
            prime_wide_index,
            prime_pair_index,
        )
        out["prime_wt_score"] = prime_wt_score
        out["prime_wt_rank"] = prime_wt_rank
    else:
        out["prime_wt_score"] = "NA"
        out["prime_wt_rank"] = "NA"

    out["bigmhc_im_score"] = _lookup_bigmhc_im(mt_peptide, allele, bigmhc_im_index)
    if wt and allele:
        out["bigmhc_im_wt_score"] = _lookup_bigmhc_im(wt, allele, bigmhc_im_index)
    else:
        out["bigmhc_im_wt_score"] = "NA"
    out["iedb_immunogenicity_score"] = _lookup_iedb_immunogenicity(
        mt_peptide,
        allele,
        iedb_immunogenicity_index,
    )
    out["iedb_immunogenicity_wt_score"] = (
        _lookup_iedb_immunogenicity(wt, allele, iedb_immunogenicity_index)
        if wt and allele else "NA"
    )

    for col in TOOL_ANNOTATION_FIELDS:
        out.setdefault(col, "")
    return out


def annotate_variant_peptide_row(
    row: dict[str, str],
    hla_alleles: list[str],
    index: dict[tuple[str, str], dict[str, str]] | None = None,
    *,
    mhcflurry_index: dict[tuple[str, str], dict[str, str]] | None = None,
    stabpan_index: dict[tuple[str, str], dict[str, str]] | None = None,
    prime_wide_index: dict[str, dict[str, str]] | None = None,
    prime_pair_index: dict[tuple[str, str], dict[str, str]] | None = None,
    bigmhc_im_index: dict[tuple[str, str], dict[str, str]] | None = None,
    iedb_immunogenicity_index: dict[tuple[str, str], dict[str, str]] | None = None,
    fetch_missing: bool = False,
    prefer_local_netmhcpan: bool = False,
    cache: dict[tuple[str, str], dict[str, str]] | None = None,
    mt_key: str = "mutant_peptide",
    wt_key: str = "wildtype_peptide",
    peptide_key: str | None = None,
) -> dict[str, str]:
    out = dict(row)
    out["peptide_source"] = resolve_peptide_source(row)
    out["sample_hla_alleles"] = format_sample_hla_alleles(hla_alleles)

    mt_peptide = first(row, [mt_key, peptide_key or "", "peptide", "Best Peptide", "MT Epitope Seq"], "")
    wt_peptide = first(row, [wt_key, "WT Epitope Seq", "WT Epitope", "wildtype_peptide"], "")

    existing = {} if prefer_local_netmhcpan else netmhcpan_columns_from_pvac_row(row)
    if existing.get("hla_allele") and existing.get("netmhcpan_mt_ic50"):
        out.update(existing)
        if not out.get("netmhcpan_wt_ic50") and wt_peptide and out.get("hla_allele"):
            wt_pred = _lookup_prediction(
                wt_peptide,
                out["hla_allele"],
                index,
                fetch_missing=fetch_missing,
                cache=cache or {},
            )
            if wt_pred:
                out["netmhcpan_wt_ic50"] = wt_pred.get("netmhcpan_ba_score", "")
                out["netmhcpan_wt_rank_ba"] = wt_pred.get("netmhcpan_ba_rank", "")
                out["netmhcpan_wt_rank_el"] = wt_pred.get("netmhcpan_el_rank", "")
        for col in NETMHCPAN_ANNOTATION_FIELDS:
            out.setdefault(col, "")
        apply_mhcflurry_annotation(
            out,
            mt_peptide,
            wt_peptide,
            hla_alleles,
            mhcflurry_index,
            bind_hla=out.get("hla_allele", ""),
        )
        return apply_extra_tool_annotations(
            out,
            mt_peptide,
            wt_peptide,
            stabpan_index=stabpan_index,
            prime_wide_index=prime_wide_index,
            prime_pair_index=prime_pair_index,
            bigmhc_im_index=bigmhc_im_index,
            iedb_immunogenicity_index=iedb_immunogenicity_index,
        )

    row_hla = first(row, ["Allele", "hla_allele", "HLA Allele"], "")
    if prefer_local_netmhcpan and row_hla and index:
        bind_hla = _normalize_hla(row_hla)
        mt_pred = _lookup_prediction(mt_peptide, bind_hla, index, fetch_missing=fetch_missing, cache=cache or {})
    else:
        bind_hla, mt_pred = pick_best_hla_for_peptide(
            mt_peptide,
            hla_alleles,
            index,
            fetch_missing=fetch_missing,
            cache=cache,
        )
        bind_hla = bind_hla or ""
    out["hla_allele"] = bind_hla or row_hla or ""
    if not out["hla_allele"] and hla_alleles and mt_peptide:
        out["hla_allele"] = _normalize_hla(hla_alleles[0])
    if mt_pred:
        out["netmhcpan_mt_ic50"] = mt_pred.get("netmhcpan_ba_score", "")
        out["netmhcpan_mt_rank_ba"] = mt_pred.get("netmhcpan_ba_rank", "")
        out["netmhcpan_mt_rank_el"] = mt_pred.get("netmhcpan_el_rank", "")

    if wt_peptide and out.get("hla_allele"):
        wt_pred = _lookup_prediction(
            wt_peptide,
            out["hla_allele"],
            index,
            fetch_missing=fetch_missing,
            cache=cache or {},
        )
        if wt_pred:
            out["netmhcpan_wt_ic50"] = wt_pred.get("netmhcpan_ba_score", "")
            out["netmhcpan_wt_rank_ba"] = wt_pred.get("netmhcpan_ba_rank", "")
            out["netmhcpan_wt_rank_el"] = wt_pred.get("netmhcpan_el_rank", "")

    for col in NETMHCPAN_ANNOTATION_FIELDS:
        out.setdefault(col, "")
    apply_mhcflurry_annotation(
        out,
        mt_peptide,
        wt_peptide,
        hla_alleles,
        mhcflurry_index,
        bind_hla=out.get("hla_allele", ""),
    )
    return apply_extra_tool_annotations(
        out,
        mt_peptide,
        wt_peptide,
        stabpan_index=stabpan_index,
        prime_wide_index=prime_wide_index,
        prime_pair_index=prime_pair_index,
        bigmhc_im_index=bigmhc_im_index,
        iedb_immunogenicity_index=iedb_immunogenicity_index,
    )


def annotate_variant_peptide_rows(
    rows: list[dict[str, str]],
    hla_alleles: list[str],
    *,
    netmhcpan_xls: str | Path | None = None,
    mhcflurry_csv: str | Path | None = None,
    netmhcstabpan_tsv: str | Path | None = None,
    prime_tsv: str | Path | None = None,
    bigmhc_im_tsv: str | Path | None = None,
    iedb_immunogenicity_tsv: str | Path | None = None,
    index: dict[tuple[str, str], dict[str, str]] | None = None,
    mhcflurry_index: dict[tuple[str, str], dict[str, str]] | None = None,
    stabpan_index: dict[tuple[str, str], dict[str, str]] | None = None,
    prime_wide_index: dict[str, dict[str, str]] | None = None,
    prime_pair_index: dict[tuple[str, str], dict[str, str]] | None = None,
    bigmhc_im_index: dict[tuple[str, str], dict[str, str]] | None = None,
    iedb_immunogenicity_index: dict[tuple[str, str], dict[str, str]] | None = None,
    fetch_missing: bool = False,
    prefer_local_netmhcpan: bool = False,
) -> list[dict[str, str]]:
    if not hla_alleles:
        return rows
    lookup_index = index
    if lookup_index is None and netmhcpan_xls:
        lookup_index = build_netmhcpan_index(netmhcpan_xls)
        prefer_local_netmhcpan = True
    lookup_mhcflurry = mhcflurry_index
    if lookup_mhcflurry is None and mhcflurry_csv:
        lookup_mhcflurry = build_mhcflurry_index(mhcflurry_csv)
    lookup_stabpan = stabpan_index
    if lookup_stabpan is None and netmhcstabpan_tsv:
        lookup_stabpan = build_netmhcstabpan_index(netmhcstabpan_tsv)
    lookup_prime_wide = prime_wide_index
    lookup_prime_pair = prime_pair_index
    if lookup_prime_wide is None and lookup_prime_pair is None and prime_tsv:
        lookup_prime_wide, lookup_prime_pair = build_prime_indexes(prime_tsv)
    lookup_bigmhc = bigmhc_im_index
    if lookup_bigmhc is None and bigmhc_im_tsv:
        lookup_bigmhc = build_bigmhc_im_index(bigmhc_im_tsv)
    lookup_iedb = iedb_immunogenicity_index
    if lookup_iedb is None and iedb_immunogenicity_tsv:
        lookup_iedb = build_iedb_immunogenicity_index(iedb_immunogenicity_tsv)
    cache: dict[tuple[str, str], dict[str, str]] = {}
    return [
        annotate_variant_peptide_row(
            row,
            hla_alleles,
            lookup_index,
            mhcflurry_index=lookup_mhcflurry,
            stabpan_index=lookup_stabpan,
            prime_wide_index=lookup_prime_wide,
            prime_pair_index=lookup_prime_pair,
            bigmhc_im_index=lookup_bigmhc,
            iedb_immunogenicity_index=lookup_iedb,
            fetch_missing=fetch_missing,
            prefer_local_netmhcpan=prefer_local_netmhcpan,
            cache=cache,
        )
        for row in rows
    ]


def annotate_raw_peptide_row(
    row: dict[str, str],
    index: dict[tuple[str, str], dict[str, str]] | None = None,
    *,
    mhcflurry_index: dict[tuple[str, str], dict[str, str]] | None = None,
    cache: dict[tuple[str, str], dict[str, str]] | None = None,
) -> dict[str, str]:
    """Fill MT/WT NetMHCpan / MHCflurry binding on a raw_peptides row."""
    out = dict(row)
    peptide = first(row, ["peptide", "mutant_peptide", "Best Peptide"], "")
    hla = first(row, ["hla_allele", "HLA Allele", "Allele"], "")
    wt_peptide = first(row, ["wildtype_peptide", "WT Epitope Seq"], "")
    lookup_cache = cache if cache is not None else {}

    if peptide and hla and index:
        mt_pred = _lookup_prediction(
            peptide, hla, index, fetch_missing=False, cache=lookup_cache
        )
        if mt_pred:
            out["netmhcpan_mt_ic50"] = mt_pred.get("netmhcpan_ba_score", "")
            out["netmhcpan_mt_rank_ba"] = mt_pred.get("netmhcpan_ba_rank", "")
            out["netmhcpan_mt_rank_el"] = mt_pred.get("netmhcpan_el_rank", "")
            out["netmhcpan_ba_rank"] = mt_pred.get("netmhcpan_ba_rank", "")
            out["netmhcpan_el_rank"] = mt_pred.get("netmhcpan_el_rank", "")
            out["binding_rank"] = mt_pred.get("netmhcpan_ba_rank", out.get("binding_rank", "99"))

    if wt_peptide and hla and index:
        wt_pred = _lookup_prediction(
            wt_peptide, hla, index, fetch_missing=False, cache=lookup_cache
        )
        if wt_pred:
            out["netmhcpan_wt_ic50"] = wt_pred.get("netmhcpan_ba_score", "")
            out["netmhcpan_wt_rank_ba"] = wt_pred.get("netmhcpan_ba_rank", "")
            out["netmhcpan_wt_rank_el"] = wt_pred.get("netmhcpan_el_rank", "")
            out["wildtype_binding_rank"] = wt_pred.get("netmhcpan_ba_rank", out.get("wildtype_binding_rank", "99"))

    if peptide and hla and mhcflurry_index:
        mt_mhc = _lookup_mhcflurry(peptide, hla, mhcflurry_index)
        if mt_mhc:
            out["mhcflurry_affinity_percentile"] = mt_mhc.get("mhcflurry_affinity_percentile", "")
            out["mhcflurry_processing_score"] = mt_mhc.get("mhcflurry_processing_score", "")
            out["mhcflurry_presentation_score"] = mt_mhc.get("mhcflurry_presentation_score", "")
            out.update(_mhcflurry_columns_from_pred(mt_mhc, prefix="mt"))
        if wt_peptide:
            wt_mhc = _lookup_mhcflurry(wt_peptide, hla, mhcflurry_index)
            out.update(_mhcflurry_columns_from_pred(wt_mhc, prefix="wt"))

    return out


def annotate_raw_peptides_tsv(
    raw_peptides_tsv: str | Path,
    netmhcpan_xls: str | Path | None = None,
    mhcflurry_csv: str | Path | None = None,
) -> dict[str, Any]:
    from ..schemas import PEPTIDE_FIELDS

    path = Path(raw_peptides_tsv)
    rows_in = read_tsv(path)
    if not rows_in:
        return {"raw_peptides": str(path), "rows": 0}
    index = build_netmhcpan_index(netmhcpan_xls) if netmhcpan_xls else None
    mhc_index = build_mhcflurry_index(mhcflurry_csv) if mhcflurry_csv else None
    cache: dict[tuple[str, str], dict[str, str]] = {}
    rows_out = [
        annotate_raw_peptide_row(row, index, mhcflurry_index=mhc_index, cache=cache)
        for row in rows_in
    ]
    fieldnames = list(rows_in[0].keys())
    for col in PEPTIDE_FIELDS:
        if col not in fieldnames:
            fieldnames.append(col)
    for col in MHCFLURRY_ANNOTATION_FIELDS:
        if col not in fieldnames:
            fieldnames.append(col)
    write_tsv(path, rows_out, fieldnames)
    return {"raw_peptides": str(path), "rows": len(rows_out)}


def annotate_variant_peptide_tsv(
    tsv_path: str | Path,
    hla_alleles: list[str],
    *,
    netmhcpan_xls: str | Path | None = None,
    mhcflurry_csv: str | Path | None = None,
    netmhcstabpan_tsv: str | Path | None = None,
    prime_tsv: str | Path | None = None,
    bigmhc_im_tsv: str | Path | None = None,
    iedb_immunogenicity_tsv: str | Path | None = None,
    output_tsv: str | Path | None = None,
    fetch_missing: bool = False,
    fieldnames: list[str] | None = None,
) -> dict[str, Any]:
    path = Path(tsv_path)
    rows_in = read_tsv(path)
    rows_out = annotate_variant_peptide_rows(
        rows_in,
        hla_alleles,
        netmhcpan_xls=netmhcpan_xls,
        mhcflurry_csv=mhcflurry_csv,
        netmhcstabpan_tsv=netmhcstabpan_tsv,
        prime_tsv=prime_tsv,
        bigmhc_im_tsv=bigmhc_im_tsv,
        iedb_immunogenicity_tsv=iedb_immunogenicity_tsv,
        fetch_missing=fetch_missing,
    )
    out_path = Path(output_tsv) if output_tsv else path
    names = list(fieldnames or (list(rows_in[0].keys()) if rows_in else []))
    if "peptide_source" not in names:
        if "generation_method" in names:
            names.insert(names.index("generation_method"), "peptide_source")
        else:
            names.append("peptide_source")
    for col in TOOL_ANNOTATION_FIELDS:
        if col not in names:
            names.append(col)
    write_tsv(out_path, rows_out, names)
    return {
        "input_tsv": str(path),
        "output_tsv": str(out_path),
        "rows": len(rows_out),
        "sample_hla_alleles": format_sample_hla_alleles(hla_alleles),
    }
