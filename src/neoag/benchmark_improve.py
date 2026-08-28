"""Benchmark presentation predictors against IMPROVE / CEDAR peptide–HLA immunogenicity labels."""

from __future__ import annotations

import math
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .adapters.mhcflurry import parse_mhcflurry, write_mhcflurry_evidence
from .adapters.netmhcpan import parse_netmhcpan, write_netmhcpan_evidence
from .adapters.netmhcstabpan import parse_netmhcstabpan, write_netmhcstabpan_evidence
from .config import load_profile
from .gates import evaluate_presentation_gate
from .immunogenicity_composite import (
    apply_immunogenicity_composite,
    apply_immunogenicity_evidence,
    immunogenicity_config,
    run_immunogenicity_predictors,
)
from .presentation import build_presentation_evidence
from .scoring import compute_peptide_efficacy
from .schemas import PEPTIDE_FIELDS
from .tools.registry import RunContext
from .tools.runner import run_mhcflurry, run_netmhcpan, run_netmhcstabpan
from .utils import norm_rank, read_tsv, safe_id, to_float, write_json, write_tsv

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CEDAR = ROOT / "data" / "improve" / "Neoepitopes_CEDAR_benchmark_data.tsv"
DEFAULT_INHOUSE = ROOT / "data" / "improve" / "In_house_neoepitope_for_CV.tsv"
DEFAULT_TOOLS_CACHE = ROOT / "results" / "benchmark_cedar_v2" / "tools"
TOOL_CACHE_FILES = {
    "netmhcpan": "netmhcpan.xls",
    "mhcflurry": "mhcflurry.csv",
    "stabpan": "netmhcstabpan.tsv",
    "prime": "prime.tsv",
    "bigmhc_im": "bigmhc_im.tsv",
}

METRIC_FIELDS = [
    "cohort",
    "predictor",
    "n",
    "n_positive",
    "auroc",
    "auprc",
    "mean_score_pos",
    "mean_score_neg",
    "direction",
]
HLA_METRIC_FIELDS = [
    "hla_allele",
    "predictor",
    "n",
    "n_positive",
    "auroc",
    "auprc",
    "direction",
]
COHORT_LABELS = {
    "all": "All peptide–HLA pairs",
    "binder_el2": "Strong binders (NetMHCpan EL rank ≤ 2%)",
    "binder_pres05": "Presentation score ≥ 0.5",
    "binder_grade_ab": "Presentation grade A or B",
    "gated_tesla": "TESLA-style gate pass (grade A/B, EL≤2%, stabpan≥1.4h)",
}

BENCHMARK_EVENT_STUB = {
    "event_score": "0.5000",
    "event_expression": "0",
    "tumor_vaf": "0",
    "safety_status": "PASS",
}
BENCHMARK_PEPTIDE_STUB = {
    "immunogenicity_score": "0.5",
    "safety_status": "PASS",
}
HLA_MIN_N = 30
HLA_MIN_POS = 5
HLA_MIN_NEG = 5


@dataclass(frozen=True)
class ImproveRecord:
    peptide: str
    hla_allele: str
    response: int


def normalize_hla(hla: str) -> str:
    """Normalize IMPROVE HLA names (e.g. HLA-A02:01, HLA-A11*01) to HLA-A*02:01 style."""
    s = hla.strip().upper()
    if s.startswith("HLA-"):
        s = s[4:]
    elif s.startswith("HLA"):
        s = s[3:]
    m = re.match(r"^([ABCD])(\d{1,2})[\*:](\d{1,2})$", s)
    if not m:
        raise ValueError(f"Unrecognized HLA allele format: {hla!r}")
    return f"HLA-{m.group(1)}*{int(m.group(2)):02d}:{int(m.group(3)):02d}"


def load_improve_tsv(path: str | Path, limit: int | None = None) -> list[ImproveRecord]:
    """Load space-separated IMPROVE/CEDAR tables (HLA peptide response)."""
    rows: list[ImproveRecord] = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if parts[0] in {"HLA_allele", "Mut_peptide"} or parts[-1] == "response":
            continue
        if len(parts) < 3:
            continue
        if parts[0].startswith("HLA"):
            hla_raw, peptide, resp = parts[0], parts[1], parts[2]
        else:
            peptide, hla_raw, resp = parts[0], parts[1], parts[2]
        rows.append(
            ImproveRecord(
                peptide=peptide,
                hla_allele=normalize_hla(hla_raw),
                response=int(resp),
            )
        )
        if limit and len(rows) >= limit:
            break
    if not rows:
        raise ValueError(f"No records parsed from {path}")
    return rows


def auroc(labels: list[int], scores: list[float], higher_is_better: bool = True) -> float:
    """Binary AUROC without external deps (Mann–Whitney / rank statistic)."""
    if len(labels) != len(scores):
        raise ValueError("labels and scores length mismatch")
    n_pos = sum(labels)
    n_neg = len(labels) - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    paired = [(s, y) for s, y in zip(scores, labels)]
    if not higher_is_better:
        paired = [(-s, y) for s, y in paired]
    paired.sort(key=lambda x: x[0])
    rank = 1
    sum_ranks_pos = 0.0
    i = 0
    while i < len(paired):
        j = i
        while j < len(paired) and paired[j][0] == paired[i][0]:
            j += 1
        avg_rank = (rank + rank + (j - i) - 1) / 2.0
        for k in range(i, j):
            if paired[k][1] == 1:
                sum_ranks_pos += avg_rank
        rank += j - i
        i = j
    u = sum_ranks_pos - n_pos * (n_pos + 1) / 2.0
    return u / (n_pos * n_neg)


def auprc(labels: list[int], scores: list[float], higher_is_better: bool = True) -> float:
    """Average precision (area under PR curve, trapezoid)."""
    if sum(labels) == 0:
        return float("nan")
    paired = sorted(zip(scores, labels), reverse=higher_is_better)
    tp = fp = 0
    precisions: list[float] = []
    recalls: list[float] = []
    total_pos = sum(labels)
    for score, label in paired:
        if label == 1:
            tp += 1
        else:
            fp += 1
        if tp + fp == 0:
            continue
        precisions.append(tp / (tp + fp))
        recalls.append(tp / total_pos)
    if not precisions:
        return 0.0
    auc = 0.0
    for i in range(1, len(recalls)):
        auc += (recalls[i] - recalls[i - 1]) * (precisions[i] + precisions[i - 1]) / 2.0
    return auc


def _records_to_raw_peptides(records: list[ImproveRecord], sample_id: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for i, rec in enumerate(records):
        pid = _peptide_id(i, rec)
        rows.append({
            "peptide_id": pid,
            "event_id": pid,
            "sample_id": sample_id,
            "event_type": "Benchmark",
            "gene": "IMPROVE",
            "peptide": rec.peptide,
            "wildtype_peptide": "",
            "hla_allele": rec.hla_allele,
            "mhc_class": "I",
            "source_tool": "IMPROVE",
            "binding_rank": "99",
            "el_rank": "99",
            "presentation_score": "0",
            "immunogenicity_score": "0.5",
            "wildtype_binding_rank": "99",
            "self_similarity_score": "0",
            "normal_hla_ligand_overlap": "no",
        })
    return rows


def _metric_row(
    name: str,
    labels: list[int],
    scores: list[float],
    *,
    higher_is_better: bool = True,
) -> dict[str, str]:
    valid = [
        (y, s)
        for y, s in zip(labels, scores)
        if not math.isnan(s) and not math.isinf(s)
    ]
    if not valid:
        return {
            "predictor": name,
            "n": str(len(labels)),
            "auroc": "nan",
            "auprc": "nan",
            "mean_score_pos": "nan",
            "mean_score_neg": "nan",
            "direction": "higher_better" if higher_is_better else "lower_better",
        }
    ys, ss = zip(*valid)
    pos_scores = [s for y, s in valid if y == 1]
    neg_scores = [s for y, s in valid if y == 0]
    return {
        "predictor": name,
        "n": str(len(valid)),
        "n_positive": str(sum(ys)),
        "auroc": f"{auroc(list(ys), list(ss), higher_is_better=higher_is_better):.4f}",
        "auprc": f"{auprc(list(ys), list(ss), higher_is_better=higher_is_better):.4f}",
        "mean_score_pos": f"{(sum(pos_scores) / len(pos_scores)) if pos_scores else float('nan'):.4f}",
        "mean_score_neg": f"{(sum(neg_scores) / len(neg_scores)) if neg_scores else float('nan'):.4f}",
        "direction": "higher_better" if higher_is_better else "lower_better",
    }


def _peptide_id(i: int, rec: ImproveRecord) -> str:
    return safe_id(f"IMPROVE_{i}_{rec.peptide}_{rec.hla_allele}")


def _predictor_specs(profile: dict[str, Any] | None = None) -> list[tuple[str, Any, bool]]:
    specs = [
        ("efficacy_score", lambda r: to_float(r.get("efficacy_score"), 0.0), True),
        ("immunogenicity_composite_score", lambda r: to_float(r.get("immunogenicity_composite_score"), 0.0), True),
        ("prime_score", lambda r: to_float(r.get("prime_score"), 0.0), True),
        ("prime_rank", lambda r: to_float(r.get("prime_rank"), 99.0), False),
        ("bigmhc_im_score", lambda r: to_float(r.get("bigmhc_im_score"), 0.0), True),
        ("deepimmuno_score", lambda r: to_float(r.get("deepimmuno_score"), 0.0), True),
        ("binding_evidence_score", lambda r: to_float(r.get("binding_evidence_score"), 0.0), True),
        ("presentation_evidence_score", lambda r: to_float(r.get("presentation_evidence_score"), 0.0), True),
        ("netmhcpan_ba_rank", lambda r: norm_rank(r.get("netmhcpan_ba_rank", 99)), True),
        ("netmhcpan_el_rank", lambda r: norm_rank(r.get("netmhcpan_el_rank", 99)), True),
        ("netmhcstabpan_score", lambda r: to_float(r.get("netmhcstabpan_score"), 0.0), True),
        ("netmhcstabpan_rank", lambda r: norm_rank(r.get("netmhcstabpan_rank", 99)), True),
        ("mhcflurry_affinity_percentile", lambda r: norm_rank(r.get("mhcflurry_affinity_percentile", 99)), True),
        ("mhcflurry_presentation_score", lambda r: to_float(r.get("mhcflurry_presentation_score"), 0.0), True),
    ]
    cfg = immunogenicity_config(profile or {})
    if cfg.get("use_iedb_fallback") or "iedb" in cfg.get("sources", []):
        specs.insert(3, ("iedb_immunogenicity_score", lambda r: to_float(r.get("iedb_immunogenicity_score"), 0.0), True))
    return specs


def _enrich_pres_with_efficacy(
    records: list[ImproveRecord],
    pres_rows: list[dict[str, str]],
    profile: dict[str, Any],
) -> None:
    pres_by_id = {r["peptide_id"]: r for r in pres_rows}
    for i, rec in enumerate(records):
        pid = _peptide_id(i, rec)
        pres = pres_by_id.get(pid, {})
        peptide = {
            **BENCHMARK_PEPTIDE_STUB,
            "peptide": rec.peptide,
            "hla_allele": rec.hla_allele,
        }
        pres.update(
            compute_peptide_efficacy(
                peptide,
                BENCHMARK_EVENT_STUB,
                pres,
                profile,
                appm=1.0,
                ccf=1.0,
            )
        )


def _cohort_masks(
    records: list[ImproveRecord],
    pres_by_id: dict[str, dict[str, str]],
    profile: dict[str, Any],
) -> dict[str, list[bool]]:
    masks: dict[str, list[bool]] = {"all": [True] * len(records)}
    binder_el2: list[bool] = []
    binder_pres: list[bool] = []
    binder_ab: list[bool] = []
    gated: list[bool] = []
    for i, rec in enumerate(records):
        pres = pres_by_id.get(_peptide_id(i, rec), {})
        el = to_float(pres.get("netmhcpan_el_rank"), 99.0)
        pe = to_float(pres.get("presentation_evidence_score"), 0.0)
        grade = pres.get("presentation_evidence_grade", "")
        binder_el2.append(0 < el <= 2.0)
        binder_pres.append(pe >= 0.5)
        binder_ab.append(grade in {"A", "B"})
        gate = evaluate_presentation_gate(
            {**BENCHMARK_PEPTIDE_STUB, "peptide": rec.peptide, "hla_allele": rec.hla_allele},
            BENCHMARK_EVENT_STUB,
            pres,
            profile,
        )
        gated.append(gate["presentation_gate_status"] == "PASS")
    masks["binder_el2"] = binder_el2
    masks["binder_pres05"] = binder_pres
    masks["binder_grade_ab"] = binder_ab
    masks["gated_tesla"] = gated
    return masks


def _subset_labels_scores(
    labels: list[int],
    scores: list[float],
    mask: list[bool],
) -> tuple[list[int], list[float]]:
    sub_labels = [y for y, keep in zip(labels, mask) if keep]
    sub_scores = [s for s, keep in zip(scores, mask) if keep]
    return sub_labels, sub_scores


def _compute_cohort_metrics(
    cohort: str,
    records: list[ImproveRecord],
    labels: list[int],
    presentation_rows: list[dict[str, str]],
    mask: list[bool],
    profile: dict[str, Any],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for name, scores, higher in _collect_predictor_scores(records, presentation_rows, profile):
        sub_labels, sub_scores = _subset_labels_scores(labels, scores, mask)
        row = _metric_row(name, sub_labels, sub_scores, higher_is_better=higher)
        row["cohort"] = cohort
        rows.append(row)
    return rows


def _compute_hla_metrics(
    records: list[ImproveRecord],
    labels: list[int],
    presentation_rows: list[dict[str, str]],
    profile: dict[str, Any],
    *,
    predictor: str = "immunogenicity_composite_score",
) -> list[dict[str, str]]:
    pres_by_id = {r["peptide_id"]: r for r in presentation_rows}
    specs = {name: (fn, higher) for name, fn, higher in _predictor_specs(profile)}
    if predictor not in specs:
        return []
    fn, higher = specs[predictor]
    by_hla: dict[str, list[tuple[int, float]]] = {}
    for i, rec in enumerate(records):
        pres = pres_by_id.get(_peptide_id(i, rec), {})
        by_hla.setdefault(rec.hla_allele, []).append((rec.response, fn(pres)))
    rows: list[dict[str, str]] = []
    for hla in sorted(by_hla):
        pairs = by_hla[hla]
        ys = [y for y, _ in pairs]
        ss = [s for _, s in pairs]
        n_pos = sum(ys)
        n_neg = len(ys) - n_pos
        if len(ys) < HLA_MIN_N or n_pos < HLA_MIN_POS or n_neg < HLA_MIN_NEG:
            continue
        row = _metric_row(predictor, ys, ss, higher_is_better=higher)
        row["hla_allele"] = hla
        rows.append(row)
    rows.sort(key=lambda r: to_float(r.get("auroc"), -1), reverse=True)
    return rows


def _merge_stabpan(
    presentation_rows: list[dict[str, str]],
    stab_by_pair: dict[tuple[str, str], dict[str, str]],
) -> None:
    for row in presentation_rows:
        key = (row.get("peptide", ""), row.get("hla_allele", ""))
        stab = stab_by_pair.get(key, {})
        row["netmhcstabpan_score"] = stab.get("netmhcstabpan_score", "")
        row["netmhcstabpan_rank"] = stab.get("netmhcstabpan_rank", "")


def _report_metrics_table(metric_rows: list[dict[str, str]]) -> list[str]:
    lines = [
        "| Predictor | AUROC | AUPRC | n | pos | mean(pos) | mean(neg) |",
        "|-----------|-------|-------|---|-----|-----------|-----------|",
    ]
    for row in metric_rows:
        lines.append(
            f"| {row['predictor']} | {row['auroc']} | {row['auprc']} | {row['n']} | "
            f"{row.get('n_positive', '')} | {row['mean_score_pos']} | {row['mean_score_neg']} |"
        )
    return lines


def _seed_tools_cache(
    tools_dir: Path,
    reuse_from: Path | None,
    *,
    skip_netmhcpan: bool,
    skip_mhcflurry: bool,
    skip_stabpan: bool,
    seed_immuno: bool = False,
) -> tuple[Path | None, list[str]]:
    """Copy cached presentation tool outputs when predictors are skipped."""
    cache_root = reuse_from
    if cache_root is None and DEFAULT_TOOLS_CACHE.is_dir():
        cache_root = DEFAULT_TOOLS_CACHE
    if not cache_root or not cache_root.is_dir():
        return None, []
    copied: list[str] = []
    seeds = [
        (skip_netmhcpan, "netmhcpan"),
        (skip_mhcflurry, "mhcflurry"),
        (skip_stabpan, "stabpan"),
    ]
    if seed_immuno:
        seeds.extend([(True, "prime"), (True, "bigmhc_im")])
    for should_seed, key in seeds:
        if not should_seed:
            continue
        name = TOOL_CACHE_FILES[key]
        src = cache_root / name
        dst = tools_dir / name
        if src.is_file() and not dst.is_file():
            shutil.copy2(src, dst)
            copied.append(name)
    return cache_root, copied


def _collect_predictor_scores(
    records: list[ImproveRecord],
    presentation_rows: list[dict[str, str]],
    profile: dict[str, Any],
) -> list[tuple[str, list[float], bool]]:
    pres_by_id = {r["peptide_id"]: r for r in presentation_rows}
    scored: list[tuple[str, list[float], bool]] = []
    for name, fn, higher in _predictor_specs(profile):
        scores = [
            fn(pres_by_id.get(_peptide_id(i, rec), {}))
            for i, rec in enumerate(records)
        ]
        scored.append((name, scores, higher))
    return scored


def run_benchmark_improve(
    *,
    input_path: str | Path,
    outdir: str | Path,
    profile_name: str = "default",
    sample_id: str = "IMPROVE_CEDAR",
    stub: bool = False,
    limit: int | None = None,
    skip_netmhcpan: bool = False,
    skip_mhcflurry: bool = False,
    skip_immunogenicity: bool = False,
    skip_stabpan: bool = False,
    reuse_tools_dir: str | Path | None = None,
) -> dict[str, Any]:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    records = load_improve_tsv(input_path, limit=limit)
    labels = [r.response for r in records]
    profile = load_profile(profile_name)

    raw_pep_path = outdir / "raw_peptides.tsv"
    raw_rows = _records_to_raw_peptides(records, sample_id)
    write_tsv(raw_pep_path, raw_rows, PEPTIDE_FIELDS)

    ctx = RunContext(
        sample_id=sample_id,
        outdir=outdir,
        stub=stub,
        hla_alleles=sorted({r.hla_allele for r in records}),
        raw_peptides=raw_pep_path,
    )
    tools_dir = outdir / "tools"
    tools_dir.mkdir(parents=True, exist_ok=True)
    cache_root, seeded_tools = _seed_tools_cache(
        tools_dir,
        Path(reuse_tools_dir) if reuse_tools_dir else None,
        skip_netmhcpan=skip_netmhcpan,
        skip_mhcflurry=skip_mhcflurry,
        skip_stabpan=skip_stabpan,
        seed_immuno=bool(reuse_tools_dir),
    )

    netmhcpan_path = None
    net_xls = tools_dir / "netmhcpan.xls"
    if not skip_netmhcpan:
        run_netmhcpan(ctx, net_xls)
    if net_xls.is_file():
        netmhcpan_path = outdir / "presentation" / "netmhcpan_evidence.tsv"
        netmhcpan_path.parent.mkdir(parents=True, exist_ok=True)
        write_netmhcpan_evidence(netmhcpan_path, parse_netmhcpan(net_xls, sample_id))

    mhcflurry_path = None
    mhc_csv = tools_dir / "mhcflurry.csv"
    if not skip_mhcflurry:
        run_mhcflurry(ctx, mhc_csv)
    if mhc_csv.is_file():
        mhcflurry_path = outdir / "presentation" / "mhcflurry_evidence.tsv"
        mhcflurry_path.parent.mkdir(parents=True, exist_ok=True)
        write_mhcflurry_evidence(mhcflurry_path, parse_mhcflurry(mhc_csv, sample_id))

    stabpan_path = None
    stab_rows: list[dict[str, str]] = []
    stab_tsv = tools_dir / "netmhcstabpan.tsv"
    if not skip_stabpan:
        run_netmhcstabpan(ctx, stab_tsv)
    if stab_tsv.is_file():
        stabpan_path = outdir / "presentation" / "netmhcstabpan_evidence.tsv"
        stabpan_path.parent.mkdir(parents=True, exist_ok=True)
        stab_rows = parse_netmhcstabpan(stab_tsv, sample_id)
        write_netmhcstabpan_evidence(stabpan_path, stab_rows)

    pres_path = outdir / "presentation" / "presentation_evidence.tsv"
    pres_rows = build_presentation_evidence(
        raw_pep_path,
        str(netmhcpan_path) if netmhcpan_path else None,
        str(mhcflurry_path) if mhcflurry_path else None,
        profile,
        pres_path,
    )

    immuno_skip_sources: set[str] = set()
    if (tools_dir / "prime.tsv").is_file() and "prime.tsv" in seeded_tools:
        immuno_skip_sources.add("prime")
    if (tools_dir / "bigmhc_im.tsv").is_file() and "bigmhc_im.tsv" in seeded_tools:
        immuno_skip_sources.add("bigmhc_im")

    immuno_paths = run_immunogenicity_predictors(
        raw_pep_path,
        outdir,
        profile,
        ctx,
        skip=skip_immunogenicity,
        skip_sources=immuno_skip_sources,
    )
    apply_immunogenicity_evidence(pres_rows, immuno_paths, profile)

    if stab_rows:
        stab_by_pair = {
            (r["peptide"], r["hla_allele"]): r
            for r in stab_rows
        }
        _merge_stabpan(pres_rows, stab_by_pair)

    _enrich_pres_with_efficacy(records, pres_rows, profile)
    pres_by_id = {r["peptide_id"]: r for r in pres_rows}
    cohort_masks = _cohort_masks(records, pres_by_id, profile)

    metric_rows: list[dict[str, str]] = []
    for cohort, mask in cohort_masks.items():
        if cohort != "all" and sum(mask) == 0:
            continue
        metric_rows.extend(
            _compute_cohort_metrics(cohort, records, labels, pres_rows, mask, profile)
        )
    all_rows = [r for r in metric_rows if r["cohort"] == "all"]
    all_rows.sort(key=lambda r: to_float(r.get("auroc"), -1), reverse=True)
    metrics_path = outdir / "benchmark_metrics.tsv"
    write_tsv(metrics_path, metric_rows, METRIC_FIELDS)

    hla_metric_rows = []
    if not skip_immunogenicity or any(immuno_paths.values()):
        hla_metric_rows = _compute_hla_metrics(records, labels, pres_rows, profile)
        write_tsv(outdir / "benchmark_metrics_by_hla.tsv", hla_metric_rows, HLA_METRIC_FIELDS)
    pred_rows: list[dict[str, str]] = []
    for i, rec in enumerate(records):
        pid = _peptide_id(i, rec)
        pres = pres_by_id.get(pid, {})
        pred_rows.append({
            "peptide": rec.peptide,
            "hla_allele": rec.hla_allele,
            "response": str(rec.response),
            "prime_score": pres.get("prime_score", ""),
            "prime_rank": pres.get("prime_rank", ""),
            "bigmhc_im_score": pres.get("bigmhc_im_score", ""),
            "deepimmuno_score": pres.get("deepimmuno_score", ""),
            "immunogenicity_composite_score": pres.get("immunogenicity_composite_score", ""),
            "immunogenicity_source": pres.get("immunogenicity_source", ""),
            "iedb_immunogenicity_score": pres.get("iedb_immunogenicity_score", ""),
            "binding_evidence_score": pres.get("binding_evidence_score", ""),
            "presentation_evidence_score": pres.get("presentation_evidence_score", ""),
            "netmhcpan_ba_rank": pres.get("netmhcpan_ba_rank", ""),
            "netmhcpan_el_rank": pres.get("netmhcpan_el_rank", ""),
            "netmhcstabpan_score": pres.get("netmhcstabpan_score", ""),
            "netmhcstabpan_rank": pres.get("netmhcstabpan_rank", ""),
            "mhcflurry_affinity_percentile": pres.get("mhcflurry_affinity_percentile", ""),
            "mhcflurry_presentation_score": pres.get("mhcflurry_presentation_score", ""),
            "presentation_evidence_grade": pres.get("presentation_evidence_grade", ""),
            "efficacy_score": pres.get("efficacy_score", ""),
            "presentation_gate_status": pres.get("presentation_gate_status", ""),
            "presentation_gate_reason": pres.get("presentation_gate_reason", ""),
            "immunogenicity_resolved": pres.get("immunogenicity_resolved", ""),
        })
    predictions_path = outdir / "benchmark_predictions.tsv"
    write_tsv(predictions_path, pred_rows)

    n_pos = sum(labels)
    report_path = outdir / "benchmark_report.md"
    lines = [
        "# IMPROVE / CEDAR benchmark report",
        "",
        f"- **Input**: `{Path(input_path).resolve()}`",
        f"- **Records**: {len(records)} (positive {n_pos}, negative {len(records) - n_pos})",
        f"- **Profile**: {profile_name}",
        f"- **Stub tools**: {stub}",
        f"- **NetMHCpan**: {'skipped' if skip_netmhcpan else 'run'}",
        f"- **MHCflurry**: {'skipped' if skip_mhcflurry else 'run'}",
        f"- **Immunogenicity**: {'skipped' if skip_immunogenicity else 'run'} "
        f"(sources: {', '.join(immunogenicity_config(profile)['sources'])})",
        f"- **NetMHCstabpan**: {'skipped' if skip_stabpan else 'run (IEDB API)'}",
    ]
    if cache_root is not None:
        if seeded_tools:
            lines.append(
                f"- **Tools cache**: seeded {', '.join(seeded_tools)} from `{cache_root.resolve()}`"
            )
        else:
            lines.append(
                f"- **Tools cache**: `{cache_root.resolve()}` (no new files copied; existing outputs kept)"
            )
    lines.extend([
        "",
        "## Metrics — all pairs (sorted by AUROC)",
        "",
    ])
    lines.extend(_report_metrics_table(all_rows))
    for cohort in ("gated_tesla", "binder_el2", "binder_pres05", "binder_grade_ab"):
        cohort_rows = [r for r in metric_rows if r["cohort"] == cohort]
        if not cohort_rows:
            continue
        cohort_rows.sort(key=lambda r: to_float(r.get("auroc"), -1), reverse=True)
        n_sub = cohort_rows[0]["n"]
        lines.extend([
            "",
            f"## Metrics — {COHORT_LABELS[cohort]} (n={n_sub})",
            "",
        ])
        lines.extend(_report_metrics_table(cohort_rows))
    lines.extend([
        "",
        f"## HLA-stratified AUROC — `immunogenicity_composite_score` "
        f"(alleles with n≥{HLA_MIN_N}, pos≥{HLA_MIN_POS}, neg≥{HLA_MIN_NEG})",
        "",
    ])
    if hla_metric_rows:
        lines.extend([
            "| HLA | AUROC | AUPRC | n | pos |",
            "|-----|-------|-------|---|-----|",
        ])
        for row in hla_metric_rows:
            lines.append(
                f"| {row['hla_allele']} | {row['auroc']} | {row['auprc']} | "
                f"{row['n']} | {row['n_positive']} |"
            )
    else:
        lines.append("_No alleles met minimum size thresholds._")
    lines.extend([
        "",
        "## Outputs",
        "",
        f"- `{metrics_path.name}` (includes `cohort` column)",
        f"- `{predictions_path.name}`",
        f"- `{pres_path.relative_to(outdir)}`",
    ])
    for key, path in immuno_paths.items():
        if path:
            lines.append(f"- `{path.relative_to(outdir)}`")
    if stabpan_path:
        lines.append(f"- `{stabpan_path.relative_to(outdir)}`")
    if hla_metric_rows:
        lines.append("- `benchmark_metrics_by_hla.tsv`")
    lines.extend([
        "",
        "## Notes",
        "",
        "- Labels are experimental immunogenicity (`response` 0/1), not binding alone.",
        "- Immunogenicity uses configured composite sources (`[immunogenicity]` in profile); IEDB v3.0 is optional fallback only.",
        "- `efficacy_score` applies v0.3 weights, presentation gates (`[gates]` in profile), and configured immunogenicity composite.",
        "- `gated_tesla` cohort applies grade A/B + EL≤2% + stabpan≥1.4h (TPM/VAF gates disabled when set to 0).",
        "- NetMHCstabpan predicts pMHC half-life/stability via IEDB API (`method=netmhcstabpan`); checkpoint resume supported.",
        "- Binder subsets filter to peptides more likely presented; immunogenicity AUROC often improves in these cohorts.",
        "- Full CEDAR benchmark (~2436 pairs) with IEDB NetMHCpan + stabpan takes several hours; use `--limit` or `--skip-stabpan` for smoke tests.",
        "- Local NetMHCpan binary is bypassed; predictions use IEDB API (same as pVACtools).",
    ])
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    summary = {
        "input": str(Path(input_path).resolve()),
        "outdir": str(outdir.resolve()),
        "n_records": len(records),
        "n_positive": n_pos,
        "metrics_tsv": str(metrics_path),
        "predictions_tsv": str(predictions_path),
        "report_md": str(report_path),
        "presentation_tsv": str(pres_path),
        "prime_tsv": str(immuno_paths["prime"]) if immuno_paths.get("prime") else "",
        "bigmhc_im_tsv": str(immuno_paths["bigmhc_im"]) if immuno_paths.get("bigmhc_im") else "",
        "iedb_immunogenicity_tsv": str(immuno_paths["iedb"]) if immuno_paths.get("iedb") else "",
        "stabpan_tsv": str(stabpan_path) if stabpan_path else "",
        "tools_cache_dir": str(cache_root.resolve()) if cache_root else "",
        "tools_cache_seeded": seeded_tools,
        "hla_metrics_tsv": str(outdir / "benchmark_metrics_by_hla.tsv") if hla_metric_rows else "",
        "top_auroc": all_rows[0] if all_rows else {},
        "cohorts": {c: sum(m) for c, m in cohort_masks.items()},
    }
    write_json(outdir / "benchmark_summary.json", summary)
    return summary
