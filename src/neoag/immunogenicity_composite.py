"""Composite immunogenicity from PRIME, BigMHC_IM, and optional IEDB fallback."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from .utils import clamp, norm_rank, to_float
from .evidence_provenance import ProvenanceRegistry

IMMUNOGENICITY_FIELDS = [
    "prime_score",
    "prime_rank",
    "bigmhc_im_score",
    "deepimmuno_score",
    "iedb_immunogenicity_score",
    "immunogenicity_composite_score",
    "immunogenicity_source",
]


def immunogenicity_config(profile: Mapping[str, Any]) -> dict[str, Any]:
    cfg = profile.get("immunogenicity", {})
    weights = cfg.get("weights") or {}
    if not isinstance(weights, dict):
        weights = {}
    hla_mult = cfg.get("hla_source_multipliers") or {}
    if not isinstance(hla_mult, dict):
        hla_mult = {}
    sources = cfg.get("sources") or ["prime", "bigmhc_im"]
    return {
        "enabled": bool(cfg.get("enabled", True)),
        "sources": list(sources),
        "composite": str(cfg.get("composite", "mean")).lower(),
        "weights": {str(k): float(v) for k, v in weights.items()},
        "hla_source_multipliers": hla_mult,
        "use_iedb_fallback": bool(cfg.get("use_iedb_fallback", False)),
    }


def _canonical_hla(hla_allele: str) -> str:
    from .adapters.peptide_input import normalize_hla_allele

    return normalize_hla_allele(hla_allele) if hla_allele.strip() else ""


def effective_source_weights(
    hla_allele: str,
    parts: Mapping[str, float],
    profile: Mapping[str, Any],
) -> dict[str, float]:
    """Base source weights with optional per-HLA multipliers (0 = exclude source)."""
    cfg = immunogenicity_config(profile)
    base = cfg["weights"]
    multipliers_table = cfg["hla_source_multipliers"]
    canon = _canonical_hla(hla_allele)
    hla_mult: Mapping[str, Any] = {}
    if canon and isinstance(multipliers_table.get(canon), Mapping):
        hla_mult = multipliers_table[canon]

    out: dict[str, float] = {}
    for name in parts:
        w = float(base.get(name, 1.0))
        mult = float(hla_mult.get(name, 1.0)) if hla_mult else 1.0
        out[name] = w * mult
    return out


def normalize_iedb_score(raw: Any) -> float:
    return clamp((to_float(raw, 0.0) + 0.5) / 1.5)


def normalize_prime_component(row: Mapping[str, Any]) -> float | None:
    rank_raw = str(row.get("prime_rank", "")).strip()
    score_raw = str(row.get("prime_score", "")).strip()
    if rank_raw:
        return norm_rank(rank_raw)
    if score_raw:
        return clamp(to_float(score_raw, 0.0))
    return None


def normalize_bigmhc_im_component(row: Mapping[str, Any]) -> float | None:
    raw = str(row.get("bigmhc_im_score", "")).strip()
    if not raw:
        return None
    return clamp(to_float(raw, 0.0))


def normalize_deepimmuno_component(row: Mapping[str, Any]) -> float | None:
    raw = str(row.get("deepimmuno_score", "")).strip()
    if not raw:
        return None
    return clamp(to_float(raw, 0.0))


def component_scores(row: Mapping[str, Any], profile: Mapping[str, Any]) -> dict[str, float]:
    if _is_stub_evidence(row):
        return {}
    cfg = immunogenicity_config(profile)
    if not cfg["enabled"]:
        return {}

    parts: dict[str, float] = {}
    for src in cfg["sources"]:
        if src == "prime":
            val = normalize_prime_component(row)
            if val is not None:
                parts["prime"] = val
        elif src == "bigmhc_im":
            val = normalize_bigmhc_im_component(row)
            if val is not None:
                parts["bigmhc_im"] = val
        elif src == "deepimmuno":
            val = normalize_deepimmuno_component(row)
            if val is not None:
                parts["deepimmuno"] = val
        elif src == "iedb":
            raw = str(row.get("iedb_immunogenicity_score", "")).strip()
            if raw not in {"", "0.5", "0.5000"}:
                parts["iedb"] = normalize_iedb_score(raw)

    if not parts and cfg["use_iedb_fallback"]:
        raw = str(row.get("iedb_immunogenicity_score", "")).strip()
        if raw not in {"", "0.5", "0.5000"}:
            parts["iedb"] = normalize_iedb_score(raw)
    return parts


def _is_stub_evidence(row: Mapping[str, Any]) -> bool:
    return str(row.get("evidence_status", "")).strip().lower() == "stub"


def _non_stub_evidence_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [row for row in rows if not _is_stub_evidence(row)]


def _immunogenicity_tool_available(name: str, ctx: Any) -> bool:
    from .tools.runner import check_tool

    try:
        return check_tool(name, ctx.exe(name)).available
    except Exception:
        return False


def combine_component_scores(
    parts: Mapping[str, float],
    profile: Mapping[str, Any],
    *,
    hla_allele: str = "",
) -> float:
    if not parts:
        return 0.0
    cfg = immunogenicity_config(profile)
    mode = cfg["composite"]
    weights = effective_source_weights(hla_allele, parts, profile)
    active = {name: val for name, val in parts.items() if weights.get(name, 0.0) > 0.0}
    if not active:
        return 0.0

    if mode == "max":
        return clamp(max(active.values()))

    if mode == "weighted":
        num = 0.0
        den = 0.0
        for name, val in active.items():
            w = weights[name]
            num += w * val
            den += w
        return clamp(num / den) if den > 0 else 0.0

    canon = _canonical_hla(hla_allele)
    hla_adj = canon and isinstance(cfg["hla_source_multipliers"].get(canon), Mapping)
    if hla_adj:
        num = 0.0
        den = 0.0
        for name, val in active.items():
            w = weights[name]
            num += w * val
            den += w
        return clamp(num / den) if den > 0 else 0.0

    return clamp(sum(active.values()) / len(active))


def resolve_immunogenicity_score(
    peptide: Mapping[str, Any],
    presentation: Mapping[str, Any] | None,
    profile: Mapping[str, Any],
) -> tuple[float, dict[str, str]]:
    """Return 0–1 immunogenicity score and metadata fields for output rows."""
    row = {**peptide, **(presentation or {})}
    parts = component_scores(row, profile)
    hla = str(row.get("hla_allele", "")).strip()
    composite = combine_component_scores(parts, profile, hla_allele=hla)
    meta = {
        "immunogenicity_composite_score": f"{composite:.4f}" if parts else "",
        "immunogenicity_source": "+".join(sorted(parts)) if parts else "none",
    }
    return composite, meta


def has_resolved_immunogenicity(
    peptide: Mapping[str, Any],
    presentation: Mapping[str, Any] | None,
    profile: Mapping[str, Any],
) -> bool:
    row = {**peptide, **(presentation or {})}
    if component_scores(row, profile):
        return True
    raw = str(peptide.get("immunogenicity_score", "")).strip()
    if raw not in {"", "0.5", "0.5000"}:
        val = to_float(raw, 0.5)
        return abs(val - 0.5) >= 1e-6 and val > 0.0
    return False


def enrich_row_immunogenicity(row: dict[str, str], profile: Mapping[str, Any]) -> None:
    composite, meta = resolve_immunogenicity_score(row, row, profile)
    row.update(meta)


def merge_prime(
    presentation_rows: list[dict[str, str]],
    prime_by_pair: dict[tuple[str, str], dict[str, str]],
) -> None:
    for row in presentation_rows:
        key = (row.get("peptide", ""), row.get("hla_allele", ""))
        prime = prime_by_pair.get(key, {})
        row["prime_score"] = prime.get("prime_score", "")
        row["prime_rank"] = prime.get("prime_rank", "")


def merge_bigmhc_im(
    presentation_rows: list[dict[str, str]],
    bigmhc_by_pair: dict[tuple[str, str], dict[str, str]],
) -> None:
    for row in presentation_rows:
        key = (row.get("peptide", ""), row.get("hla_allele", ""))
        hit = bigmhc_by_pair.get(key, {})
        row["bigmhc_im_score"] = hit.get("bigmhc_im_score", "")


def merge_deepimmuno(
    presentation_rows: list[dict[str, str]],
    deepimmuno_by_pair: dict[tuple[str, str], dict[str, str]],
) -> None:
    for row in presentation_rows:
        key = (row.get("peptide", ""), row.get("hla_allele", ""))
        hit = deepimmuno_by_pair.get(key, {})
        row["deepimmuno_score"] = hit.get("deepimmuno_score", "")


def apply_immunogenicity_composite(
    presentation_rows: list[dict[str, str]],
    profile: Mapping[str, Any],
) -> None:
    for row in presentation_rows:
        enrich_row_immunogenicity(row, profile)


def run_immunogenicity_predictors(
    raw_peptides: Path,
    outdir: Path,
    profile: Mapping[str, Any],
    ctx: Any,
    *,
    skip: bool = False,
    skip_sources: set[str] | None = None,
    provenance_registry: ProvenanceRegistry | None = None,
) -> dict[str, Path | None]:
    """Run configured immunogenicity tools and write presentation evidence TSVs."""
    from .adapters.bigmhc_im import parse_bigmhc_im, write_bigmhc_im_evidence
    from .adapters.deepimmuno import parse_deepimmuno, write_deepimmuno_evidence
    from .adapters.iedb_immunogenicity import predict_pairs as predict_iedb_pairs, write_immunogenicity_evidence
    from .adapters.prime import parse_prime, write_prime_evidence
    from .tools.runner import run_bigmhc_im, run_deepimmuno, run_prime

    cfg = immunogenicity_config(profile)
    skip_sources = skip_sources or set()
    registry = provenance_registry or ProvenanceRegistry()
    tools_dir = outdir / "tools"
    pres_dir = outdir / "presentation"
    pres_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path | None] = {
        "prime": None,
        "bigmhc_im": None,
        "deepimmuno": None,
        "iedb": None,
    }

    prime_tool = tools_dir / "prime.tsv"
    if "prime" in cfg["sources"] and "prime" not in skip_sources:
        prime_can_use_tool_output = skip or ctx.stub or _immunogenicity_tool_available("prime", ctx)
        if not skip and not prime_can_use_tool_output:
            registry.register_missing("prime")
        elif not skip:
            run_prime(ctx, prime_tool)
        if prime_can_use_tool_output and prime_tool.is_file():
            prime_path = pres_dir / "prime_evidence.tsv"
            write_prime_evidence(
                prime_path,
                parse_prime(prime_tool, ctx.sample_id),
                provenance=registry.register_real("prime", prime_tool) if not ctx.stub else registry.register_stub("prime"),
            )
            paths["prime"] = prime_path
        elif (pres_dir / "prime_evidence.tsv").is_file():
            paths["prime"] = pres_dir / "prime_evidence.tsv"
            registry.register_passthrough("prime", paths["prime"])
        elif ctx.stub:
            registry.register_stub("prime")
        else:
            registry.register_missing("prime")
    else:
        registry.register_not_used("prime")

    bigmhc_tool = tools_dir / "bigmhc_im.tsv"
    if "bigmhc_im" in cfg["sources"] and "bigmhc_im" not in skip_sources:
        bigmhc_can_use_tool_output = skip or ctx.stub or _immunogenicity_tool_available("bigmhc_im", ctx)
        if not skip and not bigmhc_can_use_tool_output:
            registry.register_missing("bigmhc_im")
        elif not skip:
            run_bigmhc_im(ctx, bigmhc_tool)
        if bigmhc_can_use_tool_output and bigmhc_tool.is_file():
            bigmhc_path = pres_dir / "bigmhc_im_evidence.tsv"
            write_bigmhc_im_evidence(
                bigmhc_path,
                parse_bigmhc_im(bigmhc_tool, ctx.sample_id),
                provenance=registry.register_real("bigmhc_im", bigmhc_tool) if not ctx.stub else registry.register_stub("bigmhc_im"),
            )
            paths["bigmhc_im"] = bigmhc_path
        elif (pres_dir / "bigmhc_im_evidence.tsv").is_file():
            paths["bigmhc_im"] = pres_dir / "bigmhc_im_evidence.tsv"
            registry.register_passthrough("bigmhc_im", paths["bigmhc_im"])
        elif ctx.stub:
            registry.register_stub("bigmhc_im")
        else:
            registry.register_missing("bigmhc_im")
    else:
        registry.register_not_used("bigmhc_im")

    deepimmuno_tool = tools_dir / "deepimmuno.tsv"
    if "deepimmuno" in cfg["sources"] and "deepimmuno" not in skip_sources:
        deepimmuno_can_use_tool_output = skip or ctx.stub or _immunogenicity_tool_available("deepimmuno", ctx)
        if not skip and not deepimmuno_can_use_tool_output:
            registry.register_missing("deepimmuno")
        elif not skip:
            run_deepimmuno(ctx, deepimmuno_tool)
        if deepimmuno_can_use_tool_output and deepimmuno_tool.is_file():
            deepimmuno_path = pres_dir / "deepimmuno_evidence.tsv"
            write_deepimmuno_evidence(
                deepimmuno_path,
                parse_deepimmuno(deepimmuno_tool, ctx.sample_id),
                provenance=registry.register_real("deepimmuno", deepimmuno_tool) if not ctx.stub else registry.register_stub("deepimmuno"),
            )
            paths["deepimmuno"] = deepimmuno_path
        elif (pres_dir / "deepimmuno_evidence.tsv").is_file():
            paths["deepimmuno"] = pres_dir / "deepimmuno_evidence.tsv"
            registry.register_passthrough("deepimmuno", paths["deepimmuno"])
        elif ctx.stub:
            registry.register_stub("deepimmuno")
        else:
            registry.register_missing("deepimmuno")
    else:
        registry.register_not_used("deepimmuno")

    iedb_path = pres_dir / "iedb_immunogenicity.tsv"
    if "iedb" in cfg["sources"] or cfg["use_iedb_fallback"]:
        if not skip and "iedb" in cfg["sources"]:
            from .tools.prep import unique_peptide_hla_pairs
            pairs = unique_peptide_hla_pairs(raw_peptides)
            write_immunogenicity_evidence(
                iedb_path,
                predict_iedb_pairs(pairs),
                provenance=registry.register_derived("iedb", iedb_path),
            )
        if iedb_path.is_file():
            paths["iedb"] = iedb_path
            if not registry.has("iedb"):
                registry.register_derived("iedb", iedb_path)
        elif "iedb" in cfg["sources"]:
            registry.register_not_used("iedb")
    else:
        registry.register_not_used("iedb")
    return paths


def apply_immunogenicity_evidence(
    pres_rows: list[dict[str, str]],
    immuno_paths: dict[str, Path | None],
    profile: Mapping[str, Any],
) -> None:
    from .adapters.bigmhc_im import bigmhc_by_pair
    from .adapters.deepimmuno import deepimmuno_by_pair
    from .adapters.prime import prime_by_pair
    from .utils import read_tsv

    if immuno_paths.get("prime"):
        merge_prime(pres_rows, prime_by_pair(_non_stub_evidence_rows(read_tsv(immuno_paths["prime"]))))
    if immuno_paths.get("bigmhc_im"):
        merge_bigmhc_im(pres_rows, bigmhc_by_pair(_non_stub_evidence_rows(read_tsv(immuno_paths["bigmhc_im"]))))
    if immuno_paths.get("deepimmuno"):
        merge_deepimmuno(pres_rows, deepimmuno_by_pair(_non_stub_evidence_rows(read_tsv(immuno_paths["deepimmuno"]))))
    if immuno_paths.get("iedb"):
        iedb_by_pair = {
            (r.get("peptide", ""), r.get("hla_allele", "")): r.get("iedb_immunogenicity_score", "")
            for r in read_tsv(immuno_paths["iedb"])
        }
        for row in pres_rows:
            key = (row.get("peptide", ""), row.get("hla_allele", ""))
            row["iedb_immunogenicity_score"] = iedb_by_pair.get(key, "")
    apply_immunogenicity_composite(pres_rows, profile)
