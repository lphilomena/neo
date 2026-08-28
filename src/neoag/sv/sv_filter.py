from __future__ import annotations

from dataclasses import dataclass
from .sv_merge import SVCluster
from .sv_callset import SVRecord


def clamp(x: float, lo: float = 0.0, hi: float = 1.0) -> float:
    return max(lo, min(hi, x))


@dataclass
class ConfidenceResult:
    tier: str
    score: float
    reason: str


def score_cluster_confidence(
    cluster: SVCluster,
    *,
    min_tumor_sr: int = 2,
    min_tumor_pe: int = 4,
    max_normal_sr: int = 0,
    max_normal_pe: int = 1,
    max_breakpoint_ci_bp: int = 500,
    rna_junction_reads: int = 0,
) -> ConfidenceResult:
    rep = cluster.representative
    caller_count = len(cluster.callers)
    tumor_sr = cluster.max_int("tumor_sr")
    tumor_pe = cluster.max_int("tumor_pe")
    sniffles_support = cluster.max_int("sniffles_support")
    normal_sr = cluster.max_int("normal_sr")
    normal_pe = cluster.max_int("normal_pe")
    precision = cluster.best_precision() or rep.breakpoint_precision_bp

    caller_score = clamp(caller_count / 3.0)
    lr_support_score = clamp(sniffles_support / max(min_tumor_pe * 2, 1))
    support_score = max(
        lr_support_score,
        clamp(0.6 * (tumor_sr / max(min_tumor_sr * 2, 1)) + 0.4 * (tumor_pe / max(min_tumor_pe * 2, 1))),
    )
    if precision <= 0:
        precision_score = 0.70
    elif precision <= 50:
        precision_score = 1.00
    elif precision <= 200:
        precision_score = 0.80
    elif precision <= max_breakpoint_ci_bp:
        precision_score = 0.50
    else:
        precision_score = 0.20
    normal_penalty = 0.0
    if normal_sr > max_normal_sr or normal_pe > max_normal_pe:
        normal_penalty = 0.30
    if normal_sr >= min_tumor_sr or normal_pe >= min_tumor_pe:
        normal_penalty = 0.60
    rna_bonus = 0.15 if rna_junction_reads >= 3 else (0.05 if rna_junction_reads > 0 else 0.0)
    score = clamp(0.30 * caller_score + 0.35 * support_score + 0.20 * precision_score + rna_bonus - normal_penalty)

    reasons: list[str] = []
    if caller_count >= 2:
        reasons.append("multi_caller")
    if tumor_sr >= min_tumor_sr:
        reasons.append("tumor_split_read_support")
    if tumor_pe >= min_tumor_pe:
        reasons.append("tumor_pair_support")
    if sniffles_support >= min_tumor_pe:
        reasons.append("long_read_support")
    if rna_junction_reads >= 3:
        reasons.append("rna_junction_supported")
    if normal_penalty:
        reasons.append("normal_support_penalty")
    if precision and precision > max_breakpoint_ci_bp:
        reasons.append("imprecise_breakpoint")

    assembly_singleton = caller_count == 1 and rep.caller.lower() in {"gridss", "gridss2", "svaba"} and (tumor_sr >= min_tumor_sr or rep.tumor_alt_support >= min_tumor_sr)
    if (caller_count >= 2 and score >= 0.55 and normal_penalty < 0.60) or (assembly_singleton and rna_junction_reads >= 3):
        tier = "Tier1"
    elif score >= 0.35 and normal_penalty < 0.60:
        tier = "Tier2"
    else:
        tier = "Tier3"
    return ConfidenceResult(tier=tier, score=score, reason=";".join(reasons) if reasons else "limited_evidence")


def passes_phase1(cluster: SVCluster, confidence: ConfidenceResult, allow_tier2: bool = True) -> bool:
    return confidence.tier == "Tier1" or (allow_tier2 and confidence.tier == "Tier2")
