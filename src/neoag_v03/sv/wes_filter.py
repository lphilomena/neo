"""WES-specific SV confidence tiers for Phase 1.5.

WES SV evidence is exome-capture-limited.  This module deliberately caps WES-only
candidates unless RNA junction or orthogonal evidence upgrades the event.
"""

from __future__ import annotations

from .sv_filter import ConfidenceResult, score_cluster_confidence
from .sv_merge import SVCluster
from .wes_capture import cluster_capture_annotation


def classify_wes_tier(event: dict) -> str:
    """Classify an SV event row into WES Phase 1.5 tiers."""
    rna = int(float(event.get("rna_junction_reads") or 0))
    base = str(event.get("event_confidence_tier") or event.get("final_sv_confidence") or "")
    capture = str(event.get("capture_interpretability") or "UNASSESSED")
    if capture in {"", "UNASSESSED", "NOT_ASSESSED", "CAPTURE_NOT_PROVIDED"}:
        if rna >= 3:
            return "WES_Tier1"
        if base == "Tier1":
            return "WES_Tier1"
        if base == "Tier2" or rna >= 1:
            return "WES_Tier2"
        return "WES_Tier3"
    cap_ok = capture in {"HIGH", "MEDIUM"}
    cap_low = capture == "LOW"

    # Tier1 requires RNA support; WES DNA support alone is not enough to become high priority.
    if cap_ok and rna >= 3 and base in {"Tier1", "Tier2", "WES_Tier1", "WES_Tier2"}:
        return "WES_Tier1"
    if cap_ok and (base in {"Tier1", "Tier2"} or rna >= 1):
        return "WES_Tier2"
    if cap_low and (base in {"Tier1", "Tier2"} or rna >= 1):
        return "WES_Tier3"
    return "WES_UNINTERPRETABLE" if capture == "UNINTERPRETABLE" else "WES_Tier3"


def priority_cap_for_wes_tier(tier: str, *, rna_reads: int = 0) -> str:
    if tier == "WES_Tier1":
        return "B" if rna_reads >= 3 else "B_CAUTION"
    if tier == "WES_Tier2":
        return "B_CAUTION"
    if tier == "WES_Tier3":
        return "C"
    return "D"


def score_wes_confidence(
    cluster: SVCluster,
    *,
    min_tumor_sr: int = 1,
    min_tumor_pe: int = 2,
    max_normal_sr: int = 0,
    max_normal_pe: int = 1,
    max_breakpoint_ci_bp: int = 750,
    rna_junction_reads: int = 0,
    capture_interpretability: str | None = None,
) -> ConfidenceResult:
    """WES-tuned confidence scoring with lower split/pair thresholds plus capture context."""
    base = score_cluster_confidence(
        cluster,
        min_tumor_sr=min_tumor_sr,
        min_tumor_pe=min_tumor_pe,
        max_normal_sr=max_normal_sr,
        max_normal_pe=max_normal_pe,
        max_breakpoint_ci_bp=max_breakpoint_ci_bp,
        rna_junction_reads=rna_junction_reads,
    )
    cap = cluster_capture_annotation(cluster)
    wes_tier = classify_wes_tier(
        {
            "rna_junction_reads": rna_junction_reads,
            "event_confidence_tier": base.tier,
            "capture_interpretability": cap.get("capture_interpretability", "UNASSESSED"),
        }
    )
    score = base.score
    if wes_tier == "WES_Tier1":
        score = min(0.92, max(score, 0.70))
    elif wes_tier == "WES_Tier2":
        score = min(score, 0.72)
    elif wes_tier == "WES_Tier3":
        score = min(score, 0.45)
    else:
        score = min(score, 0.15)
    return ConfidenceResult(tier=wes_tier, score=score, reason=f"wes:{base.reason};capture={cap.get('capture_interpretability','UNASSESSED')}")


def passes_wes_phase1_5(
    cluster: SVCluster,
    confidence: ConfidenceResult,
    allow_tier2: bool = True,
) -> bool:
    cap = cluster_capture_annotation(cluster)
    if cap.get("capture_filter_status") == "FAIL":
        return False
    if confidence.tier == "WES_Tier1":
        return True
    if confidence.tier == "WES_Tier2" and allow_tier2:
        return True
    return False
