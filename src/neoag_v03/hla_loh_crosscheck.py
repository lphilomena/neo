"""Cross-check HLA LOH calls across LOHHLA, SpecHLA, and other normalized tables."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping

from .adapters.peptide_input import normalize_hla_allele
from .utils import first, read_tsv, write_tsv

HLA_LOH_CROSSCHECK_FIELDS = [
    "hla_allele", "lohhla_status", "lohhla_confidence", "spechla_status", "spechla_confidence",
    "consensus_loh_status", "crosscheck_status", "source_tools", "reason",
]

HLA_LOH_CONSENSUS_FIELDS = ["hla_allele", "loh_status", "method", "confidence", "source"]

_LOH_POSITIVE = {"loh", "loss", "lost", "yes", "y", "true", "1"}
_LOH_NEGATIVE = {"no", "no_loh", "retained", "kept", "false", "0", "n"}


def _norm_status(raw: str) -> str:
    token = str(raw or "").strip().lower()
    if token in _LOH_POSITIVE:
        return "loh"
    if token in _LOH_NEGATIVE:
        return "no"
    return token or "unassessed"


def _load_tool(path: str | Path | None) -> dict[str, dict[str, str]]:
    if not path or not Path(path).is_file():
        return {}
    out: dict[str, dict[str, str]] = {}
    for row in read_tsv(path):
        allele = normalize_hla_allele(first(row, ["hla_allele", "allele", "HLA", "LossAllele", "loss_allele"], ""))
        if not allele:
            continue
        out[allele] = {
            "status": _norm_status(first(row, ["loh_status", "status", "LOH", "loss", "Loss"], "")),
            "confidence": first(row, ["confidence", "loh_confidence", "Pval_unique", "pval", "evidence_level"], ""),
        }
    return out


def _consensus(lohhla_status: str, spechla_status: str) -> tuple[str, str, str]:
    statuses = [s for s in (lohhla_status, spechla_status) if s != "unassessed"]
    positives = [s for s in statuses if s == "loh"]
    negatives = [s for s in statuses if s == "no"]
    if len(statuses) == 2 and len(positives) == 2:
        return "loh", "CONSENSUS_LOH", "LOHHLA and SpecHLA both call LOH"
    if len(statuses) == 2 and positives and negatives:
        return "discordant", "DISCORDANT", "LOHHLA and SpecHLA disagree"
    if len(statuses) == 2 and len(negatives) == 2:
        return "no", "CONSENSUS_NO_LOH", "Both tools retain allele"
    if positives:
        return "loh", "SINGLE_TOOL_LOH", "Only one tool reports LOH"
    if negatives:
        return "no", "SINGLE_TOOL_NO_LOH", "Only one tool reports retained allele"
    return "unassessed", "UNASSESSED", "No usable HLA LOH status"


def crosscheck_hla_loh(
    *,
    lohhla_hla_loh: str | Path | None = None,
    spechla_hla_loh: str | Path | None = None,
) -> list[dict[str, str]]:
    lohhla = _load_tool(lohhla_hla_loh)
    spechla = _load_tool(spechla_hla_loh)
    alleles = sorted(set(lohhla) | set(spechla))
    rows: list[dict[str, str]] = []
    for allele in alleles:
        l = lohhla.get(allele, {"status": "unassessed", "confidence": ""})
        s = spechla.get(allele, {"status": "unassessed", "confidence": ""})
        consensus, crosscheck, reason = _consensus(l["status"], s["status"])
        tools = []
        if allele in lohhla:
            tools.append("lohhla")
        if allele in spechla:
            tools.append("spechla")
        rows.append({
            "hla_allele": allele,
            "lohhla_status": l["status"],
            "lohhla_confidence": l["confidence"],
            "spechla_status": s["status"],
            "spechla_confidence": s["confidence"],
            "consensus_loh_status": consensus,
            "crosscheck_status": crosscheck,
            "source_tools": ";".join(tools),
            "reason": reason,
        })
    return rows


def consensus_hla_loh_rows(crosscheck_rows: list[Mapping[str, str]], *, include_single_tool: bool = True) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    for row in crosscheck_rows:
        status = row.get("consensus_loh_status", "")
        cross = row.get("crosscheck_status", "")
        if status != "loh":
            continue
        if cross == "SINGLE_TOOL_LOH" and not include_single_tool:
            continue
        conf = row.get("lohhla_confidence") or row.get("spechla_confidence") or cross
        out.append({
            "hla_allele": row.get("hla_allele", ""),
            "loh_status": "loh",
            "method": "hla_loh_crosscheck",
            "confidence": conf,
            "source": row.get("source_tools", ""),
        })
    return out


def write_hla_loh_crosscheck(
    out: str | Path,
    *,
    lohhla_hla_loh: str | Path | None = None,
    spechla_hla_loh: str | Path | None = None,
    consensus_out: str | Path | None = None,
    include_single_tool: bool = True,
) -> dict[str, str]:
    rows = crosscheck_hla_loh(lohhla_hla_loh=lohhla_hla_loh, spechla_hla_loh=spechla_hla_loh)
    write_tsv(out, rows, HLA_LOH_CROSSCHECK_FIELDS)
    result = {"hla_loh_crosscheck": str(out)}
    if consensus_out:
        consensus_rows = consensus_hla_loh_rows(rows, include_single_tool=include_single_tool)
        write_tsv(consensus_out, consensus_rows, HLA_LOH_CONSENSUS_FIELDS)
        result["hla_loh_consensus"] = str(consensus_out)
    return result
