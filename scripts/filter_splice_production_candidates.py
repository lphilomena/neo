#!/usr/bin/env python3
"""Select production-ready splice peptide-HLA candidates with an audit trail."""

from __future__ import annotations

import argparse
import json
import math
import re
from collections import Counter
from pathlib import Path

from neoag.schemas import EVENT_FIELDS, PEPTIDE_FIELDS
from neoag.utils import read_tsv, write_tsv


AA_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")
HLA_RE = re.compile(r"^HLA-[A-Z0-9]+\*\d{2,3}:\d{2,3}", re.IGNORECASE)
EXACT_STATUS = "CROSS_DOMAIN_CONFIRMED_EXACT_JUNCTION"


def valid_rank(value: str, maximum: float) -> bool:
    try:
        rank = float(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(rank) and 0.0 <= rank <= maximum


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--events", required=True, type=Path)
    ap.add_argument("--peptides", required=True, type=Path)
    ap.add_argument("--consensus", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--min-length", type=int, default=8)
    ap.add_argument("--max-length", type=int, default=12)
    ap.add_argument("--max-source-binding-rank", type=float, default=2.0)
    args = ap.parse_args()

    events = read_tsv(args.events)
    peptides = read_tsv(args.peptides)
    consensus = {row.get("event_id", ""): row for row in read_tsv(args.consensus)}
    selected: list[dict[str, str]] = []
    rejected: list[dict[str, str]] = []
    reason_counts: Counter[str] = Counter()

    for row in peptides:
        event_id = str(row.get("event_id") or "").strip()
        peptide = str(row.get("peptide") or "").strip().upper()
        hla = str(row.get("hla_allele") or "").strip()
        reasons: list[str] = []
        if not event_id:
            reasons.append("MISSING_EVENT_ID")
        if not peptide:
            reasons.append("MISSING_PEPTIDE")
        elif not AA_RE.fullmatch(peptide):
            reasons.append("INVALID_PEPTIDE_SEQUENCE")
        elif not args.min_length <= len(peptide) <= args.max_length:
            reasons.append("PEPTIDE_LENGTH_OUTSIDE_PRODUCTION_RANGE")
        if not hla:
            reasons.append("MISSING_HLA")
        elif not HLA_RE.match(hla):
            reasons.append("INVALID_HLA")
        status = str(consensus.get(event_id, {}).get("status") or "")
        if status != EXACT_STATUS:
            reasons.append("NO_EXACT_REGTOOLS_JUNCTION_SUPPORT")
        source_tools = ";".join(
            str(row.get(key) or "") for key in ("source_tool", "source_tools")
        ).upper()
        if "SNAF" in source_tools and not valid_rank(
            str(row.get("binding_rank") or row.get("netmhcpan_el_rank") or ""),
            args.max_source_binding_rank,
        ):
            reasons.append("SNAF_BINDING_RANK_MISSING_OR_ABOVE_THRESHOLD")

        normalized = dict(row)
        normalized["peptide"] = peptide
        if reasons:
            unique_reasons = list(dict.fromkeys(reasons))
            for reason in unique_reasons:
                reason_counts[reason] += 1
            rejected.append({**normalized, "production_filter_reasons": ";".join(unique_reasons)})
        else:
            selected.append(normalized)

    selected_event_ids = {row["event_id"] for row in selected}
    selected_events = [row for row in events if row.get("event_id") in selected_event_ids]
    args.outdir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.outdir / "raw_events.tsv", selected_events, EVENT_FIELDS)
    write_tsv(args.outdir / "raw_peptides.tsv", selected, PEPTIDE_FIELDS)
    write_tsv(
        args.outdir / "rejected_peptides.tsv",
        rejected,
        [*PEPTIDE_FIELDS, "production_filter_reasons"],
    )
    summary = {
        "input_events": len(events),
        "input_peptides": len(peptides),
        "selected_events": len(selected_events),
        "selected_peptides": len(selected),
        "rejected_peptides": len(rejected),
        "peptide_length_range": [args.min_length, args.max_length],
        "exact_junction_status_required": EXACT_STATUS,
        "max_snaf_source_binding_rank": args.max_source_binding_rank,
        "rejection_reasons": dict(sorted(reason_counts.items())),
    }
    (args.outdir / "production_filter_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    if not selected:
        raise SystemExit("No splice peptide-HLA candidates passed the production filter")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
