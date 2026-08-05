#!/usr/bin/env python3
"""Build an explicit HLA-I LOH consensus with retained/conflict/unassessed states."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from neoag.hla_loh_crosscheck import crosscheck_hla_loh


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def display_status(row: dict[str, str]) -> str:
    return row.get("consensus_status", "UNASSESSED")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--lohhla", help="Normalized LOHHLA hla_loh.tsv")
    ap.add_argument("--spechla", help="Normalized SpecHLA hla_loh.tsv")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = crosscheck_hla_loh(lohhla_hla_loh=args.lohhla, spechla_hla_loh=args.spechla)
    standardized: list[dict[str, str]] = []
    for row in rows:
        allele = row.get("hla_allele", "")
        locus = allele.replace("HLA-", "").split("*", 1)[0]
        if locus not in {"A", "B", "C"}:
            continue
        standardized.append({
            "sample_id": args.sample_id,
            "hla_allele": allele,
            "locus": locus,
            "consensus_status": display_status(row),
            "lohhla_status": row.get("lohhla_status", "unassessed"),
            "spechla_status": row.get("spechla_status", "unassessed"),
            "crosscheck_status": row.get("crosscheck_status", "UNASSESSED"),
            "source_tools": row.get("source_tools", ""),
            "reason": row.get("reason", ""),
            **{key: value for key, value in row.items() if key.startswith("lohhla_") or key.startswith("spechla_")},
        })
    base_columns = ["sample_id", "hla_allele", "locus", "consensus_status", "lohhla_status", "spechla_status", "crosscheck_status", "source_tools", "reason"]
    evidence_columns = sorted({
        key for row in standardized for key in row
        if (key.startswith("lohhla_") or key.startswith("spechla_")) and key not in base_columns
    })
    columns = [*base_columns, *evidence_columns]
    write_tsv(outdir / "hla_loh.standardized.tsv", standardized, columns)
    write_tsv(outdir / "hla_loh_consensus.tsv", standardized, columns)

    recommended = [{
        "hla_allele": row["hla_allele"],
        "loh_status": "loh",
        "method": "lohhla_spechla_consensus",
        "confidence": row["crosscheck_status"],
        "source": row["source_tools"],
    } for row in standardized if row["consensus_status"] == "CONSENSUS_LOST"]
    write_tsv(outdir / "recommended_hla_loh.tsv", recommended, ["hla_allele", "loh_status", "method", "confidence", "source"])

    status_names = ("CONSENSUS_LOST", "CONSENSUS_RETAINED", "DISCORDANT", "UNASSESSED")
    counts = {status: sum(row["consensus_status"] == status for row in standardized) for status in status_names}
    (outdir / "hla_loh_summary.json").write_text(json.dumps({"sample_id": args.sample_id, "counts": counts}, indent=2) + "\n", encoding="utf-8")
    review = [
        "# HLA-I LOH consensus review",
        "",
        f"- sample: `{args.sample_id}`",
        f"- CONSENSUS_LOST: {counts['CONSENSUS_LOST']}",
        f"- CONSENSUS_RETAINED: {counts['CONSENSUS_RETAINED']}",
        f"- DISCORDANT: {counts['DISCORDANT']}",
        f"- UNASSESSED: {counts['UNASSESSED']}",
        "",
        "Only HLA-A/B/C are included. HLA-II calls are intentionally excluded from the MHC-I LOH consensus.",
        "Missing coverage or one missing tool remains UNASSESSED/SINGLE_TOOL and is never interpreted as retained.",
    ]
    (outdir / "hla_loh_review.md").write_text("\n".join(review) + "\n", encoding="utf-8")
    (outdir / ".complete").write_text("hla_loh_consensus_complete\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
