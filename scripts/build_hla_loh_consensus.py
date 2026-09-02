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


def read_tsv(path: Path | None) -> list[dict[str, str]]:
    if not path or not path.is_file():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def display_status(row: dict[str, str]) -> str:
    return row.get("consensus_status", "UNASSESSED")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--lohhla", help="Normalized LOHHLA hla_loh.tsv")
    ap.add_argument("--spechla", help="Normalized SpecHLA hla_loh.tsv")
    ap.add_argument("--tool-status", help="Optional HLA LOH tool status TSV")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    tool_status_rows = read_tsv(Path(args.tool_status) if args.tool_status else None)
    tool_status = {row.get("tool", ""): row for row in tool_status_rows}

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

    recommended = []
    for row in standardized:
        if row["consensus_status"] == "CONSENSUS_LOST":
            recommended.append({
                "hla_allele": row["hla_allele"],
                "loh_status": "loh",
                "method": "lohhla_spechla_consensus",
                "confidence": row["crosscheck_status"],
                "source": row["source_tools"],
            })
        elif row["crosscheck_status"] == "SINGLE_TOOL_LOH":
            recommended.append({
                "hla_allele": row["hla_allele"],
                "loh_status": "loh",
                "method": "hla_loh_single_tool",
                "confidence": "SINGLE_TOOL_LOH",
                "source": row["source_tools"],
            })
    write_tsv(outdir / "recommended_hla_loh.tsv", recommended, ["hla_allele", "loh_status", "method", "confidence", "source"])

    status_names = (
        "CONSENSUS_LOST", "CONSENSUS_RETAINED", "DISCORDANT",
        "SINGLE_TOOL_LOST", "SINGLE_TOOL_RETAINED", "UNASSESSED",
    )
    counts = {status: sum(row["consensus_status"] == status for row in standardized) for status in status_names}
    single_tool_rows = [row for row in standardized if row["crosscheck_status"].startswith("SINGLE_TOOL")]
    dual_tool_rows = [row for row in standardized if row["source_tools"] == "lohhla;spechla"]
    found_tools = sorted({
        tool for tool in ("lohhla", "spechla")
        if (tool_status.get(tool, {}).get("output_status") == "FOUND")
        or (tool == "lohhla" and args.lohhla)
        or (tool == "spechla" and args.spechla)
    })
    if len(found_tools) == 2 and dual_tool_rows:
        evidence_mode = "dual_tool_consensus"
    elif len(found_tools) == 1:
        evidence_mode = "single_tool_result"
    else:
        evidence_mode = "no_hla_loh_tool_result"
    summary = {
        "sample_id": args.sample_id,
        "evidence_mode": evidence_mode,
        "found_tools": found_tools,
        "counts": counts,
        "single_tool_result_rows": len(single_tool_rows),
        "dual_tool_consensus_rows": len(dual_tool_rows),
        "tool_status": tool_status_rows,
    }
    (outdir / "hla_loh_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    review = [
        "# HLA-I LOH consensus review",
        "",
        f"- sample: `{args.sample_id}`",
        f"- evidence_mode: `{evidence_mode}`",
        f"- found_tools: `{','.join(found_tools) if found_tools else 'none'}`",
        f"- CONSENSUS_LOST: {counts['CONSENSUS_LOST']}",
        f"- CONSENSUS_RETAINED: {counts['CONSENSUS_RETAINED']}",
        f"- DISCORDANT: {counts['DISCORDANT']}",
        f"- SINGLE_TOOL_LOST: {counts['SINGLE_TOOL_LOST']}",
        f"- SINGLE_TOOL_RETAINED: {counts['SINGLE_TOOL_RETAINED']}",
        f"- UNASSESSED: {counts['UNASSESSED']}",
        f"- SINGLE_TOOL_RESULT_ROWS: {len(single_tool_rows)}",
        "",
        "Only HLA-A/B/C are included. HLA-II calls are intentionally excluded from the MHC-I LOH consensus.",
        "When one HLA LOH tool fails or has no output, available calls are carried forward as single-tool evidence and are labeled with method hla_loh_single_tool.",
    ]
    if tool_status_rows:
        review.extend(["", "## Tool status", ""])
        for row in tool_status_rows:
            review.append(
                f"- {row.get('tool', '')}: selected={row.get('selected', '')}, "
                f"exit_status={row.get('exit_status', '')}, output_status={row.get('output_status', '')}"
            )
    (outdir / "hla_loh_review.md").write_text("\n".join(review) + "\n", encoding="utf-8")
    (outdir / ".complete").write_text("hla_loh_consensus_complete\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
