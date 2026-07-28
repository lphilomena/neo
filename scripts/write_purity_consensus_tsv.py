#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from neoag.agent_skills.purity_cnv_review import collect_tool_results, consensus
from neoag.utils import write_tsv


def main() -> int:
    parser = argparse.ArgumentParser(description="Write a standardized purity consensus TSV")
    parser.add_argument("--result-dir", action="append", default=[])
    parser.add_argument("--file", action="append", default=[])
    parser.add_argument("--sample-id")
    parser.add_argument("--output", required=True)
    parser.add_argument("--details", required=True)
    args = parser.parse_args()
    paths = [Path(value) for value in args.result_dir + args.file]
    rows = collect_tool_results(paths, sample_id=args.sample_id)
    result = consensus(rows)
    purity = result.get("recommended_purity")
    output_rows = [{
        "sample_id": args.sample_id or "",
        "purity": "NA" if purity in {None, ""} else purity,
        "status": result.get("status", "NO_PURITY"),
        "range": result.get("range", ""),
        "n_tools": result.get("n_tools", 0),
        "tool_values": str(result.get("tool_values") or {}),
        "interpretation": result.get("interpretation", ""),
    }]
    write_tsv(args.output, output_rows)
    write_tsv(args.details, rows)
    if purity in {None, ""}:
        raise SystemExit("no usable purity estimate was produced")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
