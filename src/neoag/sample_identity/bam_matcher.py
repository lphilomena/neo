from __future__ import annotations

import csv
from pathlib import Path
from typing import Any


def parse_bam_matcher_short(path: str | Path) -> dict[str, Any]:
    source = Path(path)
    lines = [line.rstrip("\n") for line in source.open(encoding="utf-8", errors="replace") if line.strip()]
    header = next((line for line in lines if line.lstrip().startswith("# BAM1")), "")
    data = next((line for line in reversed(lines) if not line.lstrip().startswith("#") and "\t" in line), "")
    if not header or not data:
        raise ValueError(f"BAM-matcher short output is not parseable: {source}")
    names = [item.strip().lstrip("#").strip() for item in header.split("\t")]
    values = [item.strip() for item in data.split("\t")]
    row = dict(zip(names, values))
    conclusion = str(row.get("Conclusion") or "").upper()
    same = int(float(row.get("Same") or 0))
    different = int(float(row.get("Different") or 0))
    total = same + different
    fraction = float(row.get("FracCommon") or 0)
    if conclusion.startswith("SAME") or conclusion.startswith("LIKELY SAME"):
        status = "MATCH"
    elif conclusion.startswith("DIFFERENT") or conclusion.startswith("LIKELY DIFFERENT"):
        status = "MISMATCH"
    else:
        status = "INSUFFICIENT_DATA"
    confidence = "high" if conclusion in {"SAME", "DIFFERENT"} and total > 100 else "moderate" if total > 20 else "low"
    return {
        "sample_identity_status": status,
        "official_conclusion": conclusion or "INCONCLUSIVE",
        "confidence": confidence,
        "fraction_common": fraction,
        "sites_compared": total,
        "same_genotypes": same,
        "different_genotypes": different,
        "depth_threshold": row.get("DP_thresh", ""),
        "bam1": row.get("BAM1", ""),
        "bam2": row.get("BAM2", ""),
        "source_file": str(source),
    }


def write_identity_tsv(parsed: dict[str, Any], output: str | Path) -> None:
    target = Path(output)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(parsed), delimiter="\t")
        writer.writeheader()
        writer.writerow(parsed)
