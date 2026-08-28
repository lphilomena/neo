from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path
from .common import count_by, ensure_dir, markdown_table, read_tsv, safe_float, write_tsv


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Review CCF, clonality, purity, and CCF modifier evidence")
    ap.add_argument("--ranked-peptides")
    ap.add_argument("--ccf")
    ap.add_argument("--purity-table")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args(argv)
    outdir = ensure_dir(args.outdir)
    rows: list[dict[str, str]] = []
    if args.ranked_peptides and Path(args.ranked_peptides).exists():
        _, rows = read_tsv(args.ranked_peptides)
    elif args.ccf and Path(args.ccf).exists():
        _, rows = read_tsv(args.ccf)

    clonality_counts = count_by(rows, "clonality_status") if rows else Counter()
    confidence_counts = count_by(rows, "ccf_confidence") if rows else Counter()
    method_counts = count_by(rows, "ccf_method") if rows else Counter()
    low_conf = []
    for r in rows[:5000]:
        conf = (r.get("ccf_confidence", "") or "").lower()
        if conf in {"low", "unresolved", "not_mapped", "insufficient"}:
            low_conf.append({"peptide_id": r.get("peptide_id", ""), "gene": r.get("gene", ""), "peptide": r.get("peptide", ""), "hla_allele": r.get("hla_allele", ""), "clonality_status": r.get("clonality_status", ""), "ccf_estimate": r.get("ccf_estimate", r.get("ccf_best", "")), "ccf_confidence": r.get("ccf_confidence", ""), "ccf_method": r.get("ccf_method", ""), "final_priority": r.get("final_priority", "")})
    write_tsv(outdir / "ccf_confidence_flags.tsv", low_conf, ["peptide_id", "gene", "peptide", "hla_allele", "clonality_status", "ccf_estimate", "ccf_confidence", "ccf_method", "final_priority"])

    priority_mod = []
    mapping = {
        "clonal_like": ("1.00", "No clonality down-weight; do not upgrade alone if confidence is low"),
        "subclonal_like": ("0.75", "Moderate down-weight; retain but lower than clonal-like candidates"),
        "low_frequency_subclonal": ("0.45", "Strong down-weight; usually not first-line therapeutic candidate"),
        "unresolved": ("0.60", "Conservative down-weight/review; do not infer clonality"),
    }
    for status, (mult, note) in mapping.items():
        priority_mod.append({"clonality_status": status, "typical_ccf_multiplier": mult, "interpretation": note})
    write_tsv(outdir / "ccf_modifier_summary.tsv", priority_mod)

    purity_rows = []
    if args.purity_table and Path(args.purity_table).exists():
        _, purity_rows = read_tsv(args.purity_table)
    md = ["# CCF / clonality review", "", "CCF is an estimated fraction of cancer cells carrying an event. It is not VAF and should be interpreted with purity and copy-number confidence.", "", "## Clonality distribution"]
    for k, v in clonality_counts.items():
        md.append(f"- {k}: {v}")
    md.append("\n## CCF confidence distribution")
    for k, v in confidence_counts.items():
        md.append(f"- {k}: {v}")
    md.append("\n## CCF method distribution")
    for k, v in method_counts.items():
        md.append(f"- {k}: {v}")
    if purity_rows:
        md.append("\n## Purity / ploidy input")
        md.append(markdown_table(purity_rows, max_rows=20))
    md.append("\n## Suggested CCF modifier interpretation")
    md.append(markdown_table(priority_mod))
    if low_conf:
        md.append(f"\nLow-confidence CCF flags were found for {len(low_conf)} rows. In low-purity or copy-number-unstable samples, CCF should be a soft ranking modifier, not a hard clinical conclusion.")
    (outdir / "ccf_clonality_review.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
