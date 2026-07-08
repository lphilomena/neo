from __future__ import annotations

import argparse
import math
from pathlib import Path
from .common import count_by, ensure_dir, read_tsv, write_tsv, markdown_table

TOP_NS = [10, 20, 50, 100, 200, 500, 1000]


def spearman_from_ranks(rank_a: dict[str, int], rank_b: dict[str, int]) -> float | None:
    keys = [k for k in rank_a if k in rank_b]
    n = len(keys)
    if n < 2:
        return None
    xs = [rank_a[k] for k in keys]
    ys = [rank_b[k] for k in keys]
    mx = sum(xs) / n
    my = sum(ys) / n
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx <= 0 or vy <= 0:
        return None
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / math.sqrt(vx * vy)


def rank_rows(rows: list[dict[str, str]]) -> dict[str, int]:
    return {r.get("peptide_id", f"row_{i}"): i + 1 for i, r in enumerate(rows)}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Compare NetMHCpan42 ranking with recommendation ranking")
    ap.add_argument("--netmhcpan42", required=True)
    ap.add_argument("--recommendation", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args(argv)
    outdir = ensure_dir(args.outdir)
    _, net_rows = read_tsv(args.netmhcpan42)
    _, rec_rows = read_tsv(args.recommendation)
    net_rank = rank_rows(net_rows)
    rec_rank = rank_rows(rec_rows)
    common = set(net_rank) & set(rec_rank)
    corr = spearman_from_ranks(net_rank, rec_rank)
    overlap_rows = []
    for n in TOP_NS:
        a = {r.get("peptide_id", "") for r in net_rows[:n]}
        b = {r.get("peptide_id", "") for r in rec_rows[:n]}
        inter = a & b
        overlap_rows.append({"top_n": n, "overlap": len(inter), "fraction": f"{len(inter)/n:.4f}"})
    write_tsv(outdir / "topn_overlap.tsv", overlap_rows, ["top_n", "overlap", "fraction"])

    shift_rows = []
    rec_by_id = {r.get("peptide_id", ""): r for r in rec_rows}
    net_by_id = {r.get("peptide_id", ""): r for r in net_rows}
    for pid in common:
        d = abs(net_rank[pid] - rec_rank[pid])
        r = rec_by_id.get(pid) or net_by_id.get(pid) or {}
        shift_rows.append({
            "peptide_id": pid,
            "netmhcpan42_rank": net_rank[pid],
            "recommendation_rank": rec_rank[pid],
            "abs_rank_shift": d,
            "gene": r.get("gene", ""),
            "peptide": r.get("peptide", ""),
            "hla_allele": r.get("hla_allele", ""),
            "event_type": r.get("event_type", ""),
            "final_priority": r.get("final_priority", ""),
        })
    shift_rows.sort(key=lambda x: int(x["abs_rank_shift"]), reverse=True)
    write_tsv(outdir / "rank_shift.tsv", shift_rows, ["peptide_id", "netmhcpan42_rank", "recommendation_rank", "abs_rank_shift", "gene", "peptide", "hla_allele", "event_type", "final_priority"])

    shared20 = []
    net_top20 = {r.get("peptide_id", "") for r in net_rows[:20]}
    rec_top20 = {r.get("peptide_id", "") for r in rec_rows[:20]}
    for pid in sorted(net_top20 & rec_top20, key=lambda p: rec_rank.get(p, 999999)):
        r = rec_by_id[pid]
        shared20.append({"recommendation_rank": rec_rank[pid], "netmhcpan42_rank": net_rank[pid], "gene": r.get("gene", ""), "peptide": r.get("peptide", ""), "hla_allele": r.get("hla_allele", ""), "event_type": r.get("event_type", ""), "final_priority": r.get("final_priority", "")})
    write_tsv(outdir / "shared_top20_candidates.tsv", shared20)

    conflict = []
    for r in net_rows[:50]:
        pid = r.get("peptide_id", "")
        rr = rec_by_id.get(pid, r)
        if rec_rank.get(pid, 10**9) > 50 or rr.get("final_priority", "") == "D":
            conflict.append({"netmhcpan42_rank": net_rank.get(pid, ""), "recommendation_rank": rec_rank.get(pid, ""), "gene": rr.get("gene", ""), "peptide": rr.get("peptide", ""), "hla_allele": rr.get("hla_allele", ""), "event_type": rr.get("event_type", ""), "final_priority": rr.get("final_priority", ""), "recommended_use": rr.get("recommended_use", "")})
    write_tsv(outdir / "netmhc_high_reco_low.tsv", conflict)

    p_counts = count_by(rec_rows, "final_priority")
    event_counts_net20 = count_by(net_rows[:20], "event_type")
    event_counts_rec20 = count_by(rec_rows[:20], "event_type")
    md = ["# Ranking comparison", "", f"Rows: NetMHCpan42={len(net_rows)}, Recommendation={len(rec_rows)}", f"Common peptide_id: {len(common)}", f"Spearman rank correlation: {corr:.4f}" if corr is not None else "Spearman rank correlation: NA", "", "## Top-N overlap", markdown_table(overlap_rows), "", "## Recommendation final priority distribution", ""]
    for k, v in sorted(p_counts.items()):
        md.append(f"- {k}: {v}")
    md += ["", "## Top20 event-type composition", "", "Recommendation Top20: " + ", ".join(f"{k}={v}" for k, v in event_counts_rec20.items()), "NetMHCpan42 Top20: " + ", ".join(f"{k}={v}" for k, v in event_counts_net20.items()), "", "## Shared Top20 candidates", markdown_table(shared20, max_rows=20), "", "## Interpretation", "NetMHCpan42 ranking primarily reflects peptide-HLA binding/presentation strength. Recommendation ranking integrates expression, CCF/clonality, APPM, HLA LOH/immune escape, safety, event confidence, and recommended validation route. Use recommendation as the primary experimental triage ranking, and use NetMHCpan42 as a presentation-strength reference."]
    (outdir / "ranking_compare_report.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
