from __future__ import annotations

import argparse
import math

from .common import count_by, ensure_dir, read_tsv, write_tsv, markdown_table

TOP_NS = [10, 20, 50, 100, 200, 500, 1000]


def candidate_id(row: dict[str, str], idx: int) -> str:
    explicit = (row.get("peptide_id") or "").strip()
    if explicit:
        return explicit
    parts = [
        row.get("gene", ""),
        row.get("peptide", ""),
        row.get("hla_allele", ""),
        row.get("event_type", ""),
    ]
    key = "|".join(p.strip() for p in parts if p and p.strip())
    return key or f"row_{idx + 1}"


def row_ids(rows: list[dict[str, str]]) -> list[str]:
    seen: dict[str, int] = {}
    ids: list[str] = []
    for i, row in enumerate(rows):
        base = candidate_id(row, i)
        n = seen.get(base, 0)
        seen[base] = n + 1
        ids.append(base if n == 0 else f"{base}#{n + 1}")
    return ids


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


def rank_from_ids(ids: list[str]) -> dict[str, int]:
    return {pid: i + 1 for i, pid in enumerate(ids)}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Compare NetMHCpan42 ranking with recommendation ranking")
    ap.add_argument("--netmhcpan42", required=True)
    ap.add_argument("--recommendation", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args(argv)
    outdir = ensure_dir(args.outdir)
    _, net_rows = read_tsv(args.netmhcpan42)
    _, rec_rows = read_tsv(args.recommendation)
    net_ids = row_ids(net_rows)
    rec_ids = row_ids(rec_rows)
    net_rank = rank_from_ids(net_ids)
    rec_rank = rank_from_ids(rec_ids)
    common = set(net_rank) & set(rec_rank)
    corr = spearman_from_ranks(net_rank, rec_rank)

    overlap_rows = []
    for n in TOP_NS:
        denom = min(n, max(len(net_ids), len(rec_ids), 1))
        a = set(net_ids[:n])
        b = set(rec_ids[:n])
        inter = a & b
        overlap_rows.append({"top_n": n, "overlap": len(inter), "fraction": f"{len(inter) / denom:.4f}"})
    write_tsv(outdir / "topn_overlap.tsv", overlap_rows, ["top_n", "overlap", "fraction"])

    rec_by_id = dict(zip(rec_ids, rec_rows))
    net_by_id = dict(zip(net_ids, net_rows))
    shift_rows = []
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
    net_top20 = set(net_ids[:20])
    rec_top20 = set(rec_ids[:20])
    for pid in sorted(net_top20 & rec_top20, key=lambda p: rec_rank.get(p, 999999)):
        r = rec_by_id.get(pid) or net_by_id.get(pid) or {}
        shared20.append({
            "recommendation_rank": rec_rank.get(pid, ""),
            "netmhcpan42_rank": net_rank.get(pid, ""),
            "gene": r.get("gene", ""),
            "peptide": r.get("peptide", ""),
            "hla_allele": r.get("hla_allele", ""),
            "event_type": r.get("event_type", ""),
            "final_priority": r.get("final_priority", ""),
        })
    write_tsv(outdir / "shared_top20_candidates.tsv", shared20)

    conflict = []
    for pid, r in zip(net_ids[:50], net_rows[:50]):
        rr = rec_by_id.get(pid, r)
        if rec_rank.get(pid, 10**9) > 50 or rr.get("final_priority", "") == "D":
            conflict.append({
                "candidate_id": pid,
                "netmhcpan42_rank": net_rank.get(pid, ""),
                "recommendation_rank": rec_rank.get(pid, ""),
                "gene": rr.get("gene", ""),
                "peptide": rr.get("peptide", ""),
                "hla_allele": rr.get("hla_allele", ""),
                "event_type": rr.get("event_type", ""),
                "final_priority": rr.get("final_priority", ""),
                "recommended_use": rr.get("recommended_use", ""),
            })
    write_tsv(outdir / "netmhc_high_reco_low.tsv", conflict)

    p_counts = count_by(rec_rows, "final_priority")
    event_counts_net20 = count_by(net_rows[:20], "event_type")
    event_counts_rec20 = count_by(rec_rows[:20], "event_type")
    md = [
        "# Ranking comparison",
        "",
        f"Rows: NetMHCpan42={len(net_rows)}, Recommendation={len(rec_rows)}",
        f"Common candidate IDs: {len(common)}",
        f"Spearman rank correlation: {corr:.4f}" if corr is not None else "Spearman rank correlation: NA",
        "",
        "## Top-N overlap",
        markdown_table(overlap_rows),
        "",
        "## Recommendation final priority distribution",
        "",
    ]
    for k, v in sorted(p_counts.items()):
        md.append(f"- {k}: {v}")
    md += [
        "",
        "## Top20 event-type composition",
        "",
        "Recommendation Top20: " + ", ".join(f"{k}={v}" for k, v in event_counts_rec20.items()),
        "NetMHCpan42 Top20: " + ", ".join(f"{k}={v}" for k, v in event_counts_net20.items()),
        "",
        "## Shared Top20 candidates",
        markdown_table(shared20, max_rows=20),
        "",
        "## Interpretation",
        "NetMHCpan42 ranking primarily reflects peptide-HLA binding/presentation strength. Recommendation ranking integrates expression, CCF/clonality, APPM, HLA LOH/immune escape, safety, event confidence, and recommended validation route. Use recommendation as the primary experimental triage ranking, and use NetMHCpan42 as a presentation-strength reference.",
    ]
    (outdir / "ranking_compare_report.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
