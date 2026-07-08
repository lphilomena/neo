from __future__ import annotations

import argparse
import re
from pathlib import Path
from .common import ensure_dir, markdown_table, read_tsv, strip_html_text, write_tsv


def parse_evidence_html(path: str | None) -> dict[str, str]:
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        return {}
    raw = p.read_text(encoding="utf-8", errors="ignore")
    txt = strip_html_text(raw)
    out: dict[str, str] = {}
    for label in ["MHC-I", "MHC-II", "IFNG response"]:
        m = re.search(label + r"\s+Score:\s*([0-9.]+)\s+Status:\s*([A-Z0-9_]+)", txt)
        if m:
            out[label + " score"] = m.group(1)
            out[label + " status"] = m.group(2)
    m = re.search(r"APPM call confidence:\s*([A-Za-z_]+)\s+score=([0-9.]+)", txt)
    if m:
        out["appm_call_confidence"] = m.group(1)
        out["appm_call_confidence_score"] = m.group(2)
    m = re.search(r"Evidence completeness:\s*([A-Za-z_]+)\s*\(([0-9.]+)\)", txt)
    if m:
        out["evidence_completeness"] = m.group(1)
        out["evidence_completeness_score"] = m.group(2)
    m = re.search(r"Overall risk:\s*([A-Z_]+);\s*Mechanisms:\s*([^;]+);\s*Context:\s*([^\s<]+)", txt)
    if m:
        out["immune_escape_risk"] = m.group(1)
        out["immune_escape_mechanisms"] = m.group(2).strip()
        out["therapy_context"] = m.group(3).strip()
    return out


def summarize_hla_loh(path: str | None) -> list[dict[str, str]]:
    if not path or not Path(path).exists():
        return []
    header, rows = read_tsv(path)
    out = []
    for r in rows:
        allele = r.get("hla_allele") or r.get("allele") or r.get("HLA") or r.get("hla") or ""
        status = r.get("loh_status") or r.get("status") or r.get("LOH") or r.get("loss") or ""
        tool = r.get("evidence_tool") or r.get("tool") or r.get("source_tool") or ""
        if allele or status:
            out.append({"hla_allele": allele, "loh_status": status, "tool": tool, "confidence": r.get("confidence", "")})
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Review HLA LOH, APPM, and immune escape evidence")
    ap.add_argument("--evidence-report")
    ap.add_argument("--appm-gene-status")
    ap.add_argument("--appm-submodule-scores")
    ap.add_argument("--hla-loh")
    ap.add_argument("--ranked-peptides")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args(argv)
    outdir = ensure_dir(args.outdir)
    html = parse_evidence_html(args.evidence_report)

    hla_rows = summarize_hla_loh(args.hla_loh)
    if hla_rows:
        write_tsv(outdir / "hla_loh_consensus.tsv", hla_rows, ["hla_allele", "loh_status", "tool", "confidence"])
    else:
        write_tsv(outdir / "hla_loh_consensus.tsv", [], ["hla_allele", "loh_status", "tool", "confidence"])

    sub_rows: list[dict[str, str]] = []
    if args.appm_submodule_scores and Path(args.appm_submodule_scores).exists():
        _, sub_rows = read_tsv(args.appm_submodule_scores)
    gene_rows: list[dict[str, str]] = []
    if args.appm_gene_status and Path(args.appm_gene_status).exists():
        _, gene_rows = read_tsv(args.appm_gene_status)
    driver_rows = []
    for r in gene_rows:
        status = " ".join(str(r.get(k, "")) for k in ["biallelic_status", "functional_status", "gene_integrity_status", "expression_status"])
        if any(x in status.lower() for x in ["loss", "defect", "low", "caution"]):
            driver_rows.append({k: r.get(k, "") for k in ["gene", "pathway", "biallelic_status", "functional_status", "copy_number_status", "loh_status", "expression_status", "gene_integrity_status", "reason"]})
    write_tsv(outdir / "appm_driver_defects_summary.tsv", driver_rows)

    affected = []
    lost_alleles = {r["hla_allele"] for r in hla_rows if str(r.get("loh_status", "")).lower() in {"loss", "lost", "yes", "loh", "true"}}
    if args.ranked_peptides and Path(args.ranked_peptides).exists() and lost_alleles:
        _, pep_rows = read_tsv(args.ranked_peptides)
        for r in pep_rows:
            if r.get("hla_allele") in lost_alleles:
                affected.append({"peptide_id": r.get("peptide_id", ""), "gene": r.get("gene", ""), "peptide": r.get("peptide", ""), "hla_allele": r.get("hla_allele", ""), "final_priority": r.get("final_priority", ""), "recommended_action": "restricting_hla_lost_review"})
    write_tsv(outdir / "affected_peptides.tsv", affected)

    md = ["# HLA LOH / APPM / immune escape review", "", "## Sample-level APPM summary"]
    if html:
        for k, v in html.items():
            md.append(f"- {k}: **{v}**")
    else:
        md.append("- Evidence report not provided; APPM summary is incomplete.")
    md.append("\n## HLA LOH")
    if hla_rows:
        md.append(markdown_table(hla_rows, max_rows=30))
    else:
        md.append("No HLA LOH table provided. Do not interpret missing HLA LOH evidence as no LOH.")
    md.append("\n## APPM submodules")
    if sub_rows:
        md.append(markdown_table(sub_rows, max_rows=20))
    else:
        md.append("Submodule score table not provided.")
    md.append("\n## Driver defects / cautions")
    if driver_rows:
        md.append(markdown_table(driver_rows, max_rows=30))
    else:
        md.append("No driver defects/cautions detected in the provided APPM gene table, or table not provided.")
    md.append("\n## Interpretation")
    md.append("APPM/HLA LOH review is mechanism evidence for antigen presentation and immune escape. A PASS or no major signal means no strong negative evidence was detected in the provided inputs; it does not prove intact presentation when evidence completeness is partial.")
    (outdir / "appm_escape_review.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
