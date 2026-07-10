from __future__ import annotations

import argparse
from pathlib import Path
from .common import ensure_dir, find_files, read_tsv, write_json, write_tsv

FEATURE_PATTERNS = {
    "has_raw_events": ["raw_events*.tsv"],
    "has_raw_peptides": ["raw_peptides*.tsv"],
    "has_presentation": ["presentation_evidence*.tsv", "*presentation*.tsv"],
    "has_ranked_peptides_recommendation": ["ranked_peptides.recommendation*.tsv"],
    "has_ranked_peptides_netmhcpan42": ["ranked_peptides.netmhcpan42*.tsv"],
    "has_ranked_peptides": ["ranked_peptides*.tsv"],
    "has_evidence_report": ["evidence_report*.html"],
    "has_hla": ["*hla*.tsv", "*HLA*.tsv", "hla*.txt"],
    "has_hla_loh": ["*hla_loh*.tsv", "*lohhla*.tsv", "*spechla*.tsv"],
    "has_expression": ["*expression*.tsv", "*TPM*.tsv", "*tpm*.tsv"],
    "has_rna_alt": ["*rna_alt*.tsv", "*rna_vaf*.tsv"],
    "has_rna_junction": ["*junction*.tsv", "*fusion*.tsv"],
    "has_purity": ["*purity*.tsv", "*facets*.tsv", "*purple*.tsv"],
    "has_cnv": ["*cnv*.tsv", "*segments*.tsv", "*facets*.tsv", "*purple*.tsv"],
    "has_appm": ["appm_summary*.tsv", "appm_gene_status*.tsv", "appm_submodule_scores*.tsv"],
    "has_ccf": ["ccf*.tsv", "ccf_lite*.tsv", "ccf_2*.tsv"],
    "has_somatic_vcf": ["*.vcf", "*.vcf.gz"],
    "has_sv_vcf": ["*sv*.vcf", "*SV*.vcf", "*manta*.vcf.gz", "*gridss*.vcf.gz"],
}

MANIFEST_COLUMNS = {
    "has_hla": ["hla", "hla_txt", "hla_tsv", "hla_file"],
    "has_somatic_vcf": ["somatic_vcf", "vcf", "snv_vcf"],
    "has_sv_vcf": ["sv_vcf", "manta_vcf", "gridss_vcf"],
    "has_expression": ["expression_tsv", "expression", "gene_expression", "tpm"],
    "has_rna_junction": ["rna_junctions_tsv", "fusion_tsv", "splice_junction_tsv"],
    "has_hla_loh": ["hla_loh_tsv", "lohhla", "spechla"],
    "has_purity": ["purity_tsv", "facets", "purple"],
    "has_cnv": ["cnv_segments", "cnv_tsv"],
    "has_ranked_peptides": ["ranked_peptides", "ranked_peptides_tsv"],
    "has_evidence_report": ["evidence_report", "evidence_report_html"],
}


def inspect_manifest(path: Path) -> dict[str, object]:
    header, rows = read_tsv(path, limit=1000)
    cols = {c.lower(): c for c in header}
    features: dict[str, bool] = {}
    found: dict[str, list[str]] = {}
    for feat, candidates in MANIFEST_COLUMNS.items():
        present_cols = [cols[c.lower()] for c in candidates if c.lower() in cols]
        vals: list[str] = []
        for col in present_cols:
            for row in rows:
                v = str(row.get(col, "") or "").strip()
                if v:
                    vals.append(v)
        features[feat] = bool(vals)
        found[feat] = vals[:10]
    sample_ids = sorted({str(r.get("sample_id", "") or "").strip() for r in rows if str(r.get("sample_id", "") or "").strip()})
    return {"features": features, "found": found, "sample_ids": sample_ids, "rows": len(rows), "columns": header}


def inspect_result_dir(path: Path) -> dict[str, object]:
    features: dict[str, bool] = {}
    found: dict[str, list[str]] = {}
    for feat, pats in FEATURE_PATTERNS.items():
        files = find_files(path, pats)
        features[feat] = bool(files)
        found[feat] = [str(p.relative_to(path)) for p in files[:20]]
    return {"features": features, "found": found}


def recommend_workflow(features: dict[str, bool]) -> str:
    if features.get("has_ranked_peptides_recommendation") and features.get("has_ranked_peptides_netmhcpan42"):
        return "ranking_compare"
    if features.get("has_ranked_peptides") or features.get("has_ranked_peptides_recommendation"):
        return "result_review"
    if features.get("has_raw_events") and features.get("has_raw_peptides") and features.get("has_presentation"):
        return "evidence_scoring"
    if features.get("has_sv_vcf"):
        return "dna_sv_workflow"
    if features.get("has_somatic_vcf") and features.get("has_hla"):
        return "snv_indel_workflow"
    if features.get("has_rna_junction") and features.get("has_hla"):
        return "fusion_splice_review"
    if not features.get("has_hla"):
        return "hla_typing_required"
    return "input_incomplete"


def build_missing(features: dict[str, bool], workflow: str) -> list[dict[str, str]]:
    missing: list[dict[str, str]] = []
    def add(key: str, severity: str, reason: str):
        if not features.get(key):
            missing.append({"input": key, "severity": severity, "reason": reason})
    if workflow in {"snv_indel_workflow", "evidence_scoring", "result_review", "ranking_compare"}:
        add("has_hla", "required", "HLA alleles are required for peptide-HLA presentation prediction")
    if workflow == "ranking_compare":
        add("has_ranked_peptides_recommendation", "required", "recommendation ranked peptide table is required for ranking comparison")
        add("has_ranked_peptides_netmhcpan42", "required", "NetMHCpan 4.2 ranked peptide table is required for ranking comparison")
        add("has_evidence_report", "recommended", "evidence report improves interpretation of ranking shifts")
    elif workflow == "result_review":
        add("has_ranked_peptides", "required", "ranked peptides table is required for result review")
        add("has_evidence_report", "recommended", "evidence report improves APPM/escape/report interpretation")
        add("has_expression", "recommended", "missing RNA expression cannot be interpreted as not expressed")
        add("has_hla_loh", "recommended", "missing HLA LOH means immune escape is incompletely assessed")
        add("has_purity", "recommended", "missing purity reduces CCF interpretability")
    elif workflow == "evidence_scoring":
        add("has_raw_events", "required", "raw_events.tsv is needed")
        add("has_raw_peptides", "required", "raw_peptides.tsv is needed")
        add("has_presentation", "required", "presentation evidence is needed")
    elif workflow == "snv_indel_workflow":
        add("has_somatic_vcf", "required", "somatic VCF is needed")
        add("has_expression", "recommended", "RNA expression improves candidate prioritization")
    elif workflow == "dna_sv_workflow":
        add("has_sv_vcf", "required", "SV VCF is needed")
        add("has_rna_junction", "recommended", "RNA junction evidence upgrades DNA SV candidates")
    return missing


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Project B P0 skill: input QC and workflow recommendation")
    ap.add_argument("--manifest")
    ap.add_argument("--result-dir")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args(argv)
    outdir = ensure_dir(args.outdir)
    features: dict[str, bool] = {}
    found: dict[str, list[str]] = {}
    sample_ids: list[str] = []
    sources: list[str] = []
    if args.result_dir:
        rd = inspect_result_dir(Path(args.result_dir))
        features.update(rd["features"])
        found.update(rd["found"])
        sources.append(f"result_dir={args.result_dir}")
    if args.manifest:
        md = inspect_manifest(Path(args.manifest))
        for k, v in md["features"].items():
            features[k] = bool(features.get(k)) or bool(v)
        found.update({f"manifest::{k}": v for k, v in md["found"].items()})
        sample_ids = md.get("sample_ids", [])  # type: ignore[assignment]
        sources.append(f"manifest={args.manifest}")
    workflow = recommend_workflow(features)
    missing = build_missing(features, workflow)
    report_rows = [{"feature": k, "present": str(v), "examples": ";".join(found.get(k, []))} for k, v in sorted(features.items())]
    write_tsv(outdir / "input_qc_report.tsv", report_rows, ["feature", "present", "examples"])
    write_tsv(outdir / "missing_inputs.tsv", missing, ["input", "severity", "reason"])
    status = {
        "status": "PASS" if workflow != "input_incomplete" else "INCOMPLETE",
        "sources": sources,
        "sample_ids": sample_ids,
        "features": features,
        "recommended_workflow": workflow,
        "missing_inputs": missing,
    }
    write_json(outdir / "input_status.json", status)
    md = ["# NeoAg Input QC", "", f"Recommended workflow: **{workflow}**", ""]
    if missing:
        md.append("## Missing / incomplete inputs")
        for m in missing:
            md.append(f"- `{m['input']}` [{m['severity']}]: {m['reason']}")
    else:
        md.append("No required missing input detected for the recommended workflow.")
    (outdir / "recommended_workflow.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
