from __future__ import annotations

from pathlib import Path
from typing import Any

from .io import ensure_dir, normalize_hla, read_fasta_sequences, read_table, row_get, safe_float, write_json, write_tsv


def run_hla_typing_loh(args: dict[str, Any]) -> dict[str, Any]:
    outdir = ensure_dir(args["outdir"])
    hla_path = Path(args.get("hla") or args.get("input") or "")
    sample_id = args.get("sample_id") or hla_path.stem
    alleles: list[dict[str, str]] = []
    if hla_path.suffix.lower() in {".txt", ".hla"}:
        text = hla_path.read_text(encoding="utf-8", errors="replace").replace(",", "\n").replace(";", "\n")
        vals = [normalize_hla(x.strip()) for x in text.splitlines() if x.strip()]
        alleles = [{"sample_id": sample_id, "hla_allele": v, "hla_gene": v.split("*")[0].replace("HLA-", "") if "*" in v else "", "source_tool": "provided", "typing_status": "provided"} for v in vals if v]
    else:
        _, rows = read_table(hla_path)
        for row in rows:
            val = normalize_hla(row_get(row, ["hla_allele", "allele", "hla", "HLA", "result"], ""))
            if not val:
                # allow columns HLA-A,HLA-B,HLA-C
                for key, v in row.items():
                    if str(key).upper().startswith("HLA") and str(v).strip():
                        vv = normalize_hla(str(v))
                        alleles.append({"sample_id": sample_id, "hla_allele": vv, "hla_gene": vv.split("*")[0].replace("HLA-", "") if "*" in vv else key, "source_tool": row_get(row, ["source_tool", "tool"], "provided"), "typing_status": "provided"})
                continue
            alleles.append({"sample_id": sample_id, "hla_allele": val, "hla_gene": val.split("*")[0].replace("HLA-", "") if "*" in val else "", "source_tool": row_get(row, ["source_tool", "tool"], "provided"), "typing_status": row_get(row, ["typing_status", "status"], "provided")})
    seen = set()
    unique = []
    for a in alleles:
        if a["hla_allele"] not in seen:
            seen.add(a["hla_allele"])
            unique.append(a)
    loh_rows = []
    if args.get("hla_loh"):
        _, rows = read_table(args["hla_loh"])
        for row in rows:
            allele = normalize_hla(row_get(row, ["hla_allele", "allele", "hla", "HLA"], ""))
            status = row_get(row, ["loh_status", "LOH", "status", "loss"], "") or "unknown"
            loh_rows.append({"sample_id": sample_id, "hla_allele": allele, "loh_status": status.lower(), "source_tool": row_get(row, ["source_tool", "tool"], "provided"), "loh_confidence": row_get(row, ["confidence", "loh_confidence"], "")})
    else:
        loh_rows = [{"sample_id": sample_id, "hla_allele": a["hla_allele"], "loh_status": "unassessed", "source_tool": "missing", "loh_confidence": "insufficient"} for a in unique]
    flags = []
    if args.get("ranked_peptides"):
        _, pep_rows = read_table(args["ranked_peptides"])
        lost = {r["hla_allele"] for r in loh_rows if r.get("loh_status") in {"yes", "lost", "loh", "true"}}
        for row in pep_rows:
            hla = normalize_hla(row_get(row, ["hla_allele", "hla", "allele"], ""))
            if hla in lost:
                flags.append({"peptide_id": row_get(row, ["peptide_id", "id"], ""), "peptide": row_get(row, ["peptide"], ""), "hla_allele": hla, "restricting_hla_lost": "true", "recommended_action": "cap_or_reject_for_vaccine_tcr"})
    write_tsv(outdir / "hla_typing.normalized.tsv", unique)
    write_tsv(outdir / "hla_typing_consensus.tsv", unique)
    write_tsv(outdir / "hla_loh_consensus.tsv", loh_rows)
    write_tsv(outdir / "restricting_hla_peptide_flags.tsv", flags, ["peptide_id", "peptide", "hla_allele", "restricting_hla_lost", "recommended_action"])
    (outdir / "hla_review.md").write_text(f"# HLA typing / LOH review\n\nNormalized HLA alleles: {len(unique)}. HLA LOH rows: {len(loh_rows)}. Lost restricting HLA peptide flags: {len(flags)}.\n\nMissing HLA LOH is reported as unassessed, not as no LOH.\n", encoding="utf-8")
    res = {"status": "PASS", "skill": "neoag-hla-typing-loh", "summary": f"Normalized {len(unique)} HLA alleles", "outputs": {"hla_typing": str(outdir / "hla_typing.normalized.tsv"), "hla_loh": str(outdir / "hla_loh_consensus.tsv")}}
    write_json(outdir / "skill_result.json", res)
    return res


def _grade_from_rank_ic50(el_rank: float | None, ba_rank: float | None, ic50: float | None) -> str:
    best_rank = min([x for x in [el_rank, ba_rank] if x is not None], default=None)
    if best_rank is not None:
        if best_rank <= 0.5:
            return "A"
        if best_rank <= 2.0:
            return "B"
        if best_rank <= 5.0:
            return "C_BINDING_ONLY"
        return "D_WEAK"
    if ic50 is not None:
        if ic50 <= 150:
            return "A"
        if ic50 <= 500:
            return "B"
        if ic50 <= 2000:
            return "C_BINDING_ONLY"
        return "D_WEAK"
    return "UNASSESSED"


def run_presentation(args: dict[str, Any]) -> dict[str, Any]:
    outdir = ensure_dir(args["outdir"])
    path = args.get("predictions") or args.get("input") or args.get("raw_peptides")
    if not path:
        res = {"status": "FAIL", "skill": "neoag-presentation", "failure_reason": "MISSING_INPUT"}
        write_json(outdir / "skill_result.json", res)
        return res
    _, rows = read_table(path)
    pres = []
    for i, row in enumerate(rows, 1):
        pep = row_get(row, ["peptide", "mt_peptide", "sequence"], "")
        hla = normalize_hla(row_get(row, ["hla_allele", "hla", "allele"], ""))
        pid = row_get(row, ["peptide_id", "id"], "") or f"PEP_{i}_{pep}_{hla}"
        ic50 = safe_float(row_get(row, ["ic50", "mt_ic50", "netmhcpan_mt_ic50", "affinity"], ""))
        el = safe_float(row_get(row, ["el_rank", "netmhcpan_mt_rank_el", "rank_el", "presentation_rank"], ""))
        ba = safe_float(row_get(row, ["ba_rank", "netmhcpan_mt_rank_ba", "rank_ba"], ""))
        grade = row_get(row, ["presentation_evidence_grade", "grade"], "") or _grade_from_rank_ic50(el, ba, ic50)
        pres.append({"peptide_id": pid, "event_id": row_get(row, ["event_id"], ""), "gene": row_get(row, ["gene"], ""), "peptide": pep, "hla_allele": hla, "ic50": ic50 if ic50 is not None else "", "el_rank": el if el is not None else "", "ba_rank": ba if ba is not None else "", "presentation_evidence_grade": grade, "source_tool": row_get(row, ["source_tool", "tool"], "provided_or_normalized")})
    write_tsv(outdir / "presentation_evidence.tsv", pres)
    counts: dict[str, int] = {}
    for r in pres:
        counts[r["presentation_evidence_grade"]] = counts.get(r["presentation_evidence_grade"], 0) + 1
    write_tsv(outdir / "presentation_summary.tsv", [{"grade": k, "count": v} for k, v in sorted(counts.items())])
    write_tsv(outdir / "presentation_qc.tsv", [{"metric": "input_rows", "value": len(rows)}, {"metric": "presentation_rows", "value": len(pres)}])
    res = {"status": "PASS", "skill": "neoag-presentation", "summary": f"Normalized {len(pres)} presentation rows", "outputs": {"presentation": str(outdir / "presentation_evidence.tsv")}}
    write_json(outdir / "skill_result.json", res)
    return res


def run_expression(args: dict[str, Any]) -> dict[str, Any]:
    outdir = ensure_dir(args["outdir"])
    expr = Path(args.get("expression_tsv") or args.get("input") or "")
    _, rows = read_table(expr)
    out = []
    for row in rows:
        gene = row_get(row, ["gene", "gene_name", "symbol", "Gene"], "")
        tpm = safe_float(row_get(row, ["tpm", "TPM", "gene_tpm", "expression", "value"], ""), 0.0) or 0.0
        status = "expressed" if tpm >= 1.0 else ("low" if tpm > 0 else "not_detected")
        out.append({"gene": gene, "gene_tpm": f"{tpm:.4f}", "expression_status": status, "source_file": str(expr)})
    write_tsv(outdir / "expression_evidence.tsv", out)
    write_tsv(outdir / "expression_qc.tsv", [{"metric": "genes", "value": len(out)}, {"metric": "expressed_ge_1_tpm", "value": sum(1 for r in out if r["expression_status"] == "expressed")}])
    res = {"status": "PASS", "skill": "neoag-expression", "summary": f"Normalized expression for {len(out)} genes", "outputs": {"expression": str(outdir / "expression_evidence.tsv")}}
    write_json(outdir / "skill_result.json", res)
    return res


def run_rna_evidence(args: dict[str, Any]) -> dict[str, Any]:
    outdir = ensure_dir(args["outdir"])
    rna = Path(args.get("rna_tsv") or args.get("input") or "")
    _, rows = read_table(rna)
    alt = []
    junction = []
    for row in rows:
        event_id = row_get(row, ["event_id", "variant_id", "junction_id"], "")
        reads = row_get(row, ["rna_alt_reads", "alt_reads", "junction_reads", "read_count", "reads"], "")
        vaf = row_get(row, ["rna_vaf", "vaf"], "")
        typ = row_get(row, ["event_type", "source_type", "type"], "")
        rec = {"event_id": event_id, "gene": row_get(row, ["gene"], ""), "rna_reads": reads, "rna_vaf": vaf, "source_file": str(rna)}
        if "junction" in typ.lower() or row_get(row, ["junction_id"], ""):
            junction.append(rec)
        else:
            alt.append(rec)
    write_tsv(outdir / "rna_alt_evidence.tsv", alt)
    write_tsv(outdir / "rna_junction_evidence.tsv", junction)
    write_tsv(outdir / "rna_evidence_qc.tsv", [{"metric": "input_rows", "value": len(rows)}, {"metric": "rna_alt_rows", "value": len(alt)}, {"metric": "junction_rows", "value": len(junction)}])
    res = {"status": "PASS", "skill": "neoag-rna-evidence", "summary": f"Normalized {len(alt)} RNA alt rows and {len(junction)} junction rows", "outputs": {"rna_alt": str(outdir / "rna_alt_evidence.tsv"), "rna_junction": str(outdir / "rna_junction_evidence.tsv")}}
    write_json(outdir / "skill_result.json", res)
    return res


def run_ccf(args: dict[str, Any]) -> dict[str, Any]:
    outdir = ensure_dir(args["outdir"])
    table = Path(args.get("event_table_or_ranked_peptides") or args.get("input") or "")
    purity = safe_float(args.get("purity"), None)
    if purity is None and args.get("purity_tsv"):
        _, pr = read_table(args["purity_tsv"])
        if pr:
            purity = safe_float(row_get(pr[0], ["purity", "tumor_purity", "facets_purity", "purple_purity"], ""), None)
    _, rows = read_table(table)
    out = []
    for row in rows:
        vaf = safe_float(row_get(row, ["tumor_vaf", "vaf", "variant_allele_fraction"], ""), None)
        best = safe_float(row_get(row, ["ccf_best", "ccf_estimate", "ccf"], ""), None)
        method = row_get(row, ["ccf_method"], "")
        if best is None and vaf is not None and purity and purity > 0:
            best = min(1.0, max(0.0, vaf / max(purity, 0.001) * 2.0))
            method = "VAF_PURITY_APPROX"
        if best is None:
            status = "unresolved"
            mult = 0.60
            conf = "low"
        elif best >= 0.8:
            status = "clonal_like"; mult = 1.00; conf = "medium" if purity and purity >= 0.3 else "low"
        elif best >= 0.35:
            status = "subclonal_like"; mult = 0.75; conf = "medium" if purity and purity >= 0.3 else "low"
        else:
            status = "low_frequency_subclonal"; mult = 0.45; conf = "low"
        out.append({"event_id": row_get(row, ["event_id"], ""), "peptide_id": row_get(row, ["peptide_id"], ""), "gene": row_get(row, ["gene"], ""), "ccf_best": f"{best:.4f}" if best is not None else "", "clonality_status": status, "ccf_confidence": conf, "ccf_multiplier": f"{mult:.2f}", "ccf_method": method or "UNRESOLVED"})
    write_tsv(outdir / "ccf_lite.tsv", out)
    write_tsv(outdir / "ccf_modifier_summary.tsv", out)
    write_tsv(outdir / "ccf_input_qc.tsv", [{"metric": "input_rows", "value": len(rows)}, {"metric": "purity", "value": purity if purity is not None else "missing"}])
    res = {"status": "PASS", "skill": "neoag-ccf", "summary": f"Generated CCF/clonality rows: {len(out)}", "outputs": {"ccf_lite": str(outdir / "ccf_lite.tsv")}}
    write_json(outdir / "skill_result.json", res)
    return res


def run_appm_escape(args: dict[str, Any]) -> dict[str, Any]:
    outdir = ensure_dir(args["outdir"])
    src = args.get("gene_status_or_appm") or args.get("input")
    rows = []
    if src:
        _, rows = read_table(src)
    driver = []
    risk_genes = {"B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2", "TAPBP", "JAK1", "JAK2", "STAT1", "IFNGR1", "IFNGR2", "NLRC5", "CIITA"}
    for row in rows:
        gene = row_get(row, ["gene", "symbol"], "")
        status_text = " ".join(str(v) for v in row.values()).lower()
        if gene in risk_genes and any(x in status_text for x in ["loss", "low", "defect", "loh", "frameshift", "stop"]):
            driver.append({"gene": gene, "gene_integrity_status": row_get(row, ["gene_integrity_status", "status", "functional_status"], "caution"), "reason": row_get(row, ["reason"], status_text[:120])})
    mhc_i_status = "MHC_I_INTACT" if not any(d["gene"] == "B2M" and "biallelic" in d["reason"].lower() for d in driver) else "MHC_I_DEFECTIVE"
    summary = [{"module": "MHC-I", "status": mhc_i_status, "score": "1.0000" if mhc_i_status.endswith("INTACT") else "0.0000", "evidence_completeness": "PARTIAL" if rows else "UNASSESSED", "confidence": "medium" if rows else "insufficient"}]
    mods = []
    flags = []
    if args.get("ranked_peptides"):
        _, peps = read_table(args["ranked_peptides"])
        for row in peps:
            mods.append({"peptide_id": row_get(row, ["peptide_id"], ""), "appm_action": "PASS" if mhc_i_status.endswith("INTACT") else "CAP_D", "appm_multiplier": "1.0" if mhc_i_status.endswith("INTACT") else "0.2", "appm_reason": mhc_i_status})
            flags.append({"peptide_id": row_get(row, ["peptide_id"], ""), "escape_status": "ESCAPE_PASS" if mhc_i_status.endswith("INTACT") else "ESCAPE_RISK", "escape_reason": mhc_i_status})
    write_tsv(outdir / "appm_summary.tsv", summary)
    write_tsv(outdir / "appm_gene_status.tsv", driver, ["gene", "gene_integrity_status", "reason"])
    write_tsv(outdir / "appm_peptide_modifiers.tsv", mods)
    write_tsv(outdir / "immune_escape_summary.tsv", [{"overall_risk": "LOW" if mhc_i_status.endswith("INTACT") else "HIGH", "mechanisms": ";".join(d["gene"] for d in driver) or "no_major_signal"}])
    write_tsv(outdir / "peptide_escape_flags.tsv", flags)
    res = {"status": "PASS", "skill": "neoag-appm-escape", "summary": f"APPM status {mhc_i_status}; driver defects {len(driver)}", "outputs": {"appm_summary": str(outdir / "appm_summary.tsv")}}
    write_json(outdir / "skill_result.json", res)
    return res


def run_safety(args: dict[str, Any]) -> dict[str, Any]:
    outdir = ensure_dir(args["outdir"])
    pep_path = Path(args.get("raw_peptides_or_ranked_peptides") or args.get("input") or "")
    _, rows = read_table(pep_path)
    proteome = {}
    if args.get("normal_proteome"):
        proteome = read_fasta_sequences(args["normal_proteome"])
    all_seq = "\n".join(proteome.values()) if proteome else ""
    out = []
    for row in rows:
        pep = row_get(row, ["peptide"], "").upper()
        exact = bool(pep and all_seq and pep in all_seq)
        wt = row_get(row, ["wildtype_peptide", "wt_peptide"], "")
        tier = "SAFETY_FAIL" if exact else ("SAFETY_REVIEW" if wt and wt == pep else "SAFETY_PASS")
        out.append({"peptide_id": row_get(row, ["peptide_id"], ""), "event_id": row_get(row, ["event_id"], ""), "peptide": pep, "normal_proteome_exact_match": str(exact).lower(), "safety_tier": tier, "safety_status": "FAIL" if exact else ("REVIEW" if tier == "SAFETY_REVIEW" else "PASS"), "safety_reason": "normal_proteome_exact_match" if exact else ("wt_peptide_review" if tier == "SAFETY_REVIEW" else "no_major_signal")})
    write_tsv(outdir / "peptide_safety.tsv", out)
    write_tsv(outdir / "event_safety.tsv", out)
    write_tsv(outdir / "safety_review.tsv", out)
    res = {"status": "PASS", "skill": "neoag-safety", "summary": f"Safety rows {len(out)}; exact matches {sum(1 for r in out if r['normal_proteome_exact_match']=='true')}", "outputs": {"peptide_safety": str(outdir / "peptide_safety.tsv")}}
    write_json(outdir / "skill_result.json", res)
    return res


def run_ranking(args: dict[str, Any]) -> dict[str, Any]:
    outdir = ensure_dir(args["outdir"])
    raw = Path(args.get("raw_peptides") or args.get("input") or "")
    _, peps = read_table(raw)
    pres_by = {}
    if args.get("presentation"):
        _, pres = read_table(args["presentation"])
        for r in pres:
            pres_by[row_get(r, ["peptide_id"], "") or (row_get(r, ["peptide"], "") + "|" + row_get(r, ["hla_allele"], ""))] = r
    expr_by = {}
    if args.get("expression"):
        _, expr = read_table(args["expression"])
        for r in expr:
            expr_by[row_get(r, ["gene"], "")] = r
    ccf_by = {}
    if args.get("ccf"):
        _, ccf = read_table(args["ccf"])
        for r in ccf:
            ccf_by[row_get(r, ["event_id"], "") or row_get(r, ["peptide_id"], "")] = r
    safety_by = {}
    if args.get("safety"):
        _, safety = read_table(args["safety"])
        for r in safety:
            safety_by[row_get(r, ["peptide_id"], "")] = r
    ranked = []
    for row in peps:
        pid = row_get(row, ["peptide_id"], "") or (row_get(row, ["peptide"], "") + "|" + row_get(row, ["hla_allele"], ""))
        p = pres_by.get(pid) or pres_by.get(row_get(row, ["peptide"], "") + "|" + row_get(row, ["hla_allele"], "")) or {}
        grade = row_get(p, ["presentation_evidence_grade"], "UNASSESSED")
        base = {"A": 0.9, "B": 0.75, "C_BINDING_ONLY": 0.45, "D_WEAK": 0.15, "UNASSESSED": 0.2}.get(grade, 0.3)
        gene = row_get(row, ["gene"], "")
        expr = expr_by.get(gene, {})
        expr_mult = 1.0 if row_get(expr, ["expression_status"], "") == "expressed" else (0.8 if row_get(expr, ["expression_status"], "") == "low" else 0.7)
        c = ccf_by.get(row_get(row, ["event_id"], ""), {}) or ccf_by.get(pid, {})
        ccf_mult = safe_float(row_get(c, ["ccf_multiplier"], ""), 0.6) or 0.6
        s = safety_by.get(pid, {})
        safety_status = row_get(s, ["safety_status"], "PASS") or "PASS"
        safety_mult = 0.0 if safety_status == "FAIL" else (0.75 if safety_status == "REVIEW" else 1.0)
        score = base * expr_mult * ccf_mult * safety_mult
        priority = "A" if score >= 0.75 else ("B" if score >= 0.55 else ("C_CAUTION" if safety_status == "REVIEW" else ("C" if score >= 0.30 else "D")))
        ranked.append({**row, "peptide_id": pid, "presentation_evidence_grade": grade, "expression_status": row_get(expr, ["expression_status"], "unassessed"), "ccf_status": row_get(c, ["clonality_status"], "unresolved"), "ccf_multiplier": f"{ccf_mult:.2f}", "safety_status": safety_status, "final_priority": priority, "efficacy_score": f"{score:.4f}", "recommended_use": "computed_triage; validate experimentally"})
    ranked.sort(key=lambda r: float(r["efficacy_score"]), reverse=True)
    write_tsv(outdir / "ranked_peptides.recommendation.tsv", ranked)
    # NetMHCpan sort proxy: sort by EL rank if present in pres, else score.
    net = list(ranked)
    def net_key(r: dict[str, str]):
        p = pres_by.get(r.get("peptide_id", ""), {})
        return safe_float(row_get(p, ["el_rank", "ba_rank", "ic50"], ""), 999999) or 999999
    net.sort(key=net_key)
    write_tsv(outdir / "ranked_peptides.netmhcpan42.tsv", net)
    write_tsv(outdir / "ranked_events.tsv", ranked)
    write_tsv(outdir / "validation_plan.tsv", [{"peptide_id": r.get("peptide_id", ""), "recommended_validation": "short_peptide" if r.get("source_type", "").lower() in {"snv", "peptide_csv"} else "long_peptide_or_minigene"} for r in ranked])
    res = {"status": "PASS", "skill": "neoag-ranking", "summary": f"Ranked {len(ranked)} peptide candidates", "outputs": {"recommendation": str(outdir / "ranked_peptides.recommendation.tsv")}}
    write_json(outdir / "skill_result.json", res)
    return res
