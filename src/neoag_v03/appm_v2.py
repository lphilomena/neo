"""APPM 2.0 evidence layer.

This module upgrades the original APPM-lite heuristic into a gene/pathway/
peptide evidence layer. It remains a computational triage module: it reports
antigen-presentation machinery evidence and flags, not a clinical diagnosis.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from .adapters.peptide_input import normalize_hla_allele
from .utils import clamp, first, read_tsv, to_float, write_tsv

MHC_I_CORE_GENES = ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "TAPBP", "NLRC5"]
MHC_I_PROCESSING_GENES = ["PSMB8", "PSMB9", "PSMB10", "ERAP1", "ERAP2", "CALR", "CANX", "PDIA3", "B2M"]
MHC_II_CORE_GENES = [
    "HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1",
    "CIITA", "RFX5", "RFXANK", "RFXAP",
]
IFNG_RESPONSE_GENES = ["IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1", "STAT2", "IRF1", "NLRC5", "CIITA"]
IMMUNE_CONTEXT_GENES = ["CXCL9", "CXCL10", "CXCL11", "IDO1", "GZMB", "PRF1", "CD8A", "CD8B"]
MHC_II_ACCESSORY_GENES = ["CD74", "HLA-DMA", "HLA-DMB"]

MHC_I_GENES = list(dict.fromkeys(MHC_I_CORE_GENES + MHC_I_PROCESSING_GENES))
IFNG_GENES = IFNG_RESPONSE_GENES
MHC_II_GENES = list(dict.fromkeys(MHC_II_CORE_GENES + MHC_II_ACCESSORY_GENES))
APPM_INTEGRITY_GENES = list(dict.fromkeys(MHC_I_GENES + IFNG_GENES + MHC_II_GENES))
APPM_CONTEXT_GENES = IMMUNE_CONTEXT_GENES

DAMAGING_TERMS = {
    "frameshift_variant", "stop_gained", "splice_acceptor_variant", "splice_donor_variant",
    "start_lost", "stop_lost", "protein_truncating_variant", "transcript_ablation",
    "feature_truncation", "loss_of_function",
}

APPM_GENE_FIELDS = [
    "sample_id", "gene", "pathway", "module", "gene_set", "appm_integrity_role", "damaging_variant", "mutation_status",
    "mutation_consequence", "total_cn", "major_cn", "minor_cn", "copy_number_status",
    "loh_status", "expression_tpm", "expression_percentile", "expression_status",
    "expression_input_status", "rna_assay_status", "protein_status", "ligandome_support",
    "biallelic_status", "functional_status", "gene_integrity_status", "gene_integrity_score",
    "functional_validation_status", "validation_evidence_source", "evidence_level",
    "evidence_completeness", "claim_strength", "reason", "risk_reason",
]

APPM_PATHWAY_FIELDS = [
    "sample_id", "pathway", "pathway_score", "pathway_status", "driver_defects",
    "evidence_level", "functional_validation_status", "claim_strength", "reason",
]

PEPTIDE_APPM_FIELDS = [
    "peptide_id", "event_id", "hla_allele", "mhc_class", "mhc_i_integrity_status", "mhc_ii_integrity_status",
    "ifng_response_status", "restricting_locus_expression_status", "restricting_locus_loh",
    "b2m_status", "tap_status", "nlrc5_status", "ciita_status", "appm_action", "appm_multiplier",
    "priority_cap", "appm_reason",
]

APPM_MODULE_SCORE_FIELDS = [
    "sample_id", "module", "score", "status", "driver_defects",
    "evidence_level", "functional_validation_status", "claim_strength", "reason",
]

APPM_IMMUNE_CONTEXT_FIELDS = [
    "sample_id", "gene", "gene_set", "expression_tpm", "expression_percentile",
    "expression_status", "context_marker_class", "context_interpretation", "appm_integrity_role",
]

APPM_COMPLETENESS_FIELDS = [
    "sample_id", "mutation_assessed", "cnv_assessed", "rna_assessed", "hla_loh_assessed",
    "protein_assessed", "flow_assessed", "ligandome_assessed", "cnv_input_status",
    "rna_input_status", "functional_validation_status", "validation_evidence_source",
    "appm_evidence_completeness_score", "appm_evidence_completeness_status",
    "missing_evidence", "interpretation",
]

APPM_INPUT_STATUS_FIELDS = [
    "sample_id", "input_type", "input_path", "input_status", "parse_status",
    "records_total", "records_appm_relevant", "assay_scope", "notes",
]

APPM_INPUT_STATUSES = {"assessed", "not_provided", "provided_empty", "failed_parse", "not_applicable"}

APPM_CONFLICT_FIELDS = [
    "sample_id", "gene", "conflict_type", "evidence_conflict_flag", "evidence_conflict_reason",
]

APPM_PEPTIDE_MODIFIER_FIELDS = [
    "peptide_id", "event_id", "hla_allele", "mhc_class", "appm_multiplier", "appm_multiplier_reason",
    "appm_integrity_status", "appm_evidence_completeness", "functional_validation_status",
    "appm_review_required", "priority_cap", "restricting_locus_expression_status",
    "restricting_locus_loh", "appm_action",
]


def _path_exists(path: str | Path | None) -> bool:
    return bool(path and Path(path).exists())


def _safe_read_tsv(path: str | Path | None) -> tuple[list[dict[str, str]], str, str, str]:
    if not path:
        return [], "not_provided", "not_provided", "path_not_provided"
    p = Path(path)
    if not p.exists():
        return [], "failed_parse", "missing_file", "file_not_found"
    try:
        rows = read_tsv(p)
    except Exception as exc:  # pragma: no cover - defensive around external inputs
        return [], "failed_parse", "failed_parse", f"{type(exc).__name__}: {exc}"
    if not rows:
        return [], "provided_empty", "parsed_empty", "no_data_rows"
    return rows, "assessed", "parsed", "ok"


def _has_any(row: Mapping[str, Any], keys: list[str]) -> bool:
    return bool(first(row, keys, ""))


def _count_rows_with(rows: list[dict[str, str]], keys: list[str]) -> int:
    return sum(1 for r in rows if _has_any(r, keys))


def _count_appm_gene_rows(rows: list[dict[str, str]]) -> int:
    appm_genes = set(APPM_INTEGRITY_GENES + APPM_CONTEXT_GENES)
    n = 0
    for r in rows:
        gene = first(r, ["gene", "Gene", "gene_name", "SYMBOL", "symbol"], "").upper()
        if gene in appm_genes:
            n += 1
    return n


def _input_row(
    sample_id: str,
    input_type: str,
    path: str | Path | None,
    rows: list[dict[str, str]],
    status: str,
    parse_status: str,
    notes: str,
    *,
    records_appm_relevant: int,
    assay_scope: str,
) -> dict[str, str]:
    if status == "assessed" and records_appm_relevant == 0:
        status = "failed_parse"
        parse_status = "parsed_no_appm_relevant_records"
        notes = "parsed_but_no_appm_relevant_columns_or_records"
    if status not in APPM_INPUT_STATUSES:
        status = "failed_parse"
    return {
        "sample_id": sample_id,
        "input_type": input_type,
        "input_path": str(path or ""),
        "input_status": status,
        "parse_status": parse_status,
        "records_total": str(len(rows)),
        "records_appm_relevant": str(records_appm_relevant),
        "assay_scope": assay_scope,
        "notes": notes,
    }


def _build_input_status_rows(
    sample_id: str,
    *,
    vep_tsv: str | Path | None,
    expression_tsv: str | Path | None,
    cnv_tsv: str | Path | None,
    hla_loh_tsv: str | Path | None,
    hla_typing_tsv: str | Path | None,
    tumor_purity_tsv: str | Path | None,
    proteomics_tsv: str | Path | None,
    phosphoproteomics_tsv: str | Path | None,
    hla_ligandome_tsv: str | Path | None,
    raw_peptides: str | Path | None,
) -> list[dict[str, str]]:
    rows_out: list[dict[str, str]] = []

    vep_rows, status, parse_status, notes = _safe_read_tsv(vep_tsv)
    rows_out.append(_input_row(sample_id, "mutation_vep", vep_tsv, vep_rows, status, parse_status, notes, records_appm_relevant=_count_rows_with(vep_rows, ["gene", "Gene", "SYMBOL", "symbol"]), assay_scope="gene_variant"))

    cnv_rows, status, parse_status, notes = _safe_read_tsv(cnv_tsv)
    cnv_gene = _count_rows_with(cnv_rows, ["gene", "Gene", "gene_name", "SYMBOL", "symbol"])
    cnv_segment = sum(1 for r in cnv_rows if _has_any(r, ["chrom", "Chromosome", "chr"]) and _has_any(r, ["start", "Start", "loc.start"]))
    scope = "gene_level_cnv_loh" if cnv_gene else ("segment_level_unmapped" if cnv_segment else "copy_number_loh")
    rows_out.append(_input_row(sample_id, "cnv_loh", cnv_tsv, cnv_rows, status, parse_status, notes, records_appm_relevant=_count_appm_gene_rows(cnv_rows) or cnv_gene or cnv_segment, assay_scope=scope))

    expr_rows, status, parse_status, notes = _safe_read_tsv(expression_tsv)
    rows_out.append(_input_row(sample_id, "rna_expression", expression_tsv, expr_rows, status, parse_status, notes, records_appm_relevant=_count_appm_gene_rows(expr_rows) or _count_rows_with(expr_rows, ["gene", "Gene", "gene_name", "symbol", "SYMBOL"]), assay_scope="gene_expression_tpm"))

    if hla_typing_tsv:
        hla_rows, status, parse_status, notes = _safe_read_tsv(hla_typing_tsv)
        hla_path = hla_typing_tsv
        scope = "hla_typing"
        relevant = _count_rows_with(hla_rows, ["hla_allele", "allele", "HLA", "hla", "hla_a", "hla_b", "hla_c"])
    else:
        hla_rows, status, parse_status, notes = _safe_read_tsv(raw_peptides)
        hla_path = raw_peptides
        scope = "derived_from_raw_peptides" if raw_peptides else "hla_typing"
        relevant = _count_rows_with(hla_rows, ["hla_allele", "allele", "HLA", "hla"])
        if raw_peptides and status == "assessed":
            notes = "raw_peptides_used_as_hla_typing_context"
    rows_out.append(_input_row(sample_id, "hla_typing", hla_path, hla_rows, status, parse_status, notes, records_appm_relevant=relevant, assay_scope=scope))

    hla_loh_rows, status, parse_status, notes = _safe_read_tsv(hla_loh_tsv)
    rows_out.append(_input_row(sample_id, "hla_loh", hla_loh_tsv, hla_loh_rows, status, parse_status, notes, records_appm_relevant=_count_rows_with(hla_loh_rows, ["hla_allele", "allele", "HLA", "LossAllele", "loss_allele"]), assay_scope="allele_level_loh"))

    purity_rows, status, parse_status, notes = _safe_read_tsv(tumor_purity_tsv)
    rows_out.append(_input_row(sample_id, "tumor_purity", tumor_purity_tsv, purity_rows, status, parse_status, notes, records_appm_relevant=_count_rows_with(purity_rows, ["purity", "tumor_purity", "purity_estimate"]), assay_scope="sample_level_purity"))

    prot_rows, status, parse_status, notes = _safe_read_tsv(proteomics_tsv)
    rows_out.append(_input_row(sample_id, "proteomics", proteomics_tsv, prot_rows, status, parse_status, notes, records_appm_relevant=_count_rows_with(prot_rows, ["gene", "Gene", "protein", "Protein", "SYMBOL", "symbol"]), assay_scope="protein_abundance"))

    phospho_rows, status, parse_status, notes = _safe_read_tsv(phosphoproteomics_tsv)
    rows_out.append(_input_row(sample_id, "phosphoproteomics", phosphoproteomics_tsv, phospho_rows, status, parse_status, notes, records_appm_relevant=_count_rows_with(phospho_rows, ["gene", "Gene", "protein", "Protein", "phosphosite", "site", "SYMBOL", "symbol"]), assay_scope="phosphoprotein_activation"))

    ligand_rows, status, parse_status, notes = _safe_read_tsv(hla_ligandome_tsv)
    rows_out.append(_input_row(sample_id, "hla_ligandome", hla_ligandome_tsv, ligand_rows, status, parse_status, notes, records_appm_relevant=_count_rows_with(ligand_rows, ["peptide", "sequence", "hla_allele", "allele"]), assay_scope="hla_ligandome_ms"))
    return rows_out


def _input_assessed(input_rows: list[dict[str, str]], input_type: str) -> bool:
    return any(r.get("input_type") == input_type and r.get("input_status") == "assessed" and int(r.get("records_appm_relevant") or 0) > 0 for r in input_rows)


def load_expression(path: str | Path | None) -> dict[str, float]:
    if not _path_exists(path):
        return {}
    d: dict[str, float] = {}
    for r in read_tsv(path):
        gene = first(r, ["gene", "Gene", "gene_name", "symbol", "SYMBOL"], "").upper()
        if gene:
            d[gene] = to_float(first(r, ["TPM", "tpm", "expression", "expr", "gene_tpm"], "0"), 0.0)
    return d


def expression_input_status(path: str | Path | None, expr: dict[str, float]) -> str:
    if not _path_exists(path):
        return "not_provided"
    if not expr:
        return "provided_no_gene_symbols"
    return "gene_level_tpm"


def load_variants(path: str | Path | None) -> dict[str, list[str]]:
    if not _path_exists(path):
        return {}
    d: dict[str, list[str]] = {}
    for r in read_tsv(path):
        gene = first(r, ["gene", "Gene", "SYMBOL", "symbol"], "").upper()
        cons = first(r, ["mutation_consequence", "consequence", "Consequence", "variant_consequence", "effect", "IMPACT"], "")
        damaging_class = first(r, ["damaging_class", "mutation_class"], "")
        if first(r, ["is_damaging_missense"], "").lower() == "yes" and "damaging_missense" not in cons:
            cons = f"{cons}&damaging_missense" if cons else "damaging_missense"
        if first(r, ["damaging_variant"], "").lower() == "yes" and damaging_class and damaging_class != "none":
            cons = f"{cons}&{damaging_class}" if cons else damaging_class
        if gene and cons:
            bucket = d.setdefault(gene, [])
            for item in str(cons).replace(";", "&").replace(",", "&").split("&"):
                item = item.strip()
                if item and item not in bucket:
                    bucket.append(item)
    return d


def _loss_status(status: str) -> bool:
    s = str(status or "").strip().lower()
    if not s or s in {"neutral", "normal", "intact", "no_loh", "not_assessed", "unknown", "na", "n/a"}:
        return False
    if s.startswith("no_") or s.startswith("not_"):
        return False
    return s in {
        "loss",
        "loh",
        "copy_loss",
        "copy_number_loss",
        "deletion",
        "deep_deletion",
        "biallelic_deletion",
        "homdel",
        "homozygous_deletion",
        "copy_neutral_loh",
        "cnloh",
    }


def _homdel_status(status: str) -> bool:
    s = str(status or "").lower()
    return any(x in s for x in ["homdel", "homozygous", "deep_deletion", "biallelic_deletion"])


def cnv_input_status(path: str | Path | None, cnv: dict[str, dict[str, str]]) -> str:
    if not _path_exists(path):
        return "not_provided"
    rows = read_tsv(path)
    if cnv:
        return "gene_level_cnv"
    if rows and any(first(r, ["chrom", "Chromosome", "chr"], "") and first(r, ["start", "Start", "loc.start"], "") for r in rows[:50]):
        return "segment_level_unmapped"
    return "provided_no_gene_symbols"


def load_cnv_gene_status(path: str | Path | None) -> dict[str, dict[str, str]]:
    """Load gene-level CNV when available.

    Segment-only CNV tables are not gene-resolved here. Segment→gene mapping is
    intentionally left to upstream CNV adapters because it requires a transcript
    or gene interval model.
    """
    if not _path_exists(path):
        return {}
    out: dict[str, dict[str, str]] = {}
    for r in read_tsv(path):
        gene = first(r, ["gene", "Gene", "gene_name", "SYMBOL", "symbol"], "").upper()
        if not gene:
            continue
        status = first(r, ["copy_number_status", "cn_status", "loh_status", "status", "call"], "")
        total = first(r, ["total_cn", "copy_number", "total_copy_number", "tcn", "cn"], "")
        major = first(r, ["major_cn", "major_copy_number", "lcn.em", "cf.em"], "")
        minor = first(r, ["minor_cn", "minor_copy_number", "lcn", "minor"], "")
        loh = first(r, ["loh_status", "LOH", "loh", "cnloh"], "")
        out[gene] = {
            "copy_number_status": status or loh or "altered",
            "total_cn": total,
            "major_cn": major,
            "minor_cn": minor,
            "loh_status": loh or ("loh" if _loss_status(status) else ""),
        }
    return out


def load_hla_loh(path: str | Path | None) -> set[str]:
    if not _path_exists(path):
        return set()
    lost: set[str] = set()
    for r in read_tsv(path):
        allele = first(r, ["hla_allele", "allele", "HLA", "LossAllele", "loss_allele"], "")
        status = first(r, ["loh_status", "LOH", "status", "loss", "Loss"], "")
        if allele and (not status or status.lower() in {"loh", "loss", "lost", "yes", "true", "1"}):
            lost.add(normalize_hla_allele(allele))
    return lost


def _gene_group(gene: str) -> str:
    if gene in IMMUNE_CONTEXT_GENES:
        return "optional_immune_context"
    groups: list[str] = []
    if gene in MHC_I_CORE_GENES:
        groups.append("mhc_i_core")
    if gene in MHC_I_PROCESSING_GENES:
        groups.append("mhc_i_processing")
    if gene in MHC_II_CORE_GENES:
        groups.append("mhc_ii_core")
    if gene in IFNG_RESPONSE_GENES:
        groups.append("ifng_response")
    if gene in MHC_II_ACCESSORY_GENES:
        groups.append("mhc_ii_accessory")
    return ";".join(groups) if groups else "other"


def _gene_pathway(gene: str) -> str:
    if gene in IMMUNE_CONTEXT_GENES:
        return "immune context annotation"
    if gene in MHC_I_CORE_GENES or gene in MHC_I_PROCESSING_GENES:
        return "MHC-I antigen presentation"
    if gene in IFNG_RESPONSE_GENES:
        return "IFNG-JAK-STAT"
    if gene in MHC_II_GENES:
        return "MHC-II antigen presentation"
    return "other"


def _context_marker_class(gene: str) -> str:
    if gene in {"CXCL9", "CXCL10", "CXCL11", "IDO1"}:
        return "interferon_inflamed_context"
    if gene in {"GZMB", "PRF1", "CD8A", "CD8B"}:
        return "cytotoxic_t_cell_context"
    return "immune_context"


def _damaging(cons: list[str]) -> bool:
    flat = ",".join(cons).lower()
    return any(term in flat for term in DAMAGING_TERMS)


def _expression_status(gene: str, expr: dict[str, float], assessed: bool, low_tpm: float) -> tuple[str, str]:
    if not assessed:
        return "", "unassessed"
    if gene not in expr:
        return "", "missing_from_expression_matrix"
    tpm = expr.get(gene, 0.0)
    return f"{tpm:.4f}", "low" if tpm < low_tpm else "expressed"


def _biallelic_status(gene: str, variants: dict[str, list[str]], cnv: dict[str, dict[str, str]], expr: dict[str, float], expression_assessed: bool, low_tpm: float) -> tuple[str, str]:
    cons = variants.get(gene, [])
    has_dmg = _damaging(cons)
    cn = cnv.get(gene, {})
    cn_status = cn.get("copy_number_status", "") or cn.get("loh_status", "")
    loh_status = cn.get("loh_status", "")
    low_expr = expression_assessed and gene in expr and expr.get(gene, 999.0) < low_tpm
    if _homdel_status(cn_status):
        return "BIALLELIC_LOSS", "homozygous_or_deep_deletion"
    if has_dmg and (_loss_status(cn_status) or _loss_status(loh_status)):
        return "BIALLELIC_LOSS", "damaging_variant_plus_loh_or_copy_loss"
    if has_dmg and low_expr:
        return "BIALLELIC_LOSS", "damaging_variant_plus_expression_loss"
    if has_dmg:
        return "MONOALLELIC_ALTERATION", "damaging_variant_only"
    if _loss_status(cn_status) or _loss_status(loh_status):
        return "MONOALLELIC_LOSS", "copy_loss_or_loh_only"
    if low_expr:
        return "EXPRESSION_LOSS_ONLY", "low_expression_only"
    return "NO_EVIDENCE", "no_major_signal"


def _functional_status(biallelic: str) -> str:
    if biallelic == "BIALLELIC_LOSS":
        return "defective"
    if biallelic in {"MONOALLELIC_ALTERATION", "MONOALLELIC_LOSS", "EXPRESSION_LOSS_ONLY"}:
        return "caution"
    return "intact"


def _mutation_status(cons: list[str]) -> str:
    if _damaging(cons):
        return "damaging_variant"
    if cons:
        return "variant_non_damaging_or_uncertain"
    return "intact"


def _gene_integrity_status(
    *,
    biallelic_status: str,
    damaging: bool,
    copy_status: str,
    loh_status: str,
    expression_status: str,
    conflict: bool = False,
) -> str:
    if conflict:
        return "conflicting"
    if biallelic_status == "BIALLELIC_LOSS":
        return "biallelic_loss"
    if biallelic_status in {"MONOALLELIC_LOSS", "MONOALLELIC_ALTERATION"}:
        return "monoallelic_loss" if biallelic_status == "MONOALLELIC_LOSS" else "damaging_variant"
    if biallelic_status == "EXPRESSION_LOSS_ONLY" or expression_status == "low":
        return "low_expression"
    if _homdel_status(copy_status):
        return "biallelic_loss"
    if _loss_status(copy_status):
        return "copy_loss"
    if _loss_status(loh_status):
        return "loh"
    if expression_status in {"unassessed", "missing_from_expression_matrix"} and not damaging and not copy_status:
        return "not_assessed"
    return "intact"


def _gene_integrity_score(status: str) -> float:
    return {
        "intact": 1.0,
        "not_assessed": 1.0,
        "loh": 0.75,
        "copy_loss": 0.70,
        "damaging_variant": 0.65,
        "monoallelic_loss": 0.60,
        "low_expression": 0.50,
        "conflicting": 0.50,
        "biallelic_loss": 0.05,
    }.get(status, 1.0)



def _status_from_score(prefix: str, score: float, assessed_any: bool) -> str:
    if not assessed_any:
        return f"{prefix}_UNASSESSED"
    if score < 0.30:
        return f"{prefix}_DEFECTIVE"
    if score < 0.75:
        return f"{prefix}_CAUTION"
    return f"{prefix}_INTACT"


def _functional_validation_status(*, protein_assessed: bool = False, flow_assessed: bool = False, ligandome_assessed: bool = False) -> tuple[str, str, str]:
    sources: list[str] = []
    if protein_assessed:
        sources.append("protein")
    if flow_assessed:
        sources.append("flow")
    if ligandome_assessed:
        sources.append("ligandome")
    if ligandome_assessed:
        return "ligandome_supported", ";".join(sources), "validation_supported_hypothesis"
    if flow_assessed:
        return "flow_supported", ";".join(sources), "validation_supported_hypothesis"
    if protein_assessed:
        return "protein_supported", ";".join(sources), "validation_supported_hypothesis"
    return "computational_proxy", "DNA_CNV_RNA_HLA_LOH_only", "computational_triage_only"


def _evidence_completeness(
    sample_id: str,
    *,
    variants: dict[str, list[str]],
    cnv: dict[str, dict[str, str]],
    expression_assessed: bool,
    lost_hla: set[str],
    hla_loh_provided: bool,
    cnv_status: str,
    rna_status: str,
    protein_assessed: bool = False,
    flow_assessed: bool = False,
    ligandome_assessed: bool = False,
) -> dict[str, str]:
    functional_status, validation_source, claim_strength = _functional_validation_status(
        protein_assessed=protein_assessed,
        flow_assessed=flow_assessed,
        ligandome_assessed=ligandome_assessed,
    )
    checks = {
        "mutation_assessed": bool(variants),
        "cnv_assessed": bool(cnv),
        "rna_assessed": bool(expression_assessed),
        "hla_loh_assessed": bool(hla_loh_provided),
        "protein_assessed": protein_assessed,
        "flow_assessed": flow_assessed,
        "ligandome_assessed": ligandome_assessed,
    }
    score = sum(1 for v in checks.values() if v) / len(checks)
    missing = [k.replace("_assessed", "") for k, v in checks.items() if not v]
    if cnv_status == "segment_level_unmapped":
        missing.append("cnv_gene_mapping")
    if score >= 0.67:
        status = "HIGH"
    elif score >= 0.34:
        status = "PARTIAL"
    elif score > 0:
        status = "LOW"
    else:
        status = "UNASSESSED"
    return {
        "sample_id": sample_id,
        **{k: "yes" if v else "no" for k, v in checks.items()},
        "cnv_input_status": cnv_status,
        "rna_input_status": rna_status,
        "functional_validation_status": functional_status,
        "validation_evidence_source": validation_source,
        "appm_evidence_completeness_score": f"{score:.4f}",
        "appm_evidence_completeness_status": status,
        "missing_evidence": ";".join(missing) if missing else "none",
        "interpretation": "missing_evidence_limits_appm_claim_strength" if missing else "appm_evidence_inputs_available",
    }


def _detect_conflicts(sample_id: str, gene_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for r in gene_rows:
        gene = r.get("gene", "")
        expr_status = r.get("expression_status", "")
        bial = r.get("biallelic_status", "")
        cn_status = r.get("copy_number_status", "")
        damaging = r.get("damaging_variant") == "yes"
        reasons: list[str] = []
        ctype = ""
        if bial == "BIALLELIC_LOSS" and expr_status == "expressed":
            ctype = "loss_but_rna_expressed"
            reasons.append("biallelic_loss_with_detectable_rna_expression")
        if damaging and expr_status == "expressed" and bial != "BIALLELIC_LOSS":
            ctype = ctype or "damaging_variant_but_rna_expressed"
            reasons.append("damaging_variant_without_expression_loss")
        if _loss_status(cn_status) and expr_status == "expressed":
            ctype = ctype or "copy_loss_but_rna_expressed"
            reasons.append("copy_loss_with_detectable_rna_expression")
        if reasons:
            rows.append({
                "sample_id": sample_id,
                "gene": gene,
                "conflict_type": ctype,
                "evidence_conflict_flag": "yes",
                "evidence_conflict_reason": ";".join(reasons),
            })
    if not rows:
        rows.append({
            "sample_id": sample_id,
            "gene": "NA",
            "conflict_type": "none",
            "evidence_conflict_flag": "no",
            "evidence_conflict_reason": "no_conflicts_detected",
        })
    return rows

def build_appm_2(
    *,
    sample_id: str,
    outdir: str | Path,
    vep_tsv: str | Path | None = None,
    expression_tsv: str | Path | None = None,
    hla_loh_tsv: str | Path | None = None,
    cnv_tsv: str | Path | None = None,
    raw_peptides: str | Path | None = None,
    hla_typing_tsv: str | Path | None = None,
    tumor_purity_tsv: str | Path | None = None,
    proteomics_tsv: str | Path | None = None,
    phosphoproteomics_tsv: str | Path | None = None,
    hla_ligandome_tsv: str | Path | None = None,
    profile: Mapping[str, Any] | None = None,
) -> dict[str, str]:
    profile = profile or {}
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    low_tpm = float(profile.get("appm", {}).get("low_expression_tpm", 1.0))
    input_status_rows = _build_input_status_rows(
        sample_id,
        vep_tsv=vep_tsv,
        expression_tsv=expression_tsv,
        cnv_tsv=cnv_tsv,
        hla_loh_tsv=hla_loh_tsv,
        hla_typing_tsv=hla_typing_tsv,
        tumor_purity_tsv=tumor_purity_tsv,
        proteomics_tsv=proteomics_tsv,
        phosphoproteomics_tsv=phosphoproteomics_tsv,
        hla_ligandome_tsv=hla_ligandome_tsv,
        raw_peptides=raw_peptides,
    )
    protein_assessed = _input_assessed(input_status_rows, "proteomics") or _input_assessed(input_status_rows, "phosphoproteomics")
    ligandome_assessed = _input_assessed(input_status_rows, "hla_ligandome")
    variants = load_variants(vep_tsv)
    expr = load_expression(expression_tsv)
    rna_status = expression_input_status(expression_tsv, expr)
    expression_assessed = rna_status == "gene_level_tpm"
    cnv = load_cnv_gene_status(cnv_tsv)
    cnv_status = cnv_input_status(cnv_tsv, cnv)
    lost_hla = load_hla_loh(hla_loh_tsv)
    hla_loh_provided = _path_exists(hla_loh_tsv)

    genes = sorted(set(APPM_INTEGRITY_GENES + APPM_CONTEXT_GENES) | set(variants))
    expression_values = sorted(expr.values())

    def expression_percentile(gene: str) -> str:
        if not expression_assessed or gene not in expr or not expression_values:
            return ""
        val = expr[gene]
        le = sum(1 for x in expression_values if x <= val)
        return f"{(100.0 * le / len(expression_values)):.2f}"

    gene_rows: list[dict[str, str]] = []
    for g in genes:
        cons = variants.get(g, [])
        biallelic, reason = _biallelic_status(g, variants, cnv, expr, expression_assessed, low_tpm)
        functional = _functional_status(biallelic)
        expr_tpm, expr_status = _expression_status(g, expr, expression_assessed, low_tpm)
        gene_functional_validation_status, gene_validation_source, gene_claim_strength = _functional_validation_status()
        cn = cnv.get(g, {})
        pathway = _gene_pathway(g)
        gene_group = _gene_group(g)
        appm_integrity_role = "context_annotation_only" if g in IMMUNE_CONTEXT_GENES else "integrity_scored"
        damaging = _damaging(cons)
        copy_status = cn.get("copy_number_status", "unknown")
        loh_gene_status = cn.get("loh_status", "not_assessed")
        integrity_status = _gene_integrity_status(
            biallelic_status=biallelic,
            damaging=damaging,
            copy_status=copy_status,
            loh_status=loh_gene_status,
            expression_status=expr_status,
        )
        risk_reason = reason if expression_assessed or variants or cnv else "expression_not_assessed"
        if appm_integrity_role == "context_annotation_only":
            biallelic = "NOT_APPLICABLE"
            functional = "context_annotation"
            integrity_status = "not_assessed"
            reason = "immune_context_not_used_for_appm_integrity"
            risk_reason = reason
        gene_rows.append({
            "sample_id": sample_id,
            "gene": g,
            "pathway": pathway,
            "module": pathway,
            "gene_set": gene_group,
            "appm_integrity_role": appm_integrity_role,
            "damaging_variant": "yes" if damaging else "no",
            "mutation_status": _mutation_status(cons),
            "mutation_consequence": ",".join(cons),
            "total_cn": cn.get("total_cn", ""),
            "major_cn": cn.get("major_cn", ""),
            "minor_cn": cn.get("minor_cn", ""),
            "copy_number_status": copy_status,
            "loh_status": loh_gene_status,
            "expression_tpm": expr_tpm,
            "expression_percentile": expression_percentile(g),
            "expression_status": expr_status,
            "expression_input_status": "gene_present" if g in expr else ("gene_missing" if expression_assessed else "not_assessed"),
            "rna_assay_status": rna_status,
            "protein_status": "not_assessed",
            "ligandome_support": "not_assessed",
            "biallelic_status": biallelic,
            "functional_status": functional,
            "gene_integrity_status": integrity_status,
            "gene_integrity_score": f"{_gene_integrity_score(integrity_status):.4f}",
            "functional_validation_status": gene_functional_validation_status,
            "validation_evidence_source": gene_validation_source,
            "evidence_level": "HIGH" if biallelic == "BIALLELIC_LOSS" else ("MEDIUM" if functional == "caution" else "LOW"),
            "evidence_completeness": "PARTIAL" if any([variants, cnv, expression_assessed, hla_loh_provided]) else "UNASSESSED",
            "claim_strength": gene_claim_strength,
            "reason": risk_reason,
            "risk_reason": risk_reason,
        })

    context_rows = [
        {
            "sample_id": r["sample_id"],
            "gene": r["gene"],
            "gene_set": r["gene_set"],
            "expression_tpm": r["expression_tpm"],
            "expression_percentile": r["expression_percentile"],
            "expression_status": r["expression_status"],
            "context_marker_class": _context_marker_class(r["gene"]),
            "context_interpretation": "background_annotation_not_appm_integrity",
            "appm_integrity_role": r["appm_integrity_role"],
        }
        for r in gene_rows
        if r.get("appm_integrity_role") == "context_annotation_only"
    ]

    by_gene = {r["gene"]: r for r in gene_rows}

    def bial(g: str) -> bool:
        return by_gene.get(g, {}).get("biallelic_status") == "BIALLELIC_LOSS"

    def caution(g: str) -> bool:
        return by_gene.get(g, {}).get("functional_status") == "caution"

    # Rule-first pathway states.
    mhc_i_score = 1.0
    mhc_i_reasons: list[str] = []
    if bial("B2M"):
        mhc_i_score = min(mhc_i_score, 0.05); mhc_i_reasons.append("B2M_biallelic_loss")
    if bial("TAP1") or bial("TAP2"):
        mhc_i_score = min(mhc_i_score, 0.25); mhc_i_reasons.append("TAP1_TAP2_biallelic_defect")
    if lost_hla:
        mhc_i_score = min(mhc_i_score, 0.70 if len(lost_hla) == 1 else 0.45); mhc_i_reasons.append("HLA_LOH")
    if caution("NLRC5"):
        mhc_i_score = min(mhc_i_score, 0.75); mhc_i_reasons.append("NLRC5_caution")
    if bial("JAK1") or bial("JAK2") or bial("STAT1"):
        mhc_i_reasons.append("IFNG_response_defect")
        mhc_i_score = min(mhc_i_score, 0.65)

    mhc_ii_score = 1.0
    mhc_ii_reasons: list[str] = []
    if bial("CIITA") or bial("RFX5") or bial("RFXANK") or bial("RFXAP"):
        mhc_ii_score = min(mhc_ii_score, 0.15); mhc_ii_reasons.append("CIITA_RFX_biallelic_defect")
    elif caution("CIITA") or caution("RFX5") or caution("RFXANK") or caution("RFXAP"):
        mhc_ii_score = min(mhc_ii_score, 0.60); mhc_ii_reasons.append("MHC_II_regulator_caution")

    ifng_score = 1.0
    ifng_reasons: list[str] = []
    if bial("JAK1") or bial("JAK2") or bial("IFNGR1") or bial("IFNGR2") or bial("STAT1"):
        ifng_score = 0.20; ifng_reasons.append("JAK_STAT_IFNGR_biallelic_defect")
    elif any(caution(g) for g in ["JAK1", "JAK2", "IFNGR1", "IFNGR2", "STAT1", "IRF1", "PTPN2"]):
        ifng_score = 0.65; ifng_reasons.append("IFNG_pathway_caution")

    assessed_any = bool(expression_assessed or variants or cnv or lost_hla)
    completeness = _evidence_completeness(
        sample_id,
        variants=variants,
        cnv=cnv,
        expression_assessed=expression_assessed,
        lost_hla=lost_hla,
        hla_loh_provided=hla_loh_provided,
        cnv_status=cnv_status,
        rna_status=rna_status,
        protein_assessed=protein_assessed,
        ligandome_assessed=ligandome_assessed,
    )

    module_rows = [
        {"sample_id": sample_id, "module": "MHC-I", "score": f"{clamp(mhc_i_score):.4f}", "status": _status_from_score("MHC_I", mhc_i_score, assessed_any), "driver_defects": ";".join(mhc_i_reasons), "evidence_level": "HIGH" if mhc_i_score < 0.3 else ("MEDIUM" if mhc_i_score < 0.75 else "LOW"), "functional_validation_status": completeness["functional_validation_status"], "claim_strength": "computational_triage_only", "reason": ";".join(mhc_i_reasons) or ("unassessed" if not assessed_any else "no_major_signal")},
        {"sample_id": sample_id, "module": "MHC-II", "score": f"{clamp(mhc_ii_score):.4f}", "status": _status_from_score("MHC_II", mhc_ii_score, assessed_any), "driver_defects": ";".join(mhc_ii_reasons), "evidence_level": "HIGH" if mhc_ii_score < 0.3 else ("MEDIUM" if mhc_ii_score < 0.75 else "LOW"), "functional_validation_status": completeness["functional_validation_status"], "claim_strength": "computational_triage_only", "reason": ";".join(mhc_ii_reasons) or ("unassessed" if not assessed_any else "no_major_signal")},
        {"sample_id": sample_id, "module": "IFNG-JAK-STAT", "score": f"{clamp(ifng_score):.4f}", "status": _status_from_score("IFNG_RESPONSE", ifng_score, assessed_any), "driver_defects": ";".join(ifng_reasons), "evidence_level": "HIGH" if ifng_score < 0.3 else ("MEDIUM" if ifng_score < 0.75 else "LOW"), "functional_validation_status": completeness["functional_validation_status"], "claim_strength": "computational_triage_only", "reason": ";".join(ifng_reasons) or ("unassessed" if not assessed_any else "no_major_signal")},
    ]
    pathway_rows = [
        {
            "sample_id": r["sample_id"],
            "pathway": r["module"],
            "pathway_score": r["score"],
            "pathway_status": r["status"],
            "driver_defects": r["driver_defects"],
            "evidence_level": r["evidence_level"],
            "functional_validation_status": r["functional_validation_status"],
            "claim_strength": r["claim_strength"],
            "reason": r["reason"],
        }
        for r in module_rows
    ]
    conflict_rows = _detect_conflicts(sample_id, gene_rows)
    conflict_genes = {r.get("gene", "") for r in conflict_rows if r.get("evidence_conflict_flag") == "yes"}
    for row in gene_rows:
        if row.get("gene") in conflict_genes:
            row["gene_integrity_status"] = "conflicting"
            row["gene_integrity_score"] = f"{_gene_integrity_score('conflicting'):.4f}"
            row["risk_reason"] = ";".join(x for x in [row.get("risk_reason", ""), "evidence_conflict"] if x)

    peptide_flags: list[dict[str, str]] = []
    peptide_modifiers: list[dict[str, str]] = []
    if raw_peptides and Path(raw_peptides).exists():
        for p in read_tsv(raw_peptides):
            mhc = str(p.get("mhc_class", "I")).upper()
            hla = normalize_hla_allele(p.get("hla_allele", ""))
            locus = ""
            if hla.startswith("HLA-A"):
                locus = "HLA-A"
            elif hla.startswith("HLA-B"):
                locus = "HLA-B"
            elif hla.startswith("HLA-C"):
                locus = "HLA-C"
            locus_expr = by_gene.get(locus, {}).get("expression_status", "unassessed") if locus else "unassessed"
            reasons: list[str] = []
            mult = 1.0; cap = ""; action = "PASS"
            if mhc in {"I", "MHC-I", "CLASSI"}:
                if bial("B2M"):
                    action = "REJECT_OR_CAP"; mult = 0.0; cap = "D"; reasons.append("B2M_biallelic_loss")
                if hla in lost_hla:
                    # HLA-allele LOH is handled by the immune_escape peptide layer; APPM records review context
                    # without applying a second multiplier penalty.
                    action = "REVIEW" if action == "PASS" else action; reasons.append("restricting_HLA_LOH_review_in_immune_escape")
                if bial("TAP1") or bial("TAP2"):
                    if mult > 0:
                        action = "CAP"; mult = min(mult, 0.35); cap = cap or "C"; reasons.append("TAP_processing_defect")
                if caution("NLRC5") and mult > 0:
                    action = "CAUTION" if action == "PASS" else action; mult = min(mult, 0.70); cap = cap or "B_CAUTION"; reasons.append("NLRC5_caution")
            else:
                if bial("CIITA") or bial("RFX5") or bial("RFXANK") or bial("RFXAP"):
                    action = "CAP"; mult = min(mult, 0.30); cap = "C"; reasons.append("MHC_II_regulator_defect")
            if ifng_score < 0.30 and mult > 0:
                action = "CAUTION" if action == "PASS" else action; mult = min(mult, 0.65); cap = cap or "B_CAUTION"; reasons.append("IFNG_response_defect")
            if not reasons:
                reasons.append("no_major_signal")
            flag_row = {
                "peptide_id": p.get("peptide_id", ""), "event_id": p.get("event_id", ""),
                "hla_allele": p.get("hla_allele", ""), "mhc_class": p.get("mhc_class", ""),
                "mhc_i_integrity_status": pathway_rows[0]["pathway_status"],
                "mhc_ii_integrity_status": pathway_rows[1]["pathway_status"],
                "ifng_response_status": pathway_rows[2]["pathway_status"],
                "restricting_locus_expression_status": locus_expr,
                "restricting_locus_loh": "yes" if hla in lost_hla else "no",
                "b2m_status": by_gene.get("B2M", {}).get("biallelic_status", "NO_EVIDENCE"),
                "tap_status": "DEFECT" if bial("TAP1") or bial("TAP2") else "NO_HIGH_RISK_SIGNAL",
                "nlrc5_status": by_gene.get("NLRC5", {}).get("functional_status", "unassessed"),
                "ciita_status": by_gene.get("CIITA", {}).get("functional_status", "unassessed"),
                "appm_action": action,
                "appm_multiplier": f"{mult:.4f}",
                "priority_cap": cap,
                "appm_reason": ";".join(reasons),
            }
            peptide_flags.append(flag_row)
            if mhc in {"II", "MHC-II", "CLASSII"}:
                integrity = pathway_rows[1]["pathway_status"]
            else:
                integrity = pathway_rows[0]["pathway_status"]
            peptide_modifiers.append({
                "peptide_id": flag_row["peptide_id"],
                "event_id": flag_row["event_id"],
                "hla_allele": flag_row["hla_allele"],
                "mhc_class": flag_row["mhc_class"],
                "appm_multiplier": flag_row["appm_multiplier"],
                "appm_multiplier_reason": flag_row["appm_reason"],
                "appm_integrity_status": integrity,
                "appm_evidence_completeness": completeness["appm_evidence_completeness_status"],
                "functional_validation_status": completeness["functional_validation_status"],
                "appm_review_required": "yes" if action in {"REVIEW", "CAUTION", "CAP", "REJECT_OR_CAP"} or completeness["appm_evidence_completeness_status"] in {"LOW", "UNASSESSED"} else "no",
                "priority_cap": flag_row["priority_cap"],
                "restricting_locus_expression_status": flag_row["restricting_locus_expression_status"],
                "restricting_locus_loh": flag_row["restricting_locus_loh"],
                "appm_action": flag_row["appm_action"],
            })

    summary_status = "UNASSESSED" if not assessed_any else ("DEFECTIVE" if min(mhc_i_score, mhc_ii_score) < 0.30 else ("CAUTION" if min(mhc_i_score, mhc_ii_score, ifng_score) < 0.75 else "PASS"))
    summary = {
        "sample_id": sample_id,
        "mhc_i_integrity_score": f"{clamp(mhc_i_score):.4f}",
        "mhc_ii_integrity_score": f"{clamp(mhc_ii_score):.4f}",
        "ifng_response_score": f"{clamp(ifng_score):.4f}",
        "mhc_i_integrity_status": pathway_rows[0]["pathway_status"],
        "mhc_ii_integrity_status": pathway_rows[1]["pathway_status"],
        "ifng_response_status": pathway_rows[2]["pathway_status"],
        "hla_i_loh_flag": "yes" if lost_hla else "no",
        "hla_loh_alleles": ",".join(sorted(lost_hla)),
        "b2m_risk": "yes" if bial("B2M") else ("caution" if caution("B2M") else "no"),
        "tap_risk": "yes" if bial("TAP1") or bial("TAP2") else ("caution" if caution("TAP1") or caution("TAP2") else "no"),
        "nlrc5_risk": by_gene.get("NLRC5", {}).get("functional_status", "unassessed"),
        "ciita_risk": by_gene.get("CIITA", {}).get("functional_status", "unassessed"),
        "expression_assessment_status": "assessed" if expression_assessed else "unassessed",
        "rna_input_status": rna_status,
        "cnv_input_status": cnv_status,
        "functional_validation_status": completeness["functional_validation_status"],
        "validation_evidence_source": completeness["validation_evidence_source"],
        "claim_strength": "computational_triage_only",
        "appm_evidence_completeness": completeness["appm_evidence_completeness_status"],
        "appm_evidence_completeness_score": completeness["appm_evidence_completeness_score"],
        "evidence_conflict_flag": "yes" if any(r.get("evidence_conflict_flag") == "yes" for r in conflict_rows) else "no",
        "appm_overall_status": summary_status,
        "interpretation": "appm_evidence_not_clinical_resistance_diagnosis",
    }

    paths = {
        "appm_gene_status": str(out / "appm_gene_status.tsv"),
        "appm_pathway_status": str(out / "appm_pathway_status.tsv"),
        "appm_module_scores": str(out / "appm_module_scores.tsv"),
        "appm_submodule_scores": str(out / "appm_submodule_scores.tsv"),
        "appm_immune_context": str(out / "appm_immune_context.tsv"),
        "appm_evidence_completeness": str(out / "appm_evidence_completeness.tsv"),
        "appm_input_status": str(out / "appm_input_status.tsv"),
        "appm_conflicts": str(out / "appm_conflicts.tsv"),
        "peptide_appm_flags": str(out / "peptide_appm_flags.tsv"),
        "appm_peptide_modifiers": str(out / "appm_peptide_modifiers.tsv"),
        "appm_summary": str(out / "appm_summary.tsv"),
    }
    write_tsv(paths["appm_gene_status"], gene_rows, APPM_GENE_FIELDS)
    write_tsv(paths["appm_pathway_status"], pathway_rows, APPM_PATHWAY_FIELDS)
    write_tsv(paths["appm_module_scores"], module_rows, APPM_MODULE_SCORE_FIELDS)
    write_tsv(paths["appm_immune_context"], context_rows, APPM_IMMUNE_CONTEXT_FIELDS)
    write_tsv(paths["appm_evidence_completeness"], [completeness], APPM_COMPLETENESS_FIELDS)
    write_tsv(paths["appm_input_status"], input_status_rows, APPM_INPUT_STATUS_FIELDS)
    write_tsv(paths["appm_conflicts"], conflict_rows, APPM_CONFLICT_FIELDS)
    write_tsv(paths["peptide_appm_flags"], peptide_flags, PEPTIDE_APPM_FIELDS)
    write_tsv(paths["appm_peptide_modifiers"], peptide_modifiers, APPM_PEPTIDE_MODIFIER_FIELDS)
    write_tsv(paths["appm_summary"], [summary])
    # v0.4.2 P1: append APPM call confidence and submodule explainability sidecars.
    try:
        from .appm_explainability import enhance_appm_outputs_v042
        paths.update(enhance_appm_outputs_v042(out, raw_peptides=raw_peptides, profile=profile))
    except Exception as exc:  # pragma: no cover - explanatory sidecars must not break legacy APPM output
        write_tsv(out / "appm_explainability_warning.tsv", [{"warning": type(exc).__name__, "message": str(exc)}])
    return paths


def legacy_appm_rows_from_v2(gene_rows: list[dict[str, str]], expression_assessed: bool) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for r in gene_rows:
        pathway = "MHC-I" if r.get("pathway") in {"MHC-I antigen presentation", "IFNG-JAK-STAT"} else "MHC-II"
        risk_flag = "yes" if r.get("functional_status") in {"defective", "caution"} else "no"
        rows.append({
            "sample_id": r.get("sample_id", ""),
            "pathway": pathway,
            "gene": r.get("gene", ""),
            "mutation_status": "damaging" if r.get("damaging_variant") == "yes" else "none",
            "mutation_consequence": r.get("mutation_consequence", ""),
            "expression_tpm": r.get("expression_tpm", ""),
            "expression_status": r.get("expression_status", "unassessed"),
            "copy_number_status": r.get("copy_number_status", "unknown"),
            "loh_status": r.get("loh_status", "not_assessed"),
            "risk_flag": risk_flag,
            "risk_reason": r.get("reason", "") or ("expression_not_assessed" if not expression_assessed else "no_major_signal"),
        })
    return rows
