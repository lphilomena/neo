"""Extract APPM mutation/expression inputs from VCF, VEP TSV, or MAF.

The APPM mutation layer is intentionally gene-centric. It preserves the legacy
``gene``/``consequence`` columns while adding enough VEP/MAF detail to separate
protein-truncating, splice, start-lost, and damaging-missense evidence.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Mapping

from ..utils import first, open_text_maybe_gz, write_tsv

APPM_MUTATION_GENES = {
    "B2M",
    "HLA-A",
    "HLA-B",
    "HLA-C",
    "TAP1",
    "TAP2",
    "TAPBP",
    "JAK1",
    "JAK2",
    "STAT1",
    "NLRC5",
    "CIITA",
    "RFX5",
    "RFXANK",
    "RFXAP",
}

# Backwards-compatible public name used by older tests/imports.
APPM_GENES = APPM_MUTATION_GENES

PTV_TERMS = {
    "frameshift_variant",
    "stop_gained",
    "protein_truncating_variant",
    "transcript_ablation",
    "feature_truncation",
}
SPLICE_TERMS = {"splice_donor_variant", "splice_acceptor_variant"}
START_LOST_TERMS = {"start_lost"}
DAMAGING_TERMS = PTV_TERMS | SPLICE_TERMS | START_LOST_TERMS | {"stop_lost", "loss_of_function"}
MISSENSE_TERMS = {"missense_variant"}

MAF_TO_VEP_CONSEQUENCE = {
    "Frame_Shift_Del": "frameshift_variant",
    "Frame_Shift_Ins": "frameshift_variant",
    "Nonsense_Mutation": "stop_gained",
    "Nonstop_Mutation": "stop_lost",
    "Splice_Site": "splice_donor_variant",
    "Translation_Start_Site": "start_lost",
    "Missense_Mutation": "missense_variant",
    "In_Frame_Del": "inframe_deletion",
    "In_Frame_Ins": "inframe_insertion",
    "Silent": "synonymous_variant",
}


def _normalize_consequence(value: str) -> str:
    raw = str(value or "").strip()
    return MAF_TO_VEP_CONSEQUENCE.get(raw, raw)

VEP_APPM_FIELDS = [
    "gene",
    "consequence",
    "mutation_consequence",
    "damaging_variant",
    "damaging_class",
    "is_ptv",
    "is_splice",
    "is_start_lost",
    "is_damaging_missense",
    "impact",
    "sift",
    "polyphen",
    "variant_id",
    "location",
    "allele",
    "gene_id",
    "feature",
    "feature_type",
    "hgvsc",
    "hgvsp",
    "protein_position",
    "amino_acids",
    "codons",
    "canonical",
    "variant_class",
    "source_format",
]


def _open_vcf(path: Path):
    return open_text_maybe_gz(path)


def _norm_gene(value: str) -> str:
    return str(value or "").strip().upper()


def _split_terms(value: str) -> list[str]:
    terms: list[str] = []
    for chunk in str(value or "").replace(",", "&").split("&"):
        t = chunk.strip()
        if t and t not in terms:
            terms.append(t)
    return terms


def _is_damaging_missense(consequence: str, impact: str = "", sift: str = "", polyphen: str = "") -> bool:
    terms = set(_split_terms(consequence))
    if not terms & MISSENSE_TERMS:
        return False
    sift_l = str(sift or "").lower()
    poly_l = str(polyphen or "").lower()
    impact_u = str(impact or "").upper()
    return (
        "deleterious" in sift_l
        or "damaging" in poly_l
        or "probably_damaging" in poly_l
        or "possibly_damaging" in poly_l
        or impact_u == "HIGH"
    )


def _classify(consequence: str, impact: str = "", sift: str = "", polyphen: str = "") -> dict[str, str]:
    terms = set(_split_terms(consequence))
    is_ptv = bool(terms & PTV_TERMS)
    is_splice = bool(terms & SPLICE_TERMS)
    is_start_lost = bool(terms & START_LOST_TERMS)
    is_dmg_missense = _is_damaging_missense(consequence, impact, sift, polyphen)
    classes: list[str] = []
    if is_ptv:
        classes.append("protein_truncating_variant")
    if is_splice:
        classes.append("splice_disrupting_variant")
    if is_start_lost:
        classes.append("start_lost")
    if is_dmg_missense:
        classes.append("damaging_missense")
    damaging = bool(classes or (terms & DAMAGING_TERMS))
    if damaging and not classes:
        classes.append("damaging_variant")
    return {
        "damaging_variant": "yes" if damaging else "no",
        "damaging_class": ";".join(classes) if classes else "none",
        "is_ptv": "yes" if is_ptv else "no",
        "is_splice": "yes" if is_splice else "no",
        "is_start_lost": "yes" if is_start_lost else "no",
        "is_damaging_missense": "yes" if is_dmg_missense else "no",
    }


def _csq_fields(header_line: str) -> list[str] | None:
    marker = "Format: "
    if marker not in header_line:
        return None
    text = header_line.split(marker, 1)[1].strip()
    text = text.rstrip('>").')
    return text.split("|")


def _csq_indices(header_line: str) -> tuple[int, int] | None:
    fields = _csq_fields(header_line)
    if not fields:
        return None
    return fields.index("SYMBOL"), fields.index("Consequence")


def _row_from_vep_like(row: Mapping[str, str], *, source_format: str) -> dict[str, str] | None:
    gene = _norm_gene(first(row, ["SYMBOL", "Hugo_Symbol", "gene", "Gene", "symbol"], ""))
    if gene not in APPM_MUTATION_GENES:
        return None
    consequence = _normalize_consequence(first(row, ["Consequence", "Variant_Classification", "consequence", "effect", "One_Consequence"], ""))
    impact = first(row, ["IMPACT", "impact"], "")
    sift = first(row, ["SIFT", "sift"], "")
    polyphen = first(row, ["PolyPhen", "polyphen"], "")
    cls = _classify(consequence, impact, sift, polyphen)
    mutation_consequence = consequence
    if cls["is_damaging_missense"] == "yes" and "damaging_missense" not in mutation_consequence:
        mutation_consequence = f"{mutation_consequence}&damaging_missense" if mutation_consequence else "damaging_missense"
    return {
        "gene": gene,
        "consequence": consequence,
        "mutation_consequence": mutation_consequence,
        **cls,
        "impact": impact,
        "sift": sift,
        "polyphen": polyphen,
        "variant_id": first(row, ["#Uploaded_variation", "Uploaded_variation", "Mutation_ID", "variant_id", "ID"], ""),
        "location": first(row, ["Location", "Chromosome", "chromosome", "Chromosome_Start"], ""),
        "allele": first(row, ["Allele", "Tumor_Seq_Allele2", "alt", "ALT"], ""),
        "gene_id": first(row, ["Gene", "Entrez_Gene_Id", "gene_id"], ""),
        "feature": first(row, ["Feature", "Transcript_ID", "transcript_id"], ""),
        "feature_type": first(row, ["Feature_type", "feature_type"], ""),
        "hgvsc": first(row, ["HGVSc", "HGVSc_VEP", "cDNA_change"], ""),
        "hgvsp": first(row, ["HGVSp", "HGVSp_Short", "Protein_Change", "protein_change"], ""),
        "protein_position": first(row, ["Protein_position", "Protein_position_VEP"], ""),
        "amino_acids": first(row, ["Amino_acids", "Amino_Acid_Change"], ""),
        "codons": first(row, ["Codons"], ""),
        "canonical": first(row, ["CANONICAL", "canonical"], ""),
        "variant_class": first(row, ["VARIANT_CLASS", "Variant_Type", "variant_class"], ""),
        "source_format": source_format,
    }


def _dedup_rows(rows: Iterable[dict[str, str]]) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    seen: set[tuple[str, str, str, str, str]] = set()
    for r in rows:
        key = (r.get("gene", ""), r.get("variant_id", ""), r.get("feature", ""), r.get("allele", ""), r.get("consequence", ""))
        if key in seen:
            continue
        seen.add(key)
        out.append(r)
    return out


def _read_vep_table(path: Path) -> list[dict[str, str]]:
    header: list[str] | None = None
    rows: list[dict[str, str]] = []
    with open_text_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                header = line.lstrip("#").split("\t")
                continue
            if header is None:
                header = line.split("\t")
                continue
            vals = line.split("\t")
            rows.append(dict(zip(header, vals)))
    return rows


def extract_gx_expression(
    vcf_path: str | Path,
    *,
    tumor_sample: str | None = None,
) -> dict[str, float]:
    """Parse FORMAT/GX (gene|TPM pairs) from tumor sample column."""
    path = Path(vcf_path)
    expr: dict[str, float] = {}
    sample_idx: int | None = None
    with _open_vcf(path) as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                header = line.rstrip().split("\t")
                if tumor_sample and tumor_sample in header:
                    sample_idx = header.index(tumor_sample)
                elif len(header) > 10:
                    sample_idx = 10
                continue
            if line.startswith("#"):
                continue
            if sample_idx is None:
                continue
            parts = line.rstrip().split("\t")
            if len(parts) <= sample_idx:
                continue
            fmt = parts[8].split(":")
            if "GX" not in fmt:
                continue
            gx_idx = fmt.index("GX")
            vals = parts[sample_idx].split(":")
            if len(vals) <= gx_idx:
                continue
            gx = vals[gx_idx]
            if not gx or gx == ".":
                continue
            for item in gx.split(","):
                if "|" not in item:
                    continue
                gene, raw_tpm = item.split("|", 1)
                gene = gene.strip()
                if not gene:
                    continue
                try:
                    tpm = float(raw_tpm)
                except ValueError:
                    continue
                expr[gene] = max(expr.get(gene, 0.0), tpm)
    return expr


def extract_appm_variant_rows(path: str | Path) -> list[dict[str, str]]:
    """Collect APPM mutation rows from VCF CSQ, VEP tabular TSV, or MAF."""
    p = Path(path)
    first_nonempty = ""
    with open_text_maybe_gz(p) as fh:
        for line in fh:
            if line.strip():
                first_nonempty = line
                break
    is_vcf = first_nonempty.startswith("##") or first_nonempty.startswith("#CHROM")
    # VEP tabular also starts with ## comments, so detect the #Uploaded header.
    has_vep_table_header = False
    if first_nonempty.startswith("##"):
        with open_text_maybe_gz(p) as fh:
            for line in fh:
                if line.startswith("#Uploaded_variation"):
                    has_vep_table_header = True
                    break
    if not is_vcf or has_vep_table_header:
        rows = [_row_from_vep_like(r, source_format="vep_or_maf_tsv") for r in _read_vep_table(p)]
        return _dedup_rows([r for r in rows if r])

    csq_fields: list[str] | None = None
    out: list[dict[str, str]] = []
    with _open_vcf(p) as fh:
        for line in fh:
            if line.startswith("##INFO=<ID=CSQ"):
                csq_fields = _csq_fields(line)
            if line.startswith("#"):
                continue
            if not csq_fields:
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 8:
                continue
            info = parts[7]
            if "CSQ=" not in info:
                continue
            csq = info.split("CSQ=", 1)[1].split(";", 1)[0]
            for alt_entry in csq.split(","):
                fields = alt_entry.split("|")
                row = {k: fields[i] if i < len(fields) else "" for i, k in enumerate(csq_fields)}
                row["#Uploaded_variation"] = parts[2] if parts[2] != "." else f"{parts[0]}_{parts[1]}_{parts[3]}/{parts[4]}"
                row["Location"] = f"{parts[0]}:{parts[1]}"
                parsed = _row_from_vep_like(row, source_format="vcf_csq")
                if parsed:
                    out.append(parsed)
    return _dedup_rows(out)


def extract_appm_variants(vcf_path: str | Path) -> dict[str, list[str]]:
    """Collect consequences for APPM mutation genes."""
    genes: dict[str, list[str]] = {}
    for row in extract_appm_variant_rows(vcf_path):
        gene = row.get("gene", "")
        cons = row.get("mutation_consequence") or row.get("consequence") or ""
        if not gene or not cons:
            continue
        bucket = genes.setdefault(gene, [])
        for term in _split_terms(cons):
            if term not in bucket:
                bucket.append(term)
    return genes


def write_expression_tsv(expr: dict[str, float], out_path: str | Path) -> None:
    rows = [{"gene": g, "TPM": f"{expr[g]:.6f}"} for g in sorted(expr)]
    write_tsv(out_path, rows, ["gene", "TPM"])


def _no_variant_row(gene: str, source_format: str) -> dict[str, str]:
    row = {field: "" for field in VEP_APPM_FIELDS}
    row.update({
        "gene": gene,
        "consequence": "",
        "mutation_consequence": "",
        "damaging_variant": "no",
        "damaging_class": "none",
        "is_ptv": "no",
        "is_splice": "no",
        "is_start_lost": "no",
        "is_damaging_missense": "no",
        "source_format": source_format,
    })
    return row


def write_vep_appm(path: str | Path, out_path: str | Path, *, include_no_variant_genes: bool = True) -> None:
    rows = extract_appm_variant_rows(path)
    if include_no_variant_genes:
        present = {r.get("gene", "") for r in rows}
        for gene in sorted(APPM_MUTATION_GENES - present):
            rows.append(_no_variant_row(gene, "assessed_no_variant"))
    write_tsv(out_path, rows, VEP_APPM_FIELDS)


def write_vep_appm_from_vcf(vcf_path: str | Path, out_path: str | Path) -> None:
    write_vep_appm(vcf_path, out_path)


def extract_appm_inputs_from_vcf(
    vcf_path: str | Path,
    outdir: str | Path,
    *,
    tumor_sample: str | None = None,
) -> dict[str, str]:
    """Write gene_expression.tsv and vep_appm.tsv derived from annotated VCF/TSV/MAF."""
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    expr_path = out / "gene_expression.tsv"
    vep_path = out / "vep_appm.tsv"
    try:
        expr = extract_gx_expression(vcf_path, tumor_sample=tumor_sample)
    except Exception:
        expr = {}
    write_expression_tsv(expr, expr_path)
    write_vep_appm(vcf_path, vep_path)
    return {"expression": str(expr_path), "vep_appm": str(vep_path)}


def summarize_appm_risks(
    expr: dict[str, float],
    variants: dict[str, list[str]],
    *,
    low_tpm: float = 1.0,
) -> dict[str, str]:
    """Quick summary for logging / reports."""
    out: dict[str, str] = {}
    for gene in sorted(APPM_MUTATION_GENES):
        cons = variants.get(gene, [])
        damaging = any(c in DAMAGING_TERMS or c == "damaging_missense" for c in cons)
        tpm = expr.get(gene)
        if damaging:
            out[gene] = "damaging_variant"
        elif tpm is not None and tpm < low_tpm:
            out[f"{gene}_low_expr"] = f"{tpm:.2f}"
    return out
