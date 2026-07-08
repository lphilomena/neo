"""Validation assay and long-peptide / minigene design by neoantigen class."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from .utils import read_tsv

VALIDATION_PLAN_FIELDS = [
    "peptide_id",
    "event_id",
    "gene",
    "peptide",
    "wildtype_peptide",
    "hla_allele",
    "event_type",
    "peptide_consequence",
    "priority",
    "safety_status",
    "presentation_grade",
    "appm_multiplier",
    "ccf_multiplier",
    "validation_mode",
    "validation_strategy",
    "recommended_assay",
    "short_peptide_mt",
    "short_peptide_wt",
    "long_peptide",
    "minigene",
    "minigene_nt",
    "validation_notes",
    "required_followup",
]

_LONG_MODES = frozenset({
    "frameshift_long",
    "splice_junction_long",
    "fusion_junction_long",
    "insertion_long",
    "mhc_ii_long",
})


def _norm(s: Any) -> str:
    return str(s or "").strip()


def _upper_pep(s: Any) -> str:
    return _norm(s).upper()


def minigene_to_long_peptide(minigene: str) -> str:
    """Collapse minigene segments (flank|center|tail) into a single synthesis sequence."""
    mg = _norm(minigene)
    if not mg:
        return ""
    if "|" in mg:
        return "".join(part for part in mg.split("|") if part)
    return mg


def resolve_peptide_catalog(
    peptide_catalog_tsv: str | Path | None = None,
    *,
    outdir: str | Path | None = None,
) -> str | None:
    """Resolve variant/easyfuse peptide catalog for minigene lookup."""
    candidates: list[Path] = []
    if peptide_catalog_tsv:
        candidates.append(Path(peptide_catalog_tsv))
    if outdir:
        tools = Path(outdir) / "upstream" / "tools"
        candidates.extend([
            tools / "variant_peptides.annotated.tsv",
            tools / "variant_peptides.tsv",
            tools / "easyfuse_variant_peptides.tsv",
            tools / "extra_variant_peptides.tsv",
        ])
    for path in candidates:
        if path.is_file():
            return str(path)
    return None


def build_peptide_catalog_index(rows: list[Mapping[str, Any]]) -> dict[tuple[str, str], dict[str, str]]:
    """Index catalog rows by (event_key, mutant_peptide)."""
    index: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        pep = _upper_pep(row.get("mutant_peptide") or row.get("peptide"))
        if not pep:
            continue
        for key in {
            _norm(row.get("variant_key")),
            _norm(row.get("event_id")),
        }:
            if key:
                index[(key, pep)] = dict(row)
    return index


def lookup_peptide_catalog(
    peptide: Mapping[str, Any],
    catalog_index: Mapping[tuple[str, str], Mapping[str, Any]],
) -> dict[str, str]:
    if not catalog_index:
        return {}
    pep = _upper_pep(peptide.get("peptide"))
    event_id = _norm(peptide.get("event_id"))
    for key in (event_id,):
        hit = catalog_index.get((key, pep))
        if hit:
            return dict(hit)
    return {}


def classify_validation_mode(peptide: Mapping[str, Any]) -> str:
    priority = _norm(peptide.get("final_priority")).upper()
    if priority == "D":
        return "do_not_advance"
    if "CAUTION" in priority and _norm(peptide.get("safety_status")).upper() != "PASS":
        return "safety_caution"

    pc = _norm(peptide.get("peptide_consequence")).lower()
    et = _norm(peptide.get("event_type")).upper()
    mhc = _norm(peptide.get("mhc_class")).upper()
    crosses = _norm(peptide.get("crosses_junction")).lower() == "yes"

    if et == "FUSION" or pc == "fusion" or (crosses and et == "FUSION"):
        return "fusion_junction_long"
    if pc == "frameshift" or "frameshift" in pc:
        return "frameshift_long"
    if pc in {"splice_junction", "exon_deletion_junction"} or "splice" in pc:
        return "splice_junction_long"
    if pc in {"insertion", "exon_deletion_junction"}:
        return "insertion_long"
    if mhc in {"II", "MHC-II", "CLASSII"}:
        return "mhc_ii_long"
    return "missense_short_pair"


def validation_strategy_text(mode: str) -> str:
    return {
        "missense_short_pair": (
            "Mutant short peptide (8-11mer) with matched WT short peptide control; "
            "MHC-I ELISpot or tetramer"
        ),
        "frameshift_long": (
            "Novel C-terminal tail: pooled long peptides (15-27aa) and/or frameshift minigene "
            "transfection; short peptide optional secondary"
        ),
        "splice_junction_long": (
            "Abnormal splice/exon-junction long peptide (15-27aa) and/or splice minigene; "
            "do not rely on short peptide alone"
        ),
        "insertion_long": (
            "Inframe insertion/deletion junction long peptide and/or minigene centered on altered region"
        ),
        "fusion_junction_long": (
            "Fusion junction long peptide spanning breakpoint and/or fusion minigene transfection"
        ),
        "mhc_ii_long": (
            "CD4 long peptide (15-30aa) stimulation with cytokine readout"
        ),
        "safety_caution": "Safety-focused validation before efficacy assay",
        "do_not_advance": "Do not advance",
    }.get(mode, "MHC-I peptide ELISpot/tetramer")


def recommended_assay_text(mode: str) -> str:
    return {
        "missense_short_pair": "MHC-I short peptide ELISpot/tetramer (MT + WT pair)",
        "frameshift_long": "Frameshift minigene transfection or novel-tail long-peptide pool + ELISpot",
        "splice_junction_long": "Splice/exon-junction long peptide or minigene + ELISpot",
        "insertion_long": "Inframe indel junction long peptide or minigene + ELISpot",
        "fusion_junction_long": "Fusion junction long peptide or fusion minigene + ELISpot",
        "mhc_ii_long": "CD4 long peptide stimulation + cytokine assay",
        "safety_caution": "Safety-focused validation before efficacy assay",
        "do_not_advance": "Do not advance",
    }.get(mode, "MHC-I peptide ELISpot/tetramer")


def _fallback_peptide_centered_minigene(peptide: str, total_len: int = 27) -> tuple[str, str]:
    pep = _upper_pep(peptide)
    if not pep or len(pep) >= total_len:
        return (f"|{pep}|" if pep else "", "")
    pad = total_len - len(pep)
    left = pad // 2
    right = pad - left
    return f"{'X' * left}|{pep}|{'X' * right}", ""


def design_validation_row(
    peptide: Mapping[str, Any],
    *,
    catalog_row: Mapping[str, Any] | None = None,
) -> dict[str, str]:
    mode = classify_validation_mode(peptide)
    mt = _upper_pep(peptide.get("peptide"))
    wt = _upper_pep(peptide.get("wildtype_peptide"))
    catalog_row = catalog_row or {}

    minigene = _norm(catalog_row.get("minigene") or peptide.get("minigene"))
    minigene_nt = _norm(catalog_row.get("minigene_nt") or peptide.get("minigene_nt"))
    notes: list[str] = []

    if mode in _LONG_MODES and not minigene and mt:
        minigene, minigene_nt = _fallback_peptide_centered_minigene(mt)
        notes.append("minigene_placeholder: pass --variant-peptides for sequence-accurate design")

    long_peptide = minigene_to_long_peptide(minigene) if mode in _LONG_MODES else ""

    if mode == "frameshift_long":
        notes.append("prioritize novel-tail long peptide pool or frameshift minigene over short peptide alone")
    elif mode == "splice_junction_long":
        notes.append("validate abnormal junction processing; short peptide insufficient")
    elif mode == "fusion_junction_long":
        notes.append("include breakpoint-spanning long peptide or fusion minigene")
    elif mode == "missense_short_pair" and not wt:
        notes.append("WT short peptide control missing; derive from reference transcript")

    if catalog_row:
        notes.append("minigene_from_peptide_catalog")

    return {
        "peptide_id": _norm(peptide.get("peptide_id")),
        "event_id": _norm(peptide.get("event_id")),
        "gene": _norm(peptide.get("gene")),
        "peptide": mt,
        "wildtype_peptide": wt,
        "hla_allele": _norm(peptide.get("hla_allele")),
        "event_type": _norm(peptide.get("event_type")),
        "peptide_consequence": _norm(peptide.get("peptide_consequence")),
        "priority": _norm(peptide.get("final_priority")),
        "safety_status": _norm(peptide.get("safety_status")),
        "presentation_grade": _norm(peptide.get("presentation_evidence_grade")),
        "appm_multiplier": _norm(peptide.get("appm_multiplier")),
        "ccf_multiplier": _norm(peptide.get("ccf_multiplier")),
        "validation_mode": mode,
        "validation_strategy": validation_strategy_text(mode),
        "recommended_assay": recommended_assay_text(mode),
        "short_peptide_mt": mt if mode != "do_not_advance" else "",
        "short_peptide_wt": wt if mode == "missense_short_pair" else "",
        "long_peptide": long_peptide,
        "minigene": minigene,
        "minigene_nt": minigene_nt,
        "validation_notes": "; ".join(notes),
        "required_followup": _norm(peptide.get("recommended_use")),
    }


def load_peptide_catalog_index(peptide_catalog_tsv: str | Path | None) -> dict[tuple[str, str], dict[str, str]]:
    if not peptide_catalog_tsv:
        return {}
    path = Path(peptide_catalog_tsv)
    if not path.is_file():
        return {}
    return build_peptide_catalog_index(read_tsv(path))
