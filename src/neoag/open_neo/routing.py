from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from neoag.controlled_execution.io_utils import load_limited_yaml, write_json, write_tsv
from neoag.skill_taxonomy.io import detect_delimiter, open_text

PATH_ALIASES: dict[str, tuple[str, ...]] = {
    "somatic_vcf": ("somatic_vcf", "variants_vcf", "tumor_vcf", "annotated_vcf", "vcf"),
    "fusion_tsv": ("fusion_tsv", "easyfuse_tsv", "easyfuse_pass_csv", "fusion_csv", "fusions"),
    "splice_junction_tsv": ("splice_junction_tsv", "regtools_tsv", "junctions", "splice_tsv"),
    "sv_vcf": ("sv_vcf", "sv_vcfs", "structural_variant_vcf"),
    "capture_bed": ("capture_bed", "capture_regions", "target_bed"),
    "peptide_csv": ("peptide_csv", "peptide_table", "raw_peptide_table"),
    "raw_events": ("raw_events",),
    "raw_peptides": ("raw_peptides",),
    "sv_raw_events": ("sv_raw_events",),
    "sv_raw_peptides": ("sv_raw_peptides",),
    "hla_file": ("hla_file", "hla", "hla_typing"),
    "expression_tsv": ("expression_tsv", "expression", "gene_expression"),
    "transcript_expression_tsv": ("transcript_expression_tsv", "transcript_expression"),
    "rna_evidence_tsv": ("rna_evidence_tsv", "rna_vaf_tsv", "rna_vaf"),
    "purity_tsv": ("purity_tsv", "purity"),
    "cnv_tsv": ("cnv_tsv", "cnv"),
    "hla_loh_tsv": ("hla_loh_tsv", "hla_loh"),
    "normal_expression": ("normal_expression",),
    "normal_hla_ligands": ("normal_hla_ligands",),
    "reference_proteome": ("reference_proteome", "normal_proteome"),
    "normal_junctions": ("normal_junctions",),
    "reference_fasta": ("reference_fasta", "fasta"),
    "gencode_gtf": ("gencode_gtf", "gtf"),
    "vep_cache": ("vep_cache",),
    "production_manifest": ("production_manifest",),
    "result_dir": ("result_dir", "results"),
    "comprehensive_evidence": ("comprehensive_evidence", "all_tool_results"),
    "weighted_baseline": ("weighted_baseline", "ranked_peptides"),
}

LIST_KEYS = {"sv_vcf"}


@dataclass
class Route:
    route: str
    skill: str
    status: str
    reason: str
    inputs: dict[str, Any] = field(default_factory=dict)


@dataclass
class RoutingResult:
    case_id: str
    sample_id: str
    genome_build: str
    inputs: dict[str, Any]
    inventory: list[dict[str, Any]]
    routes: list[Route]
    missing: list[dict[str, str]]
    ambiguous: list[dict[str, str]]
    warnings: list[str]

    @property
    def status(self) -> str:
        if self.ambiguous:
            return "BLOCKED"
        if not self.routes:
            return "BLOCKED"
        blocking_fields = {"hla_alleles_or_hla_file", "reference_fasta", "gencode_gtf"}
        if any(item.get("field") in blocking_fields or str(item.get("reason", "")).startswith("file not found:") for item in self.missing):
            return "BLOCKED"
        if self.missing:
            return "PARTIAL"
        return "PASS"


def load_sample_manifest(path: str | Path) -> dict[str, Any]:
    return load_limited_yaml(path)


def _first(mapping: dict[str, Any], names: tuple[str, ...]) -> Any:
    for name in names:
        if name in mapping:
            value = mapping[name]
            if value is not None and value != "" and value != []:
                return value
    return None


def _resolve_path(value: Any, base: Path) -> Any:
    if isinstance(value, list):
        return [str((base / Path(v)).resolve()) if not Path(str(v)).is_absolute() else str(Path(str(v))) for v in value]
    if not isinstance(value, str) or not value.strip():
        return value
    # HLA lists and non-path literals are left unchanged.
    if value.upper().startswith("HLA-") or "," in value and "HLA-" in value.upper():
        return value
    p = Path(value)
    if p.is_absolute():
        return str(p)
    return str((base / p).resolve())


def normalize_manifest(data: dict[str, Any], *, base_dir: str | Path = ".") -> dict[str, Any]:
    base = Path(base_dir).resolve()
    source: dict[str, Any] = {}
    for section in (data, data.get("inputs") or {}, data.get("tumor") or {}, data.get("normal") or {}):
        if isinstance(section, dict):
            source.update(section)
    out: dict[str, Any] = {}
    for canonical, aliases in PATH_ALIASES.items():
        value = _first(source, aliases)
        if value is not None:
            out[canonical] = _resolve_path(value, base)
    hla = source.get("hla_alleles") or data.get("hla_alleles") or []
    if isinstance(hla, str):
        hla = [x for x in re.split(r"[,;\s]+", hla) if x]
    out["hla_alleles"] = [str(x).strip() for x in hla if str(x).strip()]
    hla_file = out.get("hla_file")
    if not out["hla_alleles"] and hla_file and Path(str(hla_file)).is_file():
        text = Path(str(hla_file)).read_text(encoding="utf-8", errors="replace")
        alleles = re.findall(r"HLA-[A-Z0-9]+\*[0-9]{2,3}:[0-9]{2,3}", text, flags=re.IGNORECASE)
        if not alleles:
            alleles = [x for x in re.split(r"[,;\s]+", text) if re.match(r"^(?:HLA-)?[A-Z0-9]+\*[0-9]{2,3}:[0-9]{2,3}$", x, flags=re.IGNORECASE)]
        out["hla_alleles"] = [a.upper() if a.upper().startswith("HLA-") else "HLA-" + a.upper() for a in dict.fromkeys(alleles)]
    out["sample_id"] = str(data.get("sample_id") or data.get("case_id") or (data.get("sample") or {}).get("id") or "SAMPLE001")
    out["case_id"] = str(data.get("case_id") or data.get("patient_id") or out["sample_id"])
    out["genome_build"] = str(data.get("genome_build") or source.get("genome_build") or "GRCh38")
    execution = data.get("execution") or {}
    out["profile"] = str(execution.get("profile") or data.get("profile") or "default")
    out["backend"] = str(execution.get("backend") or "production-run")
    out["reuse_existing"] = bool(execution.get("reuse_existing", True))
    return out


def _table_header(path: Path) -> list[str]:
    try:
        delim = detect_delimiter(path)
        with open_text(path) as fh:
            first = fh.readline().rstrip("\n\r")
        return [x.strip() for x in first.split(delim)]
    except Exception:
        return []


def _vcf_metadata(path: Path) -> dict[str, Any]:
    samples: list[str] = []
    has_csq = False
    records = 0
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("##INFO=<ID=CSQ"):
                    has_csq = True
                elif line.startswith("#CHROM"):
                    samples = line.rstrip().split("	")[9:]
                elif not line.startswith("#"):
                    records += 1
                    if records >= 1000:
                        break
    except Exception:
        pass
    return {"vcf_samples": ",".join(samples), "has_vep_csq": has_csq, "records_scanned": records}


def build_inventory(inputs: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for key, value in inputs.items():
        if key in {"sample_id", "case_id", "genome_build", "profile", "backend", "reuse_existing", "hla_alleles"}:
            continue
        values = value if isinstance(value, list) else [value]
        for index, item in enumerate(values):
            if not isinstance(item, str) or not item:
                continue
            p = Path(item)
            row: dict[str, Any] = {
                "input_key": key,
                "index": index,
                "path": str(p),
                "exists": p.exists(),
                "is_file": p.is_file(),
                "size_bytes": p.stat().st_size if p.is_file() else 0,
                "header": "",
                "status": "OK" if p.exists() else "MISSING",
            }
            if p.is_file() and p.suffix.lower() in {".tsv", ".csv", ".txt"}:
                row["header"] = ",".join(_table_header(p))
            if p.is_file() and (p.suffix.lower() == ".vcf" or str(p).lower().endswith(".vcf.gz")):
                row.update(_vcf_metadata(p))
            rows.append(row)
    if inputs.get("hla_alleles"):
        rows.append({"input_key": "hla_alleles", "index": 0, "path": "", "exists": True, "is_file": False, "size_bytes": 0, "header": "", "status": "OK", "value": ",".join(inputs["hla_alleles"])})
    return rows


def _route(name: str, skill: str, reason: str, inputs: dict[str, Any], status: str = "READY") -> Route:
    return Route(name, skill, status, reason, inputs)


def route_inputs(inputs: dict[str, Any]) -> tuple[list[Route], list[dict[str, str]], list[dict[str, str]], list[str]]:
    routes: list[Route] = []
    missing: list[dict[str, str]] = []
    ambiguous: list[dict[str, str]] = []
    warnings: list[str] = []

    if inputs.get("production_manifest"):
        routes.append(_route("production_manifest", "neoag-production-run", "Explicit production manifest takes precedence", {"production_manifest": inputs["production_manifest"]}))
        return routes, missing, ambiguous, warnings
    if inputs.get("result_dir") or (inputs.get("comprehensive_evidence") and inputs.get("weighted_baseline")):
        routes.append(_route("existing_results", "neoag-ranking", "Reuse existing result/evidence outputs", {k: inputs.get(k) for k in ["result_dir", "comprehensive_evidence", "weighted_baseline"] if inputs.get(k)}))
        return routes, missing, ambiguous, warnings

    if inputs.get("somatic_vcf"):
        routes.append(_route("snv_indel", "neoag-vcf", "Somatic VCF detected", {"vcf": inputs["somatic_vcf"]}))
    if inputs.get("fusion_tsv"):
        routes.append(_route("fusion", "neoag-fusion", "Fusion caller table detected", {"fusion": inputs["fusion_tsv"]}))
    if inputs.get("splice_junction_tsv"):
        routes.append(_route("splice", "neoag-splice", "Splice junction table detected", {"junctions": inputs["splice_junction_tsv"], "vcf": inputs.get("somatic_vcf", "")}))
    if inputs.get("peptide_csv"):
        routes.append(_route("peptide_only", "neoag-peptide-csv", "Existing peptide-HLA table detected", {"peptide_csv": inputs["peptide_csv"]}))
    if inputs.get("raw_events") and inputs.get("raw_peptides"):
        routes.append(_route("intermediates", "neoag-run-full", "Standard raw intermediates detected", {"raw_events": inputs["raw_events"], "raw_peptides": inputs["raw_peptides"]}))
    if inputs.get("sv_raw_events") and inputs.get("sv_raw_peptides"):
        routes.append(_route("sv_intermediates", "neoag-sv-wgs", "Prebuilt SV raw tables detected", {"sv_raw_events": inputs["sv_raw_events"], "sv_raw_peptides": inputs["sv_raw_peptides"]}))
    elif inputs.get("sv_vcf"):
        if inputs.get("capture_bed"):
            skill = "neoag-sv-wes"
            route = "sv_wes"
        else:
            skill = "neoag-sv-wgs"
            route = "sv_wgs"
        required = ["reference_fasta", "gencode_gtf"]
        absent = [k for k in required if not inputs.get(k)]
        status = "BLOCKED" if absent else "READY"
        routes.append(_route(route, skill, "DNA SV VCF detected", {"sv_vcf": inputs["sv_vcf"], "capture_bed": inputs.get("capture_bed", "")}, status=status))
        for key in absent:
            missing.append({"field": key, "reason": f"required for {route}"})

    if not routes:
        ambiguous.append({"field": "inputs", "reason": "No supported VCF/fusion/splice/SV/peptide/intermediate/result input detected"})
    hla_required = any(r.route in {"snv_indel", "fusion", "splice", "peptide_only", "sv_wgs", "sv_wes", "sv_intermediates"} for r in routes)
    if hla_required and not (inputs.get("hla_alleles") or inputs.get("hla_file")):
        missing.append({"field": "hla_alleles_or_hla_file", "reason": "required for peptide-HLA prediction"})
    if inputs.get("genome_build", "GRCh38").upper() not in {"GRCH38", "HG38"}:
        warnings.append("Non-GRCh38 build detected; all references and variants must use the same build")
    if inputs.get("fusion_tsv") and not inputs.get("rna_evidence_tsv"):
        warnings.append("Fusion table present without separate RNA evidence table; caller junction fields will be used when available")
    return routes, missing, ambiguous, warnings


def inspect_manifest(path: str | Path) -> RoutingResult:
    p = Path(path).resolve()
    data = load_sample_manifest(p)
    inputs = normalize_manifest(data, base_dir=p.parent)
    inventory = build_inventory(inputs)
    routes, missing, ambiguous, warnings = route_inputs(inputs)
    for row in inventory:
        if row.get("status") == "MISSING":
            missing.append({"field": str(row.get("input_key")), "reason": f"file not found: {row.get('path')}"})
    return RoutingResult(inputs["case_id"], inputs["sample_id"], inputs["genome_build"], inputs, inventory, routes, missing, ambiguous, warnings)


def write_routing_outputs(result: RoutingResult, outdir: str | Path) -> dict[str, str]:
    od = Path(outdir)
    od.mkdir(parents=True, exist_ok=True)
    route_rows = [r.__dict__ for r in result.routes]
    write_tsv(od / "input_inventory.tsv", result.inventory)
    write_tsv(od / "missing_inputs.tsv", result.missing)
    write_tsv(od / "ambiguous_inputs.tsv", result.ambiguous)
    write_json(od / "route_plan.json", {
        "status": result.status,
        "case_id": result.case_id,
        "sample_id": result.sample_id,
        "genome_build": result.genome_build,
        "routes": route_rows,
        "warnings": result.warnings,
        "missing_inputs": result.missing,
        "ambiguous_inputs": result.ambiguous,
    })
    write_json(od / "input_status.json", {
        "status": result.status,
        "inputs": result.inputs,
        "inventory": result.inventory,
        "missing_inputs": result.missing,
        "ambiguous_inputs": result.ambiguous,
    })
    return {
        "input_inventory": str(od / "input_inventory.tsv"),
        "input_status": str(od / "input_status.json"),
        "missing_inputs": str(od / "missing_inputs.tsv"),
        "ambiguous_inputs": str(od / "ambiguous_inputs.tsv"),
        "route_plan": str(od / "route_plan.json"),
    }
