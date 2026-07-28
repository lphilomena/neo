from __future__ import annotations

import gzip
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from neoag.controlled_execution.io_utils import load_limited_yaml, write_json, write_tsv
from neoag.skill_taxonomy.io import detect_delimiter, open_text

PATH_ALIASES: dict[str, tuple[str, ...]] = {
    "tumor_dna_bam": ("tumor_dna_bam", "tumor_bam", "tumor_wgs_bam", "tumor_wes_bam"),
    "normal_dna_bam": ("normal_dna_bam", "normal_bam", "blood_bam", "matched_normal_bam"),
    "tumor_rna_bam": ("tumor_rna_bam", "rna_bam", "wts_bam"),
    "tumor_dna_fastq": ("tumor_dna_fastq", "tumor_dna_fastqs"),
    "normal_dna_fastq": ("normal_dna_fastq", "normal_dna_fastqs"),
    "tumor_rna_fastq": ("tumor_rna_fastq", "tumor_rna_fastqs", "rna_fastqs"),
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
    "salmon_index": ("salmon_index",),
    "tx2gene": ("tx2gene", "transcript_to_gene"),
    "rsem_reference": ("rsem_reference",),
    "star_index": ("star_index", "star_genome_index"),
    "ctat_genome_lib": ("ctat_genome_lib", "star_fusion_reference"),
    "easyfuse_ref": ("easyfuse_ref", "easyfuse_reference"),
    "normal_readthrough": ("normal_readthrough", "normal_readthrough_db"),
    "snaf_workflow": ("snaf_workflow",),
    "snaf_db": ("snaf_db", "snaf_reference"),
    "snaf_python": ("snaf_python",),
    "altanalyze_image": ("altanalyze_image",),
    "splicemutr_workflow": ("splicemutr_workflow",),
    "production_manifest": ("production_manifest",),
    "result_dir": ("result_dir", "results"),
    "comprehensive_evidence": ("comprehensive_evidence", "all_tool_results"),
    "weighted_baseline": ("weighted_baseline", "ranked_peptides"),
    "input_dir": ("input_dir", "data_dir"),
}

LIST_KEYS = {"sv_vcf", "tumor_dna_fastq", "normal_dna_fastq", "tumor_rna_fastq"}
HLA_RE = re.compile(r"^(?:HLA-)?(?:A|B|C|DRB1|DQB1|DPB1|DQA1|DPA1)\*[0-9]{2,3}(?::[0-9A-Z]{2,3}){1,4}$", re.I)
FUSION_HEADER_GROUPS = ({"fusionname"}, {"fusion_gene"}, {"gene1", "gene2"}, {"gene5", "gene3"}, {"leftgene", "rightgene"})
SPLICE_HEADER_GROUPS = ({"chrom", "start", "end"}, {"chromosome", "intron_start", "intron_end"}, {"junction", "read_count"})
PEPTIDE_HEADERS = {"peptide", "hla_allele"}


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
    manifest_path: str = ""

    @property
    def status(self) -> str:
        if self.ambiguous:
            return "BLOCKED"
        if not self.routes:
            return "BLOCKED"
        blocking_fields = {"hla_alleles_or_hla_file", "reference_fasta", "gencode_gtf"}
        if any(
            item.get("field") in blocking_fields
            or str(item.get("reason", "")).startswith(("file not found:", "input status ", "BAM index not found", "required for capture-limited"))
            for item in self.missing
        ):
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
    out["execution_mode"] = str(execution.get("mode") or "")
    out["tumor_type"] = str(data.get("tumor_type") or data.get("disease") or "")
    out["tumor_sample_id"] = str(data.get("tumor_sample_id") or (data.get("tumor") or {}).get("sample_id") or out["sample_id"])
    out["normal_sample_id"] = str(data.get("normal_sample_id") or (data.get("normal") or {}).get("sample_id") or "")
    out["assay_type"] = str(data.get("assay_type") or source.get("assay_type") or "").upper()
    out["rna_threads"] = int(source.get("rna_threads") or execution.get("threads") or 16)
    raw_tool_results = data.get("tool_results") if isinstance(data.get("tool_results"), dict) else {}
    out["tool_results"] = {
        str(domain): {
            str(tool): _resolve_path(value, base)
            for tool, value in values.items()
            if isinstance(value, str) and value
        }
        for domain, values in raw_tool_results.items()
        if isinstance(values, dict)
    }
    return out


def _load_hla_alleles(path: str | Path) -> list[str]:
    p = Path(path)
    if not p.is_file():
        return []
    text = p.read_text(encoding="utf-8", errors="replace")
    alleles = re.findall(r"(?:HLA-)?[A-Z0-9]+\*[0-9]{2,3}(?::[0-9A-Z]{2,3}){1,4}", text, flags=re.IGNORECASE)
    normalized = [a.upper() if a.upper().startswith("HLA-") else "HLA-" + a.upper() for a in alleles]
    return list(dict.fromkeys(normalized))


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
    contigs: list[str] = []
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("##INFO=<ID=CSQ"):
                    has_csq = True
                elif line.startswith("##contig=<ID=") and len(contigs) < 30:
                    contigs.append(line.split("ID=", 1)[1].split(",", 1)[0].split(">", 1)[0])
                elif line.startswith("#CHROM"):
                    samples = line.rstrip().split("	")[9:]
                elif not line.startswith("#"):
                    records += 1
                    if records >= 1000:
                        break
    except Exception:
        pass
    build = "UNASSESSED"
    text = ",".join(contigs).lower()
    if "grch38" in text or "hg38" in text:
        build = "GRCh38"
    elif "grch37" in text or "hg19" in text:
        build = "GRCh37"
    return {"vcf_samples": ",".join(samples), "vcf_sample_count": len(samples), "has_vep_csq": has_csq, "records_scanned": records, "detected_genome_build": build}


def _header_matches(header: list[str], groups: tuple[set[str], ...] | set[str]) -> bool:
    normalized = {re.sub(r"[^a-z0-9]+", "", value.lower()) for value in header}
    if isinstance(groups, set):
        expected = {re.sub(r"[^a-z0-9]+", "", value.lower()) for value in groups}
        return expected <= normalized
    return any({re.sub(r"[^a-z0-9]+", "", value.lower()) for value in group} <= normalized for group in groups)


def _classify_file(path: Path) -> list[str]:
    low = path.name.lower()
    if low.endswith((".vcf", ".vcf.gz")):
        try:
            opener = gzip.open if low.endswith(".gz") else open
            with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
                prefix = "".join(line for _, line in zip(range(250), handle))
            return ["sv_vcf"] if "SVTYPE=" in prefix or "##ALT=<ID=BND" in prefix else ["somatic_vcf"]
        except Exception:
            return ["somatic_vcf"]
    if low.endswith(".bam"):
        if "rna" in low or "wts" in low or "transcript" in low:
            return ["tumor_rna_bam"]
        if any(token in low for token in ("normal", "blood", "germline")):
            return ["normal_dna_bam"]
        if any(token in low for token in ("tumor", "tumour", "somatic")):
            return ["tumor_dna_bam"]
        return []
    if low.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz")):
        if "rna" in low or "wts" in low:
            return ["tumor_rna_fastq"]
        if any(token in low for token in ("normal", "blood", "germline")):
            return ["normal_dna_fastq"]
        if any(token in low for token in ("tumor", "tumour", "somatic")):
            return ["tumor_dna_fastq"]
        return []
    if low.endswith((".tsv", ".csv", ".txt")):
        header = _table_header(path)
        if _header_matches(header, FUSION_HEADER_GROUPS):
            return ["fusion_tsv"]
        if _header_matches(header, SPLICE_HEADER_GROUPS):
            return ["splice_junction_tsv"]
        if _header_matches(header, PEPTIDE_HEADERS):
            return ["peptide_csv"]
        text = path.read_text(encoding="utf-8", errors="replace")[:10000]
        if re.search(r"(?:HLA-)?(?:A|B|C)\*[0-9]{2,3}:[0-9A-Z]{2,3}", text, re.I):
            return ["hla_file"]
    return []


def scan_input_directory(path: str | Path, *, max_files: int = 1000) -> tuple[dict[str, Any], list[dict[str, str]], list[str]]:
    root = Path(path).resolve()
    candidates: dict[str, list[str]] = {}
    warnings: list[str] = []
    if not root.is_dir():
        return {}, [{"field": "input_dir", "reason": f"directory not found: {root}"}], warnings
    for index, file_path in enumerate(sorted(p for p in root.rglob("*") if p.is_file())):
        if index >= max_files:
            warnings.append(f"input directory scan capped at {max_files} files")
            break
        for key in _classify_file(file_path):
            candidates.setdefault(key, []).append(str(file_path))
    detected: dict[str, Any] = {}
    ambiguous: list[dict[str, str]] = []
    for key, values in candidates.items():
        unique = sorted(set(values))
        if key in LIST_KEYS:
            detected[key] = unique
        elif len(unique) == 1:
            detected[key] = unique[0]
        else:
            ambiguous.append({"field": key, "reason": "multiple directory-scan candidates: " + ";".join(unique[:20])})
    return detected, ambiguous, warnings


def build_inventory(inputs: dict[str, Any], *, output_dir: str | Path | None = None) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for key, value in inputs.items():
        if key in {"sample_id", "case_id", "genome_build", "profile", "backend", "reuse_existing", "execution_mode", "tumor_type", "tumor_sample_id", "normal_sample_id", "assay_type", "rna_threads", "hla_alleles", "tool_results"}:
            continue
        values = value if isinstance(value, list) else [value]
        for index, item in enumerate(values):
            if not isinstance(item, str) or not item:
                continue
            p = Path(item)
            prefix_matches = list(p.parent.glob(p.name + "*")) if key == "rsem_reference" and p.parent.is_dir() else []
            exists = p.exists() or bool(prefix_matches)
            row: dict[str, Any] = {
                "input_key": key,
                "index": index,
                "path": str(p),
                "exists": exists,
                "is_file": p.is_file(),
                "size_bytes": p.stat().st_size if p.is_file() else 0,
                "header": "",
                "status": "OK" if exists and (not p.is_file() or p.stat().st_size > 0) else ("EMPTY" if p.is_file() else "MISSING"),
            }
            if key == "rsem_reference":
                row["reference_prefix_files"] = len(prefix_matches)
            if p.is_file() and p.suffix.lower() in {".tsv", ".csv", ".txt"}:
                row["header"] = ",".join(_table_header(p))
            if p.is_file() and (p.suffix.lower() == ".vcf" or str(p).lower().endswith(".vcf.gz")):
                row.update(_vcf_metadata(p))
            if p.is_file() and p.suffix.lower() == ".bam":
                index_candidates = [Path(str(p) + ".bai"), p.with_suffix(".bai")]
                row["bam_index"] = next((str(x) for x in index_candidates if x.is_file()), "")
                row["bam_index_status"] = "OK" if row["bam_index"] else "MISSING"
            rows.append(row)
    if inputs.get("hla_alleles"):
        invalid = [value for value in inputs["hla_alleles"] if not HLA_RE.match(str(value))]
        rows.append({"input_key": "hla_alleles", "index": 0, "path": "", "exists": True, "is_file": False, "size_bytes": 0, "header": "", "status": "INVALID" if invalid else "OK", "value": ",".join(inputs["hla_alleles"]), "invalid_values": ",".join(invalid)})
    for domain, tools in (inputs.get("tool_results") or {}).items():
        if not isinstance(tools, dict):
            continue
        for tool, value in tools.items():
            p = Path(str(value))
            rows.append({
                "input_key": f"tool_results.{domain}.{tool}", "index": 0, "path": str(p),
                "exists": p.exists(), "is_file": p.is_file(), "size_bytes": p.stat().st_size if p.is_file() else 0,
                "header": ",".join(_table_header(p)) if p.is_file() and p.suffix.lower() in {".tsv", ".csv", ".txt"} else "",
                "status": "OK" if p.is_file() and p.stat().st_size > 0 else ("EMPTY" if p.is_file() else "MISSING"),
            })
    if output_dir:
        target = Path(output_dir).resolve()
        parent = target if target.exists() else target.parent
        writable = parent.exists() and os.access(parent, os.W_OK)
        rows.append({"input_key": "output_dir", "index": 0, "path": str(target), "exists": target.exists(), "is_file": False, "size_bytes": 0, "header": "", "status": "OK" if writable else "NOT_WRITABLE"})
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
    production_inputs = [key for key in ("tumor_dna_bam", "normal_dna_bam", "tumor_rna_bam", "tumor_dna_fastq", "normal_dna_fastq", "tumor_rna_fastq") if inputs.get(key)]
    if production_inputs:
        routes.append(_route("production_inputs", "neoag-production-run", "BAM/FASTQ production inputs detected", {key: inputs[key] for key in production_inputs}))
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


def _merge_missing(base: dict[str, Any], supplemental: dict[str, Any]) -> dict[str, Any]:
    merged = dict(base)
    for key, value in supplemental.items():
        if value is None or value == "" or value == () or value == []:
            continue
        current = merged.get(key)
        if isinstance(current, dict) and isinstance(value, dict):
            nested = {subkey: (dict(subvalue) if isinstance(subvalue, dict) else subvalue) for subkey, subvalue in current.items()}
            for subkey, subvalue in value.items():
                if isinstance(nested.get(subkey), dict) and isinstance(subvalue, dict):
                    for leaf, leaf_value in subvalue.items():
                        nested[subkey].setdefault(leaf, leaf_value)
                else:
                    nested.setdefault(subkey, subvalue)
            merged[key] = nested
            continue
        if key not in merged or current is None or current == "" or current == () or current == []:
            merged[key] = value
    return merged


def inspect_manifest(
    path: str | Path,
    *,
    overrides: dict[str, Any] | None = None,
    input_dir: str | Path | None = None,
    output_dir: str | Path | None = None,
) -> RoutingResult:
    p = Path(path).resolve()
    data = load_sample_manifest(p)
    inputs = normalize_manifest(data, base_dir=p.parent)
    # Deterministic precedence: manifest, explicit CLI fill-ins, directory scan.
    if overrides:
        inputs = _merge_missing(inputs, overrides)
    scan_ambiguous: list[dict[str, str]] = []
    scan_warnings: list[str] = []
    scan_root = input_dir or inputs.get("input_dir")
    if scan_root:
        detected, scan_ambiguous, scan_warnings = scan_input_directory(scan_root)
        inputs = _merge_missing(inputs, detected)
    if inputs.get("hla_file") and not inputs.get("hla_alleles"):
        inputs["hla_alleles"] = _load_hla_alleles(inputs["hla_file"])
    inventory = build_inventory(inputs, output_dir=output_dir)
    routes, missing, ambiguous, warnings = route_inputs(inputs)
    ambiguous.extend(scan_ambiguous)
    warnings.extend(scan_warnings)
    for row in inventory:
        if row.get("status") in {"MISSING", "EMPTY", "INVALID", "NOT_WRITABLE"}:
            missing.append({"field": str(row.get("input_key")), "reason": f"input status {row.get('status')}: {row.get('path') or row.get('value', '')}"})
        if row.get("bam_index_status") == "MISSING":
            missing.append({"field": f"{row.get('input_key')}_index", "reason": f"BAM index not found for: {row.get('path')}"})
        detected_build = row.get("detected_genome_build")
        if detected_build not in {None, "", "UNASSESSED"} and str(detected_build).upper() != str(inputs.get("genome_build", "")).upper():
            ambiguous.append({"field": "genome_build", "reason": f"manifest={inputs.get('genome_build')}; detected={detected_build}; file={row.get('path')}"})
        if row.get("input_key") == "somatic_vcf" and row.get("vcf_sample_count") == 0:
            warnings.append("Somatic VCF has no genotype sample columns; tumor-normal identity cannot be verified")
        if row.get("input_key") == "somatic_vcf" and row.get("vcf_samples"):
            samples = {value for value in str(row.get("vcf_samples")).split(",") if value}
            for role, expected in (("tumor", inputs.get("tumor_sample_id")), ("normal", inputs.get("normal_sample_id"))):
                if expected and str(expected) not in samples:
                    ambiguous.append({"field": "tumor_normal", "reason": f"declared {role}_sample_id={expected} not found in VCF samples={','.join(sorted(samples))}"})
        header = [value for value in str(row.get("header") or "").split(",") if value]
        if row.get("input_key") == "fusion_tsv" and header and not _header_matches(header, FUSION_HEADER_GROUPS):
            ambiguous.append({"field": "fusion_tsv", "reason": f"unsupported fusion table header: {row.get('header')}"})
        if row.get("input_key") == "splice_junction_tsv" and header and not _header_matches(header, SPLICE_HEADER_GROUPS):
            ambiguous.append({"field": "splice_junction_tsv", "reason": f"unsupported splice table header: {row.get('header')}"})
        if row.get("input_key") == "peptide_csv" and header and not _header_matches(header, PEPTIDE_HEADERS):
            ambiguous.append({"field": "peptide_csv", "reason": f"unsupported peptide table header: {row.get('header')}"})
    if inputs.get("hla_file") and not inputs.get("hla_alleles"):
        missing.append({"field": "hla_alleles_or_hla_file", "reason": "HLA file contains no valid typed alleles"})
    if inputs.get("tumor_sample_id") and inputs.get("normal_sample_id") and inputs["tumor_sample_id"] == inputs["normal_sample_id"]:
        ambiguous.append({"field": "tumor_normal", "reason": "tumor_sample_id and normal_sample_id are identical"})
    assay = str(inputs.get("assay_type") or "").upper()
    if inputs.get("sv_vcf") and assay in {"WES", "PANEL", "CAPTURE"} and not inputs.get("capture_bed"):
        missing.append({"field": "capture_bed", "reason": "required for capture-limited SV input"})
    return RoutingResult(inputs["case_id"], inputs["sample_id"], inputs["genome_build"], inputs, inventory, routes, missing, ambiguous, warnings, str(p))


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
