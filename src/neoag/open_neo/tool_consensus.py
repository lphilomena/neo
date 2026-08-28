from __future__ import annotations

import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from neoag.agent_skills.fusion_rna_run import collect_results as collect_fusions
from neoag.agent_skills.fusion_rna_run import consensus as fusion_consensus
from neoag.agent_skills.hla_typing_compare import collect_typing, consensus as hla_consensus
from neoag.agent_skills.hla_typing_compare import locus_of, lowres, norm_allele
from neoag.agent_skills.purity_cnv_review import collect_tool_results as collect_purity
from neoag.agent_skills.purity_cnv_review import consensus as purity_consensus
from neoag.hla_loh_crosscheck import crosscheck_hla_loh
from neoag.splice.coordinates import iter_junction_records, peptide_metadata
from neoag.splice.registry import JunctionRegistry, unresolved_event_id
from neoag.utils import read_tsv, write_tsv


EXPECTED_TOOLS = {
    "sample_identity": ("bam-matcher",),
    "hla_typing": ("optitype", "spechla", "hla-la", "hla-hd"),
    "hla_loh": ("lohhla", "spechla", "facets_hla_cnv", "purple_hla_cnv", "hla_expression"),
    "fusion": ("easyfuse", "arriba", "star-fusion", "fusioncatcher"),
    "splice_dna": ("spliceai", "snap"),
    "splice_rna": ("regtools", "rmats", "majiq"),
    "splice_neoantigen": ("pvacsplice", "snaf", "splicemutr"),
    "presentation": ("netmhcpan", "mhcflurry", "netmhcstabpan", "prime", "bigmhc", "deepimmuno"),
    "purity_cnv": ("facets", "ascat", "purple", "sequenza"),
    "ccf": ("pyclone-vi",),
}


def _declared(inputs: dict[str, Any], domain: str) -> dict[str, str]:
    tool_results = inputs.get("tool_results") or {}
    values = tool_results.get(domain) if isinstance(tool_results, dict) else {}
    if not isinstance(values, dict):
        return {}
    return {str(k).lower(): str(v) for k, v in values.items() if isinstance(v, str) and v}


def _all_declared(inputs: dict[str, Any]) -> dict[str, dict[str, str]]:
    result = {domain: _declared(inputs, domain) for domain in EXPECTED_TOOLS}
    direct = {
        "sample_identity": inputs.get("sample_identity_tsv"),
        "hla_typing": inputs.get("hla_file"),
        "hla_loh": inputs.get("hla_loh_tsv"),
        "fusion": inputs.get("fusion_tsv"),
        "splice_rna": inputs.get("splice_junction_tsv"),
        "purity_cnv": inputs.get("purity_tsv"),
        "ccf": inputs.get("ccf_tsv"),
    }
    for domain, path in direct.items():
        if path:
            result[domain].setdefault("declared_input", str(path))
    return result


def _tool_status(declared: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for domain, expected in EXPECTED_TOOLS.items():
        seen = set()
        for tool in expected:
            path = declared.get(domain, {}).get(tool, "")
            seen.add(tool)
            rows.append({
                "evidence_domain": domain,
                "tool": tool,
                "status": "AVAILABLE" if path and Path(path).exists() else ("MISSING" if path else "NOT_DECLARED"),
                "version": "UNASSESSED",
                "input_comparable": "yes" if path and Path(path).exists() else "unassessed",
                "source_file": path,
                "reason": "declared result exists" if path and Path(path).is_file() else ("declared path missing" if path else "no result declared"),
            })
        for tool, path in sorted(declared.get(domain, {}).items()):
            if tool in seen:
                continue
            rows.append({
                "evidence_domain": domain, "tool": tool,
                "status": "AVAILABLE" if Path(path).exists() else "MISSING",
                "version": "UNASSESSED", "input_comparable": "yes" if Path(path).exists() else "no",
                "source_file": path, "reason": "additional declared result",
            })
    return rows


def _write_hla(inputs: dict[str, Any], declared: dict[str, dict[str, str]], outdir: Path) -> tuple[str, list[dict[str, str]]]:
    paths = [Path(path) for path in declared["hla_typing"].values() if Path(path).exists()]
    rows = collect_typing(paths, sample_id=None)
    if not rows and inputs.get("hla_alleles"):
        grouped: dict[str, list[str]] = defaultdict(list)
        for raw in inputs["hla_alleles"]:
            allele = norm_allele(str(raw))
            locus = locus_of(allele)
            if allele and allele not in grouped[locus]:
                grouped[locus].append(allele)
        for locus, alleles in grouped.items():
            rows.append({
                "tool": "ProvidedHLA", "locus": locus,
                "allele1": alleles[0] if alleles else "", "allele2": alleles[1] if len(alleles) > 1 else "",
                "lowres1": lowres(alleles[0]) if alleles else "", "lowres2": lowres(alleles[1]) if len(alleles) > 1 else "",
                "source_file": "manifest_or_cli",
            })
    consensus_rows = hla_consensus(rows)
    mapping = {"CONSENSUS": "CONSISTENT", "SINGLE_TOOL": "SINGLE_TOOL_ONLY", "MISSING": "UNASSESSED", "DISCORDANT": "DISCORDANT"}
    for row in consensus_rows:
        row["status"] = mapping.get(str(row.get("status")), str(row.get("status")))
    path = outdir / "hla_typing_consensus.tsv"
    write_tsv(path, consensus_rows)
    states = [str(row.get("status")) for row in consensus_rows]
    overall = "DISCORDANT" if "DISCORDANT" in states else ("PARTIAL_CONSISTENCY" if "SINGLE_TOOL_ONLY" in states else ("CONSISTENT" if "CONSISTENT" in states else "UNASSESSED"))
    return overall, consensus_rows


def _write_hla_loh(declared: dict[str, dict[str, str]], outdir: Path) -> tuple[str, list[dict[str, str]], list[dict[str, str]]]:
    paths = declared["hla_loh"]
    rows = crosscheck_hla_loh(lohhla_hla_loh=paths.get("lohhla"), spechla_hla_loh=paths.get("spechla"))
    provided = paths.get("provided_consensus") or paths.get("declared_input")
    if provided and Path(provided).is_file() and not rows:
        for record in read_tsv(provided):
            allele = str(record.get("hla_allele") or record.get("allele") or record.get("HLA") or "")
            status = str(record.get("loh_status") or record.get("status") or record.get("LOH") or "unassessed").lower()
            if allele:
                rows.append({"hla_allele": allele, "consensus_loh_status": status, "crosscheck_status": "SINGLE_TOOL_ONLY", "source_tools": "provided_consensus", "reason": "precomputed HLA LOH input"})
    conflicts: list[dict[str, str]] = []
    for row in rows:
        allele = str(row.get("hla_allele") or "")
        locus = allele.replace("HLA-", "").split("*", 1)[0]
        row["hla_class"] = "I" if locus in {"A", "B", "C"} else "II"
        if row.get("crosscheck_status") == "DISCORDANT":
            conflicts.append({"evidence_domain": "hla_loh", "record_id": allele, "conflict_type": "TOOL_DISCORDANCE", "details": str(row.get("reason") or "")})
    if not rows:
        rows = [{"hla_allele": "", "hla_class": "", "consensus_loh_status": "unassessed", "crosscheck_status": "UNASSESSED", "source_tools": "", "reason": "No comparable allele-level LOHHLA/SpecHLA result"}]
    path = outdir / "hla_loh_consensus.tsv"
    write_tsv(path, rows)
    states = {str(row.get("crosscheck_status")) for row in rows}
    overall = "DISCORDANT" if "DISCORDANT" in states else ("CONSISTENT" if states & {"CONSENSUS_LOH", "CONSENSUS_NO_LOH"} else ("SINGLE_TOOL_ONLY" if any(value.startswith("SINGLE_TOOL") for value in states) else "UNASSESSED"))
    return overall, rows, conflicts


def _write_fusion(inputs: dict[str, Any], declared: dict[str, dict[str, str]], outdir: Path) -> tuple[str, list[dict[str, Any]]]:
    paths = [Path(path) for path in declared["fusion"].values() if Path(path).exists()]
    _, calls = collect_fusions(paths, sample_id=None)
    if not calls:
        for tool, path in declared["fusion"].items():
            if not Path(path).is_file():
                continue
            for row in read_tsv(path):
                fusion = str(row.get("fusion") or row.get("FusionName") or row.get("fusion_name") or "")
                if not fusion:
                    left = str(row.get("gene1") or row.get("gene5") or row.get("LeftGene") or "")
                    right = str(row.get("gene2") or row.get("gene3") or row.get("RightGene") or "")
                    fusion = f"{left}--{right}" if left and right else ""
                if fusion:
                    calls.append({"tool": tool, "fusion": fusion, "source_file": path})
    rows = fusion_consensus(calls)
    metrics: dict[str, dict[str, Any]] = defaultdict(lambda: {"junction_reads": 0, "spanning_reads": 0, "frame": "UNASSESSED", "readthrough": False})
    for _, path in declared["fusion"].items():
        if not Path(path).is_file():
            continue
        for record in read_tsv(path):
            fusion = str(record.get("fusion") or record.get("FusionName") or record.get("fusion_name") or "")
            if not fusion:
                left = str(record.get("gene1") or record.get("gene5") or record.get("LeftGene") or "")
                right = str(record.get("gene2") or record.get("gene3") or record.get("RightGene") or "")
                fusion = f"{left}--{right}" if left and right else ""
            if not fusion:
                continue
            folded = {re.sub(r"[^a-z0-9]+", "", str(key).lower()): value for key, value in record.items()}
            for target, aliases in (("junction_reads", ("junctionreads", "junctionreadcount", "splitreads")), ("spanning_reads", ("spanningreads", "spanningpairs"))):
                for alias in aliases:
                    try:
                        metrics[fusion][target] = max(int(float(folded.get(alias, 0) or 0)), int(metrics[fusion][target]))
                    except (TypeError, ValueError):
                        pass
            frame = str(record.get("frame") or record.get("reading_frame") or record.get("Frame") or "")
            if frame:
                metrics[fusion]["frame"] = frame
            flag_text = " ".join(str(value) for value in record.values()).lower()
            metrics[fusion]["readthrough"] = metrics[fusion]["readthrough"] or "read-through" in flag_text or "readthrough" in flag_text
    for row in rows:
        evidence = metrics.get(str(row.get("fusion")), {})
        row.update(evidence)
        if evidence.get("readthrough"):
            row["status"] = "READTHROUGH_CAUTION"
        elif int(row.get("n_tools") or 0) >= 2:
            row["status"] = "MULTI_CALLER_STRONG"
        elif int(evidence.get("junction_reads") or 0) >= 3 or int(evidence.get("spanning_reads") or 0) >= 5:
            row["status"] = "SINGLE_CALLER_SUPPORTED"
        elif evidence:
            row["status"] = "LOW_SUPPORT"
        else:
            row["status"] = "SINGLE_CALLER_SUPPORTED"
    if not rows:
        rows = [{"fusion": "", "support_tools": "", "n_tools": 0, "status": "UNASSESSED"}]
    write_tsv(outdir / "fusion_consensus.tsv", rows)
    states = {str(row.get("status")) for row in rows}
    overall = "MULTI_CALLER_STRONG" if "MULTI_CALLER_STRONG" in states else ("SINGLE_CALLER_SUPPORTED" if "SINGLE_CALLER_SUPPORTED" in states else ("READTHROUGH_CAUTION" if "READTHROUGH_CAUTION" in states else ("LOW_SUPPORT" if "LOW_SUPPORT" in states else "UNASSESSED")))
    return overall, rows


def _write_splice(
    declared: dict[str, dict[str, str]],
    outdir: Path,
) -> tuple[str, list[dict[str, str]]]:
    """Create exact splice consensus without gene/locus string matching.

    RNA evidence is registered first. DNA-prediction and neoantigen rows may
    join it only through a canonical junction, exact unique source alias, or an
    explicit unique variant-to-junction relation implemented by
    :class:`JunctionRegistry`. Unresolved rows receive source-scoped IDs and can
    never produce ``CROSS_DOMAIN_CONFIRMED``.
    """

    registry = JunctionRegistry()
    evidence: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    provenance: list[dict[str, str]] = []

    def add_record(domain: str, tool: str, path: str, row_number: int, record, resolution) -> None:
        canonical_id = resolution.junction_id
        event_id = canonical_id or unresolved_event_id(record)
        evidence[event_id][domain].add(tool)
        evidence[event_id]["source_records"].add(record.source_record_id)
        provenance.append(
            {
                "event_id": event_id,
                "canonical_junction_id": canonical_id,
                "evidence_domain": domain,
                "tool": tool,
                "source_file": path,
                "source_row_number": str(row_number),
                "source_record_id": record.source_record_id,
                "source_junction_id": record.source_junction_id,
                "resolution_status": resolution.status,
                "resolution_method": resolution.method,
                "coordinate_warning": resolution.warning or record.coordinate_warning,
                "peptide_present": "yes" if peptide_metadata(record)["peptide"] else "no",
            }
        )

    # RNA evidence establishes the primary exact registry.
    for tool, path in sorted(declared["splice_rna"].items()):
        source = Path(path)
        if not source.is_file():
            continue
        for row_number, record in enumerate(
            iter_junction_records(
                source,
                sample_id="",
                source_tool=tool,
                genome_build="GRCh38",
                coordinate_system="auto",
                strict=False,
            ),
            1,
        ):
            resolution = registry.add(record)
            add_record("splice_rna", tool, str(source), row_number, record, resolution)

    # Secondary domains are resolved against the RNA registry. Exact
    # coordinate-only secondary entities are retained, but do not become RNA
    # support merely because they were normalized.
    for domain in ("splice_dna", "splice_neoantigen"):
        for tool, path in sorted(declared[domain].items()):
            source = Path(path)
            if not source.is_file():
                continue
            for row_number, record in enumerate(
                iter_junction_records(
                    source,
                    sample_id="",
                    source_tool=tool,
                    genome_build="GRCh38",
                    coordinate_system="auto",
                    strict=False,
                ),
                1,
            ):
                resolution = registry.resolve(record)
                if resolution.junction is not None:
                    registry.add(record, junction=resolution.junction)
                add_record(domain, tool, str(source), row_number, record, resolution)

    rows: list[dict[str, str]] = []
    for event_id, domains in sorted(evidence.items()):
        canonical_id = event_id if event_id.startswith("SJ|") else ""
        dna = ";".join(sorted(domains.get("splice_dna", set())))
        rna = ";".join(sorted(domains.get("splice_rna", set())))
        neo = ";".join(sorted(domains.get("splice_neoantigen", set())))
        if canonical_id and rna and neo:
            status = "CROSS_DOMAIN_CONFIRMED_EXACT_JUNCTION"
        elif canonical_id and rna and dna:
            status = "DNA_RNA_EXACT_JUNCTION_SUPPORTED"
        elif canonical_id and rna:
            status = "RNA_JUNCTION_SUPPORTED"
        elif canonical_id and neo and dna:
            status = "DNA_NEO_EXACT_COORDINATE_NO_RNA_SUPPORT"
        elif canonical_id and neo:
            status = "NEOANTIGEN_TOOL_ONLY_EXACT_COORDINATE"
        elif dna and not rna and not neo:
            status = "PREDICTION_ONLY_UNRESOLVED"
        elif neo and not rna:
            status = "NEOANTIGEN_TOOL_ONLY_UNRESOLVED"
        else:
            status = "UNRESOLVED_SOURCE_RECORD"
        rows.append(
            {
                "event_id": event_id,
                "canonical_junction_id": canonical_id,
                "junction_resolution_status": "RESOLVED" if canonical_id else "UNRESOLVED",
                "dna_prediction_tools": dna,
                "rna_junction_tools": rna,
                "neoantigen_tools": neo,
                "source_records": ";".join(sorted(domains.get("source_records", set()))),
                "status": status,
            }
        )

    if not rows:
        rows = [
            {
                "event_id": "",
                "canonical_junction_id": "",
                "junction_resolution_status": "UNASSESSED",
                "dna_prediction_tools": "",
                "rna_junction_tools": "",
                "neoantigen_tools": "",
                "source_records": "",
                "status": "UNASSESSED",
            }
        ]

    write_tsv(
        outdir / "splice_consensus.tsv",
        rows,
        [
            "event_id",
            "canonical_junction_id",
            "junction_resolution_status",
            "dna_prediction_tools",
            "rna_junction_tools",
            "neoantigen_tools",
            "source_records",
            "status",
        ],
    )
    write_tsv(
        outdir / "splice_consensus_provenance.tsv",
        provenance,
        [
            "event_id",
            "canonical_junction_id",
            "evidence_domain",
            "tool",
            "source_file",
            "source_row_number",
            "source_record_id",
            "source_junction_id",
            "resolution_status",
            "resolution_method",
            "coordinate_warning",
            "peptide_present",
        ],
    )
    write_tsv(
        outdir / "splice_consensus_conflicts.tsv",
        registry.conflicts
        or [
            {
                "evidence_domain": "splice_junction",
                "record_id": "",
                "conflict_type": "NONE",
                "details": "",
                "source_tool": "",
                "source_file": "",
                "source_row_number": "",
            }
        ],
    )
    states = {row["status"] for row in rows}
    overall = (
        "CROSS_DOMAIN_CONFIRMED"
        if "CROSS_DOMAIN_CONFIRMED_EXACT_JUNCTION" in states
        else ("PARTIAL" if states != {"UNASSESSED"} else "UNASSESSED")
    )
    return overall, rows


def _write_purity(inputs: dict[str, Any], declared: dict[str, dict[str, str]], outdir: Path) -> tuple[str, dict[str, Any]]:
    paths = [Path(path) for path in declared["purity_cnv"].values() if Path(path).exists()]
    rows = collect_purity(paths, sample_id=None)
    by_tool = {str(row.get("tool", "")).lower(): row for row in rows}
    purity_keys = ("purity", "cellularity", "tumor_purity", "tumour_purity", "aberrantcellfraction", "rho")
    ploidy_keys = ("ploidy", "tumor_ploidy", "tumour_ploidy", "psi")
    for declared_tool, path in declared["purity_cnv"].items():
        canonical = declared_tool.replace("_", "-").lower()
        existing = next((key for key in by_tool if canonical in key or key in canonical), "")
        if existing and by_tool[existing].get("status") == "FOUND":
            continue
        parsed = read_tsv(path) if Path(path).is_file() else []
        purity = ""
        ploidy = ""
        for record in parsed[:100]:
            folded = {re.sub(r"[^a-z0-9]+", "", str(key).lower()): value for key, value in record.items()}
            for key in purity_keys:
                value = folded.get(re.sub(r"[^a-z0-9]+", "", key))
                if value not in {None, "", "NA", "NaN"}:
                    try:
                        purity = float(value)
                        break
                    except (TypeError, ValueError):
                        pass
            for key in ploidy_keys:
                value = folded.get(re.sub(r"[^a-z0-9]+", "", key))
                if value not in {None, "", "NA", "NaN"}:
                    try:
                        ploidy = float(value)
                        break
                    except (TypeError, ValueError):
                        pass
            if purity != "" or ploidy != "":
                break
        if purity != "" or ploidy != "":
            label = declared_tool.upper() if declared_tool.lower() != "sequenza" else "Sequenza"
            replacement = {"tool": label, "status": "FOUND", "purity": purity, "ploidy": ploidy, "source_file": path, "parse_method": "declared generic table", "notes": ""}
            rows = [row for row in rows if str(row.get("tool", "")).lower() != declared_tool.lower()]
            rows.append(replacement)
    result = purity_consensus(rows)
    write_tsv(outdir / "purity_cnv_tool_results.tsv", rows)
    write_tsv(outdir / "purity_cnv_consensus.tsv", [{**result, "tool_values": str(result.get("tool_values") or {})}])
    ccf_confidence = "low" if result.get("status") in {"STRONG_DISCORDANCE", "SINGLE_TOOL", "NO_PURITY"} else "moderate" if result.get("status") == "MODERATE_DISCORDANCE" else "high"
    write_tsv(outdir / "ccf_consensus.tsv", [{"status": "UNASSESSED" if result.get("status") == "NO_PURITY" else "PURITY_INFORMED", "ccf_confidence": ccf_confidence, "purity_consensus_status": result.get("status", "UNASSESSED"), "reason": result.get("interpretation", "")}])
    return str(result.get("status") or "UNASSESSED"), result


def _presentation_and_ccf_from_evidence(evidence_path: str | Path | None, outdir: Path) -> tuple[str, str]:
    rows = read_tsv(evidence_path) if evidence_path and Path(evidence_path).is_file() else []
    presentation = Counter(str(row.get("presentation_consensus_state") or "PRESENTATION_UNASSESSED") for row in rows)
    presentation_rows = [{"status": key, "candidate_count": value} for key, value in sorted(presentation.items())] or [{"status": "PRESENTATION_UNASSESSED", "candidate_count": 0}]
    write_tsv(outdir / "presentation_consensus.tsv", presentation_rows)
    ccf = Counter(str(row.get("clonality_state") or row.get("ccf_confidence") or "UNASSESSED") for row in rows)
    write_tsv(outdir / "ccf_candidate_summary.tsv", [{"status": key, "candidate_count": value} for key, value in sorted(ccf.items())] or [{"status": "UNASSESSED", "candidate_count": 0}])
    restricting_rows = [{
        "peptide_id": str(row.get("peptide_id") or ""),
        "event_id": str(row.get("event_id") or ""),
        "hla_allele": str(row.get("hla_allele") or row.get("restricting_hla") or ""),
        "restricting_hla_lost": str(row.get("restricting_hla_lost") or "UNASSESSED"),
        "escape_status": str(row.get("escape_status") or "UNASSESSED"),
    } for row in rows]
    write_tsv(outdir / "restricting_hla_peptide_flags.tsv", restricting_rows or [{"peptide_id": "", "event_id": "", "hla_allele": "", "restricting_hla_lost": "UNASSESSED", "escape_status": "UNASSESSED"}])
    pstatus = "UNASSESSED" if not rows else ("DISCORDANT" if any("DISCORDANT" in key for key in presentation) else "ASSESSED")
    cstatus = "UNASSESSED" if not rows else ("LOW_CONFIDENCE" if any("LOW" in key or "UNASSESSED" in key for key in ccf) else "ASSESSED")
    return pstatus, cstatus


def _write_sample_identity(declared: dict[str, dict[str, str]], outdir: Path) -> tuple[str, list[dict[str, str]]]:
    rows: list[dict[str, str]] = []
    for tool, path in declared["sample_identity"].items():
        if not Path(path).is_file():
            continue
        for record in read_tsv(path):
            rows.append({
                "tool": tool,
                "sample_identity_status": str(record.get("sample_identity_status") or record.get("status") or "UNASSESSED"),
                "official_conclusion": str(record.get("official_conclusion") or ""),
                "confidence": str(record.get("confidence") or ""),
                "fraction_common": str(record.get("fraction_common") or ""),
                "sites_compared": str(record.get("sites_compared") or ""),
                "source_file": path,
            })
    if not rows:
        rows = [{"tool": "", "sample_identity_status": "UNASSESSED", "official_conclusion": "", "confidence": "", "fraction_common": "", "sites_compared": "", "source_file": ""}]
    write_tsv(outdir / "sample_identity_consensus.tsv", rows)
    states = {row["sample_identity_status"] for row in rows}
    status = "MISMATCH" if "MISMATCH" in states else ("MATCH" if "MATCH" in states else "INSUFFICIENT_DATA" if "INSUFFICIENT_DATA" in states else "UNASSESSED")
    return status, rows


def build_tool_consensus(inputs: dict[str, Any], outdir: str | Path, *, evidence_path: str | Path | None = None) -> dict[str, str]:
    root = Path(outdir)
    root.mkdir(parents=True, exist_ok=True)
    declared = _all_declared(inputs)
    status_rows = _tool_status(declared)
    write_tsv(root / "tool_run_status.tsv", status_rows)
    hla_status, _ = _write_hla(inputs, declared, root)
    loh_status, _, conflicts = _write_hla_loh(declared, root)
    fusion_status, _ = _write_fusion(inputs, declared, root)
    splice_status, _ = _write_splice(declared, root)
    purity_status, _ = _write_purity(inputs, declared, root)
    presentation_status, ccf_status = _presentation_and_ccf_from_evidence(evidence_path, root)
    sample_identity_status, _ = _write_sample_identity(declared, root)
    if evidence_path:
        ranked_conflicts = Path(evidence_path).parent / "evidence_conflicts.tsv"
        if ranked_conflicts.is_file():
            for row in read_tsv(ranked_conflicts):
                conflicts.append({
                    "evidence_domain": str(row.get("domain") or row.get("evidence_domain") or "evidence_normalization"),
                    "record_id": str(row.get("peptide_id") or row.get("record_id") or ""),
                    "conflict_type": str(row.get("conflict_type") or row.get("reason") or "FIELD_CONFLICT"),
                    "details": str(row.get("details") or row.get("values") or ""),
                })
    summary = [
        {"evidence_domain": "sample_identity", "consensus_status": sample_identity_status, "adopted_evidence": "genotype concordance across build-matched common SNPs", "reason": "MISMATCH blocks paired tumor-normal analyses; insufficient coverage requires review"},
        {"evidence_domain": "hla_typing", "consensus_status": hla_status, "adopted_evidence": "locus-level consensus", "reason": "allele calls compared by locus"},
        {"evidence_domain": "hla_loh", "consensus_status": loh_status, "adopted_evidence": "allele-level class-aware LOH", "reason": "HLA-I and HLA-II remain separate"},
        {"evidence_domain": "fusion", "consensus_status": fusion_status, "adopted_evidence": "caller and junction evidence", "reason": "caller count is supportive, not sufficient alone"},
        {"evidence_domain": "splice", "consensus_status": splice_status, "adopted_evidence": "DNA/RNA/neoantigen evidence chain", "reason": "three splice tool classes are not pooled"},
        {"evidence_domain": "presentation", "consensus_status": presentation_status, "adopted_evidence": "core/stability/immunogenicity evidence groups", "reason": "correlated immunogenicity models are grouped"},
        {"evidence_domain": "purity_cnv", "consensus_status": purity_status, "adopted_evidence": "median/range with discordance", "reason": "strong disagreement lowers CCF confidence"},
        {"evidence_domain": "ccf", "consensus_status": ccf_status, "adopted_evidence": "purity/CNV-informed candidate state", "reason": "missing RNA-only CCF remains unresolved"},
    ]
    write_tsv(root / "tool_consensus_summary.tsv", summary)
    write_tsv(root / "evidence_conflicts.tsv", conflicts or [{"evidence_domain": "", "record_id": "", "conflict_type": "NONE", "details": ""}])
    write_tsv(root / "evidence_source_conflicts.tsv", conflicts or [{"evidence_domain": "", "record_id": "", "conflict_type": "NONE", "details": ""}])
    long_rows = []
    for row in status_rows:
        long_rows.append({"record_type": "tool_run", "evidence_domain": row["evidence_domain"], "tool": row["tool"], "status": row["status"], "value": row["source_file"], "reason": row["reason"]})
    for row in summary:
        long_rows.append({"record_type": "domain_consensus", "evidence_domain": row["evidence_domain"], "tool": "consensus", "status": row["consensus_status"], "value": row["adopted_evidence"], "reason": row["reason"]})
    write_tsv(root / "tool_evidence.long.tsv", long_rows)
    return {name: str(root / name) for name in (
        "tool_run_status.tsv", "tool_consensus_summary.tsv", "evidence_conflicts.tsv", "evidence_source_conflicts.tsv",
        "tool_evidence.long.tsv", "hla_typing_consensus.tsv", "hla_loh_consensus.tsv", "fusion_consensus.tsv",
        "splice_consensus.tsv", "presentation_consensus.tsv", "purity_cnv_consensus.tsv", "ccf_consensus.tsv",
        "restricting_hla_peptide_flags.tsv",
        "sample_identity_consensus.tsv",
    )}


def enrich_all_tool_results(path: str | Path, consensus_summary: str | Path) -> None:
    target = Path(path)
    if not target.is_file() or not Path(consensus_summary).is_file():
        return
    rows = read_tsv(target)
    summary_rows = read_tsv(consensus_summary)
    statuses = {str(row.get("evidence_domain")): str(row.get("consensus_status")) for row in summary_rows}
    for row in rows:
        conflicts = sorted(key for key, value in statuses.items() if value in {"DISCORDANT", "STRONG_DISCORDANCE", "MISMATCH"})
        row["tool_consensus_overall"] = "BLOCKED" if statuses.get("sample_identity") == "MISMATCH" else ("REVIEW_REQUIRED" if conflicts else "ASSESSED")
        row["tool_consensus_conflicts"] = ",".join(conflicts)
        for domain, status in statuses.items():
            row[f"{domain}_consensus_status"] = status
    write_tsv(target, rows)
