from __future__ import annotations

from pathlib import Path
from typing import Any

from .io import ensure_dir, normalize_hla, read_table, read_vcf_records, row_get, write_json, write_tsv


def _event_type(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    return "InDel"


def run_vcf(args: dict[str, Any]) -> dict[str, Any]:
    outdir = ensure_dir(args["outdir"])
    vcf = Path(args.get("vcf") or args.get("input") or "")
    sample_id = args.get("sample_id") or vcf.stem.replace(".vcf", "")
    _, records = read_vcf_records(vcf)
    events = []
    for rec in records:
        eid = f"{rec.get('gene') or 'NA'}|{rec['chrom']}:{rec['pos']}{rec['ref']}>{rec['alt']}"
        events.append({
            "sample_id": sample_id,
            "event_id": eid,
            "event_type": _event_type(rec["ref"], rec["alt"]),
            "event_subtype": rec.get("consequence", ""),
            "gene": rec.get("gene", ""),
            "transcript": rec.get("transcript", ""),
            "chrom": rec["chrom"],
            "pos": rec["pos"],
            "ref": rec["ref"],
            "alt": rec["alt"],
            "protein_change": rec.get("protein_change", ""),
            "source_type": "somatic_vcf",
            "source_file": str(vcf),
            "confidence": "vcf_pass" if rec.get("filter") in {"PASS", ".", ""} else "vcf_nonpass_review",
        })
    write_tsv(outdir / "raw_events.tsv", events)
    write_tsv(outdir / "raw_peptides.tsv", [], ["sample_id", "event_id", "peptide_id", "peptide", "hla_allele", "source_type", "generation_status"])
    write_tsv(outdir / "vcf_parse_qc.tsv", [{"metric": "records", "value": len(records)}, {"metric": "events", "value": len(events)}])
    (outdir / "candidate_generation_plan.md").write_text(
        "# VCF entry skill\n\nParsed somatic VCF into `raw_events.tsv`. Peptide generation should be performed by the SNV/InDel peptide generation stage before presentation prediction.\n",
        encoding="utf-8",
    )
    res = {"status": "PASS", "skill": "neoag-vcf", "outputs": {"raw_events": str(outdir / "raw_events.tsv"), "raw_peptides": str(outdir / "raw_peptides.tsv")}, "summary": f"Parsed {len(events)} VCF events"}
    write_json(outdir / "skill_result.json", res)
    return res


def run_fusion(args: dict[str, Any]) -> dict[str, Any]:
    outdir = ensure_dir(args["outdir"])
    fusion = Path(args.get("fusion") or args.get("input") or "")
    sample_id = args.get("sample_id") or fusion.stem
    header, rows = read_table(fusion)
    normal_pairs: set[str] = set()
    normal_db = args.get("normal_readthrough_db")
    if normal_db:
        try:
            _, normal_rows = read_table(normal_db)
            for nr in normal_rows:
                ng5 = row_get(nr, ["gene5", "gene1", "left_gene", "Gene1", "FusionName"], "")
                ng3 = row_get(nr, ["gene3", "gene2", "right_gene", "Gene2"], "")
                if "::" in ng5 and not ng3:
                    ng5, ng3 = ng5.split("::", 1)
                normal_pairs.add(f"{ng5}::{ng3}".upper())
        except Exception:
            normal_pairs = set()
    events = []
    peptides = []
    evidence = []
    for idx, row in enumerate(rows, 1):
        g5 = row_get(row, ["gene5", "gene1", "left_gene", "Gene1", "#gene1", "FusionName"], "")
        g3 = row_get(row, ["gene3", "gene2", "right_gene", "Gene2"], "")
        if "::" in g5 and not g3:
            g5, g3 = g5.split("::", 1)
        if "--" in g5 and not g3:
            g5, g3 = g5.split("--", 1)
        event_id = row_get(row, ["event_id", "fusion_id", "FusionName"], "") or f"FUSION_{idx}_{g5}_{g3}"
        jr = row_get(row, ["junction_reads", "JunctionReadCount", "split_reads", "supporting_reads", "reads"], "")
        frame = row_get(row, ["frame", "Frame", "reading_frame", "in_frame"], "")
        pair = f"{g5}::{g3}"
        in_normal = pair.upper() in normal_pairs
        events.append({"sample_id": sample_id, "event_id": event_id, "event_type": "Fusion", "gene": pair, "gene5": g5, "gene3": g3, "rna_junction_reads": jr, "frame_status": frame, "source_file": str(fusion), "confidence": "normal_background_review" if in_normal else "requires_junction_frame_review", "normal_background_hit": str(in_normal).lower()})
        pep = row_get(row, ["peptide", "junction_peptide", "neoepitope", "mutant_peptide"], "")
        hla = normalize_hla(row_get(row, ["hla", "hla_allele", "HLA", "allele"], ""))
        if pep:
            peptides.append({"sample_id": sample_id, "event_id": event_id, "peptide_id": f"{event_id}_{pep}_{hla or 'HLA_NA'}", "peptide": pep, "hla_allele": hla, "source_type": "fusion", "crosses_junction": "true", "generation_status": "provided"})
        evidence.append({"event_id": event_id, "gene5": g5, "gene3": g3, "junction_reads": jr, "frame_status": frame, "review_flag": "normal_background_hit" if in_normal else "readthrough_normal_background_review"})
    write_tsv(outdir / "fusion_events.tsv", events)
    write_tsv(outdir / "raw_events.tsv", events)
    write_tsv(outdir / "raw_peptides.tsv", peptides, ["sample_id", "event_id", "peptide_id", "peptide", "hla_allele", "source_type", "crosses_junction", "generation_status"])
    write_tsv(outdir / "fusion_evidence.tsv", evidence)
    write_tsv(outdir / "fusion_qc.tsv", [{"metric": "input_rows", "value": len(rows)}, {"metric": "events", "value": len(events)}, {"metric": "provided_peptides", "value": len(peptides)}])
    res = {"status": "PASS", "skill": "neoag-fusion", "outputs": {"raw_events": str(outdir / "raw_events.tsv"), "raw_peptides": str(outdir / "raw_peptides.tsv")}, "summary": f"Normalized {len(events)} fusion events"}
    write_json(outdir / "skill_result.json", res)
    return res


def run_splice(args: dict[str, Any]) -> dict[str, Any]:
    outdir = ensure_dir(args["outdir"])
    junctions = Path(args.get("junctions") or args.get("input") or "")
    sample_id = args.get("sample_id") or junctions.stem
    header, rows = read_table(junctions)
    events = []
    peptides = []
    for idx, row in enumerate(rows, 1):
        chrom = row_get(row, ["chrom", "chr", "seqnames"], "")
        start = row_get(row, ["start", "junction_start", "donor"], "")
        end = row_get(row, ["end", "junction_end", "acceptor"], "")
        gene = row_get(row, ["gene", "gene_name", "symbol"], "")
        reads = row_get(row, ["junction_reads", "read_count", "reads", "score"], "")
        event_id = row_get(row, ["event_id", "junction_id"], "") or f"SPLICE_{idx}_{chrom}:{start}-{end}"
        events.append({"sample_id": sample_id, "event_id": event_id, "event_type": "Splice", "gene": gene, "chrom": chrom, "start": start, "end": end, "rna_junction_reads": reads, "source_file": str(junctions), "confidence": "requires_rna_junction_review"})
        pep = row_get(row, ["peptide", "junction_peptide", "mutant_peptide"], "")
        if pep:
            peptides.append({"sample_id": sample_id, "event_id": event_id, "peptide_id": f"{event_id}_{pep}", "peptide": pep, "hla_allele": normalize_hla(row_get(row, ["hla", "hla_allele"], "")), "source_type": "splice_junction", "crosses_junction": "true", "generation_status": "provided"})
    write_tsv(outdir / "splice_events.tsv", events)
    write_tsv(outdir / "raw_events.tsv", events)
    write_tsv(outdir / "raw_peptides.tsv", peptides, ["sample_id", "event_id", "peptide_id", "peptide", "hla_allele", "source_type", "crosses_junction", "generation_status"])
    write_tsv(outdir / "splice_qc.tsv", [{"metric": "junction_rows", "value": len(rows)}, {"metric": "provided_peptides", "value": len(peptides)}])
    res = {"status": "PASS", "skill": "neoag-splice", "outputs": {"raw_events": str(outdir / "raw_events.tsv"), "raw_peptides": str(outdir / "raw_peptides.tsv")}, "summary": f"Normalized {len(events)} splice/junction events"}
    write_json(outdir / "skill_result.json", res)
    return res


def _read_bed_intervals(path: str | Path | None) -> list[tuple[str, int, int]]:
    intervals: list[tuple[str, int, int]] = []
    if not path:
        return intervals
    with Path(path).open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                try:
                    intervals.append((parts[0].replace("chr", ""), int(parts[1]), int(parts[2])))
                except Exception:
                    continue
    return intervals


def _bed_overlap(chrom: str, pos: str, intervals: list[tuple[str, int, int]]) -> str:
    if not intervals:
        return "unassessed"
    try:
        c = chrom.replace("chr", "")
        p = int(pos)
    except Exception:
        return "unassessed"
    return "true" if any(c == bc and start <= p <= end for bc, start, end in intervals) else "false"


def _run_sv(args: dict[str, Any], capture_limited: bool) -> dict[str, Any]:
    outdir = ensure_dir(args["outdir"])
    sv_vcf = Path(args.get("sv_vcf") or args.get("input") or "")
    sample_id = args.get("sample_id") or sv_vcf.stem.replace(".vcf", "")
    capture_intervals = _read_bed_intervals(args.get("capture_bed")) if capture_limited else []
    _, records = read_vcf_records(sv_vcf)
    events = []
    tasks = []
    for idx, rec in enumerate(records, 1):
        svtype = rec.get("svtype") or ("BND" if "[" in rec["alt"] or "]" in rec["alt"] else "SV")
        event_id = f"SV_{idx}_{rec['chrom']}:{rec['pos']}_{svtype}"
        confidence = "capture_limited_low" if capture_limited else "wgs_sv_requires_reconstruction"
        cap = "C" if capture_limited else "B_CAUTION"
        events.append({"sample_id": sample_id, "event_id": event_id, "event_type": "SV", "event_subtype": svtype, "chrom": rec["chrom"], "pos": rec["pos"], "ref": rec["ref"], "alt": rec["alt"], "capture_limited": str(capture_limited).lower(), "capture_overlap": _bed_overlap(rec["chrom"], rec["pos"], capture_intervals) if capture_limited else "NA", "source_file": str(sv_vcf), "confidence": confidence, "priority_cap": cap})
        tasks.append({"event_id": event_id, "task": "transcript_protein_reconstruction", "capture_limited": str(capture_limited).lower(), "required_inputs": "GTF;reference_fasta;RNA_junctions_optional", "priority_cap": cap})
    write_tsv(outdir / "sv_events.tsv", events)
    write_tsv(outdir / "raw_events.tsv", events)
    write_tsv(outdir / "raw_peptides.tsv", [], ["sample_id", "event_id", "peptide_id", "peptide", "hla_allele", "source_type", "generation_status"])
    write_tsv(outdir / ("sv_wes_confidence.tsv" if capture_limited else "sv_reconstruction_tasks.tsv"), tasks)
    skill = "neoag-sv-wes" if capture_limited else "neoag-sv-wgs"
    res = {"status": "PASS", "skill": skill, "outputs": {"raw_events": str(outdir / "raw_events.tsv"), "raw_peptides": str(outdir / "raw_peptides.tsv")}, "summary": f"Parsed {len(events)} SV events; capture_limited={capture_limited}"}
    write_json(outdir / "skill_result.json", res)
    return res


def run_sv_wgs(args: dict[str, Any]) -> dict[str, Any]:
    return _run_sv(args, capture_limited=False)


def run_sv_wes(args: dict[str, Any]) -> dict[str, Any]:
    outdir = ensure_dir(args["outdir"])
    cap_bed = args.get("capture_bed")
    if not cap_bed:
        res = {"status": "FAIL", "skill": "neoag-sv-wes", "failure_reason": "MISSING_CAPTURE_BED", "summary": "WES/capture-limited SV workflow requires --capture-bed"}
        write_json(outdir / "skill_result.json", res)
        return res
    if not Path(cap_bed).exists():
        res = {"status": "FAIL", "skill": "neoag-sv-wes", "failure_reason": "CAPTURE_BED_NOT_FOUND", "summary": f"capture_bed not found: {cap_bed}"}
        write_json(outdir / "skill_result.json", res)
        return res
    return _run_sv(args, capture_limited=True)


def run_peptide_csv(args: dict[str, Any]) -> dict[str, Any]:
    outdir = ensure_dir(args["outdir"])
    peptide_csv = Path(args.get("peptide_csv") or args.get("input") or "")
    sample_id = args.get("sample_id") or peptide_csv.stem
    _, rows = read_table(peptide_csv)
    peptides = []
    pres = []
    for idx, row in enumerate(rows, 1):
        pep = row_get(row, ["peptide", "mt_peptide", "mutant_peptide", "sequence"], "")
        hla = normalize_hla(row_get(row, ["hla_allele", "hla", "allele", "HLA"], ""))
        gene = row_get(row, ["gene", "symbol", "Gene"], "")
        event_id = row_get(row, ["event_id", "variant_id", "mutation", "Mutation"], "") or f"PEPTIDE_INPUT_{idx}_{gene or 'NA'}"
        pid = row_get(row, ["peptide_id", "id"], "") or f"{event_id}_{pep}_{hla or 'HLA_NA'}"
        valid = pep and all(c in "ACDEFGHIKLMNPQRSTVWY" for c in pep.upper())
        peptides.append({"sample_id": sample_id, "event_id": event_id, "peptide_id": pid, "peptide": pep.upper(), "hla_allele": hla, "gene": gene, "source_type": "peptide_csv", "valid_peptide": str(bool(valid)).lower(), "generation_status": "provided"})
        ic50 = row_get(row, ["ic50", "mt_ic50", "netmhcpan_mt_ic50", "affinity"], "")
        el = row_get(row, ["el_rank", "netmhcpan_mt_rank_el", "rank_el", "presentation_rank"], "")
        ba = row_get(row, ["ba_rank", "netmhcpan_mt_rank_ba", "rank_ba"], "")
        if ic50 or el or ba:
            pres.append({"sample_id": sample_id, "peptide_id": pid, "event_id": event_id, "peptide": pep.upper(), "hla_allele": hla, "gene": gene, "ic50": ic50, "el_rank": el, "ba_rank": ba, "source_tool": "peptide_csv"})
    write_tsv(outdir / "raw_peptides.tsv", peptides)
    write_tsv(outdir / "presentation_evidence.tsv", pres)
    write_tsv(outdir / "peptide_input_qc.tsv", [{"metric": "rows", "value": len(rows)}, {"metric": "valid_peptides", "value": sum(1 for p in peptides if p["valid_peptide"] == "true")}, {"metric": "presentation_rows", "value": len(pres)}])
    res = {"status": "PASS", "skill": "neoag-peptide-csv", "outputs": {"raw_peptides": str(outdir / "raw_peptides.tsv"), "presentation": str(outdir / "presentation_evidence.tsv")}, "summary": f"Normalized {len(peptides)} peptide-HLA rows"}
    write_json(outdir / "skill_result.json", res)
    return res
