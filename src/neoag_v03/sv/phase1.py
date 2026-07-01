from __future__ import annotations

import datetime
from pathlib import Path
from typing import Iterable, Any

from ..schemas import EVENT_FIELDS, PEPTIDE_FIELDS
from ..config import load_profile
from ..model_layers import enrich_event_layers, enrich_peptide_layers, infer_mutation_source, infer_peptide_consequence
from ..utils import write_tsv, write_json, safe_id, clamp
from .schemas_sv import SV_EVENT_FULL_FIELDS, SV_PROTEIN_FIELDS, SV_EVENT_TO_PEPTIDE_FIELDS, SV_VALIDATION_DESIGN_FIELDS
from .sv_callset import read_sv_inputs
from .sv_merge import cluster_sv_records
from .sv_filter import score_cluster_confidence, passes_phase1
from .reference import FastaReference, GtfAnnotation
from .protein_reconstruct import reconstruct_cluster_protein, ProteinReconstruction
from .peptide_builder import build_mhc1_peptides, read_hla_alleles
from .evidence import load_expression, load_normal_expression, normal_expr_for_genes, load_normal_ligands, load_junction_reads
from .wes_capture import CaptureRegions, write_expanded_beds, annotate_cluster_capture




_WES_PRIORITY_CAP_DEFAULTS = {
    "WES_Tier1": "B",
    "WES_Tier2": "B_CAUTION",
    "WES_Tier3": "C",
    "WES_UNINTERPRETABLE": "D",
}

_WES_PRIORITY_CAP_KEYS = {
    "WES_Tier1": "wes_tier1_priority_cap",
    "WES_Tier2": "wes_tier2_priority_cap",
    "WES_Tier3": "wes_tier3_priority_cap",
    "WES_UNINTERPRETABLE": "wes_uninterpretable_cap",
}


def _wes_priority_cap_from_profile(profile: dict[str, Any], tier: str) -> str:
    caps = profile.get("wes_confidence_caps", {}) or {}
    key = _WES_PRIORITY_CAP_KEYS.get(tier, "wes_uninterpretable_cap")
    fallback = _WES_PRIORITY_CAP_DEFAULTS.get(tier, _WES_PRIORITY_CAP_DEFAULTS["WES_UNINTERPRETABLE"])
    return str(caps.get(key) or fallback)

def _event_id(sample_id: str, svtype: str, gene1: str, gene2: str, chrom1: str, pos1: int, chrom2: str, pos2: int) -> str:
    gene_part = f"{gene1}_{gene2}" if gene2 and gene2 != gene1 else gene1 or "NA"
    return safe_id(f"SV_{sample_id}_{svtype}_{gene_part}_{chrom1}_{pos1}_{chrom2}_{pos2}")


def _event_expression(effect_class: str, gene1: str, gene2: str, expr: dict[str, float], rna_reads: int) -> float:
    if rna_reads > 0:
        # Junction evidence is not TPM; use a bounded proxy only for event-level scoring.
        return min(20.0, max(1.0, float(rna_reads)))
    if effect_class == "SV_Fusion" and gene1 and gene2:
        vals = [expr.get(g, 0.0) for g in (gene1, gene2) if g]
        return min(vals) if vals else 0.0
    return expr.get(gene1, 0.0)


def _sv_event_row(cluster, event_id: str, sample_id: str, meta: dict[str, Any], conf, protein: ProteinReconstruction | None, rna_reads: int, extra: dict[str, Any] | None = None) -> dict[str, Any]:
    rep = cluster.representative
    t_alt = max(cluster.max_int("tumor_alt_support"), rep.tumor_alt_support)
    t_depth = max(cluster.max_int("tumor_local_depth"), rep.tumor_local_depth)
    n_alt = max(cluster.max_int("normal_alt_support"), rep.normal_alt_support)
    n_depth = max(cluster.max_int("normal_local_depth"), rep.normal_local_depth)
    vaf = (t_alt / t_depth) if t_depth else rep.vaf_like
    extra = extra or {}
    return {
        "sv_event_id": event_id,
        "event_id": event_id,
        "sample_id": sample_id,
        "svtype": rep.svtype,
        "chrom1": rep.chrom1,
        "pos1": rep.pos1,
        "strand1": rep.strand1,
        "chrom2": rep.chrom2,
        "pos2": rep.pos2,
        "strand2": rep.strand2,
        "bnd_alt": rep.alt,
        "cipos": rep.cipos,
        "ciend": rep.ciend,
        "svlen": rep.svlen,
        "inserted_sequence": rep.inserted_sequence,
        "callers": ",".join(cluster.callers),
        "caller_count": len(cluster.callers),
        "record_ids": ",".join(cluster.record_ids),
        "tumor_sr": cluster.max_int("tumor_sr"),
        "tumor_pe": cluster.max_int("tumor_pe"),
        "tumor_alt_support": t_alt,
        "tumor_local_depth": t_depth,
        "normal_sr": cluster.max_int("normal_sr"),
        "normal_pe": cluster.max_int("normal_pe"),
        "normal_alt_support": n_alt,
        "normal_local_depth": n_depth,
        "sniffles_support": cluster.max_int("sniffles_support"),
        "sniffles_coverage": cluster.max_int("sniffles_coverage"),
        "sniffles_precise": rep.sniffles_precise,
        "sniffles_rnames_count": cluster.max_int("sniffles_rnames_count"),
        "sv_vaf_like": f"{clamp(vaf):.4f}",
        "pon_overlap": "not_assessed",
        "population_sv_overlap": "not_assessed",
        "blacklist_overlap": "not_assessed",
        "breakpoint_precision_bp": cluster.best_precision() or rep.breakpoint_precision_bp,
        "event_confidence_tier": conf.tier,
        "event_confidence_score": f"{conf.score:.4f}",
        "gene1": meta.get("gene1", ""),
        "gene2": meta.get("gene2", ""),
        "transcript1": meta.get("transcript1", ""),
        "transcript2": meta.get("transcript2", ""),
        "exon1": "",
        "exon2": "",
        "cds_phase1": "",
        "cds_phase2": "",
        "effect_class": meta.get("effect_class", "SV_Noncoding"),
        "fusion_in_frame": meta.get("fusion_in_frame", "unknown"),
        "frameshift": meta.get("frameshift", "unknown"),
        "protein_sequence_id": protein.protein_sequence_id if protein else "",
        "junction_aa_position": meta.get("junction_aa_position", ""),
        "rna_junction_reads": rna_reads,
        "rna_support_status": "RNA_JUNCTION_SUPPORTED" if rna_reads >= 3 else ("RNA_WEAK_JUNCTION" if rna_reads > 0 else "RNA_NOT_AVAILABLE_OR_NOT_DETECTED"),
        "final_sv_confidence": conf.tier,
        "reconstruction_status": meta.get("reconstruction_status", "failed"),
        "reconstruction_reason": meta.get("reconstruction_reason", ""),
        "evidence_scope": extra.get("evidence_scope", "GENOME_WIDE"),
        "breakend1_capture_status": extra.get("breakend1_capture_status", ""),
        "breakend2_capture_status": extra.get("breakend2_capture_status", ""),
        "breakend1_capture_distance_bp": extra.get("breakend1_capture_distance_bp", ""),
        "breakend2_capture_distance_bp": extra.get("breakend2_capture_distance_bp", ""),
        "capture_interpretability": extra.get("capture_interpretability", ""),
        "wes_confidence_tier": extra.get("wes_confidence_tier", conf.tier if str(conf.tier).startswith("WES_") else ""),
        "priority_cap": extra.get("priority_cap", ""),
        "filter_status": extra.get("filter_status", "PASS"),
        "filter_reason": extra.get("filter_reason", ""),
    }


def _raw_event_from_sv(svrow: dict[str, Any], protein: ProteinReconstruction | None, expr: dict[str, float], normal_expr: dict[str, dict[str, Any]], profile_name: str) -> dict[str, Any]:
    gene1 = str(svrow.get("gene1") or "")
    gene2 = str(svrow.get("gene2") or "")
    genes = [g for g in [gene1, gene2] if g]
    effect_class = str(svrow.get("effect_class") or "SV_Noncoding")
    rna_reads = int(float(svrow.get("rna_junction_reads") or 0))
    exp = _event_expression(effect_class, gene1, gene2, expr, rna_reads)
    ne = normal_expr_for_genes(genes, normal_expr)
    t_alt = float(svrow.get("tumor_alt_support") or 0)
    t_depth = float(svrow.get("tumor_local_depth") or 0)
    vaf = float(svrow.get("sv_vaf_like") or 0)
    specificity = 1.0
    if float(svrow.get("normal_alt_support") or 0) > 0:
        specificity -= 0.3
    if str(svrow.get("population_sv_overlap", "not_assessed")) == "yes":
        specificity -= 0.4
    gene_label = f"{gene1}::{gene2}" if effect_class == "SV_Fusion" and gene2 and gene2 != gene1 else gene1
    base = {
        "event_id": svrow["event_id"],
        "sample_id": svrow["sample_id"],
        "disease_profile": profile_name,
        "event_type": effect_class,
        "mutation_source": "SV",
        "peptide_consequence": infer_peptide_consequence(event_type=effect_class),
        "evidence_scope": svrow.get("evidence_scope", "GENOME_WIDE"),
        "priority_cap": svrow.get("priority_cap", ""),
        "wes_confidence_tier": svrow.get("wes_confidence_tier", ""),
        "gene": gene_label,
        "event_name": gene_label or svrow["event_id"],
        "chrom": svrow.get("chrom1", ""),
        "pos": svrow.get("pos1", ""),
        "ref": "N",
        "alt": svrow.get("bnd_alt") or f"<{effect_class}>",
        "transcript_id": protein.transcript_id if protein else svrow.get("transcript1", ""),
        "consequence": "fusion_transcript" if effect_class == "SV_Fusion" else ("frameshift_variant" if effect_class == "SV_Frameshift" else "protein_altering_variant"),
        "rna_junction_reads": str(rna_reads),
        "event_confidence": svrow.get("event_confidence_score", "0"),
        "event_expression": f"{exp:.4f}",
        "driver_relevance": "0.5000",
        "tumor_vaf": f"{vaf:.4f}",
        "tumor_depth": f"{t_depth:.0f}",
        "tumor_alt_count": f"{t_alt:.0f}",
        "clonality": f"{min(1.0, vaf / 0.35 if vaf > 0 else 0.0):.4f}",
        "persistence": "0.7000" if effect_class in {"SV_Fusion", "SV_Frameshift"} else "0.5000",
        "tumor_specificity": f"{clamp(specificity):.4f}",
        "ccf_estimate": "",
        "ccf_status": "",
        "clonality_multiplier": "1.0",
        "normal_tissue_max_tpm": f"{ne['normal_tissue_max_tpm']:.4f}",
        "normal_hspc_tpm": f"{ne['normal_hspc_tpm']:.4f}",
        "critical_tissue_hit": ne["critical_tissue_hit"],
        "safety_status": "",
        "safety_reason": "",
        "appm_mhc_i_integrity": "",
        "appm_mhc_ii_integrity": "",
        "event_score": "",
        "evidence_scope": svrow.get("evidence_scope", "GENOME_WIDE"),
        "priority_cap": svrow.get("priority_cap", ""),
        "wes_confidence_tier": svrow.get("wes_confidence_tier", ""),
        "source": ("WES_SV_Phase1_5:" if svrow.get("evidence_scope") == "EXOME_CAPTURE_LIMITED" else "SVPhase1:") + str(svrow.get('callers','')),
    }
    return enrich_event_layers(base)


def _raw_peptide_from_candidate(cand, effect_class: str, gene_label: str, rna_reads: int = 0, evidence_scope: str = "", priority_cap: str = "", wes_confidence_tier: str = "") -> dict[str, Any]:
    base = {
        "peptide_id": cand.peptide_id,
        "event_id": cand.event_id,
        "sample_id": cand.sample_id,
        "event_type": effect_class,
        "mutation_source": "SV",
        "peptide_consequence": infer_peptide_consequence(event_type=effect_class),
        "evidence_scope": evidence_scope,
        "priority_cap": priority_cap,
        "wes_confidence_tier": wes_confidence_tier,
        "gene": gene_label,
        "peptide": cand.peptide,
        "wildtype_peptide": cand.wildtype_peptide,
        "crosses_junction": cand.crosses_junction,
        "contains_novel_aa": cand.contains_novel_aa,
        "rna_junction_reads": str(rna_reads),
        "hla_allele": cand.hla_allele,
        "mhc_class": cand.mhc_class,
        "source_tool": "SVPhase1",
        "binding_rank": "99",
        "el_rank": "99",
        "presentation_score": "0",
        "immunogenicity_score": "",
        "wildtype_binding_rank": "99",
        "self_similarity_score": "0",
        "normal_hla_ligand_overlap": cand.normal_hla_ligand_overlap,
        "netmhcpan_ba_rank": "",
        "netmhcpan_el_rank": "",
        "netmhcstabpan_score": "",
        "netmhcstabpan_rank": "",
        "mhcflurry_affinity_percentile": "",
        "mhcflurry_processing_score": "",
        "mhcflurry_presentation_score": "",
        "binding_evidence_score": "",
        "presentation_evidence_score": "",
        "presentation_evidence_grade": "",
        "iedb_immunogenicity_score": "",
        "immunogenicity_resolved": "no",
        "prime_score": "",
        "prime_rank": "",
        "bigmhc_im_score": "",
        "immunogenicity_composite_score": "",
        "immunogenicity_source": "",
        "presentation_gate_status": "",
        "presentation_gate_reason": "",
        "presentation_gate_multiplier": "",
        "appm_multiplier": "",
        "ccf_multiplier": "",
        "safety_status": "",
        "safety_reason": "",
        "efficacy_score": "",
        "final_priority": "",
        "recommended_use": "",
        "safety_tier": "",
        "safety_multiplier": "",
        "review_required": "",
        "reference_proteome_exact_match": "",
        "normal_ligand_tissue": "",
        "mutation_anchor_only": "",
        "escape_status": "",
        "escape_reason": "",
        "escape_multiplier": "",
        "restricting_hla_lost": "",
    }
    return enrich_peptide_layers(base)


def write_fasta(path: str | Path, rows: list[ProteinReconstruction]) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", encoding="utf-8") as fh:
        for r in rows:
            seq = r.protein_sequence or ""
            if not seq:
                continue
            fh.write(f">{r.protein_sequence_id}|{r.event_id}|{r.gene}|{r.protein_type}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + "\n")


def _int_or_zero(value: Any) -> int:
    try:
        return int(float(value))
    except Exception:
        return 0


def _long_design_from_protein(protein: ProteinReconstruction) -> dict[str, Any]:
    seq = (protein.protein_sequence or "").replace("*", "").upper()
    if not seq:
        return {}
    junction = _int_or_zero(protein.junction_aa_position)
    novel_start = _int_or_zero(protein.novel_start_aa)
    center = junction if protein.protein_type == "SV_Fusion" and junction > 0 else novel_start
    flank = 15
    start = max(0, center - flank)
    end = min(len(seq), center + flank)
    if end - start < 15:
        start = max(0, min(start, len(seq) - 31))
        end = min(len(seq), max(end, start + 31))
    long_peptide = seq[start:end]
    if protein.protein_type == "SV_Fusion":
        design_type = "fusion_junction_long_peptide_or_minigene"
        reason = "Fusion/BND candidates should validate the processed junction region, not only short exogenous peptides."
    elif protein.protein_type == "SV_Frameshift":
        design_type = "frameshift_novel_tail_long_peptide_or_minigene"
        reason = "Frameshift candidates should validate the novel tail sequence and processing context."
    elif protein.protein_type == "SV_Insertion":
        design_type = "insertion_long_peptide_or_minigene"
        reason = "Insertion candidates should validate the inserted/novel sequence in local context."
    else:
        design_type = "sv_junction_long_peptide_or_minigene"
        reason = "SV junction candidates benefit from long peptide or minigene validation."
    return {
        "event_id": protein.event_id,
        "sample_id": protein.sample_id,
        "protein_sequence_id": protein.protein_sequence_id,
        "gene": protein.gene,
        "protein_type": protein.protein_type,
        "design_type": design_type,
        "long_peptide": long_peptide,
        "long_peptide_length": len(long_peptide),
        "minigene_aa": long_peptide,
        "junction_aa_position": protein.junction_aa_position,
        "novel_start_aa": protein.novel_start_aa,
        "design_reason": reason,
    }


def build_sv_phase1_raw(
    *,
    sample_id: str,
    sv_vcfs: Iterable[str | Path],
    reference_fasta: str | Path,
    gencode_gtf: str | Path,
    hla: str | Path | list[str],
    outdir: str | Path,
    profile_name: str = "sv_wgs_phase1",
    callers: Iterable[str] | None = None,
    tumor_sample_name: str | None = None,
    normal_sample_name: str | None = None,
    expression_tsv: str | Path | None = None,
    rna_junction_tsv: str | Path | None = None,
    normal_expression_tsv: str | Path | None = None,
    normal_hla_ligands_tsv: str | Path | None = None,
    merge_distance_bp: int = 200,
    allow_tier2: bool = True,
    wes_mode: bool = False,
    capture_bed: str | Path | None = None,
    capture_near_bp: int = 250,
    capture_slop_bp: int = 1000,
) -> dict[str, str]:
    outdir = Path(outdir)
    parsed_dir = outdir / "parsed"
    sv_dir = outdir / "sv"
    parsed_dir.mkdir(parents=True, exist_ok=True)
    sv_dir.mkdir(parents=True, exist_ok=True)

    profile = load_profile(profile_name)
    hla_alleles = read_hla_alleles(hla)
    if not hla_alleles:
        raise ValueError("At least one HLA allele is required for SV peptide generation.")
    ref = FastaReference(reference_fasta)
    ann = GtfAnnotation(gencode_gtf)
    expr = load_expression(expression_tsv)
    normal_expr = load_normal_expression(normal_expression_tsv)
    normal_ligands = load_normal_ligands(normal_hla_ligands_tsv)
    junction_reads = load_junction_reads(rna_junction_tsv)
    capture_regions = CaptureRegions.from_bed(capture_bed) if (wes_mode and capture_bed) else None
    capture_paths = write_expanded_beds(capture_bed, outdir / "capture", near_bp=capture_near_bp, slop_bp=capture_slop_bp) if (wes_mode and capture_bed) else {}

    if wes_mode:
        from . import wes_filter as filt
    else:
        from . import sv_filter as filt

    records = read_sv_inputs(sv_vcfs, callers, tumor_sample_name=tumor_sample_name, normal_sample_name=normal_sample_name)
    clusters = cluster_sv_records(records, distance=merge_distance_bp)

    sv_rows: list[dict[str, Any]] = []
    proteins: list[ProteinReconstruction] = []
    validation_designs: list[dict[str, Any]] = []
    raw_events: list[dict[str, Any]] = []
    sidecar_peptides = []
    raw_peptides: list[dict[str, Any]] = []

    for cluster in clusters:
        rep = cluster.representative
        provisional_gene1 = ann.gene_at(rep.chrom1, rep.pos1)
        provisional_gene2 = ann.gene_at(rep.chrom2, rep.pos2)
        gene_pair = f"{provisional_gene1}::{provisional_gene2}" if provisional_gene1 and provisional_gene2 else ""
        rna_reads = max(
            junction_reads.get(gene_pair, 0),
            junction_reads.get(f"{provisional_gene2}::{provisional_gene1}", 0),
        )
        cap_meta = annotate_cluster_capture(cluster, capture_regions, near_bp=capture_near_bp, slop_bp=capture_slop_bp) if wes_mode else {}
        conf = filt.score_wes_confidence(cluster, rna_junction_reads=rna_reads, capture_interpretability=cap_meta.get("capture_interpretability")) if wes_mode else score_cluster_confidence(cluster, rna_junction_reads=rna_reads)
        if not (filt.passes_wes_phase1_5(cluster, conf, allow_tier2=allow_tier2) if wes_mode else passes_phase1(cluster, conf, allow_tier2=allow_tier2)):
            continue
        eid = _event_id(sample_id, rep.svtype, provisional_gene1, provisional_gene2, rep.chrom1, rep.pos1, rep.chrom2, rep.pos2)
        # Event-specific lookup can override provisional RNA evidence.
        rna_reads = max(rna_reads, junction_reads.get(eid, 0))
        cap_meta = annotate_cluster_capture(cluster, capture_regions, near_bp=capture_near_bp, slop_bp=capture_slop_bp) if wes_mode else {}
        conf = filt.score_wes_confidence(cluster, rna_junction_reads=rna_reads, capture_interpretability=cap_meta.get("capture_interpretability")) if wes_mode else score_cluster_confidence(cluster, rna_junction_reads=rna_reads)
        protein, meta = reconstruct_cluster_protein(cluster, eid, sample_id, ann, ref)
        priority_cap = _wes_priority_cap_from_profile(profile, conf.tier) if wes_mode else ""
        extra = {**cap_meta, "evidence_scope": "EXOME_CAPTURE_LIMITED" if wes_mode else "GENOME_WIDE", "wes_confidence_tier": conf.tier if wes_mode else "", "priority_cap": priority_cap}
        svrow = _sv_event_row(cluster, eid, sample_id, meta, conf, protein, rna_reads, extra)
        sv_rows.append(svrow)
        if not protein or meta.get("reconstruction_status") != "ok":
            continue
        if meta.get("effect_class") not in {"SV_Fusion", "SV_Frameshift", "SV_Junction", "SV_Insertion"}:
            continue
        proteins.append(protein)
        design = _long_design_from_protein(protein)
        if design:
            validation_designs.append(design)
        raw_event = _raw_event_from_sv(svrow, protein, expr, normal_expr, profile_name)
        raw_events.append(raw_event)
        cands = build_mhc1_peptides([protein], hla_alleles, normal_ligands=normal_ligands)
        sidecar_peptides.extend(cands)
        gene_label = raw_event["gene"]
        for cand in cands:
            raw_peptides.append(_raw_peptide_from_candidate(cand, str(svrow.get("effect_class")), gene_label, rna_reads, raw_event.get("evidence_scope", ""), raw_event.get("priority_cap", ""), svrow.get("wes_confidence_tier", "")))

    raw_events_path = parsed_dir / "raw_events.tsv"
    raw_peptides_path = parsed_dir / "raw_peptides.tsv"
    sv_events_path = sv_dir / "sv_events.full.tsv"
    proteins_path = sv_dir / "sv_protein_reconstruction.tsv"
    proteins_fa_path = sv_dir / "sv_mutant_proteins.fa"
    sidecar_peptides_path = sv_dir / "sv_event_to_peptide.tsv"
    validation_design_path = sv_dir / "sv_validation_design.tsv"
    write_tsv(raw_events_path, raw_events, EVENT_FIELDS)
    write_tsv(raw_peptides_path, raw_peptides, PEPTIDE_FIELDS)
    write_tsv(sv_events_path, sv_rows, SV_EVENT_FULL_FIELDS)
    write_tsv(proteins_path, [p.to_row() for p in proteins], SV_PROTEIN_FIELDS)
    write_fasta(proteins_fa_path, proteins)
    write_tsv(sidecar_peptides_path, [p.to_sidecar_row() for p in sidecar_peptides], SV_EVENT_TO_PEPTIDE_FIELDS)
    write_tsv(validation_design_path, validation_designs, SV_VALIDATION_DESIGN_FIELDS)
    prov_name = "provenance.sv_wes_phase1_5.json" if wes_mode else "provenance.sv_phase1.json"
    prov_path = outdir / prov_name
    prov_payload = {
        "created_at": datetime.datetime.now(datetime.UTC).isoformat().replace("+00:00", "Z"),
        "sample_id": sample_id,
        "profile": profile_name,
        "sv_vcfs": [str(p) for p in sv_vcfs],
        "reference_fasta": str(reference_fasta),
        "gencode_gtf": str(gencode_gtf),
        "hla_alleles": hla_alleles,
        "records_parsed": len(records),
        "clusters": len(clusters),
        "sv_events_written": len(sv_rows),
        "raw_events_written": len(raw_events),
        "raw_peptides_written": len(raw_peptides),
        "validation_designs_written": len(validation_designs),
        "warning": "SV Phase 1 is a computational triage adapter; validate SVs and RNA junctions experimentally.",
        "capture_bed": str(capture_bed) if capture_bed else None,
        "capture_sidecars": capture_paths,
    }
    if wes_mode:
        prov_payload["evidence_scope"] = "EXOME_CAPTURE_LIMITED"
        prov_payload["mode"] = "wes_phase1_5"
        prov_payload["warning"] = (
            "WES SV Phase 1.5 is exome-capture-limited triage; RNA junction support strongly "
            "recommended. Validate SVs experimentally."
        )
    write_json(prov_path, prov_payload)
    return {
        "raw_events": str(raw_events_path),
        "raw_peptides": str(raw_peptides_path),
        "sv_events_full": str(sv_events_path),
        "sv_protein_reconstruction": str(proteins_path),
        "sv_mutant_proteins": str(proteins_fa_path),
        "sv_event_to_peptide": str(sidecar_peptides_path),
        "sv_validation_design": str(validation_design_path),
        "provenance": str(prov_path),
    }
