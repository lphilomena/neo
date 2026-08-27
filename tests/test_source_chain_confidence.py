from __future__ import annotations

from pathlib import Path

from neoag.evidence_consensus import build_evidence_consensus, load_consensus_rules
from neoag.source_chain import (
    CONFLICT,
    INDETERMINATE_LOW_POWER,
    NOT_APPLICABLE,
    SUPPORTED,
    UNASSESSED,
    build_source_chain_table,
    derive_source_chain_confidence,
    source_chain_track,
)
from neoag.utils import read_tsv, write_tsv


def _base_peptide(pid: str, track: str) -> dict[str, str]:
    row = {
        "sample_id": "S1",
        "peptide_id": pid,
        "event_id": f"E_{pid}",
        "gene": "GENE1",
        "peptide": "ABCDEFGH",
        "hla_allele": "HLA-A*02:01",
        "presentation_gate_status": "PASS",
        "netmhcpan_mt_rank_el": "0.2",
        "mhcflurry_presentation_score": "0.8",
        "safety_status": "PASS",
        "safety_evidence_completeness": "COMPLETE",
        "restricting_hla_lost": "false",
        "hla_loh_status": "RETAINED",
        "appm_integrity_status": "INTACT",
        "appm_evidence_completeness": "HIGH",
        "l3_apm_integrity_score": "0.9",
        "ccf_estimate": "0.9",
        "ccf_status": "clonal_like",
        "ccf_confidence": "high_confidence",
        "mutant_specificity_status": "MT_SPECIFIC",
        "mutant_specificity_gate_status": "MT_SPECIFIC",
        "mutation_positions_in_peptide": "4",
        "filter_status": "PASS",
    }
    if track == "SNV":
        row.update({
            "mutation_source": "SNV",
            "event_type": "SNV",
            "cross_platform_status": "CROSS_PLATFORM_PASS_CONCORDANT",
            "tumor_depth": "100",
            "tumor_alt_count": "20",
            "normal_depth": "80",
            "normal_alt_count": "0",
            "rna_alt_reads": "12",
            "rna_vaf": "0.25",
            "rna_depth": "48",
            "transcript_id": "ENST000001",
            "protein_change": "p.Ala10Val",
            "mean_alt_base_quality": "35",
            "mean_alt_mapping_quality": "60",
            "strand_bias_status": "PASS",
            "ffpe_artifact_status": "PASS",
            "low_complexity_status": "PASS",
            "paralogous_region_status": "PASS",
        })
    elif track == "INDEL":
        row.update({
            "mutation_source": "InDel",
            "event_type": "frameshift_variant",
            "cross_platform_status": "CROSS_PLATFORM_PASS_CONCORDANT",
            "tumor_depth": "120",
            "tumor_alt_count": "15",
            "normal_depth": "90",
            "normal_alt_count": "0",
            "rna_alt_reads": "8",
            "rna_vaf": "0.20",
            "rna_depth": "40",
            "normalized_variant_key": "chr1:10:A:AT",
            "local_realign_status": "PASS",
            "repeat_context_status": "PASS",
            "frame_status": "VALID",
            "orf_id": "ORF1",
            "contains_novel_aa": "true",
            "novel_tail_status": "VALID",
            "stop_position": "25",
            "nmd_risk_status": "ASSESSED_LOW_RISK",
            "transcript_id": "ENST000001",
        })
    elif track == "FUSION":
        row.update({
            "mutation_source": "RNA_ONLY",
            "event_type": "Fusion",
            "fusion_gene": "GENE1::GENE2",
            "gene5": "GENE1",
            "gene3": "GENE2",
            "breakpoint1": "chr1:100:+",
            "breakpoint2": "chr2:200:-",
            "rna_junction_reads": "22",
            "rna_spanning_reads": "4",
            "tools_detected": "EasyFuse,Arriba",
            "normal_junction_assessment_status": "NOT_DETECTED_ADEQUATE_COVERAGE",
            "normal_junction_depth": "40",
            "readthrough_status": "EXCLUDED",
            "duplicate_removed_status": "PASS",
            "junction_sequence_uniqueness_status": "PASS",
            "independent_start_sites": "8",
            "frame_status": "IN_FRAME",
            "rna_frame_status": "IN_FRAME",
            "fusion_protein_sequence": "MABCDEFGHIK",
            "orf_id": "FUSION_ORF1",
            "transcript_hypothesis_id": "FT1",
            "crosses_junction": "true",
            "contains_novel_aa": "true",
            "junction_position_in_peptide_1based": "4",
            "fusion_left_peptide": "ABCD",
            "fusion_right_peptide": "EFGH",
            "ccf_status": "RNA_ONLY_UNRESOLVED",
            "ccf_estimate": "NA",
        })
    elif track == "SPLICE":
        row.update({
            "mutation_source": "Splice",
            "event_type": "splice_junction",
            "genome_build": "GRCh38",
            "canonical_junction_id": "GRCh38|chr1|100|200|+",
            "junction_chrom": "chr1",
            "junction_start": "100",
            "junction_end": "200",
            "junction_strand": "+",
            "rna_junction_reads": "30",
            "unique_split_reads": "30",
            "anchor_size": "20",
            "mapping_quality": "60",
            "junction_support_status": "PASS",
            "splice_event_type": "SE",
            "reference_path_status": "ANNOTATED",
            "alternative_path_status": "RESOLVED",
            "junction_match_status": "MATCHED",
            "transcript_hypothesis_id": "ST1",
            "orf_id": "SORF1",
            "orf_evidence_grade": "PASS",
            "frame_status": "IN_FRAME",
            "translation_direction_status": "SENSE",
            "nmd_risk_status": "ASSESSED_LOW_RISK",
            "psi": "0.25",
            "junction_usage_status": "SUPPORTED",
            "known_normal_isoform_status": "EXCLUDED",
            "crosses_junction": "true",
            "contains_novel_aa": "true",
            "normal_junction_assessment_status": "NOT_DETECTED_ADEQUATE_COVERAGE",
            "normal_junction_depth": "35",
            "ccf_status": "RNA_ONLY_UNRESOLVED",
            "ccf_estimate": "NA",
        })
    return row


def _requirements_by_name(result):
    return {req.name: req for req in result.requirements}


def test_source_chain_track_uses_upstream_event_source():
    row = _base_peptide("P1", "INDEL")
    row["peptide_consequence"] = "splice_junction"
    assert source_chain_track(row) == "INDEL"


def test_source_chain_track_prefers_declared_snv_over_mixed_source_label():
    row = {
        "event_type": "SNV",
        "mutation_source": "SNV_INDEL",
        "peptide_consequence": "splice_junction",
    }
    assert source_chain_track(row) == "SNV"


def test_source_chain_track_keeps_sv_fusion_on_dna_sv_track():
    row = {
        "event_type": "SV_Fusion",
        "mutation_source": "SV",
    }
    assert source_chain_track(row) == "DNA_SV"


def test_snv_complete_cross_modal_chain_is_c1():
    result = derive_source_chain_confidence(_base_peptide("P1", "SNV"), {})
    assert result.track == "SNV"
    assert result.tier == "C1"
    assert "DNA/RNA-or-protein-cross-modal" in result.orthogonal_sources
    assert _requirements_by_name(result)["read_backed_phasing"].status == NOT_APPLICABLE


def test_snv_dna_only_is_c3_not_c4():
    row = _base_peptide("P2", "SNV")
    row.update({"rna_alt_reads": "0", "rna_vaf": "0", "rna_depth": "60"})
    result = derive_source_chain_confidence(row, {})
    assert result.tier == "C3"
    assert "rna_or_direct_evidence" in result.negative_requirements
    assert result.hard_failure is False


def test_snv_phasing_that_disproves_candidate_is_c4():
    row = _base_peptide("P3", "SNV")
    row.update({
        "phasing_required": "true",
        "proximal_variant_count": "1",
        "haplotype_status": "DISPROVED_TRANS_INCOMPATIBLE",
    })
    result = derive_source_chain_confidence(row, {})
    assert result.tier == "C4"
    assert "SC_PHASING_DISPROVES_PEPTIDE" in result.hard_failure_codes


def test_indel_complete_chain_is_c1_and_local_context_is_explicit():
    result = derive_source_chain_confidence(_base_peptide("P4", "INDEL"), {})
    assert result.tier == "C1"
    reqs = _requirements_by_name(result)
    assert reqs["normalized_representation"].status == SUPPORTED
    assert reqs["local_context_qc"].status == SUPPORTED
    assert reqs["orf_reconstruction"].status == SUPPORTED


def test_indel_unresolved_local_realign_is_c3():
    row = _base_peptide("P5", "INDEL")
    row["local_realign_status"] = "LOW_CONFIDENCE_REPEAT_RISK"
    result = derive_source_chain_confidence(row, {})
    assert result.tier == "C3"
    assert "local_context_qc" in result.low_power_requirements


def test_indel_wrong_frame_is_c4():
    row = _base_peptide("P6", "INDEL")
    row["frame_status"] = "WRONG_FRAME"
    result = derive_source_chain_confidence(row, {})
    assert result.tier == "C4"
    assert "SC_INDEL_ORF_INVALID" in result.hard_failure_codes


def test_fusion_complete_computational_chain_without_orthogonal_is_c2():
    result = derive_source_chain_confidence(_base_peptide("P7", "FUSION"), {})
    assert result.track == "FUSION"
    assert result.tier == "C2"
    assert result.orthogonal_status == UNASSESSED


def test_fusion_label_without_residue_boundary_mapping_is_not_supported():
    row = _base_peptide("P7_missing_boundary", "FUSION")
    row["junction_position_in_peptide_1based"] = ""
    row["fusion_left_peptide"] = ""
    row["fusion_right_peptide"] = ""
    result = derive_source_chain_confidence(row, {})
    assert result.tier == "C4"
    assert "SC_FUSION_BOUNDARY_MAPPING_INCOMPLETE" in result.hard_failure_codes


def test_fusion_rt_pcr_confirmation_promotes_to_c1():
    row = _base_peptide("P8", "FUSION")
    row["rt_pcr_status"] = "CONFIRMED"
    row["sanger_status"] = "CONFIRMED"
    result = derive_source_chain_confidence(row, {})
    assert result.tier == "C1"
    assert "RT-PCR" in result.orthogonal_sources
    assert "Sanger" in result.orthogonal_sources


def test_fusion_gene_pair_without_breakpoint_is_c3():
    row = _base_peptide("P9", "FUSION")
    row["breakpoint1"] = ""
    row["breakpoint2"] = ""
    result = derive_source_chain_confidence(row, {})
    assert result.tier == "C3"
    assert "breakpoint_definition" in result.missing_requirements


def test_fusion_readthrough_artifact_is_c4():
    row = _base_peptide("P10", "FUSION")
    row["readthrough_status"] = "READTHROUGH_CONFIRMED"
    result = derive_source_chain_confidence(row, {})
    assert result.tier == "C4"
    assert "SC_FUSION_NORMAL_BACKGROUND_REFUTES" in result.hard_failure_codes


def test_splice_complete_chain_without_orthogonal_is_c2():
    result = derive_source_chain_confidence(_base_peptide("P11", "SPLICE"), {})
    assert result.track == "SPLICE"
    assert result.tier == "C2"
    reqs = _requirements_by_name(result)
    assert reqs["junction_definition"].status == SUPPORTED
    assert reqs["junction_mapping_qc"].status == SUPPORTED
    assert reqs["reference_alternative_paths"].status == SUPPORTED


def test_splice_low_anchor_is_c3_low_power():
    row = _base_peptide("P12", "SPLICE")
    row["anchor_size"] = "5"
    result = derive_source_chain_confidence(row, {})
    assert result.tier == "C3"
    assert "junction_mapping_qc" in result.low_power_requirements


def test_splice_wrong_junction_fill_is_c4():
    row = _base_peptide("P13", "SPLICE")
    row["junction_match_status"] = "WRONG_JUNCTION_FILL"
    result = derive_source_chain_confidence(row, {})
    assert result.tier == "C4"
    assert "SC_SPLICE_PATH_LINK_INVALID" in result.hard_failure_codes


def test_statuses_do_not_conflate_not_applicable_unassessed_and_negative():
    snv = derive_source_chain_confidence(_base_peptide("P14", "SNV"), {})
    assert _requirements_by_name(snv)["read_backed_phasing"].status == NOT_APPLICABLE

    fusion_row = _base_peptide("P15", "FUSION")
    fusion_row["rna_spanning_reads"] = ""
    fusion_row["tools_detected"] = ""
    fusion_row["independent_start_sites"] = ""
    fusion_row["rna_junction_reads"] = "3"
    fusion = derive_source_chain_confidence(fusion_row, {})
    assert _requirements_by_name(fusion)["supplementary_structure_support"].status in {UNASSESSED, INDETERMINATE_LOW_POWER}

    snv_negative = _base_peptide("P16", "SNV")
    snv_negative.update({"rna_alt_reads": "0", "rna_vaf": "0", "rna_depth": "100"})
    negative = derive_source_chain_confidence(snv_negative, {})
    assert _requirements_by_name(negative)["rna_or_direct_evidence"].status == "NEGATIVE"


def _legacy_complete_row(pid: str) -> dict[str, str]:
    # Existing EC-v2 complete row intentionally lacks several source-chain audit fields.
    return {
        "peptide_id": pid,
        "event_id": f"E{pid}",
        "event_type": "SNV",
        "efficacy_score": "0.5",
        "presentation_gate_status": "PASS",
        "netmhcpan_mt_rank_el": "0.2",
        "mhcflurry_presentation_score": "0.8",
        "rna_alt_reads": "8",
        "rna_vaf": "0.12",
        "ccf_status": "clonal_like",
        "ccf_estimate": "0.9",
        "ccf_confidence": "high_confidence",
        "cross_platform_status": "CROSS_PLATFORM_PASS_CONCORDANT",
        "mutant_specificity_gate_status": "PASS",
        "safety_evidence_completeness": "COMPLETE",
        "safety_status": "PASS",
        "l3_apm_integrity_score": "0.9",
        "appm_evidence_completeness": "HIGH",
        "restricting_hla_lost": "false",
        "hla_loh_status": "RETAINED",
    }


def test_compatibility_profile_emits_source_chain_without_changing_r_grade(tmp_path: Path):
    source = tmp_path / "input.tsv"
    write_tsv(source, [_legacy_complete_row("P_COMPAT")])
    rules = load_consensus_rules(Path(__file__).parents[1] / "configs/ranking/sarcoma_evidence_consensus_v2_1_source_chain.toml")
    result = build_evidence_consensus(source, tmp_path / "compat", rules)
    row = read_tsv(result["ranked_peptides"])[0]
    assert row["source_chain_confidence_tier"] == "C3"
    assert row["source_chain_integration_mode"] == "compatibility"
    assert row["evidence_grade"] == "R1"
    assert Path(result["source_chain_confidence"]).is_file()
    assert Path(result["source_chain_requirements"]).is_file()


def test_compatibility_profile_preserves_rna_only_fusion_cap(tmp_path: Path):
    source = tmp_path / "fusion.tsv"
    row = _base_peptide("P_COMPAT_FUSION", "FUSION")
    write_tsv(source, [row])
    rules = load_consensus_rules(
        Path(__file__).parents[1]
        / "configs/ranking/sarcoma_evidence_consensus_v2_1_source_chain.toml"
    )
    result = build_evidence_consensus(source, tmp_path / "compat_fusion", rules)
    ranked = read_tsv(result["ranked_peptides"])[0]
    assert ranked["source_chain_confidence_tier"] == "C2"
    assert ranked["source_chain_integration_mode"] == "compatibility"
    assert "CAP_RNA_ONLY_FUSION" in ranked["evidence_grade_cap_reasons"]


def test_integrated_profile_caps_c3_at_r3_and_adds_pareto_dimension(tmp_path: Path):
    source = tmp_path / "input.tsv"
    write_tsv(source, [_legacy_complete_row("P_INTEGRATED")])
    rules = load_consensus_rules(Path(__file__).parents[1] / "configs/ranking/sarcoma_evidence_consensus_v3_source_chain.toml")
    result = build_evidence_consensus(source, tmp_path / "integrated", rules)
    row = read_tsv(result["ranked_peptides"])[0]
    assert row["source_chain_confidence_tier"] == "C3"
    assert row["source_chain_integration_mode"] == "integrated"
    assert row["evidence_grade"] == "R3"
    assert "CAP_SOURCE_CHAIN_C3" in row["evidence_grade_cap_reasons"]
    assert "source_chain_confidence_grade" in row["pareto_dimensions"]


def test_integrated_complete_rna_only_fusion_is_not_automatically_capped_r3(tmp_path: Path):
    source = tmp_path / "fusion.tsv"
    write_tsv(source, [_base_peptide("P_FUSION_C2", "FUSION")])
    rules = load_consensus_rules(
        Path(__file__).parents[1] / "configs/ranking/sarcoma_evidence_consensus_v3_source_chain.toml"
    )
    result = build_evidence_consensus(source, tmp_path / "fusion_integrated", rules)
    row = read_tsv(result["ranked_peptides"])[0]
    assert row["source_chain_confidence_tier"] == "C2"
    assert "CAP_RNA_ONLY_FUSION" not in row["evidence_grade_cap_reasons"]
    assert row["evidence_grade"] in {"R1", "R2"}


def test_source_chain_cli_table_builder_writes_long_audit(tmp_path: Path):
    source = tmp_path / "source.tsv"
    output = tmp_path / "source_chain.tsv"
    requirements = tmp_path / "requirements.tsv"
    write_tsv(source, [_base_peptide("P17", "FUSION")])
    result = build_source_chain_table(source, output, requirements)
    assert result["rows"] == 1
    assert read_tsv(output)[0]["source_chain_confidence_tier"] == "C2"
    long_rows = read_tsv(requirements)
    assert any(row["requirement_name"] == "breakpoint_definition" for row in long_rows)
    assert all(row["requirement_applicability"] in {"APPLICABLE", "NOT_APPLICABLE"} for row in long_rows)
    assert all("requirement_conflict" in row for row in long_rows)


def test_normal_junction_not_detected_without_coverage_is_low_power():
    row = _base_peptide("P17B", "FUSION")
    row["normal_junction_assessment_status"] = "NOT_DETECTED"
    row["normal_junction_depth"] = ""
    result = derive_source_chain_confidence(row, {})
    assert result.tier == "C3"
    req = _requirements_by_name(result)["normal_background"]
    assert req.status == INDETERMINATE_LOW_POWER
    assert req.reason_code == "SC_FUSION_NORMAL_BACKGROUND_COVERAGE_UNASSESSED"


def test_snv_explicit_ffpe_artifact_is_c4():
    row = _base_peptide("P17C", "SNV")
    row["ffpe_artifact_status"] = "FFPE_ARTIFACT"
    result = derive_source_chain_confidence(row, {})
    assert result.tier == "C4"
    assert "SC_SNV_SEQUENCE_QC_ARTIFACT" in result.hard_failure_codes


def test_splice_known_normal_isoform_is_c4():
    row = _base_peptide("P17D", "SPLICE")
    row["known_normal_isoform_status"] = "KNOWN_NORMAL_ISOFORM_CONFIRMED"
    result = derive_source_chain_confidence(row, {})
    assert result.tier == "C4"
    assert "SC_SPLICE_KNOWN_NORMAL_ISOFORM" in result.hard_failure_codes


def test_orthogonal_not_confirmed_is_unassessed_not_c4():
    row = _base_peptide("P18", "FUSION")
    row["orthogonal_confirmation_status"] = "NOT_CONFIRMED"
    result = derive_source_chain_confidence(row, {})
    assert result.tier == "C2"
    assert result.orthogonal_status == UNASSESSED
    assert result.hard_failure is False


def test_same_bam_multi_caller_is_not_orthogonal_confirmation():
    row = _base_peptide("P19", "FUSION")
    row["tools_detected"] = "EasyFuse,Arriba,STAR-Fusion"
    result = derive_source_chain_confidence(row, {})
    assert result.tier == "C2"
    assert result.orthogonal_status == UNASSESSED
    assert not result.orthogonal_sources
