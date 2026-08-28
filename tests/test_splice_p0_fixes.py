from __future__ import annotations

from pathlib import Path

from neoag.splice.adapters.pvacsplice import parse_pvacsplice
from neoag.splice.adapters.regtools import parse_junction_source
from neoag.splice.consensus import build_consensus
from neoag.splice.evidence_chains import build_evidence_chains
from neoag.splice.pipeline import build_splice_provenance_layer
from neoag.utils import read_tsv


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_complete_junction_qc_rejects_one_read_and_missing_quality(tmp_path: Path):
    source = _write(
        tmp_path / "junctions.tsv",
        "chrom\tstart\tend\tstrand\tunique_split_reads\n"
        "chr1\t150\t200\t+\t1\n",
    )
    bundle = parse_junction_source(source, sample_id="S1", coordinate_system="regtools_annotated")
    qc = bundle["junction_read_qc"][0]
    assert qc["qc_status"] == "FAIL"
    assert "UNIQUE_SPLIT_READS_BELOW_THRESHOLD" in qc["failed_checks"]
    assert "MAX_OVERHANG" in qc["missing_checks"]


def test_complete_junction_qc_passes_only_with_all_required_metrics(tmp_path: Path):
    source = _write(
        tmp_path / "junctions.tsv",
        "chrom\tstart\tend\tstrand\tunique_split_reads\tmulti_split_reads\tunique_fragment_starts\tmax_overhang\tmedian_mapq\tmultimapping_fraction\ttumor_psi\tcaller_filter\n"
        "chr1\t150\t200\t+\t7\t1\t4\t18\t45\t0.125\t0.12\tPASS\n",
    )
    bundle = parse_junction_source(source, sample_id="S1", coordinate_system="regtools_annotated")
    assert bundle["junction_read_qc"][0]["qc_status"] == "PASS"


def test_peptide_cannot_borrow_support_from_another_junction_in_event():
    j1, j2 = "SJ|GRCh38|chr1|101|200|+", "SJ|GRCh38|chr1|101|300|+"
    tables = {
        "events": [{"splice_event_id": "SEV|1", "alternative_junction_ids": f"{j1};{j2}"}],
        "event_junction_links": [
            {"splice_event_id": "SEV|1", "junction_id": j1},
            {"splice_event_id": "SEV|1", "junction_id": j2},
        ],
        "orfs": [{"orf_id": "ORF|1", "splice_event_id": "SEV|1", "protein_sequence_sha256": "abc",
                  "frame_status": "IN_FRAME", "orf_validity_status": "VALID_FULL_LENGTH", "source_generator": "G1"}],
        "peptide_origins": [{"origin_peptide_id": "POR|1", "peptide_id": "PEP|1", "orf_id": "ORF|1",
                             "splice_event_id": "SEV|1", "junction_ids": j2, "required_junction_ids": j2,
                             "crosses_junction": "true", "contains_novel_aa": "true"}],
        "junction_read_qc": [{"junction_id": j1, "qc_status": "PASS", "resolution_status": "RESOLVED_EXACT"}],
        "tool_evidence": [{"entity_id": j1, "evidence_group": "RNA_JUNCTION", "source_assay_id": "RNA1",
                           "verified_value": "10", "resolution_status": "RESOLVED_EXACT"}],
        "presentation": [{"origin_peptide_id": "POR|1"}],
        "normal_background": [],
    }
    row = build_consensus(tables, sample_id="S1")[0]
    assert row["event_evidence_grade"] == "E0"
    assert row["required_junction_ids"] == j2
    assert row["required_junction_qc_status"] == "MISSING"
    assert "CAP_REQUIRED_JUNCTION_QC_NOT_PASS_R3" in row["cap_codes"]


def _pvac_entities() -> dict[str, list[dict[str, str]]]:
    jid = "SJ|GRCh38|chr1|151|200|+"
    event = "SEV|1"
    variant = "VAR|GRCh38|chr1|1000|A|G"
    return {
        "events": [{"splice_event_id": event, "junction_ids": jid, "alternative_junction_ids": jid, "gene": "G"}],
        "event_junction_links": [{"splice_event_id": event, "junction_id": jid}],
        "causal_links": [{"causal_link_id": "DCL|1", "variant_id": variant, "junction_id": jid,
                          "splice_event_id": event, "causal_status": "DNA_RNA_CIS_SUPPORTED"}],
    }


def test_pvacsplice_does_not_invent_crossing_or_novel_residues(tmp_path: Path):
    report = _write(
        tmp_path / "pvac.tsv",
        "Chromosome\tStart\tStop\tReference\tVariant\tJunction\tGene Name\tHLA Allele\tEpitope Seq\n"
        "chr1\t999\t1000\tA\tG\tchr1:151-200:+\tG\tHLA-A*02:01\tARNDCEQGH\n",
    )
    parsed = parse_pvacsplice(report, sample_id="S1", entity_bundle=_pvac_entities(), strict=True)
    origin = parsed["peptide_origins"][0]
    assert origin["crosses_junction"] == "UNASSESSED"
    assert origin["contains_novel_aa"] == "UNASSESSED"
    assert origin["novel_aa_positions"] == ""
    assert origin["required_junction_ids"] == "SJ|GRCh38|chr1|151|200|+"


def test_pvacsplice_uses_explicit_residue_and_boundary_positions(tmp_path: Path):
    report = _write(
        tmp_path / "pvac.tsv",
        "Chromosome\tStart\tStop\tReference\tVariant\tJunction\tGene Name\tHLA Allele\tEpitope Seq\tPos\tJunction AA Position\n"
        "chr1\t999\t1000\tA\tG\tchr1:151-200:+\tG\tHLA-A*02:01\tARNDCEQGH\t4\t5\n",
    )
    parsed = parse_pvacsplice(report, sample_id="S1", entity_bundle=_pvac_entities(), strict=True)
    origin = parsed["peptide_origins"][0]
    assert origin["contains_novel_aa"] == "true"
    assert origin["novel_aa_positions"] == "4"
    assert origin["crosses_junction"] == "true"
    assert origin["junction_offset_in_peptide"] == "5-6"


def test_normal_detected_state_is_consistent_between_chain_and_consensus():
    jid = "SJ|GRCh38|chr1|151|200|+"
    tables = {
        "events": [{"splice_event_id": "SEV|1", "alternative_junction_ids": jid}],
        "event_junction_links": [{"splice_event_id": "SEV|1", "junction_id": jid}],
        "orfs": [{"orf_id": "ORF|1", "splice_event_id": "SEV|1", "protein_sequence_sha256": "abc",
                  "frame_status": "IN_FRAME", "orf_validity_status": "VALID_FULL_LENGTH", "source_generator": "G1"}],
        "peptide_origins": [{"origin_peptide_id": "POR|1", "peptide_id": "PEP|1", "orf_id": "ORF|1",
                             "splice_event_id": "SEV|1", "required_junction_ids": jid,
                             "crosses_junction": "true", "contains_novel_aa": "true"}],
        "normal_background": [{"normal_background_id": "NBG|1", "junction_id": jid,
                               "assessment_status": "DETECTED_CRITICAL_TISSUE", "critical_tissue": "true"}],
        "tool_evidence": [], "junction_read_qc": [], "presentation": [],
    }
    normal_chain = next(row for row in build_evidence_chains(tables, sample_id="S1") if row["chain_type"] == "NORMAL_BACKGROUND")
    consensus = build_consensus({**tables, "evidence_chains": [normal_chain]}, sample_id="S1")[0]
    assert normal_chain["chain_status"] == "NORMAL_DETECTED_CRITICAL"
    assert consensus["normal_safety_grade"] == "N0"
    assert "HARD_NORMAL_BACKGROUND_DETECTED" in consensus["hard_fail_codes"]


def test_fallback_event_is_unclassified_not_novel(tmp_path: Path):
    source = _write(
        tmp_path / "known.tsv",
        "chrom\tstart\tend\tstrand\tunique_split_reads\tunique_fragment_starts\tmax_overhang\tmedian_mapq\tmultimapping_fraction\ttumor_psi\tknown_junction\n"
        "chr1\t150\t200\t+\t7\t4\t18\t45\t0.1\t0.2\tKNOWN\n",
    )
    outputs = build_splice_provenance_layer(
        sample_id="S1", outdir=tmp_path / "out", junctions=source,
        junction_coordinate_system="regtools_annotated",
        tool_versions={"RegTools": "fixture"}, strict=True,
    )
    event = read_tsv(outputs["events"])[0]
    assert event["event_type"] == "JUNCTION_ONLY_UNCLASSIFIED"
    assert event["alternative_junction_ids"] == ""


def test_epitope_only_orf_is_capped_at_r3():
    jid = "SJ|GRCh38|chr1|151|200|+"
    tables = {
        "events": [{"splice_event_id": "SEV|1", "alternative_junction_ids": jid}],
        "event_junction_links": [{"splice_event_id": "SEV|1", "junction_id": jid}],
        "orfs": [{"orf_id": "ORF|1", "splice_event_id": "SEV|1", "protein_sequence_sha256": "abc",
                  "frame_status": "IN_FRAME", "orf_validity_status": "VALID_EPITOPE_PRODUCT_ONLY",
                  "source_generator": "pVACsplice"}],
        "peptide_origins": [{"origin_peptide_id": "POR|1", "peptide_id": "PEP|1", "orf_id": "ORF|1",
                             "splice_event_id": "SEV|1", "required_junction_ids": jid,
                             "crosses_junction": "true", "contains_novel_aa": "true"}],
        "junction_read_qc": [{"junction_id": jid, "qc_status": "PASS", "resolution_status": "RESOLVED_EXACT"}],
        "tool_evidence": [{"entity_id": jid, "evidence_group": "RNA_JUNCTION", "source_assay_id": "RNA1",
                           "verified_value": "10", "resolution_status": "RESOLVED_EXACT"}],
        "normal_background": [{"normal_background_id": "NBG|1", "junction_id": jid,
                               "assessment_status": "NOT_DETECTED_ADEQUATE_COVERAGE",
                               "normal_source_type": "MATCHED_NORMAL", "normal_source": "N1"}],
        "presentation": [{"origin_peptide_id": "POR|1"}],
    }
    row = build_consensus(tables, sample_id="S1")[0]
    assert row["orf_evidence_grade"] == "O1"
    assert row["final_evidence_tier"] == "R3"
    assert "CAP_PARTIAL_OR_EPITOPE_ONLY_R3" in row["cap_codes"]
