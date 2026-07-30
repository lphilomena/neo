from __future__ import annotations

from pathlib import Path

import pytest

from neoag.splice.adapters.immunopepper import parse_immunopepper_meta
from neoag.splice.adapters.irfinder import parse_irfinder
from neoag.splice.adapters.pvacbind import parse_pvacbind
from neoag.splice.adapters.spladder import parse_spladder_gff3, parse_spladder_txt
from neoag.splice.identifiers import splice_event_id
from neoag.splice.pipeline import build_splice_provenance_layer
from neoag.utils import read_tsv, write_tsv


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def _inputs(tmp_path: Path) -> dict[str, Path]:
    reg = _write(
        tmp_path / "regtools.tsv",
        "chrom\tstart\tend\tname\tscore\tstrand\tgene_names\n"
        "chr1\t150\t200\tJ1\t17\t+\tGENE1\n"
        "chr1\t500\t600\tHIGH_OTHER\t9000\t+\tGENE1\n",
    )
    gff = _write(
        tmp_path / "merge_graphs_exon_skip_C3.confirmed.gff3",
        "##gff-version 3\n"
        "chr1\tSplAdder\tgene\t100\t250\t.\t+\t.\tID=EV1;gene_name=GENE1;gene_id=G1\n"
        "chr1\tSplAdder\tmRNA\t100\t250\t.\t+\t.\tID=EV1.iso1;Parent=EV1;role=alternative\n"
        "chr1\tSplAdder\texon\t100\t150\t.\t+\t.\tID=e1;Parent=EV1.iso1\n"
        "chr1\tSplAdder\texon\t201\t250\t.\t+\t.\tID=e2;Parent=EV1.iso1\n"
        "chr1\tSplAdder\tmRNA\t100\t250\t.\t+\t.\tID=EV1.iso2;Parent=EV1;role=retained\n"
        "chr1\tSplAdder\texon\t100\t250\t.\t+\t.\tID=e3;Parent=EV1.iso2\n",
    )
    meta = _write(
        tmp_path / "sample_peptides_meta.tsv",
        "peptide\tid\treadFrame\treadFrameAnnotated\tgeneName\tgeneChr\tgeneStrand\tmutationMode\thasStopCodon\tisJunctionList\tvariantSegExpr\tmodifiedExonsCoord\n"
        "MARNDCEQGHILKFPSTWYV\tTX1\t0\tTrue\tGENE1\tchr1\t+\tsplice\tFalse\tTrue\t12.5\tchr1:100-150;chr1:201-250\n",
    )
    normal_cov = tmp_path / "normal_coverage.tsv"
    write_tsv(normal_cov, [{
        "junction_id": "SJ|GRCh38|chr1|151|200|+",
        "detection_status": "NOT_DETECTED",
        "coverage_status": "ADEQUATE",
        "normal_source_type": "MATCHED_NORMAL",
        "normal_tissue": "blood",
    }])
    return {"reg": reg, "gff": gff, "meta": meta, "normal_cov": normal_cov}


def test_stable_event_id_is_order_invariant():
    a = splice_event_id(genome_build="GRCh38", event_type="SE", strand="+", junction_ids=["J2", "J1"], gene="G")
    b = splice_event_id(genome_build="GRCh38", event_type="SE", strand="+", junction_ids=["J1", "J2"], gene="G")
    assert a == b
    assert a.startswith("SEV|")


def test_spladder_gff3_builds_event_paths_and_exact_junction(tmp_path: Path):
    bundle = parse_spladder_gff3(_inputs(tmp_path)["gff"], sample_id="S1")
    assert len(bundle["events"]) == 1
    assert bundle["events"][0]["event_type"] == "SE"
    assert bundle["events"][0]["junction_ids"] == "SJ|GRCh38|chr1|151|200|+"
    assert len(bundle["transcripts"]) == 2
    assert {row["path_role"] for row in bundle["transcripts"]} == {"alternative", "retained"}
    assert any(row["entity_id"] == "SJ|GRCh38|chr1|151|200|+" and row["evidence_group"] == "SPLICE_GRAPH" for row in bundle["tool_evidence"])


def test_irfinder_registers_retained_and_spliced_hypotheses(tmp_path: Path):
    ir = _write(
        tmp_path / "IRFinder-IR-nondir.txt",
        "Chr\tStart\tEnd\tName\tStrand\tIntronDepth\tSpliceLeft\tSpliceRight\tSpliceExact\tIRratio\tWarnings\tCNN\n"
        "chr2\t301\t400\tGENE2\t-\t8.2\t3\t4\t7\t0.25\t\tPASS\n",
    )
    bundle = parse_irfinder(ir, sample_id="S1")
    assert bundle["events"][0]["event_type"] == "RI"
    assert bundle["junctions"][0]["junction_id"] == "SJ|GRCh38|chr2|301|400|-"
    assert {row["path_role"] for row in bundle["transcripts"]} == {"SPLICED", "RETAINED"}
    assert bundle["events"][0]["psi"] == "0.25"


def test_immunopepper_registers_partial_transcript_orf_and_origin(tmp_path: Path):
    bundle = parse_immunopepper_meta(_inputs(tmp_path)["meta"], sample_id="S1")
    assert len(bundle["events"]) == 1
    assert len(bundle["transcripts"]) == 1
    assert len(bundle["orfs"]) == 1
    assert len(bundle["peptide_origins"]) == 1
    assert bundle["transcripts"][0]["full_length_status"] == "PARTIAL_JUNCTION_TRANSLATION"
    assert bundle["orfs"][0]["orf_validity_status"] == "PARTIAL_TRANSLATED_SEGMENT"
    assert bundle["peptide_origins"][0]["crosses_junction"] == "true"


def test_pvacbind_requires_exact_fasta_index(tmp_path: Path):
    mapping = tmp_path / "map.tsv"
    write_tsv(mapping, [{
        "index": "SPORF_1", "orf_id": "ORF|1", "transcript_hypothesis_id": "STH|1",
        "splice_event_id": "SEV|1", "sample_id": "S1", "gene": "G",
    }])
    pvac = _write(
        tmp_path / "all_epitopes.tsv",
        "Index\tHLA Allele\tEpitope Seq\tSub-peptide Position\tBest IC50 Score\tBest Percentile\n"
        "SPORF_1\tHLA-A*02:01\tARNDCEQGH\t2\t25\t0.3\n"
        "UNKNOWN\tHLA-A*02:01\tRNDCEQGHI\t3\t40\t0.5\n",
    )
    entity_bundle = {
        "orfs": [{"orf_id": "ORF|1", "transcript_hypothesis_id": "STH|1", "splice_event_id": "SEV|1", "gene": "G", "protein_sequence": "MARNDCEQGHJK"}],
        "transcripts": [{"transcript_hypothesis_id": "STH|1", "junction_chain": "SJ|GRCh38|chr1|1|2|+"}],
        "peptide_origins": [],
    }
    bundle = parse_pvacbind(pvac, sample_id="S1", fasta_map=mapping, entity_bundle=entity_bundle)
    assert len(bundle["presentation"]) == 1
    assert bundle["presentation"][0]["mapping_status"] == "MAPPED_EXACT_FASTA_INDEX"
    assert any(row["conflict_type"] == "PVACBIND_INDEX_UNRESOLVED" for row in bundle["conflicts"])


def test_full_layer_has_five_level_referential_integrity_and_no_gene_read_leakage(tmp_path: Path):
    inputs = _inputs(tmp_path)
    out1 = tmp_path / "layer1"
    outputs1 = build_splice_provenance_layer(
        sample_id="S1", outdir=out1, junctions=inputs["reg"],
        spladder_gff3=[inputs["gff"]], immunopepper_meta=[inputs["meta"]],
        normal_coverage=[inputs["normal_cov"]], strict=True,
    )
    fasta_map = read_tsv(outputs1["pvacbind_fasta_map"])
    assert len(fasta_map) == 1
    pvac = _write(
        tmp_path / "pvacbind.tsv",
        "Index\tHLA Allele\tEpitope Seq\tSub-peptide Position\tBest IC50 Score\tBest Percentile\tImmunogenicity Score\n"
        f"{fasta_map[0]['index']}\tHLA-A*02:01\tARNDCEQGH\t2\t20\t0.2\t0.7\n",
    )
    out2 = tmp_path / "layer2"
    outputs = build_splice_provenance_layer(
        sample_id="S1", outdir=out2, junctions=inputs["reg"],
        spladder_gff3=[inputs["gff"]], immunopepper_meta=[inputs["meta"]],
        pvacbind=[pvac], normal_coverage=[inputs["normal_cov"]], strict=True,
    )
    junctions = {row["junction_id"]: row for row in read_tsv(outputs["junctions"])}
    assert junctions["SJ|GRCh38|chr1|151|200|+"]["total_split_reads"] == "17"
    assert junctions["SJ|GRCh38|chr1|501|600|+"]["total_split_reads"] == "9000"
    origins = read_tsv(outputs["peptide_origins"])
    presentations = read_tsv(outputs["presentation"])
    assert origins and presentations
    event_ids = {row["splice_event_id"] for row in read_tsv(outputs["events"])}
    transcript_ids = {row["transcript_hypothesis_id"] for row in read_tsv(outputs["transcripts"])}
    orf_ids = {row["orf_id"] for row in read_tsv(outputs["orfs"])}
    for origin in origins:
        assert origin["splice_event_id"] in event_ids
        assert origin["transcript_hypothesis_id"] in transcript_ids
        assert origin["orf_id"] in orf_ids
    qc = read_tsv(outputs["qc"])
    assert not [row for row in qc if row["status"] == "FAIL"]
    consensus = read_tsv(outputs["consensus"])
    assert any(row["event_evidence_grade"] == "E2" for row in consensus)
    assert all(row["orf_evidence_grade"] == "O1" for row in consensus)
    assert any("CAP_SINGLE_PEPTIDE_GENERATOR_R2" in row["cap_codes"] for row in consensus)
    raw = read_tsv(outputs["raw_peptides"])
    assert any(row["origin_peptide_id"] and row["orf_id"] and row["transcript_hypothesis_id"] for row in raw)


def test_strict_irfinder_rejects_unsupported_coordinate_system(tmp_path: Path):
    ir = _write(tmp_path / "ir.tsv", "Chr\tStart\tEnd\tName\tStrand\tIRratio\nchr1\t10\t20\tG\t+\t0.2\n")
    with pytest.raises(Exception):
        parse_irfinder(ir, sample_id="S1", coordinate_system="unknown", strict=True)


def test_spladder_headerless_build_txt_is_parsed_without_feature_columns_becoming_exons(tmp_path: Path):
    txt = _write(
        tmp_path / "merge_graphs_exon_skip_C3.confirmed.txt",
        "chr3\t+\tEV42\tGENE3\t100\t150\t201\t250\tS1:1\tS1:8\n",
    )
    bundle = parse_spladder_txt(txt, sample_id="S1")
    assert len(bundle["events"]) == 1
    assert bundle["events"][0]["event_type"] == "SE"
    assert bundle["events"][0]["junction_ids"] == "SJ|GRCh38|chr3|151|200|+"
    assert bundle["transcripts"][0]["exon_chain"] == "chr3:100-150:+;chr3:201-250:+"


def test_immunopepper_official_semicolon_coordinates_and_isolated_semantics(tmp_path: Path):
    meta = _write(
        tmp_path / "ref_sample_peptides_meta.tsv",
        "peptide\tid\treadFrame\treadFrameAnnotated\tgeneName\tgeneChr\tgeneStrand\tmutationMode\thasStopCodon\tisJunctionList\tisIsolated\tmodifiedExonsCoord\n"
        "MARNDCEQGHILKFPSTWYV\tTX2\t0\tTrue\tGENE4\tchr4\t+\tref\t0\t1\t0\t100;150;201;250\n",
    )
    bundle = parse_immunopepper_meta(meta, sample_id="S1")
    assert bundle["events"][0]["junction_ids"] == "SJ|GRCh38|chr4|151|200|+"
    origin = bundle["peptide_origins"][0]
    assert origin["crosses_junction"] == "true"
    assert origin["junction_offset_in_peptide"] == "17-18"
    assert origin["contains_novel_aa"] == "UNASSESSED"


def test_pvacbind_rejects_epitope_that_is_not_in_exact_mapped_orf(tmp_path: Path):
    mapping = tmp_path / "map.tsv"
    write_tsv(mapping, [{
        "index": "SPORF_1", "orf_id": "ORF|1", "transcript_hypothesis_id": "STH|1",
        "splice_event_id": "SEV|1", "sample_id": "S1", "gene": "G",
    }])
    pvac = _write(
        tmp_path / "all_epitopes.tsv",
        "Index\tHLA Allele\tEpitope Seq\tSub-peptide Position\tBest IC50 Score\n"
        "SPORF_1\tHLA-A*02:01\tYYYYYYYYY\t2\t25\n",
    )
    entity_bundle = {
        "orfs": [{"orf_id": "ORF|1", "transcript_hypothesis_id": "STH|1", "splice_event_id": "SEV|1", "gene": "G", "protein_sequence": "MARNDCEQGHILKFPSTWYV"}],
        "transcripts": [{"transcript_hypothesis_id": "STH|1", "junction_chain": "SJ|GRCh38|chr1|1|2|+"}],
        "peptide_origins": [],
    }
    bundle = parse_pvacbind(pvac, sample_id="S1", fasta_map=mapping, entity_bundle=entity_bundle)
    assert not bundle["presentation"]
    assert any(row["conflict_type"] == "PVACBIND_EPITOPE_ORF_MISMATCH" for row in bundle["conflicts"])
