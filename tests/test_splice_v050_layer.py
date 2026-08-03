from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from neoag.splice.adapters.immunopepper import parse_immunopepper_kmers, parse_immunopepper_meta
from neoag.splice.adapters.high_order import parse_high_order_evidence
from neoag.splice.adapters.irfinder import parse_irfinder
from neoag.splice.adapters.pvacbind import parse_pvacbind
from neoag.splice.adapters.regtools import parse_junction_source
from neoag.splice.adapters.spladder import infer_event_type, parse_spladder_gff3, parse_spladder_txt
from neoag.splice.coordinates import JunctionNormalizationError, read_source_rows
from neoag.splice.identifiers import sequence_sha256, splice_event_id
from neoag.splice.pipeline import build_splice_provenance_layer
from neoag.splice.consensus import build_consensus, consensus_reason_conflicts
from neoag.splice.normal_background import parse_normal_coverage, parse_normal_junctions
from neoag.splice.junction_queries import build_canonical_junction_queries
from neoag.utils import read_tsv, write_tsv


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_gzip_junction_rows_are_streamed(tmp_path: Path):
    source = tmp_path / "normal.tsv.gz"
    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write("junction_id\tchromosome\tstart\tend\tstrand\tnormal_reads\n")
        handle.write("chr1:151-200:+\tchr1\t151\t200\t+\t7\n")
        handle.write("chr2:301-400:-\tchr2\t301\t400\t-\t2\n")
    rows = read_source_rows(source)
    assert not isinstance(rows, list)
    assert [row["junction_id"] for row in rows] == ["chr1:151-200:+", "chr2:301-400:-"]


def test_normal_junction_panel_keeps_only_exact_targets(tmp_path: Path):
    source = tmp_path / "normal.tsv.gz"
    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write("junction_id\tchromosome\tstart\tend\tstrand\tnormal_reads\n")
        handle.write("chr1:151-200:+\tchr1\t151\t200\t+\t7\n")
        handle.write("chr2:301-400:-\tchr2\t301\t400\t-\t2\n")
    parsed = parse_normal_junctions(
        source,
        sample_id="S1",
        allowed_junction_ids={"SJ|GRCh38|chr1|151|200|+"},
    )
    assert [row["junction_id"] for row in parsed["normal_background"]] == ["SJ|GRCh38|chr1|151|200|+"]
    assert parsed["manifest"][0]["rows_scanned"] == "2"
    assert parsed["manifest"][0]["rows_retained"] == "1"


def test_canonical_junction_queries_use_exact_peptide_origins_and_strand(tmp_path: Path):
    sequence = ("ACGT" * 30)[:120]
    fasta = tmp_path / "ref.fa"
    fasta.write_text(f">1\n{sequence}\n", encoding="ascii")
    (tmp_path / "ref.fa.fai").write_text(f"1\t120\t3\t120\t121\n", encoding="ascii")
    junctions = [
        {"junction_id": "SJ|GRCh38|chr1|41|60|+", "chrom": "chr1", "intron_start_1based": "41", "intron_end_1based": "60", "strand": "+"},
        {"junction_id": "SJ|GRCh38|chr1|41|60|-", "chrom": "chr1", "intron_start_1based": "41", "intron_end_1based": "60", "strand": "-"},
        {"junction_id": "SJ|GRCh38|chr1|70|80|+", "chrom": "chr1", "intron_start_1based": "70", "intron_end_1based": "80", "strand": "+"},
    ]
    origins = [
        {"splice_event_id": "SEV|plus", "crosses_junction": "true", "junction_ids": junctions[0]["junction_id"]},
        {"splice_event_id": "SEV|minus", "crosses_junction": "true", "junction_ids": junctions[1]["junction_id"]},
    ]
    parsed = build_canonical_junction_queries(
        {"junctions": junctions, "peptide_origins": origins},
        sample_id="S1", reference_fasta=fasta, flank_bases=10,
    )
    queries = {row["splice_event_id"]: row for row in parsed["sequence_queries"]}
    assert set(queries) == {"SEV|plus", "SEV|minus"}
    assert queries["SEV|plus"]["nucleotide_sequence"] == sequence[30:40] + sequence[60:70]
    expected_minus = (sequence[60:70].translate(str.maketrans("ACGT", "TGCA"))[::-1]
                      + sequence[30:40].translate(str.maketrans("ACGT", "TGCA"))[::-1])
    assert queries["SEV|minus"]["nucleotide_sequence"] == expected_minus
    assert all(row["position_1based"] == "10" for row in queries.values())


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
    assert bundle["events"][0]["alternative_junction_ids"] == "SJ|GRCh38|chr1|151|200|+"
    assert bundle["events"][0]["reference_junction_ids"] == ""
    assert bundle["events"][0]["reference_path_status"] == "ALTERNATIVE_ONLY_REFERENCE_UNRESOLVED"


def test_spladder_uses_only_explicit_reference_and_alternative_roles(tmp_path: Path):
    gff = _write(
        tmp_path / "events.gff3",
        "##gff-version 3\n"
        "chr1\tSplAdder\tevent\t100\t350\t.\t+\t.\tID=EV2;gene_name=G\n"
        "chr1\tSplAdder\tmRNA\t100\t350\t.\t+\t.\tID=EV2.ref;Parent=EV2;role=reference\n"
        "chr1\tSplAdder\texon\t100\t150\t.\t+\t.\tID=r1;Parent=EV2.ref\n"
        "chr1\tSplAdder\texon\t201\t250\t.\t+\t.\tID=r2;Parent=EV2.ref\n"
        "chr1\tSplAdder\tmRNA\t100\t350\t.\t+\t.\tID=EV2.alt;Parent=EV2;role=alternative\n"
        "chr1\tSplAdder\texon\t100\t150\t.\t+\t.\tID=a1;Parent=EV2.alt\n"
        "chr1\tSplAdder\texon\t301\t350\t.\t+\t.\tID=a2;Parent=EV2.alt\n",
    )
    event = parse_spladder_gff3(gff, sample_id="S1")["events"][0]
    assert event["reference_junction_ids"] == "SJ|GRCh38|chr1|151|200|+"
    assert event["alternative_junction_ids"] == "SJ|GRCh38|chr1|151|300|+"
    assert event["reference_path_status"] == "RESOLVED_EXPLICIT_SOURCE_ROLE"


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


def test_immunopepper_stream_filter_keeps_only_exact_target_junction(tmp_path: Path):
    meta = _write(
        tmp_path / "meta.tsv",
        "peptide\tid\tgeneName\tgeneChr\tgeneStrand\tmutationMode\tmodifiedExonsCoord\n"
        "MARNDCEQGHILK\tTX1\tG1\tchr1\t+\tsplice\tchr1:100-150;chr1:201-250\n"
        "ACDEFGHIKLMNP\tTX2\tG2\tchr2\t+\tsplice\tchr2:300-350;chr2:401-450\n",
    )
    bundle = parse_immunopepper_meta(
        meta, sample_id="S1",
        allowed_junction_ids={"SJ|GRCh38|chr1|151|200|+"},
    )
    assert len(bundle["events"]) == 1
    assert bundle["events"][0]["gene"] == "G1"
    assert bundle["manifest"][0]["rows_scanned"] == "2"
    assert bundle["manifest"][0]["rows_retained"] == "1"
    assert bundle["manifest"][0]["filter_policy"] == "EXACT_CANONICAL_JUNCTION_INTERSECTION"


def test_immunopepper_kmer_uses_exact_retained_orf_index(tmp_path: Path):
    meta_bundle = parse_immunopepper_meta(_inputs(tmp_path)["meta"], sample_id="S1")
    kmers = _write(
        tmp_path / "kmers.tsv",
        "kmer\tcoord\tisCrossJunction\tjunctionAnnotated\n"
        "NDCEQGHI\t1:2:3:4\tTrue\tTrue\n"
        "AAAAAAAA\t1:2:3:4\tTrue\tTrue\n",
    )
    bundle = parse_immunopepper_kmers(
        kmers, sample_id="S1", meta_bundle=meta_bundle,
        record_unmapped_conflicts=False,
    )
    assert len(bundle["peptide_origins"]) == 1
    assert bundle["manifest"][0]["rows_scanned"] == "2"
    assert bundle["manifest"][0]["rows_mapped"] == "1"
    assert bundle["manifest"][0]["rows_unmapped"] == "1"
    assert not bundle["conflicts"]


def test_pvacbind_requires_exact_fasta_index(tmp_path: Path):
    protein = "MARNDCEQGHILK"
    mapping = tmp_path / "map.tsv"
    write_tsv(mapping, [{
        "index": "SPORF_1", "orf_id": "ORF|1", "transcript_hypothesis_id": "STH|1",
        "splice_event_id": "SEV|1", "sample_id": "S1", "gene": "G",
        "sequence_sha256": sequence_sha256(protein),
    }])
    pvac = _write(
        tmp_path / "all_epitopes.tsv",
        "Index\tHLA Allele\tEpitope Seq\tSub-peptide Position\tBest IC50 Score\tBest Percentile\n"
        "SPORF_1\tHLA-A*02:01\tARNDCEQGH\t2\t25\t0.3\n"
        "UNKNOWN\tHLA-A*02:01\tRNDCEQGHI\t3\t40\t0.5\n",
    )
    entity_bundle = {
        "events": [{"splice_event_id": "SEV|1"}],
        "orfs": [{"orf_id": "ORF|1", "transcript_hypothesis_id": "STH|1", "splice_event_id": "SEV|1", "gene": "G", "protein_sequence": protein, "protein_sequence_sha256": sequence_sha256(protein)}],
        "transcripts": [{"transcript_hypothesis_id": "STH|1", "splice_event_id": "SEV|1", "junction_chain": "SJ|GRCh38|chr1|1|2|+"}],
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
        normal_coverage=[inputs["normal_cov"]],
        tool_versions={"RegTools": "fixture", "SplAdder": "fixture", "ImmunoPepper": "fixture"},
        strict=True,
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
        sample_id="S1", outdir=out2, base_layer_dir=out1,
        pvacbind=[pvac], tool_versions={"pVACbind": "fixture"},
        strict=True,
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


def test_production_layer_requires_explicit_irfinder_coordinate_system(tmp_path: Path):
    ir = _write(tmp_path / "ir.tsv", "Chr\tStart\tEnd\tName\tStrand\tIRratio\nchr1\t10\t20\tG\t+\t0.2\n")
    with pytest.raises(ValueError, match="explicit --irfinder-coordinate-system"):
        build_splice_provenance_layer(sample_id="S1", outdir=tmp_path / "out", irfinder=[ir])


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
    protein = "MARNDCEQGHILKFPSTWYV"
    mapping = tmp_path / "map.tsv"
    write_tsv(mapping, [{
        "index": "SPORF_1", "orf_id": "ORF|1", "transcript_hypothesis_id": "STH|1",
        "splice_event_id": "SEV|1", "sample_id": "S1", "gene": "G",
        "sequence_sha256": sequence_sha256(protein),
    }])
    pvac = _write(
        tmp_path / "all_epitopes.tsv",
        "Index\tHLA Allele\tEpitope Seq\tSub-peptide Position\tBest IC50 Score\n"
        "SPORF_1\tHLA-A*02:01\tYYYYYYYYY\t2\t25\n",
    )
    entity_bundle = {
        "events": [{"splice_event_id": "SEV|1"}],
        "orfs": [{"orf_id": "ORF|1", "transcript_hypothesis_id": "STH|1", "splice_event_id": "SEV|1", "gene": "G", "protein_sequence": protein, "protein_sequence_sha256": sequence_sha256(protein)}],
        "transcripts": [{"transcript_hypothesis_id": "STH|1", "splice_event_id": "SEV|1", "junction_chain": "SJ|GRCh38|chr1|1|2|+"}],
        "peptide_origins": [],
    }
    bundle = parse_pvacbind(pvac, sample_id="S1", fasta_map=mapping, entity_bundle=entity_bundle)
    assert not bundle["presentation"]
    assert any(row["conflict_type"] == "PVACBIND_EPITOPE_ORF_MISMATCH" for row in bundle["conflicts"])


def test_unstranded_junction_is_not_exact_and_strict_mode_rejects_it(tmp_path: Path):
    junctions = _write(
        tmp_path / "unstranded.tsv",
        "chrom\tstart\tend\tname\tscore\nchr1\t150\t200\tJ1\t17\n",
    )
    bundle = parse_junction_source(
        junctions, sample_id="S1", coordinate_system="regtools_annotated", strict=False,
    )
    assert bundle["junctions"][0]["junction_resolution_status"] == "RESOLVED_UNSTRANDED"
    assert bundle["tool_evidence"][0]["verified_value"] == ""
    assert bundle["tool_evidence"][0]["evidence_type"] == "UNSTRANDED_SPLIT_READ_SUPPORT_UNVERIFIED"
    assert any(row["conflict_type"] == "JUNCTION_STRAND_UNRESOLVED" for row in bundle["conflicts"])
    with pytest.raises(JunctionNormalizationError, match="strand is required"):
        parse_junction_source(
            junctions, sample_id="S1", coordinate_system="regtools_annotated", strict=True,
        )


def test_pvacbind_rejects_missing_orf_and_sequence_hash_mismatch(tmp_path: Path):
    protein = "MARNDCEQGHILK"
    mapping = tmp_path / "map.tsv"
    write_tsv(mapping, [{
        "index": "SPORF_1", "orf_id": "ORF|missing", "transcript_hypothesis_id": "STH|1",
        "splice_event_id": "SEV|1", "sample_id": "S1", "gene": "G",
        "sequence_sha256": "0" * 64,
    }])
    pvac = _write(
        tmp_path / "all_epitopes.tsv",
        "Index\tHLA Allele\tEpitope Seq\tSub-peptide Position\n"
        "SPORF_1\tHLA-A*02:01\tARNDCEQGH\t2\n",
    )
    entity_bundle = {
        "events": [{"splice_event_id": "SEV|1"}],
        "transcripts": [{"transcript_hypothesis_id": "STH|1", "splice_event_id": "SEV|1"}],
        "orfs": [], "peptide_origins": [],
    }
    bundle = parse_pvacbind(pvac, sample_id="S1", fasta_map=mapping, entity_bundle=entity_bundle)
    assert not bundle["presentation"]
    assert any(row["conflict_type"] == "PVACBIND_PROVENANCE_CHAIN_INVALID" for row in bundle["conflicts"])

    entity_bundle["orfs"] = [{
        "orf_id": "ORF|missing", "transcript_hypothesis_id": "STH|1", "splice_event_id": "SEV|1",
        "protein_sequence": protein, "protein_sequence_sha256": sequence_sha256(protein),
    }]
    bundle = parse_pvacbind(pvac, sample_id="S1", fasta_map=mapping, entity_bundle=entity_bundle)
    assert not bundle["presentation"]
    assert "does not match ORF sequence" in bundle["conflicts"][0]["resolution_reason"]


def test_pvacbind_rejects_epitope_with_multiple_unresolved_positions(tmp_path: Path):
    protein = "AAAAAAAAAAA"
    mapping = tmp_path / "map.tsv"
    write_tsv(mapping, [{
        "index": "SPORF_1", "orf_id": "ORF|1", "transcript_hypothesis_id": "STH|1",
        "splice_event_id": "SEV|1", "sample_id": "S1", "gene": "G",
        "sequence_sha256": sequence_sha256(protein),
    }])
    pvac = _write(
        tmp_path / "all_epitopes.tsv",
        "Index\tHLA Allele\tEpitope Seq\tSub-peptide Position\n"
        "SPORF_1\tHLA-A*02:01\tAAAAAAAAA\t\n",
    )
    entity_bundle = {
        "events": [{"splice_event_id": "SEV|1"}],
        "transcripts": [{"transcript_hypothesis_id": "STH|1", "splice_event_id": "SEV|1"}],
        "orfs": [{
            "orf_id": "ORF|1", "transcript_hypothesis_id": "STH|1", "splice_event_id": "SEV|1",
            "protein_sequence": protein, "protein_sequence_sha256": sequence_sha256(protein),
        }],
        "peptide_origins": [],
    }
    bundle = parse_pvacbind(pvac, sample_id="S1", fasta_map=mapping, entity_bundle=entity_bundle)
    assert not bundle["presentation"]
    assert any(row["conflict_type"] == "PVACBIND_EPITOPE_POSITION_UNRESOLVED" for row in bundle["conflicts"])


@pytest.mark.parametrize(
    ("label", "expected"),
    [("cryptic exon", "CRYPTIC_EXON"), ("exitron", "EXITRON"), ("novel splice junction", "NOVEL_JUNCTION")],
)
def test_extended_splice_event_type_aliases(label: str, expected: str):
    assert infer_event_type("events.tsv", explicit=label) == expected


def test_immunopepper_marks_cryptic_exon_and_exitron(tmp_path: Path):
    meta = _write(
        tmp_path / "events.tsv",
        "peptide\tid\tgeneName\tgeneChr\tgeneStrand\tmutationMode\tmodifiedExonsCoord\n"
        "MARNDCEQGHILK\tTX1\tG1\tchr1\t+\tcryptic_exon\tchr1:100-150;chr1:201-250\n"
        "MARNDCEQGHILK\tTX2\tG2\tchr2\t+\texitron\tchr2:100-150;chr2:201-250\n",
    )
    events = parse_immunopepper_meta(meta, sample_id="S1")["events"]
    assert [row["event_type"] for row in events] == ["CRYPTIC_EXON", "EXITRON"]
    assert events[0]["cryptic_exon_status"] == "PRESENT"


def test_strict_layer_requires_versions_for_executed_tools(tmp_path: Path):
    junctions = _write(
        tmp_path / "junctions.tsv",
        "chrom\tstart\tend\tname\tscore\tstrand\nchr1\t150\t200\tJ1\t17\t+\n",
    )
    with pytest.raises(ValueError, match="RegTools"):
        build_splice_provenance_layer(
            sample_id="S1", outdir=tmp_path / "blocked", junctions=junctions,
            junction_coordinate_system="regtools_annotated", strict=True,
        )
    outputs = build_splice_provenance_layer(
        sample_id="S1", outdir=tmp_path / "ready", junctions=junctions,
        junction_coordinate_system="regtools_annotated",
        tool_versions={"RegTools": "1.0-fixture"}, strict=True,
    )
    assert outputs["junctions"].is_file()


def test_normal_background_emits_formal_seven_state_machine(tmp_path: Path):
    table = tmp_path / "normal.tsv"
    cases = [
        ("DETECTED", "ADEQUATE", "MATCHED_NORMAL", "false", "DETECTED_MATCHED_NORMAL"),
        ("DETECTED", "ADEQUATE", "GTEX_PANEL", "false", "DETECTED_BROAD_NORMAL"),
        ("DETECTED", "ADEQUATE", "TISSUE", "true", "DETECTED_CRITICAL_TISSUE"),
        ("LOW_LEVEL", "ADEQUATE", "TISSUE", "false", "LOW_LEVEL_NONCRITICAL_NORMAL"),
        ("NOT_DETECTED", "ADEQUATE", "PANEL", "false", "NOT_DETECTED_ADEQUATE_COVERAGE"),
        ("NOT_DETECTED", "LOW_COVERAGE", "PANEL", "false", "NOT_DETECTED_LOW_COVERAGE"),
        ("NOT_DETECTED", "UNASSESSED", "PANEL", "false", "UNASSESSED"),
    ]
    write_tsv(table, [
        {"junction_id": f"SJ|GRCh38|chr1|{i}|{i + 1}|+", "detection_status": detection,
         "coverage_status": coverage, "normal_source_type": source, "critical_tissue": critical}
        for i, (detection, coverage, source, critical, _) in enumerate(cases, start=10)
    ])
    observed = parse_normal_coverage(table, sample_id="S1")["normal_background"]
    assert [row["assessment_status"] for row in observed] == [case[-1] for case in cases]


def test_high_order_evidence_requires_exact_entity_and_upgrades_e3_o3(tmp_path: Path):
    entities = {
        "junctions": [{"junction_id": "SJ|GRCh38|chr1|151|200|+"}],
        "events": [{"splice_event_id": "SEV|1"}],
        "transcripts": [{"transcript_hypothesis_id": "STH|1", "splice_event_id": "SEV|1"}],
        "orfs": [{"orf_id": "ORF|1", "splice_event_id": "SEV|1", "transcript_hypothesis_id": "STH|1",
                  "protein_sequence_sha256": "abc", "frame_status": "IN_FRAME", "orf_validity_status": "VALID",
                  "source_generator": "Generator1"}],
        "peptide_origins": [{"origin_peptide_id": "POR|1", "peptide_id": "PEP|1", "orf_id": "ORF|1",
                             "splice_event_id": "SEV|1", "junction_ids": "SJ|GRCh38|chr1|151|200|+",
                             "crosses_junction": "true", "contains_novel_aa": "true"}],
        "event_junction_links": [{"splice_event_id": "SEV|1", "junction_id": "SJ|GRCh38|chr1|151|200|+"}],
        "presentation": [{"origin_peptide_id": "POR|1"}], "normal_background": [],
        "tool_evidence": [{"entity_id": "SJ|GRCh38|chr1|151|200|+", "evidence_group": "RNA_JUNCTION",
                           "source_assay_id": "RNA_BAM_SHA256_1", "verified_value": "12", "resolution_status": "RESOLVED_EXACT"}],
    }
    high = _write(
        tmp_path / "high.tsv",
        "entity_type\tentity_id\tevidence_group\tevidence_status\tsource_tool\tsource_tool_version\n"
        "SPLICE_EVENT\tSEV|1\tLONG_READ\tCONFIRMED\tIsoSeq\t1.0\n"
        "ORF\tORF|1\tPROTEIN_VALIDATION\tVALIDATED\tProteomics\t2.0\n"
        "ORF\tORF|missing\tLIGANDOME\tDETECTED\tMS\t1.0\n",
    )
    parsed = parse_high_order_evidence(high, sample_id="S1", entity_bundle=entities)
    assert len(parsed["tool_evidence"]) == 2
    assert len(parsed["conflicts"]) == 1
    entities["tool_evidence"].extend(parsed["tool_evidence"])
    consensus = build_consensus(entities, sample_id="S1")[0]
    assert consensus["event_evidence_grade"] == "E3"
    assert consensus["orf_evidence_grade"] == "O3"


def test_rna_assay_identity_collapses_same_source_and_separates_independent_sources():
    base = {
        "events": [{"splice_event_id": "SEV|1"}],
        "event_junction_links": [{"splice_event_id": "SEV|1", "junction_id": "SJ|1"}],
        "orfs": [{"orf_id": "ORF|1", "splice_event_id": "SEV|1", "protein_sequence_sha256": "abc",
                  "frame_status": "IN_FRAME", "orf_validity_status": "VALID", "source_generator": "G1"}],
        "peptide_origins": [{"origin_peptide_id": "POR|1", "peptide_id": "PEP|1", "orf_id": "ORF|1",
                             "splice_event_id": "SEV|1", "crosses_junction": "true", "contains_novel_aa": "true"}],
        "presentation": [], "normal_background": [],
        "tool_evidence": [
            {"entity_id": "SJ|1", "evidence_group": "RNA_JUNCTION", "source_assay_id": "BAM_SHA_1",
             "verified_value": "5", "resolution_status": "RESOLVED_EXACT"},
            {"entity_id": "SJ|1", "evidence_group": "RNA_JUNCTION", "source_assay_id": "BAM_SHA_1",
             "verified_value": "5", "resolution_status": "RESOLVED_EXACT"},
        ],
    }
    same = build_consensus(base, sample_id="S1")[0]
    assert same["independent_rna_sources"] == "BAM_SHA_1"
    base["tool_evidence"].append({"entity_id": "SJ|1", "evidence_group": "RNA_JUNCTION", "source_assay_id": "BAM_SHA_2",
                                  "verified_value": "4", "resolution_status": "RESOLVED_EXACT"})
    separate = build_consensus(base, sample_id="S1")[0]
    assert separate["independent_rna_sources"] == "BAM_SHA_1;BAM_SHA_2"


def test_consensus_reason_codes_are_materialized_as_conflicts():
    rows = [{"consensus_id": "CON|1", "origin_peptide_id": "POR|1", "final_evidence_tier": "R4",
             "hard_fail_codes": "HARD_ORF_INVALID", "cap_codes": "CAP_PRESENTATION_UNASSESSED_R3"}]
    conflicts = consensus_reason_conflicts(rows, sample_id="S1")
    assert {row["conflict_type"] for row in conflicts} == {"CONSENSUS_HARD_FAIL", "CONSENSUS_PRIORITY_CAP"}
    assert {row["observed_values"] for row in conflicts} == {"HARD_ORF_INVALID", "CAP_PRESENTATION_UNASSESSED_R3"}
