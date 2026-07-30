from __future__ import annotations

from pathlib import Path

import pytest

from neoag.splice.adapters.easyquant import parse_easyquant
from neoag.splice.adapters.immunopepper import parse_immunopepper_meta
from neoag.splice.adapters.k4neo import parse_k4neo
from neoag.splice.adapters.mopepgen import parse_mopepgen
from neoag.splice.adapters.pvacsplice import parse_pvacsplice
from neoag.splice.adapters.regtools import parse_junction_source
from neoag.splice.adapters.splice2neo import parse_splice2neo
from neoag.splice.evidence_chains import build_evidence_chains
from neoag.splice.pipeline import build_splice_provenance_layer
from neoag.splice.sequence_queries import write_external_query_files
from neoag.utils import read_tsv, write_tsv


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def _immunopepper(tmp_path: Path) -> Path:
    return _write(
        tmp_path / "immunopepper.tsv",
        "peptide\tid\treadFrame\treadFrameAnnotated\tgeneName\tgeneChr\tgeneStrand\tmutationMode\thasStopCodon\tisJunctionList\tvariantSegExpr\tmodifiedExonsCoord\n"
        "MARNDCEQGHILKFPSTWYV\tTX1\t0\tTrue\tGENE1\tchr1\t+\tsplice\tFalse\tTrue\t12.5\tchr1:100-150;chr1:201-250\n",
    )


def _splice2neo(tmp_path: Path) -> Path:
    return _write(
        tmp_path / "splice2neo.tsv",
        "junc_id\tevent_type\tgene\tgene_id\tchrom\tpos\tref\talt\ttx_id\ttx_mod_id\tjunction_reads\tspliceai_score\tpangolin_score\tcts_seq\tcts_junc_pos\tcts_size\tcts_id\tpeptide_context\tneo_peptide\tjunc_in_orf\n"
        "chr1:151-200:+\tSE\tGENE1\tG1\tchr1\t1000\tA\tG\tTX1\tTX1.mod\t9\t0.82\t0.71\tAAAAAAAAAACCCCCCCCCCGGGGGGGGGG\t11\t30\tCTX1\tMARNDCEQGHILKFPSTWYV\tARNDCEQGH\ttrue\n",
    )


def _regtools(tmp_path: Path) -> Path:
    return _write(
        tmp_path / "regtools.tsv",
        "chrom\tstart\tend\tname\tscore\tstrand\tgene_names\n"
        "chr1\t150\t200\tJ1\t17\t+\tGENE1\n",
    )


def test_splice2neo_creates_exact_variant_junction_causal_link_and_queries(tmp_path: Path):
    reg = parse_junction_source(_regtools(tmp_path), sample_id="S1", source_tool="RegTools")
    bundle = parse_splice2neo(_splice2neo(tmp_path), sample_id="S1", entity_bundle=reg)
    assert bundle["variants"][0]["variant_id"] == "VAR|GRCh38|chr1|1000|A|G"
    assert bundle["causal_links"][0]["junction_id"] == "SJ|GRCh38|chr1|151|200|+"
    assert bundle["causal_links"][0]["causal_status"] == "DNA_RNA_CIS_SUPPORTED"
    assert bundle["sequence_queries"][0]["query_status"] == "READY"
    assert bundle["sequence_queries"][0]["position_1based"] == "11"
    assert bundle["orfs"][0]["orf_validity_status"] == "VALID_TRANSLATED_CONTEXT"
    assert bundle["peptide_origins"][0]["peptide_sequence"] == "ARNDCEQGH"


def test_easyquant_accepts_only_exact_project_query_id_and_upgrades_causal_link(tmp_path: Path):
    bundle = parse_splice2neo(_splice2neo(tmp_path), sample_id="S1")
    maps = write_external_query_files(tmp_path / "queries", bundle["sequence_queries"])
    qid = bundle["sequence_queries"][0]["query_id"]
    quant = _write(
        tmp_path / "quantification.tsv",
        "name\tpos\tjunc\tspan\tanch\ta\tb\n"
        f"{qid}\t11\t6\t2\t12\t20\t22\n"
        "UNKNOWN\t11\t999\t0\t0\t0\t0\n",
    )
    parsed = parse_easyquant(quant, sample_id="S1", query_map=maps["easyquant_query_map"])
    assert len(parsed["targeted_quantification"]) == 1
    assert parsed["targeted_quantification"][0]["support_status"] == "TARGETED_REQUANT_SUPPORTED"
    assert parsed["causal_links"][0]["causal_status"] == "TARGETED_REQUANT_SUPPORTED"
    assert any(row["conflict_type"] == "EASYQUANT_QUERY_ID_UNRESOLVED" for row in parsed["conflicts"])


def test_mopepgen_forms_peptide_level_dual_generator_without_false_full_orf_claim(tmp_path: Path):
    reg = parse_junction_source(_regtools(tmp_path), sample_id="S1", source_tool="RegTools")
    immuno = parse_immunopepper_meta(_immunopepper(tmp_path), sample_id="S1")
    event_id = immuno["events"][0]["splice_event_id"]
    mapping = tmp_path / "mopepgen_map.tsv"
    write_tsv(mapping, [{
        "mopepgen_header": "TX1|AS-1|1", "variant_id": "AS-1",
        "splice_event_id": event_id, "junction_id": "SJ|GRCh38|chr1|151|200|+",
        "gene": "GENE1", "transcript_id": "TX1", "crosses_junction": "true",
        "contains_novel_aa": "true",
    }])
    fasta = _write(tmp_path / "mopepgen.fasta", ">TX1|AS-1|1\nMARNDCEQGHILKFPSTWYV\n")
    entity = {key: [*reg.get(key, []), *immuno.get(key, [])] for key in set(reg) | set(immuno)}
    mo = parse_mopepgen(fasta, sample_id="S1", provenance_maps=[mapping], entity_bundle=entity)
    assert mo["peptide_origins"][0]["origin_status"] == "RESOLVED_EXACT_EVENT_AND_PEPTIDE"
    assert mo["orfs"][0]["orf_validity_status"] == "VALID_PEPTIDE_PRODUCT_ONLY"
    combined = {key: [*entity.get(key, []), *mo.get(key, [])] for key in set(entity) | set(mo)}
    chains = build_evidence_chains(combined, sample_id="S1")
    rna = [row for row in chains if row["chain_type"] == "RNA_DRIVEN"]
    assert any(row["chain_status"] == "DUAL_GENERATOR_EXACT_PEPTIDE" for row in rna)
    assert not any(row["chain_status"] == "DUAL_GENERATOR_EXACT_ORF" for row in rna)


def test_pvacsplice_requires_strand_aware_junction_and_exact_event(tmp_path: Path):
    reg = parse_junction_source(_regtools(tmp_path), sample_id="S1", source_tool="RegTools")
    s2n = parse_splice2neo(_splice2neo(tmp_path), sample_id="S1", entity_bundle=reg)
    entity = {key: [*reg.get(key, []), *s2n.get(key, [])] for key in set(reg) | set(s2n)}
    pvac = _write(
        tmp_path / "pvacsplice.tsv",
        "Chromosome\tStart\tStop\tReference\tVariant\tJunction\tJunction Score\tTranscript\tEnsembl Gene ID\tGene Name\tHLA Allele\tEpitope Seq\tBest IC50 Score\tBest Percentile\n"
        "chr1\t999\t1000\tA\tG\tchr1:151-200:+\t9\tTX1\tG1\tGENE1\tHLA-A*02:01\tARNDCEQGH\t22\t0.3\n",
    )
    parsed = parse_pvacsplice(pvac, sample_id="S1", entity_bundle=entity, strict=True)
    assert len(parsed["pvacsplice_predictions"]) == 1
    assert parsed["causal_links"][0]["causal_status"] == "PVACSPLICE_SUPPORTED"
    assert parsed["presentation"][0]["mapping_status"] == "MAPPED_EXACT_VARIANT_JUNCTION_EVENT"
    bad = _write(
        tmp_path / "pvacsplice_bad.tsv",
        "Chromosome\tStart\tStop\tReference\tVariant\tJunction Start\tJunction Stop\tGene Name\tHLA Allele\tEpitope Seq\n"
        "chr1\t999\t1000\tA\tG\t150\t200\tGENE1\tHLA-A*02:01\tARNDCEQGH\n",
    )
    rejected = parse_pvacsplice(bad, sample_id="S1", entity_bundle=entity)
    assert not rejected["pvacsplice_predictions"]
    assert any(row["conflict_type"] == "PVACSPLICE_EXACT_PROVENANCE_UNRESOLVED" for row in rejected["conflicts"])


def test_k4neo_maps_exact_cts_id_and_keeps_kmer_negative_distinct_from_locus_negative(tmp_path: Path):
    bundle = parse_splice2neo(_splice2neo(tmp_path), sample_id="S1")
    maps = write_external_query_files(tmp_path / "queries", bundle["sequence_queries"])
    qid = bundle["sequence_queries"][0]["query_id"]
    rates = _write(
        tmp_path / "healthy.tsv",
        "cts_id\tdevelopmental_stage\ttissue\tsample_rate\n"
        f"{qid}\tadult\tbrain\t0\n",
    )
    parsed = parse_k4neo(
        sample_id="S1", query_map=maps["k4neo_query_map"], healthy_sample_rate=[rates],
        critical_tissues=["brain"],
    )
    row = parsed["normal_background"][0]
    assert row["assessment_status"] == "NOT_DETECTED_KMER_SCREEN"
    assert row["coverage_status"] == "SEQUENCE_INDEX_QUERIED"
    detected = _write(
        tmp_path / "healthy_detected.tsv",
        "cts_id\tdevelopmental_stage\ttissue\tsample_rate\n"
        f"{qid}\tadult\tbrain\t0.25\n",
    )
    parsed2 = parse_k4neo(
        sample_id="S1", query_map=maps["k4neo_query_map"], healthy_sample_rate=[detected],
        critical_tissues=["brain"],
    )
    assert parsed2["normal_background"][0]["assessment_status"] == "NORMAL_DETECTED"
    assert parsed2["normal_background"][0]["critical_tissue"] == "true"


def test_end_to_end_build_forms_three_independent_chains(tmp_path: Path):
    reg = _regtools(tmp_path)
    immuno = _immunopepper(tmp_path)
    s2n_path = _splice2neo(tmp_path)
    immuno_bundle = parse_immunopepper_meta(immuno, sample_id="S1")
    event_id = immuno_bundle["events"][0]["splice_event_id"]
    mo_map = tmp_path / "mo_map.tsv"
    write_tsv(mo_map, [{
        "mopepgen_header": "TX1|AS-1|1", "variant_id": "AS-1", "splice_event_id": event_id,
        "junction_id": "SJ|GRCh38|chr1|151|200|+", "gene": "GENE1", "crosses_junction": "true",
        "contains_novel_aa": "true",
    }])
    mo_fasta = _write(tmp_path / "mo.fasta", ">TX1|AS-1|1\nMARNDCEQGHILKFPSTWYV\n")
    # Phase 1: generate the exact query registry from the same canonical event
    # graph that will be used for final import.  External EasyQuant/k4neo output
    # is only valid against this project-generated map.
    phase1 = build_splice_provenance_layer(
        sample_id="S1", outdir=tmp_path / "phase1", junctions=reg,
        immunopepper_meta=[immuno], splice2neo=[s2n_path],
        mopepgen_fasta=[mo_fasta], mopepgen_provenance_map=[mo_map],
        tool_versions={"RegTools": "fixture", "ImmunoPepper": "fixture"},
        strict=True,
    )
    query_files = {
        "easyquant_query_map": phase1["easyquant_query_map"],
        "k4neo_query_map": phase1["k4neo_query_map"],
    }
    phase1_queries = read_tsv(phase1["sequence_queries"])
    qid = phase1_queries[0]["query_id"]
    easy = _write(tmp_path / "easy.tsv", f"name\tpos\tjunc\tspan\tanch\ta\tb\n{qid}\t11\t5\t1\t10\t20\t20\n")
    k4 = _write(tmp_path / "k4.tsv", f"cts_id\tdevelopmental_stage\ttissue\tsample_rate\n{qid}\tadult\tblood\t0\n")
    normal = tmp_path / "normal.tsv"
    write_tsv(normal, [{
        "junction_id": "SJ|GRCh38|chr1|151|200|+", "detection_status": "NOT_DETECTED",
        "coverage_status": "ADEQUATE", "normal_source_type": "MATCHED_NORMAL", "normal_tissue": "blood",
    }])
    out = tmp_path / "layer"
    outputs = build_splice_provenance_layer(
        sample_id="S1", outdir=out, junctions=reg, immunopepper_meta=[immuno],
        splice2neo=[s2n_path], mopepgen_fasta=[mo_fasta], mopepgen_provenance_map=[mo_map],
        easyquant=[easy], easyquant_query_map=query_files["easyquant_query_map"],
        normal_coverage=[normal], k4neo_healthy_sample_rate=[k4],
        k4neo_query_map=query_files["k4neo_query_map"], k4neo_license_accepted=True,
        tool_versions={"RegTools": "fixture", "ImmunoPepper": "fixture"},
        strict=True,
    )
    chains = read_tsv(outputs["evidence_chains"])
    statuses = {(r["chain_type"], r["chain_status"]) for r in chains}
    assert ("RNA_DRIVEN", "DUAL_GENERATOR_EXACT_PEPTIDE") in statuses
    assert ("DNA_CAUSAL", "TARGETED_REQUANT_SUPPORTED") in statuses
    assert ("NORMAL_BACKGROUND", "LOCUS_AND_KMER_NEGATIVE") in statuses
    rna_chain = next(r for r in chains if r["chain_type"] == "RNA_DRIVEN" and r["chain_status"] == "DUAL_GENERATOR_EXACT_PEPTIDE")
    assert {"ImmunoPepper", "moPepGen"}.issubset(set(rna_chain["source_tools"].split(";")))
    consensus = read_tsv(outputs["consensus"])
    exact = next(r for r in consensus if r["translation_consensus_level"] == "EXACT_PEPTIDE")
    assert set(exact["independent_peptide_generators"].split(";")) == {"ImmunoPepper", "moPepGen"}
    assert exact["orf_consensus_status"] == "MULTI_GENERATOR_EXACT_PEPTIDE"
    qc = read_tsv(outputs["qc"])
    assert not [row for row in qc if row["status"] == "FAIL"]


def test_k4neo_build_requires_explicit_license_acceptance(tmp_path: Path):
    bundle = parse_splice2neo(_splice2neo(tmp_path), sample_id="S1")
    maps = write_external_query_files(tmp_path / "queries", bundle["sequence_queries"])
    qid = bundle["sequence_queries"][0]["query_id"]
    k4 = _write(tmp_path / "k4.tsv", f"cts_id\tdevelopmental_stage\ttissue\tsample_rate\n{qid}\tadult\tblood\t0\n")
    with pytest.raises(ValueError, match="license"):
        build_splice_provenance_layer(
            sample_id="S1", outdir=tmp_path / "layer", splice2neo=[_splice2neo(tmp_path)],
            k4neo_healthy_sample_rate=[k4], k4neo_query_map=maps["k4neo_query_map"],
        )
