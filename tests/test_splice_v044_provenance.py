"""Regression tests for the v0.4.4 exact-junction/provenance repair."""

from __future__ import annotations

from pathlib import Path

from neoag.adapters.pvactools_parser import parse_pvactools_outputs
from neoag.adapters.splice_junction_adapter import (
    build_junction_support_index,
    merge_splice_into_catalog,
)
from neoag.open_neo.tool_consensus import _write_splice
from neoag.open_neo.rna_fusion_splice_profile import generate_rna_fusion_splice_manifest
from neoag.production_runner import run_production
from neoag.splice.coordinates import iter_junction_records, read_source_rows
from neoag.splice.normalization import normalize_splice_sources
from neoag.splice.registry import resolve_junction_support
from neoag.utils import read_tsv, write_tsv

ROOT = Path(__file__).resolve().parents[1]
REGTOOLS = ROOT / "data/fixtures/regtools_splice_junctions.tsv"
PVACSPLICE = ROOT / "data/fixtures/pvacsplice_aggregated.tsv"


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_regtools_annotated_coordinates_are_canonical_1based_closed():
    records = list(
        iter_junction_records(
            REGTOOLS,
            sample_id="S1",
            source_tool="RegTools",
            genome_build="GRCh38",
        )
    )
    by_alias = {record.source_junction_id: record for record in records}
    plus = by_alias["JUNC00000003"].junction
    minus = by_alias["JUNC00000017"].junction
    assert plus is not None
    assert plus.junction_id == "SJ|GRCh38|chr1|93221744|93221991|+"
    assert plus.donor_1based == 93221744
    assert plus.acceptor_1based == 93221991
    assert minus is not None
    assert minus.junction_id == "SJ|GRCh38|chr1|154589929|154590184|-"
    assert minus.donor_1based == 154590184
    assert minus.acceptor_1based == 154589929


def test_regtools_bed12_uses_blocks_to_derive_intron(tmp_path: Path):
    bed = _write(
        tmp_path / "junctions.bed",
        "chr2\t100\t250\tJX\t12\t+\t100\t250\t0\t2\t50,40,\t0,110,\n",
    )
    record = next(
        iter_junction_records(
            bed,
            sample_id="S1",
            source_tool="RegTools",
            genome_build="GRCh38",
        )
    )
    assert record.junction is not None
    # left exon is [100,150), right exon begins at 210; intron is 151..210.
    assert record.junction.junction_id == "SJ|GRCh38|chr2|151|210|+"


def test_snaf_separate_boundary_columns_match_regtools_bed12(tmp_path: Path):
    regtools = _write(
        tmp_path / "junctions.bed",
        "chr11\t124872999\t124873391\tJUNC1\t54\t?\t124872999\t124873391\t0\t2\t90,82\t0,310\n",
    )
    snaf = _write(
        tmp_path / "snaf.tsv",
        "sample_id\tevent_id\tgene\tchrom\tstart\tend\tjunction_reads\tpeptide\thla_allele\tbinding_rank\n"
        "S1\tENSG00000154134:E9.1-E10.1\tROBO3\tchr11\t124873089\t124873310\t54\tIANVQEMDM\tHLA-C*02:02\t1.19\n",
    )
    outputs = normalize_splice_sources(
        sample_id="S1",
        junctions=regtools,
        snaf=snaf,
        outdir=tmp_path / "normalized",
    )
    consensus = read_tsv(outputs["splice_consensus"])
    assert len(consensus) == 1
    assert consensus[0]["event_id"] == "SJ|GRCh38|chr11|124873090|124873309|."
    assert consensus[0]["status"] == "CROSS_DOMAIN_CONFIRMED_EXACT_JUNCTION"


def test_same_gene_read_leakage_is_forbidden(tmp_path: Path):
    primary = _write(
        tmp_path / "regtools.tsv",
        "chrom\tstart\tend\tname\tscore\tstrand\tgene_names\n"
        "chr1\t99\t200\tHIGH\t5000\t+\tGENE1\n"
        "chr1\t299\t400\tLOW\t2\t+\tGENE1\n",
    )
    index = build_junction_support_index(primary, sample_id="S1")

    low = resolve_junction_support(
        {"source_junction_id": "LOW", "gene": "GENE1", "provided_rna_junction_reads": "2"},
        index,
    )
    assert low.selected_reads == 2
    assert low.support_status == "SUPPORTED_EXACT_JUNCTION"

    unknown = resolve_junction_support(
        {
            "source_junction_id": "UNKNOWN",
            "gene": "GENE1",
            "provided_rna_junction_reads": "5000",
        },
        index,
    )
    assert unknown.selected_reads == 0
    assert unknown.provided_reads == 5000
    assert unknown.support_status == "PROVIDED_UNVERIFIED"


def test_ambiguous_source_alias_never_transfers_reads(tmp_path: Path):
    primary = _write(
        tmp_path / "ambiguous.tsv",
        "chrom\tstart\tend\tname\tscore\tstrand\tgene_names\n"
        "chr1\t99\t200\tDUP\t100\t+\tGENE1\n"
        "chr1\t299\t400\tDUP\t900\t+\tGENE1\n",
    )
    index = build_junction_support_index(primary, sample_id="S1")
    match = resolve_junction_support(
        {"source_junction_id": "DUP", "provided_rna_junction_reads": "900"},
        index,
    )
    assert match.selected_reads == 0
    assert match.match_status == "AMBIGUOUS"
    assert match.conflict == "AMBIGUOUS_LINK"


def test_normalizer_keeps_exact_schema_and_multitool_provenance(tmp_path: Path):
    primary = _write(
        tmp_path / "regtools.tsv",
        "chrom\tstart\tend\tname\tscore\tstrand\tgene_names\n"
        "chr1\t99\t200\tJ1\t7\t+\tGENE1\n"
        "chr1\t299\t400\tJ2\t5000\t+\tGENE1\n",
    )
    header = (
        "source_junction_id\tgene\tpeptide\thla_allele\trna_junction_reads\t"
        "crosses_junction\n"
    )
    snaf = _write(tmp_path / "snaf.tsv", header + "J1\tGENE1\tPEPTIDEAA\tHLA-A*02:01\t7\tyes\n")
    splicemutr = _write(
        tmp_path / "splicemutr.tsv",
        header
        + "J1\tGENE1\tPEPTIDEAA\tHLA-A*02:01\t7\tyes\n"
        + "UNKNOWN\tGENE1\tUNRESOLV\tHLA-A*02:01\t5000\tyes\n",
    )
    outputs = normalize_splice_sources(
        sample_id="S1",
        junctions=primary,
        snaf=snaf,
        splicemutr=splicemutr,
        outdir=tmp_path / "out",
    )

    events = read_tsv(outputs["raw_events"])
    peptides = read_tsv(outputs["raw_peptides"])
    resolved_event = next(row for row in events if row["canonical_junction_id"].endswith("|100|200|+"))
    assert resolved_event["junction_start"] == "100"
    assert resolved_event["junction_end"] == "200"
    assert resolved_event["junction_strand"] == "+"
    assert resolved_event["rna_junction_reads"] == "7"
    assert resolved_event["source_tools"] == "RegTools;SNAF;SpliceMutr"
    assert resolved_event["provenance_record_count"] == "3"

    merged = next(row for row in peptides if row["peptide"] == "PEPTIDEAA")
    assert merged["source_tool"] in {"SNAF", "SpliceMutr"}
    assert merged["source_tools"] == "SNAF;SpliceMutr"
    assert merged["provenance_record_count"] == "2"
    assert merged["canonical_junction_id"].endswith("|100|200|+")
    assert merged["rna_junction_reads"] == "7"

    unresolved = next(row for row in peptides if row["peptide"] == "UNRESOLV")
    assert unresolved["canonical_junction_id"] == ""
    assert unresolved["rna_junction_reads"] == "0"
    assert unresolved["provided_rna_junction_reads"] == "5000"
    assert unresolved["junction_support_status"] == "PROVIDED_UNVERIFIED"

    provenance = read_tsv(outputs["peptide_merge_provenance"])
    merged_key_rows = [row for row in provenance if "PEPTIDEAA" in row["merge_key"]]
    assert {row["source_tool"] for row in merged_key_rows} == {"SNAF", "SpliceMutr"}
    assert read_tsv(outputs["splice_qc"])[-1] == {
        "metric": "gene_or_nearest_locus_fallbacks_used",
        "value": "0",
    }
    assert Path(outputs["junction_aliases"]).is_file()
    consensus_provenance = read_tsv(outputs["splice_consensus_provenance"])
    assert {row["tool"] for row in consensus_provenance} == {
        "RegTools",
        "SNAF",
        "SpliceMutr",
    }
    assert Path(outputs["splice_consensus_conflicts"]).is_file()
    assert Path(outputs["evidence_conflicts"]).is_file()


def test_candidate_only_emits_linked_primary_records_not_unrelated_junctions(tmp_path: Path):
    primary = _write(
        tmp_path / "regtools.tsv",
        "chrom\tstart\tend\tname\tscore\tstrand\tgene_names\n"
        "chr1\t99\t200\tJ1\t7\t+\tGENE1\n"
        "chr1\t299\t400\tJ2\t50\t+\tGENE2\n",
    )
    snaf = _write(
        tmp_path / "snaf.tsv",
        "source_junction_id\tgene\tpeptide\thla_allele\tbinding_rank\n"
        "J1\tGENE1\tACDEFGHIK\tHLA-A*02:01\t0.8\n",
    )
    outputs = normalize_splice_sources(
        sample_id="S1",
        junctions=primary,
        snaf=snaf,
        outdir=tmp_path / "candidate_only",
        candidate_only=True,
    )
    events = read_tsv(outputs["raw_events"])
    evidence = read_tsv(outputs["splice_tool_evidence"])
    assert {row["source_junction_id"] for row in events} == {"J1"}
    assert {row["source_junction_id"] for row in evidence} == {"J1"}
    qc = {row["metric"]: row["value"] for row in read_tsv(outputs["splice_qc"])}
    assert qc["primary_junction_records"] == "2"
    assert qc["candidate_only_output"] == "true"


def test_tool_consensus_does_not_confirm_same_gene_different_junction(tmp_path: Path):
    rna = _write(
        tmp_path / "rna.tsv",
        "chrom\tstart\tend\tname\tscore\tstrand\tgene_names\n"
        "chr1\t99\t200\tJ1\t10\t+\tGENE1\n"
        "chr1\t299\t400\tJ2\t20\t+\tGENE1\n",
    )
    neo = _write(
        tmp_path / "neo.tsv",
        "source_junction_id\tgene\tpeptide\n"
        "J2\tGENE1\tPEPTIDEAA\n"
        "UNKNOWN\tGENE1\tUNRESOLV\n",
    )
    declared = {
        "splice_dna": {},
        "splice_rna": {"regtools": str(rna)},
        "splice_neoantigen": {"snaf": str(neo)},
    }
    overall, rows = _write_splice(declared, tmp_path / "consensus")
    assert overall == "CROSS_DOMAIN_CONFIRMED"
    by_event = {row["event_id"]: row for row in rows}
    j1 = next(row for key, row in by_event.items() if key.endswith("|100|200|+"))
    j2 = next(row for key, row in by_event.items() if key.endswith("|300|400|+"))
    unresolved = next(row for key, row in by_event.items() if key.startswith("UNRESOLVED_SJ|"))
    assert j1["status"] == "RNA_JUNCTION_SUPPORTED"
    assert j2["status"] == "CROSS_DOMAIN_CONFIRMED_EXACT_JUNCTION"
    assert unresolved["status"] == "NEOANTIGEN_TOOL_ONLY_UNRESOLVED"


def test_pvacsplice_fixture_gets_only_exact_regtools_support(tmp_path: Path):
    events, peptides = parse_pvactools_outputs([PVACSPLICE], "S1", "default")
    merged_events, merged_peptides = merge_splice_into_catalog(
        REGTOOLS,
        "S1",
        "default",
        events,
        peptides,
        tools_dir=tmp_path / "tools",
    )
    adar = next(row for row in merged_peptides if row["gene"] == "ADAR")
    assert adar["rna_junction_reads"] == "5919"
    assert adar["canonical_junction_id"] == "SJ|GRCh38|chr1|154589929|154590184|-"
    assert adar["junction_support_status"] == "SUPPORTED_EXACT_JUNCTION"
    assert any(row["gene"] == "ADAR" for row in merged_events)


def test_production_runner_merges_without_discarding_source_rows(tmp_path: Path):
    hla = _write(tmp_path / "hla.txt", "HLA-A*02:01\n")
    event_a = tmp_path / "event_a.tsv"
    event_b = tmp_path / "event_b.tsv"
    peptide_a = tmp_path / "peptide_a.tsv"
    peptide_b = tmp_path / "peptide_b.tsv"
    write_tsv(event_a, [{"event_id": "E1", "sample_id": "S1", "event_type": "Splice", "gene": "GENE1", "source_tools": "ToolA", "source_record_id": "A1"}])
    write_tsv(event_b, [{"event_id": "E1", "sample_id": "S1", "event_type": "Splice", "gene": "GENE1", "source_tools": "ToolB", "source_record_id": "B1"}])
    write_tsv(peptide_a, [{"peptide_id": "PA", "event_id": "E1", "sample_id": "S1", "gene": "GENE1", "peptide": "PEPTIDEAA", "hla_allele": "HLA-A*02:01", "source_tool": "ToolA", "source_record_id": "A2"}])
    write_tsv(peptide_b, [{"peptide_id": "PB", "event_id": "E1", "sample_id": "S1", "gene": "GENE1", "peptide": "PEPTIDEAA", "hla_allele": "HLA-A*02:01", "source_tool": "ToolB", "source_record_id": "B2"}])
    manifest = _write(
        tmp_path / "production.toml",
        f'''[run]
sample_id = "S1"
profile = "default"
hla_file = "{hla}"
expected_peptide_sources = ["ToolA", "ToolB"]

[stages.a]
required = true
source = "ToolA"
[stages.a.outputs]
raw_events = "{event_a}"
raw_peptides = "{peptide_a}"

[stages.b]
required = true
source = "ToolB"
[stages.b.outputs]
raw_events = "{event_b}"
raw_peptides = "{peptide_b}"
''',
    )
    result = run_production(
        manifest,
        outdir=tmp_path / "run",
        project_root=ROOT,
        execute=True,
        skip_ranking=True,
    )
    assert result.status == "PARTIAL"
    merged = read_tsv(tmp_path / "run/merged/raw_peptides.tsv")
    assert len(merged) == 1
    assert merged[0]["source_tools"] == "ToolA;ToolB"
    assert merged[0]["source_records"] == "A2;B2"
    provenance = read_tsv(tmp_path / "run/merged/peptide_provenance.tsv")
    assert len(provenance) == 2
    assert {row["stage_source"] for row in provenance} == {"ToolA", "ToolB"}
    assert result.provenance_outputs["peptide_provenance"].endswith("peptide_provenance.tsv")


def test_skill2_manifest_tracks_all_splice_audit_outputs(tmp_path: Path):
    fastq1 = _write(tmp_path / "tumor_R1.fastq.gz", "fixture\n")
    fastq2 = _write(tmp_path / "tumor_R2.fastq.gz", "fixture\n")
    manifest = tmp_path / "rna_fusion_splice.production.toml"
    generate_rna_fusion_splice_manifest(
        {
            "sample_id": "S1",
            "tumor_rna_fastq": [str(fastq1), str(fastq2)],
            "hla_alleles": ["HLA-A*02:01"],
        },
        manifest,
        project_root=ROOT,
        outdir=tmp_path / "run",
    )
    text = manifest.read_text(encoding="utf-8")
    for filename in (
        "junction_aliases.tsv",
        "evidence_conflicts.tsv",
        "splice_consensus_provenance.tsv",
        "splice_consensus_conflicts.tsv",
    ):
        assert filename in text
