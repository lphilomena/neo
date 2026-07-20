"""Tests for Mode C — RegTools / splice junction catalog ingestion and peptide e2e."""

from pathlib import Path

from neoag.adapters.event_catalog import parse_splice_catalog
from neoag.adapters.splice_junction_adapter import (
    build_junction_support_index,
    build_splice_peptides_from_vcf,
    merge_splice_into_catalog,
    run_splice_junction_upstream,
)
from neoag.evidence_layer import build_standard_evidence_layer
from neoag.input_router import build_raw_intermediates, resolve_entry_mode
from neoag.pipeline import run
from neoag.utils import read_tsv

ROOT = Path(__file__).resolve().parents[1]
REGTOOLS_FIXTURE = ROOT / "data/fixtures/regtools_splice_junctions.tsv"
GENERIC_FIXTURE = ROOT / "data/fixtures/splice_junctions_generic.tsv"
PVACSPLICE_FIXTURE = ROOT / "data/fixtures/pvacsplice_aggregated.tsv"
HCC1395_FIXTURE = ROOT / "data/examples/HCC1395/HCC1395_inputs/HCC1395.splice_junctions.tsv"
HCC1395_VCF = ROOT / "data/examples/HCC1395/HCC1395_inputs/annotated.expression.vcf.gz"
HCC1395_HLA = [
    "HLA-A*29:02",
    "HLA-B*45:01",
    "HLA-B*82:02",
    "HLA-C*06:02",
]


def test_resolve_entry_mode_splice_aliases():
    assert resolve_entry_mode({"inputs": {"entry_mode": "c"}}) == "splice_junction"
    assert resolve_entry_mode({"inputs": {"entry_mode": "splice_junction"}}) == "splice_junction"


def test_parse_splice_catalog_regtools_columns():
    events = parse_splice_catalog(REGTOOLS_FIXTURE, "S1", "default")
    by_gene = {e["gene"]: e for e in events}
    assert len(events) == 4
    assert by_gene["CCDC18"]["event_name"] == "JUNC00000003"
    assert by_gene["CCDC18"]["rna_junction_reads"] == "8"
    assert by_gene["CCDC18"]["chrom"] == "chr1"
    assert by_gene["CCDC18"]["pos"] == "93221743"
    assert by_gene["CCDC18"]["transcript_id"] == "ENST00000343253"
    assert by_gene["ADAR"]["rna_junction_reads"] == "5919"
    assert by_gene["chr1:93221463"]["rna_junction_reads"] == "2"
    assert all(e["peptide_consequence"] == "splice_junction" for e in events)
    assert all(e["event_type"] == "Splice" for e in events)
    assert all(e["source"].startswith("splice_catalog:") for e in events)


def test_parse_splice_catalog_generic_columns():
    events = parse_splice_catalog(GENERIC_FIXTURE, "S1", "default")
    by_gene = {e["gene"]: e for e in events}
    assert by_gene["TP53"]["rna_junction_reads"] == "42"
    assert by_gene["BRCA1"]["event_name"] == "J2"


def test_parse_hcc1395_splice_fixture_gene_names_and_reads():
    if not HCC1395_FIXTURE.is_file():
        return
    events = parse_splice_catalog(HCC1395_FIXTURE, "HCC1395", "default")
    genes = {e["gene"] for e in events}
    assert "ADAR" in genes
    assert "CCDC18" in genes
    adar = next(e for e in events if e["gene"] == "ADAR")
    assert int(adar["rna_junction_reads"]) >= 100
    named = [e for e in events if not e["gene"].startswith("chr")]
    assert len(named) >= 150


def test_build_intermediates_splice_junction_mode(tmp_path):
    cfg = {
        "sample": {"id": "SJ1", "profile": "default"},
        "inputs": {
            "entry_mode": "splice_junction",
            "splice_junction_tsv": str(REGTOOLS_FIXTURE),
        },
    }
    paths = build_raw_intermediates(cfg, tmp_path / "layer", root=ROOT)
    assert paths["entry_mode"] == "splice_junction"
    events = read_tsv(paths["raw_events"])
    peptides = read_tsv(paths["raw_peptides"])
    assert len(events) == 4
    assert peptides == []


def test_build_intermediates_splice_with_pvac_merge(tmp_path):
    cfg = {
        "sample": {"id": "SJ2", "profile": "default"},
        "inputs": {
            "entry_mode": "e2e",
            "pvac_files": [str(ROOT / "data/fixtures/pvacseq_aggregated.tsv")],
            "splice_junction_tsv": str(REGTOOLS_FIXTURE),
        },
    }
    paths = build_raw_intermediates(cfg, tmp_path / "layer", root=ROOT)
    events = read_tsv(paths["raw_events"])
    peptides = read_tsv(paths["raw_peptides"])
    assert len(peptides) > 0
    genes = {e["gene"] for e in events}
    assert "ADAR" in genes
    assert any(e["event_type"] == "SNV" for e in events)


def test_evidence_layer_rna_junction_from_splice_events(tmp_path):
    layer = tmp_path / "layer"
    build_raw_intermediates(
        {
            "sample": {"id": "SJ3", "profile": "default"},
            "inputs": {
                "entry_mode": "splice_junction",
                "regtools_tsv": str(REGTOOLS_FIXTURE),
            },
        },
        layer,
        root=ROOT,
    )
    evidence = build_standard_evidence_layer(layer, "default", sample_id="SJ3")
    rna_rows = read_tsv(evidence["rna_junction_evidence"])
    adar = [r for r in rna_rows if r.get("gene") == "ADAR"]
    assert adar and int(adar[0]["junction_reads"]) == 5919


def test_merge_splice_with_pvacsplice_fixture():
    from neoag.adapters.pvactools_parser import parse_pvactools_outputs

    _, pvac_peptides = parse_pvactools_outputs([PVACSPLICE_FIXTURE], "S1", "default")
    events, peptides = merge_splice_into_catalog(
        REGTOOLS_FIXTURE,
        "S1",
        "default",
        [],
        pvac_peptides,
    )
    assert len(events) >= 4
    assert len(peptides) == 3
    assert all(p["peptide_consequence"] == "splice_junction" for p in peptides)
    adar = next(p for p in peptides if p["gene"] == "ADAR")
    assert adar["crosses_junction"] == "yes"
    assert int(adar["rna_junction_reads"]) == 5919


def test_build_intermediates_splice_pvacsplice_mode(tmp_path):
    cfg = {
        "sample": {"id": "SJ4", "profile": "default"},
        "inputs": {
            "entry_mode": "splice_junction",
            "splice_junction_tsv": str(REGTOOLS_FIXTURE),
            "pvac_files": [str(PVACSPLICE_FIXTURE)],
        },
    }
    paths = build_raw_intermediates(cfg, tmp_path / "layer", root=ROOT)
    peptides = read_tsv(paths["raw_peptides"])
    assert len(peptides) == 3
    assert all(p["source_tool"] == "pVACsplice" for p in peptides)


def test_build_splice_peptides_from_hcc1395_vcf(tmp_path):
    if not (HCC1395_VCF.is_file() and HCC1395_FIXTURE.is_file()):
        return
    built = build_splice_peptides_from_vcf(
        HCC1395_VCF,
        HCC1395_FIXTURE,
        sample_id="HCC1395",
        profile_name="default",
        hla_alleles=HCC1395_HLA,
        cfg={"inputs": {"tumor_sample_name": "HCC1395_TUMOR_DNA"}},
        tools_dir=tmp_path / "tools",
        min_junction_reads=0,
    )
    assert built["splice_variant_rows"] > 0
    assert len(built["peptides"]) > 0
    assert all(p["peptide_consequence"] == "splice_junction" for p in built["peptides"])


def test_run_splice_junction_upstream_pvacsplice_stub(tmp_path):
    outputs = run_splice_junction_upstream(
        {"inputs": {"tumor_sample_name": "HCC1395_TUMOR_DNA"}},
        splice_path=REGTOOLS_FIXTURE,
        variants_vcf=HCC1395_VCF if HCC1395_VCF.is_file() else None,
        parsed_dir=tmp_path / "parsed",
        tools_dir=tmp_path / "tools",
        sample_id="SJ5",
        profile_name="default",
        hla_alleles=HCC1395_HLA,
        pvacsplice_tsv=PVACSPLICE_FIXTURE,
    )
    peptides = read_tsv(outputs["raw_peptides"])
    events = read_tsv(outputs["raw_events"])
    assert len(peptides) == 3
    assert len(events) >= 4
    assert outputs["peptide_source"] == "pvacsplice"


def test_e2e_run_splice_junction(tmp_path):
    layer = tmp_path / "layer"
    build_raw_intermediates(
        {
            "sample": {"id": "SJ6", "profile": "default"},
            "inputs": {
                "entry_mode": "splice_junction",
                "splice_junction_tsv": str(REGTOOLS_FIXTURE),
                "pvac_files": [str(PVACSPLICE_FIXTURE)],
            },
        },
        layer,
        root=ROOT,
    )
    evidence = build_standard_evidence_layer(
        layer,
        "default",
        normal_expression=ROOT / "resources/normal_expression.example.tsv",
        normal_hla_ligands=ROOT / "resources/normal_hla_ligands.example.tsv",
        sample_id="SJ6",
    )
    out = run(
        tmp_path / "score",
        "default",
        "SJ6",
        pvac_paths=[],
        raw_events=layer / "parsed/raw_events.tsv",
        raw_peptides=layer / "parsed/raw_peptides.tsv",
        netmhcpan=ROOT / "data/fixtures/netmhcpan_example.xls",
        mhcflurry=ROOT / "data/fixtures/mhcflurry_predictions.csv",
        vep_appm=ROOT / "data/fixtures/vep_appm.tsv",
        expression=ROOT / "data/fixtures/gene_expression.tsv",
        normal_expression=ROOT / "resources/normal_expression.example.tsv",
        normal_hla_ligands=ROOT / "resources/normal_hla_ligands.example.tsv",
        entry_mode="splice_junction",
    )
    ranked = read_tsv(out["ranked_peptides"])
    splice_rows = [r for r in ranked if r.get("peptide_consequence") == "splice_junction"]
    assert splice_rows
    assert Path(out["ranked_peptides"]).is_file()
    assert Path(evidence["rna_junction_evidence"]).is_file()



def test_evidence_layer_targeted_junction_validation_status(tmp_path):
    from neoag.evidence_layer import build_standard_evidence_layer
    from neoag.utils import read_tsv

    events = tmp_path / "raw_events.tsv"
    events.write_text(
        "event_id\tsample_id\tgene\tmutation_source\tpeptide_consequence\trna_junction_reads\n"
        "FUS1\tS1\tETV6::NTRK3\tfusion\tfusion\t0\n"
        "SJ1\tS1\tADAR\tsplice\tsplice_junction\t0\n",
        encoding="utf-8",
    )
    peptides = tmp_path / "raw_peptides.tsv"
    peptides.write_text("peptide_id\tevent_id\tsample_id\tgene\tpeptide\trna_junction_reads\n", encoding="utf-8")
    fusion = tmp_path / "fusion_evidence.tsv"
    fusion.write_text(
        "event_id\tfilter_status\trna_junction_reads\n"
        "FUS1\tpass\t12\n",
        encoding="utf-8",
    )
    junction = tmp_path / "junction.tsv"
    junction.write_text("event_id\tjunction_reads\nSJ1\t9\n", encoding="utf-8")

    outdir = tmp_path / "layer"
    paths = build_standard_evidence_layer(
        outdir,
        "default",
        raw_events=events,
        raw_peptides=peptides,
        rna_junction=junction,
        fusion_evidence=fusion,
        sample_id="S1",
    )
    by_event = {r["event_id"]: r for r in read_tsv(paths["rna_junction_evidence"])}
    assert by_event["FUS1"]["junction_reads"] == "12"
    assert by_event["FUS1"]["targeted_validation_status"] == "SUPPORTED"
    assert by_event["FUS1"]["targeted_validation_method"] == "fusion_targeted_rna"
    assert by_event["SJ1"]["junction_reads"] == "9"
    assert by_event["SJ1"]["targeted_validation_status"] == "SUPPORTED"
    assert by_event["SJ1"]["targeted_validation_method"] == "splice_junction_targeted_rna"
