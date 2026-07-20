from neoag.adapters.pvactools_parser import event_from_row
from neoag.adapters.variant_peptide_adapter import event_from_variant_row


def test_pvactools_event_from_row_parses_rna_metrics():
    row = {
        "Gene": "KRAS",
        "DNA VAF": "0.42",
        "Tumor DNA Depth": "110",
        "Tumor DNA Alt Count": "46",
        "RNA VAF": "0.18",
        "RNA Depth": "88",
        "RNA Alt Count": "16",
        "Transcript Expression": "15.0",
    }
    ev = event_from_row(row, "S1", "default", "pVACseq")
    assert ev["tumor_vaf"] == "0.42"
    assert ev["tumor_depth"] == "110"
    assert ev["tumor_alt_count"] == "46"
    assert ev["rna_vaf"] == "0.18"
    assert ev["rna_depth"] == "88"
    assert ev["rna_alt_reads"] == "16"


def test_variant_peptide_adapter_maps_rna_metrics():
    row = {
        "gene": "TP53",
        "chrom": "chr17",
        "pos": "100",
        "ref": "G",
        "alt": "A",
        "consequence": "missense_variant",
        "vaf": "0.35",
        "tumor_depth": "80",
        "tumor_alt_count": "28",
        "rna_vaf": "0.12",
        "rna_alt_reads": "9",
        "rna_depth": "75",
        "variant_key": "TP53|chr17:100G>A",
    }
    ev = event_from_variant_row(row, "S1", "default")
    assert ev["tumor_vaf"] == "0.35"
    assert ev["tumor_depth"] == "80"
    assert ev["tumor_alt_count"] == "28"
    assert ev["rna_vaf"] == "0.12"
    assert ev["rna_alt_reads"] == "9"
    assert ev["rna_depth"] == "75"



def test_parse_generic_rna_vaf_table_and_evidence_layer(tmp_path):
    from neoag.adapters.rna_vaf import parse_rna_vaf_table
    from neoag.evidence_layer import build_standard_evidence_layer
    from neoag.utils import read_tsv

    events = tmp_path / "raw_events.tsv"
    events.write_text(
        "event_id\tsample_id\tgene\tmutation_source\tpeptide_consequence\trna_junction_reads\n"
        "E1\tS1\tTP53\tSNV\tmissense\t0\n",
        encoding="utf-8",
    )
    peptides = tmp_path / "raw_peptides.tsv"
    peptides.write_text(
        "peptide_id\tevent_id\tsample_id\tgene\tpeptide\trna_junction_reads\n"
        "P1\tE1\tS1\tTP53\tAAAAAAAAA\t0\n",
        encoding="utf-8",
    )
    rna_vaf = tmp_path / "rna_vaf.tsv"
    rna_vaf.write_text(
        "event_id\tgene\trna_ref_reads\trna_alt_reads\trna_depth\trna_vaf\n"
        "E1\tTP53\t30\t10\t40\t0.25\n",
        encoding="utf-8",
    )
    parsed = parse_rna_vaf_table(rna_vaf)
    assert parsed[0].rna_alt_reads == "10"
    assert parsed[0].rna_vaf == "0.2500"

    outdir = tmp_path / "layer"
    paths = build_standard_evidence_layer(
        outdir,
        "default",
        raw_events=events,
        raw_peptides=peptides,
        rna_vaf=rna_vaf,
        sample_id="S1",
    )
    rows = read_tsv(paths["rna_junction_evidence"])
    evt = next(r for r in rows if r["event_id"] == "E1" and not r["peptide_id"])
    assert evt["rna_alt_reads"] == "10"
    assert evt["rna_ref_reads"] == "30"
    assert evt["rna_depth"] == "40"
    assert evt["rna_vaf"] == "0.2500"
    assert evt["rna_support_status"] == "RNA_ALT_SUPPORTED"
