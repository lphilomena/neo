from neoag.evidence_layer import build_standard_evidence_layer
from neoag.model_layers import compute_l3_dimension_scores
from neoag.utils import read_tsv, write_tsv


def test_wts_evidence_is_hydrated_and_changes_l3_scores(tmp_path):
    parsed = tmp_path / "parsed"
    parsed.mkdir()
    events = parsed / "raw_events.tsv"
    peptides = parsed / "raw_peptides.tsv"
    write_tsv(events, [
        {
            "event_id": "SNV1", "sample_id": "S1", "event_type": "SNV",
            "mutation_source": "SNV", "peptide_consequence": "missense",
            "gene": "GENE1", "transcript_id": "ENST1.4", "chrom": "chr1",
            "pos": "100", "ref": "A", "alt": "T",
        },
        {
            "event_id": "FUS1", "sample_id": "S1", "event_type": "Fusion",
            "mutation_source": "SV", "peptide_consequence": "fusion",
            "gene": "GENE2::GENE3", "rna_junction_reads": "8",
        },
    ], [
        "event_id", "sample_id", "event_type", "mutation_source",
        "peptide_consequence", "gene", "transcript_id", "chrom", "pos",
        "ref", "alt", "rna_junction_reads",
    ])
    write_tsv(peptides, [
        {"peptide_id": "P1", "event_id": "SNV1", "sample_id": "S1", "peptide": "AAAAAAAAA"},
        {"peptide_id": "P2", "event_id": "FUS1", "sample_id": "S1", "peptide": "CCCCCCCCC"},
    ], ["peptide_id", "event_id", "sample_id", "peptide"])

    gene_tpm = tmp_path / "gene_tpm.tsv"
    write_tsv(gene_tpm, [
        {"gene": "GENE1", "TPM": "10"},
        {"gene": "GENE2::GENE3", "TPM": "4"},
    ], ["gene", "TPM"])
    tx_tpm = tmp_path / "transcript_tpm.tsv"
    write_tsv(tx_tpm, [{"transcript_id": "ENST1", "TPM": "5"}], ["transcript_id", "TPM"])
    rna_vaf = tmp_path / "rna_vaf.tsv"
    write_tsv(rna_vaf, [{
        "chrom": "chr1", "pos": "100", "ref": "A", "alt": "T",
        "rna_ref_reads": "16", "rna_alt_reads": "4", "rna_depth": "20", "rna_vaf": "0.2",
    }], ["chrom", "pos", "ref", "alt", "rna_ref_reads", "rna_alt_reads", "rna_depth", "rna_vaf"])

    profile = {
        "_profile_name": "test", "l3_weights": {"expression": 1, "rna_junction_support": 1},
        "safety": {},
    }
    build_standard_evidence_layer(
        tmp_path, profile, raw_events=events, raw_peptides=peptides,
        expression=gene_tpm, transcript_expression=tx_tpm, rna_vaf=rna_vaf,
        sample_id="S1",
    )

    by_id = {row["event_id"]: row for row in read_tsv(events)}
    snv = by_id["SNV1"]
    fusion = by_id["FUS1"]
    assert snv["gene_expression_tpm"] == "10.0000"
    assert snv["transcript_expression_tpm"] == "5.0000"
    assert snv["expression_evidence_status"] == "GENE_AND_TRANSCRIPT_SUPPORTED"
    assert snv["rna_alt_reads"] == "4"
    assert snv["rna_support_status"] == "RNA_ALT_SUPPORTED"
    assert float(snv["rna_evidence_score"]) > 0
    assert fusion["rna_support_status"] == "RNA_JUNCTION_SUPPORTED"
    assert float(fusion["rna_evidence_score"]) > 0
    peptide_by_id = {row["peptide_id"]: row for row in read_tsv(peptides)}
    assert peptide_by_id["P1"]["rna_support_status"] == "RNA_ALT_SUPPORTED"
    assert float(peptide_by_id["P1"]["rna_evidence_score"]) > 0

    scores = compute_l3_dimension_scores({}, snv, {}, profile, high_expression_tpm=20)
    assert float(scores["l3_expression_score"]) > 0
    assert scores["l3_rna_support_score"] == scores["l3_rna_junction_support_score"]
    assert float(scores["l3_rna_junction_support_score"]) > 0


def test_coordinate_rna_count_does_not_leak_to_another_variant_in_same_gene(tmp_path):
    from neoag.adapters.rna_vaf import choose_rna_vaf_support, load_rna_vaf_support

    counts = tmp_path / "counts.tsv"
    write_tsv(counts, [{
        "gene": "GENE1", "chrom": "chr1", "pos": "100", "ref": "A", "alt": "T",
        "rna_alt_reads": "9", "rna_depth": "20",
    }], ["gene", "chrom", "pos", "ref", "alt", "rna_alt_reads", "rna_depth"])
    index = load_rna_vaf_support(counts)
    other = choose_rna_vaf_support(
        {"gene": "GENE1", "chrom": "chr1", "pos": "200", "ref": "G", "alt": "C"},
        index,
        "raw_events",
    )
    assert other.rna_alt_reads == ""
