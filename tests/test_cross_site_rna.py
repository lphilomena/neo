from neoag.cross_site_rna import build_cross_site_rna_evidence
from neoag.utils import read_tsv, write_tsv


def test_cross_site_requires_identity_for_exact_support(tmp_path):
    primary = tmp_path / "primary.tsv"
    secondary = tmp_path / "secondary.tsv"
    output = tmp_path / "cross.tsv"
    write_tsv(primary, [{"event_id": "E1", "event_type": "SNV", "chrom": "1", "pos": "10", "ref": "A", "alt": "T"}], ["event_id", "event_type", "chrom", "pos", "ref", "alt"])
    write_tsv(secondary, [{"event_id": "S1", "event_type": "SNV", "chrom": "chr1", "pos": "10", "ref": "A", "alt": "T", "rna_alt_reads": "4", "rna_depth": "20"}], ["event_id", "event_type", "chrom", "pos", "ref", "alt", "rna_alt_reads", "rna_depth"])
    build_cross_site_rna_evidence(primary, secondary, output, secondary_sample_id="ASCITES", identity_status="UNASSESSED")
    row = read_tsv(output)[0]
    assert row["cross_site_status"] == "EXACT_MATCH_IDENTITY_UNASSESSED"
    assert row["cross_site_exact_support"] == "no"
    assert row["cross_site_review_status"] == "REVIEW_REQUIRED"
    assert "指纹" in row["cross_site_review_reason"]

    build_cross_site_rna_evidence(primary, secondary, output, secondary_sample_id="ASCITES", identity_status="CONFIRMED")
    row = read_tsv(output)[0]
    assert row["cross_site_status"] == "EXACT_SHARED"
    assert row["cross_site_exact_support"] == "yes"


def test_fusion_gene_pair_is_not_exact_without_breakpoints(tmp_path):
    primary = tmp_path / "primary.tsv"
    secondary = tmp_path / "secondary.tsv"
    output = tmp_path / "cross.tsv"
    fields = ["event_id", "event_type", "gene", "rna_junction_reads"]
    write_tsv(primary, [{"event_id": "F1", "event_type": "Fusion", "gene": "EWSR1::WT1"}], fields)
    write_tsv(secondary, [{"event_id": "F2", "event_type": "Fusion", "gene": "EWSR1::WT1", "rna_junction_reads": "12"}], fields)
    build_cross_site_rna_evidence(primary, secondary, output, secondary_sample_id="ASCITES", identity_status="CONFIRMED")
    row = read_tsv(output)[0]
    assert row["cross_site_status"] == "GENE_PAIR_SHARED_BREAKPOINT_UNASSESSED"
    assert row["cross_site_exact_support"] == "no"
    assert "精确断点" in row["cross_site_review_reason"]


def test_same_variant_without_secondary_reads_is_low_power(tmp_path):
    primary = tmp_path / "primary.tsv"
    secondary = tmp_path / "secondary.tsv"
    output = tmp_path / "cross.tsv"
    fields = ["event_id", "event_type", "chrom", "pos", "ref", "alt"]
    write_tsv(primary, [{"event_id": "E1", "event_type": "SNV", "chrom": "chr1", "pos": "10", "ref": "A", "alt": "T"}], fields)
    write_tsv(secondary, [{"event_id": "S1", "event_type": "SNV", "chrom": "chr1", "pos": "10", "ref": "A", "alt": "T"}], fields)
    build_cross_site_rna_evidence(primary, secondary, output, secondary_sample_id="ASCITES", identity_status="UNASSESSED")
    assert read_tsv(output)[0]["cross_site_status"] == "SECONDARY_MATCH_LOW_POWER"
