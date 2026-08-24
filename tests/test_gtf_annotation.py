from pathlib import Path

from neoag.splice.coordinates import CanonicalJunction, JunctionSourceRecord
from neoag.splice.gtf_annotation import resolve_gtf_junction_strands


def _record() -> JunctionSourceRecord:
    return JunctionSourceRecord(
        row={}, junction=CanonicalJunction("GRCh38", "chr1", 151, 200, "."),
        sample_id="S1", source_tool="RegTools", source_tool_version="1",
        source_file="junctions.bed", source_row_number=1, source_record_id="r1",
        source_junction_id="J1", source_coordinate_system="bed12",
        source_chrom="chr1", source_start="100", source_end="250", gene="", gene_id="",
        transcript_ids="", unique_split_reads=7, multi_split_reads=0,
        total_split_reads=7, splice_motif="GT-AG", known_donor="", known_acceptor="",
        known_junction="", resolution_status="UNSTRANDED", resolution_method="BED12",
        coordinate_warning="strand unavailable", record_sha256="abc",
    )


def test_resolves_unique_same_transcript_exon_boundaries(tmp_path: Path):
    gtf = tmp_path / "annotation.gtf"
    gtf.write_text(
        'chr1\tt\texon\t100\t150\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "GENE1";\n'
        'chr1\tt\texon\t201\t250\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "GENE1";\n'
    )
    record = resolve_gtf_junction_strands([_record()], gtf)[0]
    assert record.junction.strand == "+"
    assert record.junction.junction_id == "SJ|GRCh38|chr1|151|200|+"
    assert record.transcript_ids == "T1"
    assert record.gene == "GENE1"
    assert record.resolution_method == "GTF_EXACT_TRANSCRIPT_EXON_BOUNDARIES"


def test_keeps_ambiguous_opposite_strands_unresolved(tmp_path: Path):
    gtf = tmp_path / "annotation.gtf"
    gtf.write_text(
        'chr1\tt\texon\t100\t150\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
        'chr1\tt\texon\t201\t250\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
        'chr1\tt\texon\t100\t150\t.\t-\t.\tgene_id "G2"; transcript_id "T2";\n'
        'chr1\tt\texon\t201\t250\t.\t-\t.\tgene_id "G2"; transcript_id "T2";\n'
    )
    record = resolve_gtf_junction_strands([_record()], gtf)[0]
    assert record.junction.strand == "."
    assert "unambiguous" in record.coordinate_warning
