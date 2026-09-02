from __future__ import annotations

import csv
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "build_fusion_caller_union", ROOT / "scripts/build_fusion_caller_union.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_event_junction_matching_is_breakpoint_specific_and_deduplicated(tmp_path: Path):
    junctions = tmp_path / "Chimeric.out.junction"
    junctions.write_text(
        "chr22\t101\t+\tchr11\t201\t-\t1\t0\t0\treadA\n"
        "chr22\t101\t+\tchr11\t201\t-\t1\t0\t0\treadA\n"
        "chr22\t101\t+\tchr11\t501\t-\t1\t0\t0\treadB\n",
        encoding="utf-8",
    )
    audit = [
        {"event_id": "E1", "left_breakpoint": "chr22:100", "right_breakpoint": "chr11:200"},
        {"event_id": "E2", "left_breakpoint": "11:500", "right_breakpoint": "22:100"},
    ]

    measurements, sidecar = MODULE.verify_event_junction_reads(
        audit, [junctions], None, tolerance=2
    )

    assert measurements["E1"]["verified_count"] == 1
    assert measurements["E2"]["verified_count"] == 1
    assert measurements["E1"]["status"] == "STAR_JUNCTION_VERIFIED"
    assert {row["event_id"] for row in sidecar} == {"E1", "E2"}


def test_verified_count_replaces_caller_count_without_losing_it():
    rows = [{"event_id": "E1", "rna_junction_reads": "7"}]
    MODULE.apply_junction_measurements(
        rows,
        {
            "E1": {
                "verified_count": 1,
                "status": "BAM_VERIFIED",
                "method": "star_chimeric_breakpoint_plus_bam_qname",
                "source": "junctions;rna.bam",
            }
        },
    )

    assert rows[0]["provided_rna_junction_reads"] == "7"
    assert rows[0]["rna_junction_reads"] == "1"
    assert rows[0]["junction_match_status"] == "BAM_VERIFIED"


def test_no_exact_match_clears_verified_count_but_preserves_caller_total():
    rows = [{"event_id": "E1", "rna_junction_reads": "1062"}]
    MODULE.apply_junction_measurements(
        rows,
        {
            "E1": {
                "verified_count": 0,
                "status": "NO_EXACT_JUNCTION_MATCH",
                "method": "star_chimeric_breakpoint",
                "source": "",
            }
        },
    )

    assert rows[0]["provided_rna_junction_reads"] == "1062"
    assert rows[0]["verified_rna_junction_reads"] == "0"
    assert rows[0]["rna_junction_reads"] == "0"
    assert rows[0]["junction_match_status"] == "NO_EXACT_JUNCTION_MATCH"


def test_targeted_support_count_is_exact_not_gene_region_total(tmp_path: Path):
    junctions = tmp_path / "Chimeric.out.junction"
    junctions.write_text(
        "chr22\t101\t+\tchr11\t201\t-\t1\t0\t0\texactA\n"
        "chr22\t101\t+\tchr11\t201\t-\t1\t0\t0\texactA\n"
        "chr22\t101\t+\tchr11\t501\t-\t1\t0\t0\totherBreakpoint\n",
        encoding="utf-8",
    )

    assert MODULE.star_chimeric_support_count(
        junctions, "EWSR1::WT1", "chr22:100:+", "chr11:200:-", tolerance=2
    ) == 1


def test_breakpoint_windows_prove_residue_level_fusion_boundary():
    rows = MODULE.breakpoint_window_records(
        "ACDEFGHIKLMNPQRST", 5, left_gene="EWSR1", right_gene="WT1", lengths=(9,)
    )
    assert rows
    assert all(row["crosses_junction"] == "yes" for row in rows)
    assert all(row["fusion_left_peptide"] and row["fusion_right_peptide"] for row in rows)
    assert all(
        row["fusion_left_peptide"] + row["fusion_right_peptide"] == row["peptide"]
        for row in rows
    )
    assert all("|" in row["fusion_junction_display"] for row in rows)


def test_external_confirmed_product_closes_only_matching_adjacency(tmp_path: Path):
    products = tmp_path / "expressed_products.tsv"
    fieldnames = [
        "genome_build", "chrom1", "pos1", "strand1", "chrom2", "pos2", "strand2",
        "gene1", "gene2", "transcript1", "transcript2", "protein_sequence",
        "junction_aa_position", "in_frame", "orf_status", "source_tool", "source_record_id",
    ]
    with products.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerow({
            "genome_build": "GRCh38", "chrom1": "chr22", "pos1": "100", "strand1": "+",
            "chrom2": "chr11", "pos2": "200", "strand2": "-", "gene1": "EWSR1",
            "gene2": "WT1", "transcript1": "ENST0001", "transcript2": "ENST0002",
            "protein_sequence": "ACDEFGHIKLMNPQRST", "junction_aa_position": "8",
            "in_frame": "yes", "orf_status": "CONFIRMED", "source_tool": "AGFusion",
            "source_record_id": "record-1",
        })
    event = {"event_id": "E1", "sample_id": "S1", "gene": "EWSR1::WT1"}
    audit = [{
        "event_id": "E1", "left_breakpoint": "chr22:100:+", "right_breakpoint": "chr11:200:-",
    }]
    peptides: list[dict[str, str]] = []
    origin = MODULE.apply_confirmed_expressed_products(
        [event], peptides, audit, products, ["HLA-A*02:01"], genome_build="GRCh38"
    )
    assert origin
    assert event["orf_status"] == "CONFIRMED"
    assert all(row["crosses_junction"] == "yes" for row in peptides)
    assert all(row["fusion_left_peptide"] and row["fusion_right_peptide"] for row in peptides)


def test_closed_source_chain_requires_exact_junction_support():
    rows = MODULE.fusion_peptide_origin_chain_rows([{
        "event_id": "E1", "adjacency_key": "GRCh38|chr22:100:+|chr11:200:-",
        "fusion_transcript_id": "ENST1::ENST2", "orf_id": "ORF1",
        "fusion_protein_sequence": "ACDEFGHIKLMNPQRST", "orf_status": "CONFIRMED",
        "junction_position_in_peptide_1based": "5", "fusion_left_peptide": "ACDEF",
        "fusion_right_peptide": "GHIK", "verified_rna_junction_reads": "2",
    }])
    assert rows[0]["source_chain_status"] == "ORF_CLOSED_EXACT_JUNCTION_INCOMPLETE"
