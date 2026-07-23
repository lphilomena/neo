import csv
import gzip
import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "build_hla_ligand_atlas_reference.py"


def load_builder():
    spec = importlib.util.spec_from_file_location("hla_ligand_builder", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


def write_gzip_tsv(path, fields, rows):
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def test_build_aggregates_tissues_and_hla_classes(tmp_path):
    builder = load_builder()
    peptides = tmp_path / "peptides.tsv.gz"
    hits = tmp_path / "sample_hits.tsv.gz"
    output = tmp_path / "normal_ms_ligands.tsv"
    write_gzip_tsv(peptides, ["peptide_sequence_id", "peptide_sequence"], [
        {"peptide_sequence_id": "1", "peptide_sequence": "AAAAAAAAA"},
        {"peptide_sequence_id": "2", "peptide_sequence": "CCCCCCCCC"},
    ])
    write_gzip_tsv(hits, ["peptide_sequence_id", "donor", "tissue", "hla_class"], [
        {"peptide_sequence_id": "1", "donor": "D1", "tissue": "Lung", "hla_class": "HLA-II"},
        {"peptide_sequence_id": "1", "donor": "D2", "tissue": "Liver", "hla_class": "HLA-I"},
        {"peptide_sequence_id": "2", "donor": "D1", "tissue": "Skin", "hla_class": "HLA-I"},
    ])

    stats = builder.build(peptides, hits, output, "2020.12")
    with output.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert stats["rows"] == 2
    assert stats["class_counts"] == {"HLA-I": 1, "HLA-I+II": 1}
    assert rows[0]["source_tissue"] == "Liver,Lung"
    assert rows[0]["hla_class"] == "HLA-I+II"
    assert rows[0]["hla_allele"] == ""
