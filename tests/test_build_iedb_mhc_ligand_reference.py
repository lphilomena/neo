import csv
import importlib.util
import json
import zipfile
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "build_iedb_mhc_ligand_reference.py"


def load_builder():
    spec = importlib.util.spec_from_file_location("iedb_builder", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


def write_fixture(path):
    groups = ["Assay ID", "Epitope", "Epitope", "Epitope", "Host", "in vivo Process", "Assay", "Assay", "Antigen Presenting Cell", "MHC Restriction", "MHC Restriction"]
    fields = ["IEDB IRI", "Epitope IRI", "Object Type", "Name", "Name", "Disease", "Method", "Qualitative Measurement", "Culture Condition", "Name", "Class"]
    rows = [
        ["a1", "e1", "Linear peptide", "AAAAAAAAA", "Homo sapiens (human)", "healthy", "cellular MHC/mass spectrometry", "Positive", "Direct Ex Vivo", "HLA-A*02:01", "I"],
        ["a2", "e2", "Linear peptide", "CCCCCCCCC", "Homo sapiens (human)", "cancer", "mass spectrometry", "Positive", "Cell Line / Clone", "HLA-B*07:02", "I"],
        ["a3", "e3", "Linear peptide", "DDDDDDDDD", "Homo sapiens (human)", "", "mass spectrometry", "Negative", "Direct Ex Vivo", "HLA-DR", "II"],
        ["a4", "e4", "Linear peptide", "EEEEEEEEE", "Homo sapiens (human)", "", "mass spectrometry", "Positive", "Direct Ex Vivo", "HLA-DR", "II"],
        ["a5", "e5", "Linear peptide", "PEP*IDE", "Homo sapiens (human)", "", "mass spectrometry", "Positive", "Direct Ex Vivo", "HLA-A*01:01", "I"],
    ]
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        lines = []
        lines.append("\t".join(groups))
        lines.append("\t".join(fields))
        for row in rows:
            lines.append("\t".join(row))
        archive.writestr("mhc_ligand_full.tsv", "\n".join(lines) + "\n")


def test_build_records_filters_counts_and_checksums(tmp_path):
    builder = load_builder()
    archive = tmp_path / "mhc_ligand_full_tsv.zip"
    write_fixture(archive)
    outdir = tmp_path / "build"
    args = builder.parse_args(["--archive", str(archive), "--outdir", str(outdir), "--source-release", "fixture"])
    builder.build(args)

    manifest = json.loads((outdir / "manifest.json").read_text())
    assert manifest["raw_archive"]["sha256"] == builder.sha256(archive)
    assert manifest["counts"]["raw_data_rows"] == 5
    assert manifest["counts"]["accepted_assay_rows"] == 3
    assert manifest["counts"]["strict_normal_assay_rows"] == 2
    assert manifest["outputs"]["all_human_ms_peptides"]["rows"] == 3
    assert manifest["outputs"]["strict_normal_direct_ex_vivo"]["rows"] == 2
    assert manifest["outputs"]["detail"]["sha256"] == builder.sha256(outdir / "build" / "iedb_human_ms_ligands_detail.tsv")
    assert "not_positive" in manifest["counts"]["rejection_reasons"]
    assert "noncanonical_or_empty_peptide" in manifest["counts"]["rejection_reasons"]

    with (outdir / "build" / "iedb_human_ms_ligands.tsv").open() as handle:
        peptides = [row["peptide"] for row in csv.DictReader(handle, delimiter="\t")]
    assert peptides == sorted(peptides)
    assert (outdir / "SHA256SUMS").read_text().count("\n") == 4
