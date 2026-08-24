from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "prepare_snaf_splice_branch", ROOT / "scripts/prepare_snaf_splice_branch.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_load_tpm_maps_accepts_lowercase_rsem_headers(tmp_path: Path):
    genes = tmp_path / "gene_tpm.tsv"
    transcripts = tmp_path / "transcript_tpm.tsv"
    genes.write_text("gene_id\ttpm\nENSG00000139304.16\t18.85\n", encoding="utf-8")
    transcripts.write_text("transcript_id\ttpm\nENST00000376621.8\t4.25\n", encoding="utf-8")

    gene_tpm, transcript_tpm = MODULE.load_tpm_maps(genes, transcripts)

    assert gene_tpm["ENSG00000139304"] == 18.85
    assert transcript_tpm["ENST00000376621"] == 4.25


def test_snaf_uid_ensembl_id_is_not_hidden_by_symbol():
    ensembl, gene = MODULE.resolve_candidate_gene(
        {"uid": "ENSG00000139304:E45.1-E46.1", "symbol": "PTPRQ"},
        {"ENSG00000139304": "PTPRQ"},
    )

    assert ensembl == "ENSG00000139304"
    assert gene == "PTPRQ"
