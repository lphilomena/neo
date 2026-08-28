from pathlib import Path

import shutil

from neoag.evidence_layer import build_standard_evidence_layer
from neoag.input_router import build_raw_intermediates
from neoag.pipeline import run

ROOT = Path(__file__).resolve().parents[1]


def test_build_intermediates_pvac(tmp_path):
    cfg = {
        "sample": {"id": "S1", "profile": "default"},
        "inputs": {
            "entry_mode": "snv_indel",
            "pvac_files": [
                str(ROOT / "data/fixtures/pvacseq_aggregated.tsv"),
                str(ROOT / "data/fixtures/pvacfuse_aggregated.tsv"),
            ],
        },
    }
    paths = build_raw_intermediates(cfg, tmp_path / "parsed_layer", root=ROOT)
    assert Path(paths["raw_events"]).is_file()
    assert Path(paths["raw_peptides"]).is_file()


def test_build_intermediates_passthrough(tmp_path):
    src = tmp_path / "src"
    src.mkdir()
    cfg = {
        "sample": {"id": "S1", "profile": "default"},
        "inputs": {"entry_mode": "intermediates"},
    }
    build_raw_intermediates(
        {
            "sample": {"id": "S1", "profile": "default"},
            "inputs": {
                "entry_mode": "snv_indel",
                "pvac_files": [str(ROOT / "data/fixtures/pvacseq_aggregated.tsv")],
            },
        },
        src,
        root=ROOT,
    )
    pre_e = src / "parsed/raw_events.tsv"
    pre_p = src / "parsed/raw_peptides.tsv"
    dst = tmp_path / "dst"
    paths = build_raw_intermediates(
        {
            **cfg,
            "inputs": {
                "entry_mode": "intermediates",
                "raw_events": str(pre_e),
                "raw_peptides": str(pre_p),
            },
        },
        dst,
        root=ROOT,
    )
    assert Path(paths["raw_events"]).read_text() == pre_e.read_text()


def test_evidence_layer_and_run_from_intermediates(tmp_path):
    src = tmp_path / "src"
    build_raw_intermediates(
        {
            "sample": {"id": "DEMO", "profile": "leukemia"},
            "inputs": {
                "entry_mode": "snv_indel",
                "pvac_files": [
                    str(ROOT / "data/fixtures/pvacseq_aggregated.tsv"),
                    str(ROOT / "data/fixtures/pvacfuse_aggregated.tsv"),
                ],
            },
        },
        src,
        root=ROOT,
    )
    evidence = build_standard_evidence_layer(
        src,
        "leukemia",
        normal_expression=ROOT / "resources/normal_expression.example.tsv",
        normal_hla_ligands=ROOT / "resources/normal_hla_ligands.example.tsv",
        sample_id="DEMO",
    )
    for key in ("expression_evidence", "rna_junction_evidence", "safety_evidence"):
        assert Path(evidence[key]).is_file()

    out = run(
        tmp_path / "score",
        "leukemia",
        "DEMO",
        pvac_paths=[],
        raw_events=src / "parsed/raw_events.tsv",
        raw_peptides=src / "parsed/raw_peptides.tsv",
        netmhcpan=ROOT / "data/fixtures/netmhcpan_example.xls",
        mhcflurry=ROOT / "data/fixtures/mhcflurry_predictions.csv",
        vep_appm=ROOT / "data/fixtures/vep_appm.tsv",
        expression=ROOT / "data/fixtures/gene_expression.tsv",
        hla_loh=ROOT / "data/fixtures/hla_loh.tsv",
        purity=ROOT / "data/fixtures/purity.tsv",
        cnv=ROOT / "data/fixtures/cnv_segments.tsv",
        normal_expression=ROOT / "resources/normal_expression.example.tsv",
        normal_hla_ligands=ROOT / "resources/normal_hla_ligands.example.tsv",
        entry_mode="intermediates",
    )
    assert Path(out["expression_evidence"]).exists()
    assert Path(out["safety_evidence"]).exists()
    txt = Path(out["ranked_peptides"]).read_text()
    assert "immunology_composite_score" in txt
