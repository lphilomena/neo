from pathlib import Path

from neoag.adapters.iedb_immunogenicity import predict_immunogenicity
from neoag.adapters.netmhcstabpan import parse_netmhcstabpan
from neoag.benchmark_improve import (
    DEFAULT_CEDAR,
    TOOL_CACHE_FILES,
    _seed_tools_cache,
    auroc,
    load_improve_tsv,
    normalize_hla,
    run_benchmark_improve,
)

ROOT = Path(__file__).resolve().parents[1]


def test_normalize_hla():
    assert normalize_hla("HLA-A02:01") == "HLA-A*02:01"
    assert normalize_hla("HLA-A11*01") == "HLA-A*11:01"
    assert normalize_hla("HLA-B35*01") == "HLA-B*35:01"


def test_load_cedar_sample():
    rows = load_improve_tsv(DEFAULT_CEDAR, limit=50)
    assert len(rows) == 50
    assert sum(r.response for r in rows) >= 0


def test_auroc_perfect():
    labels = [0, 0, 1, 1]
    scores = [0.1, 0.2, 0.8, 0.9]
    assert auroc(labels, scores) == 1.0


def test_iedb_immunogenicity_known_allele():
    score = predict_immunogenicity("YMDGTMSQV", "HLA-A*02:01")
    assert isinstance(score, float)
    assert -1.0 < score < 1.0


def test_parse_netmhcstabpan_fixture():
    rows = parse_netmhcstabpan(ROOT / "data" / "fixtures" / "netmhcstabpan_example.tsv")
    assert len(rows) >= 1
    assert rows[0]["netmhcstabpan_score"]
    assert rows[0]["netmhcstabpan_rank"]


def test_seed_tools_cache(tmp_path):
    cache = tmp_path / "cache" / "tools"
    cache.mkdir(parents=True)
    for name in TOOL_CACHE_FILES.values():
        (cache / name).write_text(f"fixture:{name}\n", encoding="utf-8")
    dest = tmp_path / "dest" / "tools"
    dest.mkdir(parents=True)

    cache_root, copied = _seed_tools_cache(
        dest,
        cache,
        skip_netmhcpan=True,
        skip_mhcflurry=True,
        skip_stabpan=False,
    )
    assert cache_root == cache
    assert copied == [TOOL_CACHE_FILES["netmhcpan"], TOOL_CACHE_FILES["mhcflurry"]]
    assert (dest / TOOL_CACHE_FILES["netmhcpan"]).is_file()
    assert (dest / TOOL_CACHE_FILES["mhcflurry"]).is_file()
    assert not (dest / TOOL_CACHE_FILES["stabpan"]).exists()


def test_benchmark_improve_stub(tmp_path):
    out = run_benchmark_improve(
        input_path=DEFAULT_CEDAR,
        outdir=tmp_path / "bench",
        stub=True,
        limit=25,
        skip_netmhcpan=False,
        skip_mhcflurry=False,
    )
    assert out["n_records"] == 25
    assert (tmp_path / "bench" / "benchmark_metrics.tsv").is_file()
    assert (tmp_path / "bench" / "benchmark_report.md").is_file()
    assert (tmp_path / "bench" / "presentation" / "prime_evidence.tsv").is_file()
    assert (tmp_path / "bench" / "presentation" / "bigmhc_im_evidence.tsv").is_file()
    assert (tmp_path / "bench" / "presentation" / "netmhcstabpan_evidence.tsv").is_file()
    metrics = (tmp_path / "bench" / "benchmark_metrics.tsv").read_text()
    assert "presentation_evidence_score" in metrics
    assert "immunogenicity_composite_score" in metrics
    assert "prime_score" in metrics
    assert "bigmhc_im_score" in metrics
    assert "netmhcstabpan_score" in metrics
    assert "netmhcstabpan_rank" in metrics
    assert "efficacy_score" in metrics
    assert "cohort" in metrics
    assert out["cohorts"]["all"] == 25
    assert "gated_tesla" in out["cohorts"]
    report = (tmp_path / "bench" / "benchmark_report.md").read_text()
    assert "HLA-stratified AUROC" in report
    assert "immunogenicity_composite_score" in report
    assert "Metrics — all pairs" in report
    assert "efficacy_score" in report
    preds = (tmp_path / "bench" / "benchmark_predictions.tsv").read_text()
    assert "efficacy_score" in preds
    assert "presentation_gate_status" in preds
