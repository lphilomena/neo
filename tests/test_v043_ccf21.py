from pathlib import Path

from neoag_v03.ccf_v2 import build_ccf_2
from neoag_v03.config import load_profile
from neoag_v03.utils import read_tsv, write_tsv


def test_ccf21_writes_input_qc_and_confidence_fields(tmp_path):
    events = tmp_path / "events.tsv"
    purity = tmp_path / "purity.tsv"
    cnv = tmp_path / "cnv.tsv"
    out = tmp_path / "ccf" / "ccf_lite.tsv"
    write_tsv(events, [{
        "event_id": "E1", "sample_id": "S1", "mutation_source": "SNV", "chrom": "chr1", "pos": "150",
        "tumor_vaf": "0.35", "tumor_depth": "100", "tumor_alt_count": "35",
    }])
    write_tsv(purity, [{"purity": "0.70", "ploidy": "2.4", "source": "FACETS", "confidence": "high"}])
    write_tsv(cnv, [{"chrom": "chr1", "start": "1", "end": "1000", "total_cn": "3", "major_cn": "2", "minor_cn": "1"}])
    rows = build_ccf_2(events, purity, cnv, load_profile("default"), out)
    assert out.exists()
    assert (out.parent / "ccf_input_qc.tsv").exists()
    assert (out.parent / "ccf_conflicts.tsv").exists()
    assert (out.parent / "ccf_cluster.tsv").exists()
    row = rows[0]
    assert row["multiplicity_candidates"]
    assert row["multiplicity_confidence"] in {"high", "medium", "low", "ambiguous"}
    assert row["clonality_confidence"] in {"high", "medium", "low", "unresolved"}
    assert row["probability_ccf_gt_0_8"] != ""
    qc = read_tsv(out.parent / "ccf_input_qc.tsv")[0]
    assert qc["ccf_ready_status"] in {"ready", "approximate"}
    assert qc["major_minor_cn_available"] == "yes"


def test_ccf21_rna_only_and_wes_sv_are_not_overconfident(tmp_path):
    events = tmp_path / "events.tsv"
    purity = tmp_path / "purity.tsv"
    cnv = tmp_path / "cnv.tsv"
    write_tsv(events, [
        {"event_id": "R1", "sample_id": "S1", "event_type": "RNA fusion", "source": "rna", "chrom": "chr1", "pos": "10"},
        {"event_id": "W1", "sample_id": "S1", "mutation_source": "WES SV", "event_type": "SV_Junction", "evidence_scope": "EXOME_CAPTURE_LIMITED", "chrom": "chr1", "pos": "20", "tumor_vaf": "0.20", "tumor_depth": "40", "tumor_alt_count": "8"},
    ])
    write_tsv(purity, [{"purity": "0.60", "ploidy": "2"}])
    write_tsv(cnv, [{"chrom": "chr1", "start": "1", "end": "1000", "total_cn": "2", "major_cn": "1", "minor_cn": "1"}])
    rows = build_ccf_2(events, purity, cnv, load_profile("default"), tmp_path / "ccf.tsv")
    by_id = {r["event_id"]: r for r in rows}
    assert by_id["R1"]["ccf_method"] == "RNA_ONLY_UNRESOLVED"
    assert by_id["R1"]["clonality_status"] == "unresolved"
    assert by_id["R1"]["ccf_best"] == ""
    assert by_id["W1"]["ccf_method"] == "WES_SV_CAPTURE_LIMITED_APPROX"
    assert by_id["W1"]["ccf_confidence"] == "low"
    assert "wes_sv_capture_limited" in by_id["W1"]["ccf_warning"]


def test_ccf21_external_clonality_conflict_and_clusters(tmp_path):
    events = tmp_path / "events.tsv"
    purity = tmp_path / "purity.tsv"
    cnv = tmp_path / "cnv.tsv"
    ext = tmp_path / "pyclone.tsv"
    out = tmp_path / "ccf.tsv"
    write_tsv(events, [{
        "event_id": "E1", "sample_id": "S1", "mutation_source": "SNV", "chrom": "chr1", "pos": "100",
        "tumor_vaf": "0.45", "tumor_depth": "100", "tumor_alt_count": "45",
    }])
    write_tsv(purity, [{"purity": "0.80", "ploidy": "2"}])
    write_tsv(cnv, [{"chrom": "chr1", "start": "1", "end": "1000", "total_cn": "2", "major_cn": "1", "minor_cn": "1"}])
    write_tsv(ext, [{"event_id": "E1", "cluster_id": "C2", "cellular_prevalence": "0.20", "cellular_prevalence_low": "0.15", "cellular_prevalence_high": "0.25", "cluster_assignment_probability": "0.93", "tool": "PyClone-VI"}])
    rows = build_ccf_2(events, purity, cnv, load_profile("default"), out, external_clonality_tsv=ext)
    assert rows[0]["external_clonality_tool"] == "PyClone-VI"
    assert rows[0]["external_cluster_id"] == "C2"
    assert rows[0]["ccf_resolution"] in {"external_conflict_review", "external_supported"}
    conflicts = read_tsv(out.parent / "ccf_conflicts.tsv")
    assert conflicts
    clusters = read_tsv(out.parent / "ccf_cluster.tsv")
    assert clusters[0]["cluster_id"] == "C2"
    assert clusters[0]["external_tool"] == "PyClone-VI"


def test_ccf21_svclone_preferred_for_sv_event(tmp_path):
    events = tmp_path / "events.tsv"
    purity = tmp_path / "purity.tsv"
    cnv = tmp_path / "cnv.tsv"
    svclone = tmp_path / "svclone.tsv"
    write_tsv(events, [{
        "event_id": "SV1", "sample_id": "S1", "mutation_source": "SV", "event_type": "SV_Fusion", "chrom": "chr2", "pos": "200",
        "tumor_vaf": "0.10", "tumor_depth": "60", "tumor_alt_count": "6",
    }])
    write_tsv(purity, [{"purity": "0.55", "ploidy": "2"}])
    write_tsv(cnv, [{"chrom": "chr2", "start": "1", "end": "1000", "total_cn": "2", "major_cn": "1", "minor_cn": "1"}])
    write_tsv(svclone, [{"sv_event_id": "SV1", "sv_ccf": "0.88", "sv_copy_number": "1", "confidence": "high"}])
    rows = build_ccf_2(events, purity, cnv, load_profile("default"), tmp_path / "ccf.tsv", svclone_tsv=svclone)
    assert rows[0]["external_clonality_tool"] == "SVclone"
    assert rows[0]["svclone_ccf"] == "0.88"
    assert rows[0]["ccf_resolution"] == "external_svclone_preferred"
