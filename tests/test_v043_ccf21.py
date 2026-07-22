import json
from pathlib import Path

from neoag.ccf_v2 import build_ccf_2
from neoag.config import load_profile
from neoag.utils import read_tsv, write_tsv


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
    assert by_id["R1"]["clonality_status"] == "RNA_ONLY_UNRESOLVED"
    assert by_id["R1"]["ccf_best"] == ""
    assert by_id["R1"]["ccf_estimate"] == "NA"
    assert by_id["R1"]["ccf_confidence"] == "unresolved"
    assert by_id["R1"]["clonality_multiplier"] == "1.0000"
    assert by_id["W1"]["ccf_method"] == "WES_SV_CAPTURE_LIMITED_APPROX"
    assert by_id["W1"]["ccf_confidence"] == "low"
    assert "wes_sv_capture_limited" in by_id["W1"]["ccf_warning"]


def test_ccf21_clamps_estimate_and_preserves_raw_ccf(tmp_path):
    events = tmp_path / "events.tsv"
    purity = tmp_path / "purity.tsv"
    cnv = tmp_path / "cnv.tsv"
    write_tsv(events, [{
        "event_id": "E1", "sample_id": "S1", "mutation_source": "SNV",
        "chrom": "chr1", "pos": "100", "tumor_vaf": "0.127",
        "tumor_depth": "100", "tumor_alt_count": "13",
    }])
    write_tsv(purity, [{"purity": "0.12", "ploidy": "2"}])
    write_tsv(cnv, [{
        "chrom": "chr1", "start": "1", "end": "1000", "total_cn": "2",
        "major_cn": "1", "minor_cn": "1",
    }])
    row = build_ccf_2(events, purity, cnv, load_profile("default"), tmp_path / "ccf.tsv")[0]
    assert float(row["raw_ccf"]) > 1.0
    assert row["ccf_estimate"] == "1.0000"
    assert row["ccf_best"] == "1.0000"
    assert "raw_ccf_above_1_clamped" in row["ccf_warning"]


def test_ccf21_easyfuse_without_dna_vaf_is_rna_only_unresolved(tmp_path):
    events = tmp_path / "events.tsv"
    purity = tmp_path / "purity.tsv"
    write_tsv(events, [{
        "event_id": "EF1", "sample_id": "S1", "event_type": "Fusion",
        "mutation_source": "SV", "source": "easyfuse:fusions.pass.csv",
        "rna_junction_reads": "25", "tumor_vaf": "0.0", "chrom": "chr1", "pos": "100",
    }])
    write_tsv(purity, [{"purity": "0.20", "ploidy": "2"}])
    row = build_ccf_2(events, purity, None, load_profile("default"), tmp_path / "ccf.tsv")[0]
    assert row["ccf_method"] == "RNA_ONLY_UNRESOLVED"
    assert row["ccf_status"] == "RNA_ONLY_UNRESOLVED"
    assert row["ccf_estimate"] == "NA"
    assert row["clonality_multiplier"] == "1.0000"


def test_ccf21_accepts_sequenza_segment_columns(tmp_path):
    events = tmp_path / "events.tsv"
    purity = tmp_path / "purity.tsv"
    cnv = tmp_path / "sequenza_segments.tsv"
    write_tsv(events, [{
        "event_id": "E1", "sample_id": "S1", "mutation_source": "SNV",
        "chrom": "chr1", "pos": "150", "tumor_vaf": "0.10",
        "tumor_depth": "100", "tumor_alt_count": "10",
    }])
    write_tsv(purity, [{"purity": "0.30", "ploidy": "2.2", "source": "sequenza"}])
    write_tsv(cnv, [{
        "chromosome": "chr1", "start.pos": "100", "end.pos": "200",
        "CNt": "5", "A": "4", "B": "1",
    }])
    row = build_ccf_2(events, purity, cnv, load_profile("default"), tmp_path / "ccf.tsv")[0]
    assert row["total_cn"] == "5.0000"
    assert row["major_cn"] == "4.0000"
    assert row["minor_cn"] == "1.0000"
    assert row["cnv_confidence"] == "high"


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


def _write_ccf_inputs(tmp_path):
    events = tmp_path / "events.tsv"
    cnv = tmp_path / "cnv.tsv"
    write_tsv(events, [{
        "event_id": "E1", "sample_id": "S1", "mutation_source": "SNV",
        "chrom": "chr1", "pos": "150", "tumor_vaf": "0.25",
        "tumor_depth": "100", "tumor_alt_count": "25",
    }])
    write_tsv(cnv, [{
        "chrom": "chr1", "start": "1", "end": "1000", "total_cn": "2",
        "major_cn": "1", "minor_cn": "1", "cnv_confidence": "high",
    }])
    return events, cnv


def test_ccf21_uses_concordant_multi_tool_median(tmp_path):
    events, cnv = _write_ccf_inputs(tmp_path)
    recommendation = tmp_path / "purity_recommendation.json"
    recommendation.write_text(json.dumps({
        "status": "CONCORDANT", "recommended_purity": 0.7, "range": "0.6800-0.7200",
        "n_tools": 3, "tool_values": {"FACETS": 0.68, "PURPLE": 0.7, "Sequenza": 0.72},
    }), encoding="utf-8")

    rows = build_ccf_2(events, recommendation, cnv, load_profile("default"), tmp_path / "ccf.tsv")

    row = rows[0]
    assert row["purity"] == "0.7000"
    assert row["purity_source"] == "multi_tool_median"
    assert row["purity_confidence"] == "high"
    assert row["purity_consensus_status"] == "CONCORDANT"
    assert json.loads(row["purity_tool_values"])["FACETS"] == 0.68
    assert row["ccf_confidence"] in {"high", "medium"}
    assert "purity_" not in row["ccf_warning"]


def test_ccf21_strong_purity_discordance_forces_low_confidence(tmp_path):
    events, cnv = _write_ccf_inputs(tmp_path)
    recommendation = tmp_path / "purity_recommendation.json"
    recommendation.write_text(json.dumps({
        "status": "STRONG_DISCORDANCE", "recommended_purity": 0.52,
        "range": "0.2500-0.8200", "n_tools": 3,
        "tool_values": {"FACETS": 0.25, "PURPLE": 0.52, "ASCAT": 0.82},
    }), encoding="utf-8")

    out = tmp_path / "ccf.tsv"
    rows = build_ccf_2(events, recommendation, cnv, load_profile("default"), out)

    row = rows[0]
    assert row["purity"] == "0.5200"
    assert row["purity_consensus_status"] == "STRONG_DISCORDANCE"
    assert row["ccf_confidence"] == "low"
    assert "purity_strong_discordance" in row["ccf_warning"]
    qc = read_tsv(out.parent / "ccf_input_qc.tsv")[0]
    assert qc["ccf_ready_status"] == "limited"
    assert json.loads(qc["purity_tool_values"])["ASCAT"] == 0.82
