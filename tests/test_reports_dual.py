import hashlib
import json
from pathlib import Path

from neoag.reports_dual import ReportBundle, _patient_event_grade_counts, _patient_evidence_summary, _patient_event_change, _patient_limitation, _patient_rna_measurements, _replace_gene_ids, load_report_bundle, make_dual_reports, make_patient_report, make_technical_report
from neoag.utils import write_tsv


def _bundle():
    events = [{
        "event_id": "B2M|chr15:1A>G", "event_name": "B2M variant", "event_type": "SNV",
        "gene": "B2M", "event_score": "0.6", "ccf_status": "clonal", "safety_status": "PASS",
    }]
    peptides = [{
        "peptide_id": "P1", "event_id": "B2M|chr15:1A>G", "sample_id": "S1",
        "gene": "B2M", "peptide": "AAAAAAAAA", "wildtype_peptide": "AAABAAAAA",
        "hla_allele": "HLA-A*02:01", "event_type": "SNV", "peptide_consequence": "missense",
        "source_chain_confidence_tier": "C2", "netmhcpan_el_rank": "0.5",
        "gene_expression_tpm": "12.34", "transcript_expression_tpm": "3.5",
        "rna_depth": "20", "rna_alt_reads": "4", "rna_vaf": "0.2",
        "final_priority": "B", "efficacy_score": "0.72", "safety_status": "PASS",
        "combined_protein_change": "B2M:p.A1V",
        "presentation_evidence_grade": "A", "recommended_use": "MHC-I short peptide ELISpot/tetramer (MT + WT pair)",
    }]
    validation = [{
        "peptide_id": "P1", "validation_mode": "missense_short_pair",
        "recommended_assay": "MHC-I short peptide ELISpot/tetramer (MT + WT pair)",
        "validation_strategy": "Mutant short peptide with WT control",
    }]
    return ReportBundle(
        profile={"_profile_name": "test", "gates": {"max_el_rank": 2.0}},
        events=events,
        peptides=peptides,
        appm_summary={"mhc_i_integrity_status": "MHC_I_INTACT", "mhc_i_integrity_score": "0.95"},
        validation_rows=validation,
        sample_id="S1",
        entry_mode="snv_indel",
        provenance={"sample_id": "S1", "tools": {"netmhcpan": {"status": "real", "version": "4.2", "file": "/data/netmhcpan.xls"}}},
    )


def test_patient_report_is_plain_language(tmp_path):
    out = tmp_path / "patient.html"
    make_patient_report(out, _bundle())
    text = out.read_text(encoding="utf-8")
    assert "患者沟通版" in text or "新抗原计算分析报告" in text
    assert "不能替代临床诊断" in text
    assert "netmhcpan.xls" not in text
    assert "fastq" not in text.lower()
    assert "移码变异" not in text or "点突变" in text
    assert "重点变异事件（按类型、按事件去重）" in text
    assert "关键人工审阅事件" in text
    assert "进入本表不等于自动升级为R1/R2" in text
    assert "关键证据与下一步" in text
    section6 = text.split("6. Top候选综合证据与实验建议", 1)[1].split("7. 分析方法与工具状态", 1)[0]
    assert "<th>候选</th>" in section6
    assert "<th>为什么值得关注</th>" in section6
    assert "<th>当前不确定性</th>" in section6
    assert "<th>建议下一步</th>" in section6
    assert "<th>呈递工具</th>" not in section6
    assert "基因表达 12.3400 TPM" in text
    assert "转录本表达 3.5000 TPM" in text
    assert "RNA位点深度 20" in text
    assert "RNA alt reads 4" in text
    assert "RNA VAF 0.2000" in text
    assert "5. 候选肽段Top 50（跨赛道、按事件去重）" in text
    assert "重点候选 Top 50" in text
    candidate_section = text.split("5. 候选肽段Top 50（跨赛道、按事件去重）", 1)[1].split("6. Top候选综合证据与实验建议", 1)[0]
    assert "<th>关键证据与下一步</th>" in candidate_section
    assert "<th>关键证据</th>" not in candidate_section
    assert "<th>主要限制</th>" not in candidate_section
    assert "<th>建议实验</th>" not in candidate_section
    assert "<th>审阅状态</th>" not in candidate_section
    assert "本表保留R3-GAP候选" in text
    assert "补证据后再决定" not in text
    assert "此类候选存在关键证据缺口" not in text
    assert "附录A：R1–R4证据分层" in text
    assert "附录B：术语说明" in text
    assert "附件与可追溯文件" in text
    assert "缺失证据统一视为未评估" in text
    assert "DQA1/DQB1" not in text
    assert "EWSR1::WT1" not in text


def test_patient_summary_separates_event_and_peptide_hla_counts(tmp_path):
    bundle = _bundle()
    bundle.events = [
        {"event_group_id": "E1", "best_evidence_grade": "R1"},
        {"event_group_id": "E2", "best_evidence_grade": "R3", "evidence_missing_layers": "rna"},
        {"event_group_id": "E3", "best_evidence_grade": "R3", "evidence_conflict_layers": "presentation"},
        {"event_group_id": "E4", "best_evidence_grade": "R4"},
    ]
    bundle.peptides.append({**bundle.peptides[0], "peptide_id": "P2", "peptide": "BBBBBBBBB"})
    counts = _patient_event_grade_counts(bundle.events)
    assert counts["R1"] == 1
    assert counts["R3-GAP"] == 1
    assert counts["R3-REVIEW"] == 1
    assert counts["R4"] == 1
    out = tmp_path / "patient_counts.html"
    make_patient_report(out, bundle)
    text = out.read_text(encoding="utf-8")
    assert "筛选规模" in text
    assert "事件级结论" in text
    assert "直接读取ranked_events.evidence_consensus.tsv" in text
    assert "同一事件可产生多个重叠肽段和多个HLA组合" in text
    assert "<td>R3-GAP</td><td>1</td>" in text
    assert "<td>R3-REVIEW</td><td>1</td>" in text


def test_patient_event_top_table_uses_event_level_r3_subgrade(tmp_path):
    bundle = _bundle()
    bundle.events[0].update({
        "best_evidence_grade": "R3",
        "evidence_missing_layers": "rna_variant_support",
    })
    out = tmp_path / "patient_r3_gap.html"
    make_patient_report(out, bundle)
    text = out.read_text(encoding="utf-8")
    assert "SNV Top 10" in text
    assert "<td>R3-GAP</td>" in text


def test_patient_event_evidence_is_compact_and_omits_rna_source():
    row = {
        "event_type": "Fusion",
        "gene_expression_tpm": "1.25",
        "rna_junction_reads": "18",
        "rna_junction_source": "raw_events",
    }
    text = _patient_rna_measurements(row)
    assert "基因表达 1.2500 TPM" in text
    assert "junction reads 18" in text
    assert "数据来源" not in text
    assert "raw_events" not in text


def test_patient_report_top_limits_are_configurable(tmp_path):
    out = tmp_path / "patient_custom_limits.html"
    make_patient_report(out, _bundle(), event_top_n=1, candidate_top_n=3)
    text = out.read_text(encoding="utf-8")
    assert "SNV Top 1" in text
    assert "5. 候选肽段Top 3（跨赛道、按事件去重）" in text


def test_patient_report_has_track_top5_when_present(tmp_path):
    bundle = _bundle()
    bundle.peptides.extend([
        {
            "peptide_id": "P2", "event_id": "FUS1", "gene": "GENE1::GENE2",
            "peptide": "BBBBBBBBB", "hla_allele": "HLA-B*07:02", "event_type": "Fusion",
            "pipeline_r_grade": "R3", "rna_support_status": "RNA_JUNCTION_SUPPORTED",
            "source_chain_confidence_tier": "C3", "netmhcpan_el_rank": "0.7",
        },
        {
            "peptide_id": "P3", "event_id": "SPL1", "gene": "GENE3",
            "peptide": "CCCCCCCCC", "hla_allele": "HLA-C*07:02", "event_type": "Splice",
            "pipeline_r_grade": "R2", "rna_support_status": "RNA_JUNCTION_SUPPORTED",
            "source_chain_confidence_tier": "C3", "mhcflurry_presentation_score": "0.8",
        },
    ])
    out = tmp_path / "patient_tracks.html"
    make_patient_report(out, bundle)
    text = out.read_text(encoding="utf-8")
    assert "Fusion Top 10" in text
    assert "Splice Top 10" in text
    assert "GENE1::GENE2" in text
    assert "<th>来源链</th>" not in text
    assert "<th>肽段-HLA</th>" in text
    assert "BBBBBBBBB / HLA-B*07:02" in text
    assert "CCCCCCCCC / HLA-C*07:02" in text


def test_event_top_uses_ranked_events_before_candidate_integrity_filter(tmp_path):
    bundle = _bundle()
    for index in range(10):
        event_id = f"FUS{index}"
        bundle.events.append({
            "event_id": event_id,
            "event_type": "Fusion",
            "gene": f"GENE{index}::PARTNER{index}",
            "best_evidence_grade": "R3",
            "evidence_missing_layers": "presentation" if index >= 2 else "",
        })
        peptide = {
            "peptide_id": f"FP{index}",
            "event_id": event_id,
            "event_type": "Fusion",
            "gene": f"GENE{index}::PARTNER{index}",
            "peptide": f"FUSION{index}AA",
            "hla_allele": "HLA-A*02:01",
            "source_chain_confidence_tier": "C3",
        }
        if index < 2:
            peptide["netmhcpan_el_rank"] = "0.5"
        bundle.peptides.append(peptide)

    out = tmp_path / "patient_fusion_events.html"
    make_patient_report(out, bundle, event_top_n=10)
    text = out.read_text(encoding="utf-8")
    assert "Fusion Top 10" in text
    assert "GENE0::PARTNER0" in text
    assert "GENE9::PARTNER9" in text
    assert "FUSION9AA / HLA-A*02:01" in text


def test_technical_report_has_provenance_and_thresholds(tmp_path):
    out = tmp_path / "technical.html"
    make_technical_report(out, _bundle())
    text = out.read_text(encoding="utf-8")
    assert "Tool provenance" in text
    assert "netmhcpan.xls" in text
    assert "Profile Thresholds" in text
    assert "Field Glossary" in text
    assert "Ranked Peptides (full)" in text


def test_dual_reports_writes_three_files(tmp_path):
    paths = make_dual_reports(tmp_path, _bundle())
    assert Path(paths["evidence_report_patient"]).exists()
    assert Path(paths["evidence_report_technical"]).exists()
    assert Path(paths["evidence_report"]).exists()
    assert "Technical Report" in Path(paths["evidence_report_technical"]).read_text(encoding="utf-8")


def test_technical_report_loads_independent_wes_qc(tmp_path):
    qc_dir = tmp_path / "qc" / "wes"
    qc_dir.mkdir(parents=True)
    write_tsv(qc_dir / "wes_qc.tsv", [{
        "sample_id": "WES1", "qc_status": "PASS_WITH_CAPTURE_RATE_UNASSESSED",
        "total_reads": "1000", "primary_mapping_rate_pct": "99.9",
        "properly_paired_rate_pct": "98.0", "duplicate_rate_pct": "20.0",
        "target_definition": "GENCODE_CDS_PROXY_NOT_ASSAY_CAPTURE_BED",
        "mean_target_coverage": "80", "pct_target_bases_20x": "95",
        "pct_target_bases_30x": "90", "on_target_rate_pct": "70",
        "capture_rate_status": "UNASSESSED_CAPTURE_BED_MISSING",
        "formal_capture_rate_pct": "",
    }])
    base = _bundle()
    bundle = load_report_bundle(
        profile=base.profile, events=base.events, peptides=base.peptides,
        appm_summary=base.appm_summary, validation_rows=base.validation_rows,
        outdir=tmp_path, sample_id="S1",
    )
    out = tmp_path / "technical_wes.html"
    make_technical_report(out, bundle)
    text = out.read_text(encoding="utf-8")
    assert "Independent WES QC" in text
    assert "PASS_WITH_CAPTURE_RATE_UNASSESSED" in text
    assert "assay-specific capture BED" in text


def test_event_fields_are_enriched_from_representative_peptide():
    base = _bundle()
    base.events = [{
        "event_id": "B2M|chr15:1A>G", "event_type": "SNV", "gene": "B2M",
        "best_peptide_id": "P1",
    }]
    bundle = load_report_bundle(profile=base.profile, events=base.events, peptides=base.peptides)
    assert bundle.events[0]["combined_protein_change"] == "B2M:p.A1V"


def test_event_fields_are_enriched_from_raw_event_table(tmp_path):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    write_tsv(inputs / "combined_raw_events.tsv", [{
        "event_id": "B2M|chr15:1A>G", "gene": "B2M",
        "event_name": "ENSP000001:p.Ser10Phe", "peptide_consequence": "missense",
    }])
    base = _bundle()
    base.events = [{"event_id": "B2M|chr15:1A>G", "event_type": "SNV", "gene": "B2M", "best_peptide_id": "P1"}]
    bundle = load_report_bundle(profile=base.profile, events=base.events, peptides=base.peptides, outdir=tmp_path)
    assert bundle.events[0]["event_name"] == "ENSP000001:p.Ser10Phe"


def test_rna_fields_are_enriched_from_scoring_all_tool_results(tmp_path):
    source_dir = tmp_path / "scoring" / "evidence_consensus"
    source_dir.mkdir(parents=True)
    write_tsv(source_dir / "all_tool_results.tsv", [{
        "peptide_id": "P1", "rna_support_status": "RNA_ALT_SUPPORTED",
        "rna_alt_reads": "12", "rna_depth": "80", "rna_vaf": "0.15",
    }])
    base = _bundle()
    base.peptides[0]["rna_support_status"] = "UNASSESSED"
    bundle = load_report_bundle(profile=base.profile, events=base.events, peptides=base.peptides, outdir=tmp_path)
    assert bundle.peptides[0]["rna_support_status"] == "RNA_ALT_SUPPORTED"
    assert bundle.peptides[0]["rna_alt_reads"] == "12"


def test_all_tool_results_match_uses_peptide_and_hla(tmp_path):
    source_dir = tmp_path / "scoring" / "evidence_consensus"
    source_dir.mkdir(parents=True)
    write_tsv(source_dir / "all_tool_results.tsv", [
        {"peptide_id": "P1", "hla_allele": "HLA-A*01:01", "netmhcpan_el_rank": "4.0"},
        {"peptide_id": "P1", "hla_allele": "HLA-A*02:01", "netmhcpan_el_rank": "0.2"},
    ])
    base = _bundle()
    bundle = load_report_bundle(profile=base.profile, events=base.events, peptides=base.peptides, outdir=tmp_path)
    assert bundle.peptides[0]["netmhcpan_el_rank"] == "0.2"
    assert bundle.peptides[0]["_report_evidence_matched"] == "YES"
    assert bundle.evidence_source_status == "CANONICAL_ALL_TOOL_RESULTS"


def test_all_tool_results_manifest_sha256_is_verified(tmp_path):
    source_dir = tmp_path / "scoring" / "evidence_consensus"
    source_dir.mkdir(parents=True)
    evidence_path = source_dir / "all_tool_results.tsv"
    write_tsv(evidence_path, [{"peptide_id": "P1", "hla_allele": "HLA-A*02:01"}])
    digest = hashlib.sha256(evidence_path.read_bytes()).hexdigest()
    (source_dir / "all_tool_results.manifest.json").write_text(
        json.dumps({"output": {"sha256": digest}}), encoding="utf-8"
    )
    base = _bundle()
    bundle = load_report_bundle(profile=base.profile, events=base.events, peptides=base.peptides, outdir=tmp_path)
    assert bundle.evidence_integrity["status"] == "PASS"


def test_patient_evidence_summary_includes_all_major_dimensions():
    summary = _patient_evidence_summary({
        "event_authenticity_state": "SUPPORTED",
        "rna_support_state": "RNA_CONFIRMED",
        "presentation_consensus_state": "PRESENTATION_CONSISTENT_STRONG",
        "mutant_specificity_status": "MT_SPECIFIC",
        "clonality_state": "SUPPORTED",
        "hla_appm_state": "HLA_APPM_INTACT",
        "safety_state": "SAFETY_PARTIAL",
        "source_chain_confidence_tier": "C2",
    })
    for label in ("事件=", "RNA=", "呈递=", "MT/WT=", "克隆性=", "限制性HLA=", "APPM=", "安全性=", "来源链="):
        assert label in summary


def test_candidate_hla_and_appm_use_separate_sample_level_sources():
    bundle = _bundle()
    bundle.appm_summary = {
        "mhc_i_integrity_status": "MHC_I_INTACT",
        "appm_evidence_completeness": "PARTIAL",
        "tap_risk": "caution",
    }
    bundle.hla_loh_tool_results = [{
        "_report_tool": "SpecHLA", "hla_allele": "HLA-A*02:01", "loh_status": "no",
    }]
    summary = _patient_evidence_summary(bundle.peptides[0], bundle)
    limitation = _patient_limitation(bundle.peptides[0], bundle)
    assert "限制性HLA=单工具提示保留；多工具LOH确认未完成" in summary
    assert "APPM=证据部分完整；未见整体失活，但加工环节仍有谨慎信号" in summary
    assert "限制性HLA仅单工具提示保留；多工具LOH确认未完成" in limitation
    assert "APPM证据部分完整；加工环节仍需谨慎解释" in limitation
    assert "HLA_LOH_UNASSESSED" not in summary


def test_patient_report_does_not_present_mhc_i_intact_as_confirmed_function(tmp_path):
    bundle = _bundle()
    bundle.appm_summary = {
        "mhc_i_integrity_status": "MHC_I_INTACT",
        "appm_evidence_completeness": "PARTIAL",
        "ifng_response_status": "IFNG_RESPONSE_CAUTION",
        "tap_risk": "caution",
    }
    out = tmp_path / "patient_mhc_i_language.html"
    make_patient_report(out, bundle)
    text = out.read_text(encoding="utf-8")
    assert "MHC_I_INTACT" not in text
    assert "现有结果未发现HLA-I呈递系统整体完全丧失" in text
    assert "肿瘤可能仍保留一定呈递能力" in text
    assert "抗原加工环节和HLA-LOH证据尚不完整" in text
    assert "不表示HLA-I呈递功能已被实验确认完整" in text


def test_patient_ccf_does_not_render_supported_without_reliable_estimate():
    bundle = _bundle()
    bundle.purity_consensus = {"status": "LOW_PURITY_REVIEW"}
    row = {**bundle.peptides[0], "clonality_state": "SUPPORTED", "ccf_estimate": "UNASSESSED", "ccf_confidence_state": "UNSPECIFIED"}
    summary = _patient_evidence_summary(row, bundle)
    limitation = _patient_limitation(row, bundle)
    assert "克隆性=未形成可靠估计" in summary
    assert "样本纯度低且缺少可用CCF结果" in limitation
    assert "不作为阴性，也不作为正向加分" in limitation
    assert "克隆性=SUPPORTED" not in summary


def test_patient_candidate_tables_exclude_incomplete_and_do_not_advance_rows(tmp_path):
    bundle = _bundle()
    bundle.peptides.extend([
        {
            **bundle.peptides[0], "peptide_id": "P_NO_HLA", "event_id": "E_NO_HLA",
            "gene": "MISSINGHLA", "peptide": "AAAARLVD", "hla_allele": "",
        },
        {
            **bundle.peptides[0], "peptide_id": "P_STOP", "event_id": "E_STOP",
            "gene": "STOPGENE", "peptide": "STOPPEPTI", "recommended_use": "Do not advance",
        },
    ])
    out = tmp_path / "patient_filtered.html"
    make_patient_report(out, bundle)
    text = out.read_text(encoding="utf-8")
    assert "MISSINGHLA" not in text
    assert "AAAARLVD /" not in text
    assert "STOPGENE" not in text
    assert "Do not advance" not in text
    assert "科研技术版审阅池" in text


def test_patient_status_codes_are_translated_to_fixed_chinese():
    bundle = _bundle()
    row = {
        **bundle.peptides[0],
        "event_authenticity_state": "EVENT_PARTIAL",
        "rna_support_state": "RNA_CONFIRMED",
        "presentation_consensus_state": "PRESENTATION_CONSISTENT_STRONG",
        "safety_state": "SAFETY_PARTIAL",
        "source_chain_confidence_tier": "C3",
    }
    text = _patient_evidence_summary(row, bundle)
    for translated in (
        "事件获得部分支持", "RNA中检测到直接支持", "两个核心工具呈递预测一致且较强",
        "正常组织安全性证据不完整", "候选来源链基本合理，但关键环节尚未闭合",
    ):
        assert translated in text
    for code in ("EVENT_PARTIAL", "RNA_CONFIRMED", "PRESENTATION_CONSISTENT_STRONG", "SAFETY_PARTIAL", "C3"):
        assert code not in text


def test_tool_version_manifest_is_used_in_patient_tool_table(tmp_path):
    metadata = tmp_path / "metadata"
    metadata.mkdir()
    (metadata / "tool_versions.json").write_text(json.dumps({"tools": {
        "netmhcpan": {"version": "4.2c", "evidence": "运行版本清单"},
        "mhcflurry": {"version": "2.0.6", "evidence": "原始运行来源清单"},
    }}), encoding="utf-8")
    base = _bundle()
    bundle = load_report_bundle(profile=base.profile, events=base.events, peptides=base.peptides, outdir=tmp_path)
    out = tmp_path / "patient_versions.html"
    make_patient_report(out, bundle)
    text = out.read_text(encoding="utf-8")
    assert "4.2c" in text
    assert "2.0.6" in text
    assert "版本依据" in text


def test_patient_report_writes_machine_readable_release_audit(tmp_path):
    out = tmp_path / "patient.html"
    make_patient_report(out, _bundle())
    audit = out.with_suffix(".release_audit.json")
    assert audit.is_file()
    payload = json.loads(audit.read_text(encoding="utf-8"))
    assert len(payload["checks"]) == 10
    assert {row["编号"] for row in payload["checks"]} == {str(i) for i in range(1, 11)}


def test_patient_limitation_keeps_hard_fail_and_missing_dimensions():
    limitation = _patient_limitation({"hard_failure_codes": "FAIL_EVENT", "safety_state": "SAFETY_PARTIAL"})
    assert "FAIL_EVENT" in limitation
    assert "RNA证据未评估" in limitation
    assert "安全性证据不完整" in limitation
    assert "NetMHCstabpan未评估" in limitation


def test_patient_report_contains_evidence_source_audit(tmp_path):
    bundle = _bundle()
    bundle.evidence_source_status = "CANONICAL_ALL_TOOL_RESULTS"
    bundle.peptides[0]["_report_evidence_matched"] = "YES"
    bundle.peptides[0]["_report_evidence_source"] = "scoring/evidence_consensus/all_tool_results.tsv"
    bundle.evidence_manifest = {"validation": {"status": "PASS", "errors": []}}
    out = tmp_path / "patient_audit.html"
    make_patient_report(out, bundle)
    text = out.read_text(encoding="utf-8")
    assert "证据来源审计" in text
    assert "Top候选逐行匹配" in text
    assert "Top 100证据维度完整性" in text
    assert "NetMHCstabpan" in text


def test_patient_report_lists_input_files_and_purity_consensus(tmp_path):
    base = _bundle()
    provenance = {
        "input_files": {
            "somatic_vcf": "/patient/data/tumor.somatic.vcf.gz",
            "tumor_dna_bam": "/patient/data/tumor.bam",
            "normal_dna_bam": "/patient/data/normal.bam",
            "tumor_short_rna_fastq": ["/patient/data/rna_R1.fq.gz", "/patient/data/rna_R2.fq.gz"],
        },
        "pairing_status": "已使用肿瘤和配对正常样本；指纹未评估",
        "tumor_dna_depth": "候选位点中位深度 100x",
        "normal_dna_depth": "候选位点中位深度 80x",
        "rna_qc_status": "RNA覆盖已评估",
        "genome_build": "GRCh38",
        "purity_cnv_consensus": {
            "recommended_purity": "0.12", "recommended_ploidy": "2.20",
            "selected_tool": "Sequenza", "status": "LOW_PURITY_REVIEW",
            "basis": "Sequenza与PURPLE一致；FACETS结果偏高并保留冲突。",
        },
        "purity_cnv_tools": [
            {"tool": "Sequenza", "purity": "0.12", "ploidy": "2.20", "status": "LOW_PURITY_REVIEW", "note": "推荐值"},
            {"tool": "PURPLE", "purity": "0.13", "ploidy": "2.12", "status": "WARN_LOW_PURITY", "note": "支持低纯度"},
            {"tool": "FACETS", "purity": "0.28", "ploidy": "NA", "status": "PASS_WITH_REVIEW", "note": "与其他工具冲突"},
        ],
    }
    bundle = load_report_bundle(
        profile=base.profile, events=base.events, peptides=base.peptides,
        provenance=provenance, patient_inputs={"input_files": provenance["input_files"]},
    )
    out = tmp_path / "patient_inputs.html"
    make_patient_report(out, bundle)
    text = out.read_text(encoding="utf-8")
    for name in ("tumor.somatic.vcf.gz", "tumor.bam", "normal.bam", "rna_R1.fq.gz", "rna_R2.fq.gz"):
        assert name in text
    assert "/patient/data/" not in text
    assert "推荐纯度 0.12、倍性 2.20（Sequenza；LOW_PURITY_REVIEW）" in text
    assert "Sequenza与PURPLE一致；FACETS结果偏高并保留冲突。" in text
    assert "与其他工具冲突" in text


def test_explicit_patient_inputs_override_provenance_inputs():
    base = _bundle()
    bundle = load_report_bundle(
        profile=base.profile, events=base.events, peptides=base.peptides,
        provenance={"input_files": {"somatic_vcf": "/old/old.vcf.gz"}},
        patient_inputs={"somatic_vcf": "/new/new.vcf.gz"},
    )
    assert bundle.provenance["input_files"]["somatic_vcf"] == "/new/new.vcf.gz"


def test_patient_report_shows_hla_loh_tools_and_uses_consensus_in_appm(tmp_path):
    root = tmp_path / "result"
    loh_dir = root / "hla_loh_consensus"
    loh_dir.mkdir(parents=True)
    write_tsv(loh_dir / "spechla_hla_loh.tsv", [
        {"hla_allele": "HLA-A*02:01", "loh_status": "no", "evidence_tool": "spechla"},
    ])
    base = _bundle()
    bundle = load_report_bundle(
        profile=base.profile, events=base.events, peptides=base.peptides, outdir=root,
    )
    out = tmp_path / "patient_hla_loh.html"
    make_patient_report(out, bundle)
    text = out.read_text(encoding="utf-8")
    assert "HLA-I LOH多工具结果" in text
    assert "HLA-A*02:01" in text
    assert "<td>未提供</td><td>未见LOH（保留）</td>" in text
    assert "未见限制性HLA-I LOH（仅SpecHLA，证据有限）" in text
    assert "仅SpecHLA报告未见LOH（保留），证据有限" in text
    assert "限制性HLA=单工具提示保留；多工具LOH确认未完成" in text
    assert "HLA/APPM未评估" not in text


def test_patient_hla_loh_consensus_preserves_tool_conflict(tmp_path):
    root = tmp_path / "result"
    lohhla = root / "hla_loh" / "lohhla"
    spechla = root / "hla_loh" / "spechla"
    lohhla.mkdir(parents=True)
    spechla.mkdir(parents=True)
    write_tsv(lohhla / "hla_loh.tsv", [{"hla_allele": "HLA-A*02:01", "loh_status": "loh"}])
    write_tsv(spechla / "hla_loh.tsv", [{"hla_allele": "HLA-A*02:01", "loh_status": "no"}])
    base = _bundle()
    bundle = load_report_bundle(
        profile=base.profile, events=base.events, peptides=base.peptides, outdir=root,
    )
    out = tmp_path / "patient_hla_loh_conflict.html"
    make_patient_report(out, bundle)
    text = out.read_text(encoding="utf-8")
    assert "HLA-I LOH工具结果冲突：HLA-A*02:01" in text
    assert "工具结果冲突，暂不判定" in text


def test_event_change_uses_consequence_when_exact_protein_change_is_unavailable():
    assert _patient_event_change({"peptide_consequence": "frameshift"}) == "移码变异产生的新肽段"
    assert _patient_event_change({"event_type": "Splice", "peptide": "RSSTFPKWVTK"}) == "异常剪接肽段 RSSTFPKWVTK"
    assert _patient_event_change({}) == "蛋白改变待确认"


def test_ensembl_gene_ids_are_replaced_with_symbols():
    assert _replace_gene_ids("ENSG00000047648::ENSG00000101871", {
        "GENE|ENSG00000047648": "ARHGAP6",
        "GENE|ENSG00000101871": "MID1",
    }) == "ARHGAP6::MID1"
    assert _replace_gene_ids("CPSF6::ENSG00000170846", {}) == "CPSF6::MRFAP1L2"
