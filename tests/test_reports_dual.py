import hashlib
import gzip
import json
from pathlib import Path

from neoag.reports_dual import ReportBundle, _apply_patient_gene_expression, _augment_runtime_tool_provenance, _find_bam_input, _patient_analysis_context, _patient_conflict_summary, _patient_disease_background, _patient_dna_evidence, _patient_dna_rna_interpretation, _patient_event_grade_counts, _patient_event_representatives, _patient_evidence_audit_rows, _patient_evidence_summary, _patient_event_change, _patient_expression_tpm_map, _patient_fusion_artifact_review, _patient_fusion_boundary_evidence, _patient_hla_loh_consensus, _patient_key_gaps, _patient_limitation, _patient_manual_review_rows, _patient_metric, _patient_presentation_metric, _patient_presentation_quantitative_row, _patient_rna_measurements, _patient_rna_metric, _patient_safety_dimensions, _patient_safety_gap, _patient_tool_rows, _patient_track, _patient_validation, _replace_gene_ids, load_report_bundle, make_dual_reports, make_patient_report, make_technical_report
from neoag.reports_dual import _patient_ccf_coverage_rows, _patient_clonality_boundary, _patient_coding_variant_rna_rows, _patient_experiment_entry_gate_rows, _patient_splice_funnel_rows
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
    section6 = text.split("6. 人工复核候选的综合证据与实验建议（1个）", 1)[1].split("7. 分析方法与工具状态", 1)[0]
    assert "<th>候选</th>" in section6
    assert "<th>为什么值得关注</th>" in section6
    assert "<th>当前不确定性</th>" in section6
    assert "<th>建议下一步</th>" in section6
    assert "Safety-focused validation before efficacy assay" not in section6
    assert "<th>呈递工具</th>" not in section6
    assert "基因表达 12.3400 TPM" in text
    assert "转录本表达 3.5000 TPM" in text
    assert "RNA位点深度 20" in text
    assert "RNA alt reads 4" in text
    assert "RNA VAF 0.2000" in text
    assert "5. 当前进入人工复核的候选Peptide–HLA组合（去重后1个）" in text
    assert "当前展示1个去重候选组合" in text
    candidate_section = text.split("5. 当前进入人工复核的候选Peptide–HLA组合（去重后1个）", 1)[1].split("6. 人工复核候选的综合证据与实验建议（1个）", 1)[0]
    assert "<th>关键证据与下一步</th>" in candidate_section
    assert "<th>关键证据</th>" not in candidate_section
    assert "<th>主要限制</th>" not in candidate_section
    assert "<th>建议实验</th>" not in candidate_section
    assert "<th>审阅状态</th>" not in candidate_section
    assert "本表按事件级证据等级展示候选" in text
    assert "表内不再使用未细分的R3" in text
    assert "R3-READY表示候选基本合理" in text
    assert "R3-GAP表示关键资料缺失" in text
    assert "R3-REVIEW表示证据冲突或伪影风险需人工复核" in text
    assert "补证据后再决定" not in text
    assert "此类候选存在关键证据缺口" not in text
    assert "附录A：R1–R4证据分层" in text
    assert "附录B：术语说明" in text
    assert "附件与可追溯文件" in text
    assert "缺失证据统一视为未评估" in text
    assert "异常剪接逐级筛选漏斗" not in text
    assert "DQA1/DQB1" not in text
    assert "EWSR1::WT1" not in text


def test_splice_funnel_and_ccf_coverage_are_event_level_and_explicit():
    bundle = _bundle()
    bundle.events = [
        {"event_id": "S1", "event_type": "Splice", "mutation_source": "RNA_ONLY", "ccf_status": "RNA_ONLY_UNRESOLVED"},
        {"event_id": "V1", "event_type": "SNV", "ccf_status": "UNASSESSED"},
    ]
    bundle.peptides = [{
        "event_id": "S1", "event_type": "Splice", "peptide": "ABCDEFGHI",
        "rna_junction_reads": "12", "presentation_consensus_state": "PRESENTATION_CONSISTENT_STRONG",
    }]
    bundle.provenance["splice_filter_funnel_rows"] = [{
        "stage": "UNIQUE_JUNCTION_READS", "entered_events": "10", "assessed_events": "6",
        "passed_events": "4", "failed_events": "2", "unassessed_events": "4",
        "possible_remaining_range": "4-8", "criterion": "explicit unique reads >= 3",
    }]
    funnel = _patient_splice_funnel_rows(bundle)
    assert funnel[0]["筛选阶段"] == "unique junction reads门槛"
    assert funnel[0]["阶段后可能剩余"] == "4-8"
    assert funnel[-1]["筛选阶段"] == "探索性HLA预测（不计为正式通过）"
    assert funnel[-2]["筛选阶段"] == "正式前置门控后的HLA呈递支持"
    assert funnel[-2]["进入事件"] == "0"
    assert funnel[-1]["明确通过"] == "0"
    ccf = {row["事件类型"]: row for row in _patient_ccf_coverage_rows(bundle)}
    assert ccf["Splice"]["RNA-only不适用"] == "1"
    assert ccf["SNV"]["缺失/未解析"] == "1"


def test_stale_splice_funnel_is_hidden_without_splice_results():
    bundle = _bundle()
    bundle.provenance["splice_filter_funnel_rows"] = [{
        "stage": "UNIQUE_JUNCTION_READS",
        "entered_events": "10",
        "assessed_events": "10",
        "passed_events": "8",
    }]
    assert _patient_splice_funnel_rows(bundle) == []


def test_splice_rna_measurements_do_not_substitute_tpm_for_psi_or_unique_reads():
    text = _patient_rna_measurements({
        "event_type": "Splice", "gene_expression_tpm": "20", "rna_junction_reads": "12",
    })
    assert "unique junction reads未记录" in text
    assert "junction总覆盖未记录" in text
    assert "PSI未记录（不能用基因TPM替代）" in text


def test_splice_exact_support_reads_are_not_mislabeled_as_unique_reads():
    text = _patient_rna_measurements({
        "event_type": "Splice", "rna_junction_reads": "182",
        "junction_match_status": "EXACT", "junction_match_method": "bam_recount",
    })
    assert "精确junction支持reads 182" in text
    assert "unique junction reads 182" not in text
    assert "unique junction reads未记录" in text


def test_patient_top_candidates_include_non_r4_technical_review_rows(tmp_path):
    bundle = _bundle()
    bundle.events = [
        {
            "event_id": f"E{i}",
            "gene": f"G{i}",
            "event_type": "SNV",
            "best_evidence_grade": "R3",
            "manual_review_required": "yes",
        }
        for i in range(120)
    ]
    bundle.peptides = [
        {
            "peptide_id": f"P{i}",
            "event_id": f"E{i}",
            "gene": f"G{i}",
            "event_type": "SNV",
            # This test exercises Top-N capacity, so each candidate must have
            # a distinct peptide-HLA identity under the report contract.
            "peptide": f"AAAAAAA{i:03d}",
            "hla_allele": "HLA-A*02:01",
            "evidence_grade": "R3",
            "netmhcpan_el_rank": "0.5",
            # Missing source-chain tier intentionally routes this row to
            # TECHNICAL_REVIEW without making it R4.
            "source_chain_confidence_tier": "",
        }
        for i in range(120)
    ]
    bundle.validation_rows = []
    out = tmp_path / "patient_top100.html"
    make_patient_report(out, bundle, event_top_n=20, candidate_top_n=100)
    text = out.read_text(encoding="utf-8")
    section = text.split("<h3>当前展示100个去重候选组合</h3>", 1)[1].split("</table>", 1)[0]
    assert section.count("<tr>") - 1 == 100
    assert section.count("<td>R3-REVIEW</td>") == 100
    assert "G99" in section
    assert "G100" not in section

    interpretation = text.split("6. 人工复核候选的综合证据与实验建议（20个）", 1)[1].split("7. 分析方法与工具状态", 1)[0]
    assert interpretation.count("<tr>") - 1 == 20
    assert "G19 |" in interpretation
    assert "G20 |" not in interpretation


def test_splice_dna_evidence_does_not_render_placeholder_vcf_zero():
    row = {"event_type": "Splice", "tumor_vaf": "0", "tumor_alt_count": "0"}
    text = _patient_dna_evidence(row)
    assert "点突变VCF口径不适用" in text
    assert "VAF 0.0000" not in text


def test_splice_expression_uses_exact_patient_ensembl_id(tmp_path):
    expression = tmp_path / "gene_tpm.tsv"
    expression.write_text("gene_id\ttpm\nENSG00000139304.16\t18.85\n", encoding="utf-8")
    expression_map, source = _patient_expression_tpm_map({"expression": str(expression)})
    rows = [{
        "event_type": "Splice",
        "gene": "PTPRQ",
        "source_record_id": "ENSG00000139304:E45.1-E46.1",
        "gene_expression_tpm": "0",
    }]
    _apply_patient_gene_expression(rows, expression_map, source)
    assert rows[0]["gene_expression_tpm"] == "18.85"
    assert rows[0]["expression_evidence_status"] == "GENE_EXPRESSION_MATCHED_BY_ENSEMBL_ID"
    assert rows[0]["expression_source"] == str(expression)


def test_expression_overlay_does_not_collapse_multi_gene_fusion(tmp_path):
    expression = tmp_path / "gene_tpm.tsv"
    expression.write_text(
        "gene_id\ttpm\nENSG00000182944.17\t50\nENSG00000184937.13\t25\n",
        encoding="utf-8",
    )
    expression_map, source = _patient_expression_tpm_map({"gene_expression": str(expression)})
    rows = [{
        "event_type": "Fusion",
        "gene": "EWSR1::WT1",
        "source_record_id": "ENSG00000182944::ENSG00000184937",
        "gene_expression_tpm": "7",
    }]
    _apply_patient_gene_expression(rows, expression_map, source)
    assert rows[0]["gene_expression_tpm"] == "7"


def test_splice_ambiguous_gene_ids_do_not_render_placeholder_zero(tmp_path):
    expression = tmp_path / "gene_tpm.tsv"
    expression.write_text(
        "gene_id\ttpm\nENSG00000244731.4\t3.2\nENSG00000224389.9\t5.1\n",
        encoding="utf-8",
    )
    expression_map, source = _patient_expression_tpm_map({"expression": str(expression)})
    rows = [{
        "event_type": "Splice",
        "gene": "C4A / C4B",
        "source_records": "ENSG00000244731;ENSG00000224389",
        "gene_expression_tpm": "0",
    }]
    _apply_patient_gene_expression(rows, expression_map, source)
    assert rows[0]["gene_expression_tpm"] == ""
    assert rows[0]["expression_evidence_status"] == "UNASSESSED_AMBIGUOUS_GENE_ID"


def test_patient_report_track_uses_explicit_event_type_before_vep_consequence():
    snv = {
        "event_type": "SNV",
        "peptide_consequence": "splice_junction",
        "event_name": "ENSP000001:p.Glu284Gly",
    }
    indel = {
        "event_type": "InDel",
        "peptide_consequence": "splice_junction",
        "event_name": "ENSP000002:p.Gly12fs",
    }
    splice = {"event_type": "Splice", "peptide_consequence": "splice_junction", "peptide": "ABCDEFGHI"}
    assert _patient_track(snv) == "SNV"
    assert _patient_track(indel) == "InDel"
    assert _patient_track(splice) == "Splice"
    assert _patient_event_change(snv) == "ENSP000001:p.Glu284Gly"
    assert _patient_event_change(indel) == "ENSP000002:p.Gly12fs"
    assert _patient_event_change(splice) == "异常剪接肽段 ABCDEFGHI"


def test_patient_validation_translates_safety_recommendation():
    row = {
        "peptide_id": "P1",
        "event_type": "SNV",
        "recommended_use": "Safety-focused validation before efficacy assay",
    }
    text = _patient_validation(row, {})
    assert text == "先完成正常组织数据库复核及实验性脱靶/交叉反应验证，再考虑有效性实验"
    assert "Safety-focused" not in text


def test_patient_metric_explains_junction_novel_sequence_instead_of_unassessed():
    fusion = {
        "event_type": "Fusion",
        "peptide_consequence": "fusion",
        "crosses_junction": "yes",
        "mutant_specificity_state": "NOVEL_SEQUENCE",
    }
    assert _patient_metric("MT/WT", fusion, "mutant_specificity_status", "mutant_specificity_state") == (
        "MT/WT=异常连接新序列；应使用正常连接或正常异构体肽作为对照"
    )


def test_dsrct_defining_fusion_is_retained_once_in_manual_review_top_five():
    bundle = _bundle()
    events = []
    peptides = []
    for index in range(6):
        event_id = f"E{index}"
        events.append({
            "event_id": event_id, "gene": f"GENE{index}", "event_type": "SNV",
            "manual_review_required": "yes", "cancer_driver_context": "DRIVER_CONTEXT",
            "source_tools": "caller1;caller2", "evidence_conflict_fields": "presentation",
        })
        peptides.append({
            "event_id": event_id, "gene": f"GENE{index}", "event_type": "SNV",
            "peptide": "AAAAAAAAA", "hla_allele": "HLA-A*02:01",
        })
    for suffix in ("", "_ALT"):
        events.append({
            "event_id": f"FUSION_EWSR1_WT1{suffix}", "gene": "EWSR1::WT1", "event_type": "Fusion",
            "manual_review_required": "yes",
        })
    peptides.append({
        "event_id": "FUSION_EWSR1_WT1", "gene": "EWSR1::WT1", "event_type": "Fusion",
        "peptide": "SSYGQQSEK", "hla_allele": "HLA-A*03:01",
    })
    bundle.events = events
    bundle.peptides = peptides
    bundle.disease_knowledge = {
        "anchors": [{"event": "EWSR1::WT1", "label": "DSRCT disease anchor"}],
    }
    rows = _patient_manual_review_rows(events, peptides, bundle, {}, limit=5)
    assert sum(row["事件"] == "EWSR1::WT1" for row in rows) == 1


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
    assert "1个独立事件处于R3-REVIEW状态" in text


def test_patient_rna_interpretation_distinguishes_zero_low_and_adequate_coverage():
    base = {"event_type": "SNV", "tumor_alt_count": "20", "rna_alt_reads": "0"}
    no_coverage = _patient_dna_rna_interpretation({**base, "rna_depth": "0"})
    low_coverage = _patient_dna_rna_interpretation({**base, "rna_depth": "2"})
    adequate_negative = _patient_dna_rna_interpretation({**base, "rna_depth": "1471"})
    assert "无有效覆盖" in no_coverage and "不能作为阴性证据" in no_coverage
    assert "覆盖仅2" in low_coverage and "证据不足" in low_coverage
    assert "充分覆盖（1471" in adequate_negative
    assert "显著负证据并限制候选等级" in adequate_negative


def test_patient_home_conclusion_is_direct_and_derived_from_current_results(tmp_path):
    bundle = _bundle()
    bundle.events = [
        {"event_id": "F1", "gene": "GENEA::GENEB", "event_type": "Fusion", "best_evidence_grade": "R3", "manual_review_required": "yes"},
        {"event_id": "S1", "gene": "GENEC", "event_type": "SNV", "best_evidence_grade": "R3", "evidence_missing_layers": "rna"},
        {"event_id": "S2", "gene": "GENED", "event_type": "SNV", "best_evidence_grade": "R3", "evidence_missing_layers": "safety"},
    ]
    bundle.peptides = [
        {"peptide_id": "PF", "event_id": "F1", "gene": "GENEA::GENEB", "event_type": "Fusion", "peptide": "ABCDEFGHI", "hla_allele": "HLA-A*02:01"},
        {"peptide_id": "PS1", "event_id": "S1", "gene": "GENEC", "event_type": "SNV", "peptide": "BCDEFGHIK", "hla_allele": "HLA-A*02:01"},
        {"peptide_id": "PS2", "event_id": "S2", "gene": "GENED", "event_type": "SNV", "peptide": "CDEFGHIKL", "hla_allele": "HLA-B*07:02"},
    ]
    bundle.disease_knowledge = {"anchors": [{"event": "GENEA::GENEB", "label": "disease anchor"}]}
    out = tmp_path / "patient_direct_conclusion.html"
    make_patient_report(out, bundle, candidate_top_n=5)
    text = out.read_text(encoding="utf-8")
    assert "首页结论" in text
    assert "本次分析共发现3个独立候选事件，并产生3个肽段–HLA预测组合" in text
    assert "当前没有任何独立事件达到R1或R2计算证据等级" in text
    assert "1个独立事件处于R3-REVIEW状态" in text
    assert "为方便技术审阅，后文展示3个去重的Peptide–HLA组合" in text
    assert "这3个展示组合不是3个独立候选事件，也不是3个已经确认的新抗原" in text
    assert "<td>技术审阅展示组合</td><td>3</td>" in text
    assert "技术审阅清单覆盖疾病相关锚定事件GENEA::GENEB相关候选" in text
    assert "2个SNV候选" in text
    assert "尚未证明相关肽段在肿瘤细胞表面真实呈递或能够诱导T细胞反应" in text
    assert "本次重点审阅：" not in text


def test_patient_candidate_section_counts_final_peptide_hla_deduplication(tmp_path):
    bundle = _bundle()
    bundle.events = [
        {"event_id": "E1", "gene": "GENEA", "event_type": "SNV", "best_evidence_grade": "R3", "evidence_missing_layers": "rna"},
        {"event_id": "E2", "gene": "GENEA", "event_type": "SNV", "best_evidence_grade": "R3", "evidence_missing_layers": "safety"},
    ]
    common = {"gene": "GENEA", "event_type": "SNV", "peptide": "ABCDEFGHI", "hla_allele": "HLA-A*02:01"}
    bundle.peptides = [
        {**common, "peptide_id": "P1", "event_id": "E1"},
        {**common, "peptide_id": "P2", "event_id": "E2"},
    ]
    out = tmp_path / "patient_deduplicated_count.html"
    make_patient_report(out, bundle, candidate_top_n=100)
    text = out.read_text(encoding="utf-8")
    assert "候选Peptide–HLA组合（去重后1个）" in text
    assert "人工复核候选的综合证据与实验建议（1个）" in text
    assert "Top 100" not in text
    assert "前20项" not in text
    section = text.split("<h3>当前展示1个去重候选组合</h3>", 1)[1].split("</table>", 1)[0]
    assert section.count("<tr>") - 1 == 1


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


def test_patient_fusion_table_collapses_same_peptide_hla_across_hypotheses():
    events = [
        {
            "event_id": "FUSION_ENOPH1_EWSR1_chr4_82430884_chr22_29292137",
            "event_type": "Fusion", "gene": "ENOPH1::EWSR1",
            "best_peptide": "AEVTVILL", "best_hla_allele": "HLA-B*40:01",
        },
        {
            "event_id": "FUSION_ENOPH1_EWSR1_chr4_82430913_chr22_29292137",
            "event_type": "Fusion", "gene": "ENOPH1::EWSR1",
            "best_peptide": "AEVTVILL", "best_hla_allele": "HLA-B*40:01",
        },
    ]
    selected = _patient_event_representatives(events, [], limit=10, track="Fusion")
    assert len(selected) == 1
    assert selected[0]["patient_display_hypothesis_count"] == "2"
    assert "82430884" in selected[0]["patient_display_member_event_ids"]
    assert "82430913" in selected[0]["patient_display_member_event_ids"]


def test_patient_report_top_limits_are_configurable(tmp_path):
    out = tmp_path / "patient_custom_limits.html"
    make_patient_report(out, _bundle(), event_top_n=1, candidate_top_n=3)
    text = out.read_text(encoding="utf-8")
    assert "SNV Top 1" in text
    assert "5. 当前进入人工复核的候选Peptide–HLA组合（去重后1个）" in text


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


def test_patient_track_does_not_promote_vcf_consequence_to_splice():
    assert _patient_track({
        "event_type": "InDel",
        "mutation_source": "SNV_INDEL",
        "peptide_consequence": "splice_region_variant,frameshift_variant",
    }) == "InDel"
    assert _patient_track({
        "event_type": "SNV",
        "mutation_source": "VCF",
        "peptide_consequence": "splice_region_variant,missense_variant",
    }) == "SNV"


def test_patient_tool_coverage_discovers_event_only_splicemutr(tmp_path):
    root = tmp_path / "run"
    parsed = root / "parsed"
    parsed.mkdir(parents=True)
    write_tsv(parsed / "raw_events.tsv", [{
        "event_id": "SJ|GRCh38|chr1|101|200|+",
        "event_type": "Splice",
        "source_tools": "RegTools;SNAF;SpliceMutr",
    }])
    provenance = {}
    _augment_runtime_tool_provenance(root, provenance)
    assert provenance["tools"]["splicemutr"]["evidence_event_rows"] == 1

    bundle = ReportBundle(profile={}, events=[], peptides=[], provenance=provenance)
    rows = {row["流程/工具"]: row for row in _patient_tool_rows(bundle)}
    assert rows["SpliceMutr"]["状态"] == "事件来源记录已确认（1个事件）"
    assert rows["SpliceMutr"]["作用"] == "异常剪接交叉验证"


def test_patient_disease_background_prefers_structured_clinical_context():
    bundle = ReportBundle(
        profile={"_profile_name": "sarcoma_profile"}, events=[], peptides=[],
        provenance={
            "profile": "default",
            "clinical_context": {"diagnosis": "DSRCT"},
            "rules_name": "sarcoma_rules",
        },
    )
    result, basis = _patient_disease_background(bundle)
    assert result == "DSRCT"
    assert "clinical_context.diagnosis" in basis


def test_patient_disease_background_does_not_treat_profile_or_rules_as_diagnosis():
    profile_bundle = ReportBundle(
        profile={"_profile_name": "/profiles/sarcoma_rna_supported.toml"},
        events=[], peptides=[], provenance={"profile": "default", "rules_name": "fallback_rules"},
    )
    result, basis = _patient_disease_background(profile_bundle)
    assert result == "未提供"
    assert "分析配置和分子知识库关联不会代替临床诊断" in basis
    assert _patient_analysis_context(profile_bundle)[0] == "sarcoma_rna_supported"

    rules_bundle = ReportBundle(
        profile={"_profile_name": "default"}, events=[], peptides=[],
        provenance={"profile": "default", "rules": {"path": "/rules/sarcoma_consensus.toml"}},
    )
    assert _patient_disease_background(rules_bundle)[0] == "未提供"
    result, basis = _patient_analysis_context(rules_bundle)
    assert result == "sarcoma_consensus"
    assert "不代表临床诊断" in basis


def test_patient_disease_background_does_not_guess_from_paths():
    bundle = ReportBundle(
        profile={"_profile_name": "default"}, events=[], peptides=[],
        provenance={"input_files": {"tumor_bam": "/data/dsrct/patient/tumor.bam"}},
    )
    assert _patient_disease_background(bundle)[0] == "未提供"


def test_patient_manual_review_keeps_configured_key_fusion_with_junction_guidance():
    routine = [
        {
            "event_id": f"F{i}", "gene": f"GENE{i}::PARTNER{i}", "event_type": "Fusion",
            "manual_review_required": "yes", "source_chain_orthogonal_status": "SUPPORTED",
        }
        for i in range(8)
    ]
    key = {
        "event_id": "KEY", "gene": "KEY1::KEY2", "event_type": "Fusion",
        "manual_review_required": "yes", "best_evidence_grade": "R4",
    }
    bundle = ReportBundle(
        profile={"manual_review": {"events": ["KEY1::KEY2"]}},
        events=routine + [key], peptides=[],
    )
    rows = _patient_manual_review_rows(bundle.events, [], bundle, {}, limit=5)
    assert rows[0]["事件"] == "KEY1::KEY2"
    assert "疾病/证据规则明确指定" in rows[0]["为什么重要"]
    assert "精确跨断点" in rows[0]["当前建议"]
    assert "普通融合伙伴蛋白" in rows[0]["当前建议"]


def test_patient_conflicts_ignore_none_status_and_provenance_only_fields():
    row = {
        "evidence_conflict_status": "NONE",
        "evidence_conflict_fields": "source_file,source_row_number,source_tools,safety_status",
    }
    assert _patient_conflict_summary(row) == ""


def test_patient_presentation_disagreement_names_tools_and_values():
    text = _patient_presentation_metric({
        "presentation_consensus_state": "PRESENTATION_DISCORDANT",
        "netmhcpan_el_rank": "0.261",
        "mhcflurry_presentation_score": "0.5847",
    })
    assert "NetMHCpan EL rank=0.261" in text
    assert "MHCflurry呈递分=0.5847" in text


def test_patient_rna_only_track_does_not_report_ccf_as_universal_gap():
    bundle = _bundle()
    row = {
        "event_type": "Fusion", "event_id": "F1", "peptide": "AAAAAAAAA",
        "hla_allele": "HLA-A*02:01", "source_chain_confidence_tier": "C2",
        "netmhcpan_el_rank": "0.5", "ccf_status": "RNA_ONLY_UNRESOLVED",
        "safety_status": "PASS",
    }
    assert "CCF未形成可靠估计" not in _patient_key_gaps(row, bundle)


def test_patient_splice_reports_provided_reads_separately_from_verified_reads():
    row = {
        "event_type": "Splice",
        "provided_rna_junction_reads": "218",
        "rna_junction_reads": "0",
        "rna_support_status": "UNASSESSED",
    }
    assert "上游工具报告junction reads 218" in _patient_rna_metric(row)
    assert "已核实reads为0" in _patient_rna_metric(row)
    measurements = _patient_rna_measurements(row)
    assert "junction reads 0（核实状态未确认）" in measurements
    assert "上游工具汇总junction reads 218" in measurements
    assert "差额 218 条尚未归属" in measurements


def test_patient_splice_reports_only_the_unresolved_read_difference():
    measurements = _patient_rna_measurements({
        "event_type": "Splice",
        "canonical_junction_id": "SJ|GRCh38|chr2|232337162|232344092|-",
        "junction_match_status": "EXACT",
        "junction_support_status": "SUPPORTED_EXACT_JUNCTION",
        "provided_rna_junction_reads": "126",
        "rna_junction_reads": "108",
    })
    assert "caller原始记录 126" in measurements
    assert "主比对表unique reads 108" in measurements
    assert "两者均已坐标回链" in measurements
    assert "差异不代表未归属reads" in measurements
    assert "差额 18 条尚未归属" not in measurements


def test_patient_fusion_reports_source_record_backlink_separately_from_alignment_verification():
    measurements = _patient_rna_measurements({
        "event_type": "Fusion",
        "rna_junction_reads": "11",
        "source_tool": "Arriba",
        "source_file": "/analysis/arriba/fusions.tsv",
        "source_record_id": "FUSION_A_B_chr1_100_chr2_200",
    })
    assert "融合caller原始记录已回链，junction reads 11" in measurements
    assert "尚无独立主比对表核实" in measurements
    assert "尚未独立回链核实" not in measurements


def test_patient_splice_labels_exact_cross_tool_junction_link():
    measurements = _patient_rna_measurements({
        "event_type": "Splice",
        "rna_junction_reads": "41",
        "rna_junction_source": "SNAF+SpliceMutr",
        "strict_cross_validated": "yes",
        "splicemutr_structure_exact": "yes",
        "junction_key": "chr19:8056237:8056366:+",
    })
    assert "已按event/junction精确回链的剪接工具支持reads 41" in measurements
    assert "尚未独立回链核实" not in measurements
    assert "126（尚未精确回链）" not in measurements


def test_patient_splice_recognizes_v051_exact_canonical_statuses():
    measurements = _patient_rna_measurements({
        "event_type": "Splice",
        "canonical_junction_id": "SJ|GRCh38|chr19|56177196|56182171|+",
        "junction_match_status": "MATCHED_EXACT_CANONICAL",
        "junction_resolution_status": "RESOLVED_EXACT_CANONICAL",
        "junction_support_status": "SUPPORTED_ALL_EVENT_JUNCTIONS_EXACT",
        "rna_junction_reads": "2589",
        "provided_rna_junction_reads": "2357",
    })
    assert "已按标准化精确junction坐标回链，junction reads 2589" in measurements
    assert "尚无独立主比对表核实" not in measurements


def test_patient_safety_gap_distinguishes_gene_absent_from_reference():
    gap = _patient_safety_gap({
        "safety_status": "SAFETY_PARTIAL",
        "safety_missing_layers": "normal_expression;normal_hspc",
        "safety_reason": (
            "normal_expression_gene_not_in_reference;"
            "normal_hspc_gene_not_in_reference;safety_evidence_incomplete"
        ),
    })
    assert "GTEx/HPA正常组织表达参考未收录该基因" in gap
    assert "HSPC正常造血参考未收录该基因" in gap
    assert "不能按0表达解释" in gap


def test_patient_source_chain_c4_lists_specific_traceability_gaps():
    bundle = _bundle()
    row = {
        "event_type": "Splice", "event_id": "SJ1", "peptide": "AAAAAAAAA",
        "hla_allele": "HLA-A*02:01", "source_chain_confidence_tier": "C4",
        "canonical_junction_id": "SJ|GRCh38|chr1|100|200|.",
        "junction_strand": ".", "provided_rna_junction_reads": "218",
        "rna_junction_reads": "0", "netmhcpan_el_rank": "0.5",
        "source_record_id": "JUNC_TEST_1",
        "source_chain_reason_codes": "SC_PEPTIDE_HLA_TRACEABILITY_INCOMPLETE",
        "safety_status": "PASS",
    }
    gaps = _patient_key_gaps(row, bundle)
    source_gap = next(item for item in gaps if item.startswith("完整性缺口：来源链C4"))
    assert "canonical junction缺少可用strand" in source_gap
    assert "上游同坐标报告218条reads；因strand未解析，未计入严格verified reads" in source_gap
    assert "上游caller事件已回溯，但transcript hypothesis尚未建立" in source_gap
    assert "正式ORF尚未建立" in source_gap
    assert "肽段-HLA可回溯至上游caller事件" in source_gap


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
    for label in (
        "事件=", "RNA=", "呈递=", "MT/WT=", "限制性HLA=", "APPM=",
        "自身相似性/正常组织风险=", "来源链=",
    ):
        assert label in summary
    assert "克隆性=" not in summary
    assert "CCF=" not in summary


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
    assert "限制性HLA=单工具未提示LOH；不足以确认完整保留" in summary
    assert "APPM=证据部分完整；未见整体失活，但加工环节仍有谨慎信号" in summary
    assert "限制性HLA仅单工具未提示LOH；不足以确认肿瘤中完整保留" in limitation
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


def test_patient_candidate_text_hides_ccf_details_regardless_of_confidence():
    bundle = _bundle()
    bundle.purity_consensus = {"status": "LOW_PURITY_REVIEW"}
    row = {**bundle.peptides[0], "clonality_state": "SUPPORTED", "ccf_estimate": "UNASSESSED", "ccf_confidence_state": "UNSPECIFIED"}
    summary = _patient_evidence_summary(row, bundle)
    limitation = _patient_limitation(row, bundle)
    assert "克隆性" not in summary
    assert "CCF" not in summary
    assert "CCF" not in limitation


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
        "自身相似性与正常组织风险筛查仅部分完成（具体缺口见候选说明）", "候选来源链基本合理，但关键环节尚未闭合",
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


def test_patient_tool_table_overrides_stale_deepimmuno_not_used_status(tmp_path):
    bundle = _bundle()
    bundle.peptides[0]["deepimmuno_score"] = "0.82"
    bundle.provenance["tools"]["deepimmuno"] = {
        "status": "not_used",
        "version": "deepimmuno-cnn.py",
        "mode": "derived",
    }
    out = tmp_path / "patient_deepimmuno.html"
    make_patient_report(out, bundle)
    text = out.read_text(encoding="utf-8")
    assert "综合证据表已载入结果值（1/1 个适用组合）" in text
    assert "结果值优先于可能过期的运行清单状态" in text
    deep_row = text.split("<td>DeepImmuno</td>", 1)[1].split("</tr>", 1)[0]
    assert "not_used" not in deep_row


def test_patient_tool_table_recognizes_appm_candidate_modifiers(tmp_path):
    bundle = _bundle()
    bundle.peptides[0]["appm_multiplier"] = "0.4500"
    out = tmp_path / "patient_appm.html"
    make_patient_report(out, bundle)
    text = out.read_text(encoding="utf-8")
    appm_row = text.split("<td>TAP/APPM</td>", 1)[1].split("</tr>", 1)[0]
    assert "NeoAg APPM 2.0" in appm_row
    assert "APPM已评估（候选修饰值 1/1）" in appm_row
    assert "独立TAP候选级状态未单列" in appm_row


def test_patient_splice_missing_transcript_tpm_explains_mapping_boundary():
    text = _patient_rna_measurements({
        "event_type": "Splice",
        "gene": "GENE1",
        "gene_expression_tpm": "12.5",
        "transcript_expression_tpm": "",
    })
    assert "基因表达 12.5000 TPM" in text
    assert "异常剪接异构体无法唯一对应定量矩阵中的标准转录本" in text
    assert "表达证据未接入" not in text


def test_patient_event_top_moves_track_common_not_applicable_evidence_to_intro(tmp_path):
    bundle = _bundle()
    bundle.events = [{
        "event_id": "SJ|GRCh38|chr1|101|200|+",
        "event_type": "Splice",
        "gene": "GENE1",
        "evidence_grade": "R3-REVIEW",
    }]
    bundle.peptides = [{
        "peptide_id": "SP1",
        "event_id": "SJ|GRCh38|chr1|101|200|+",
        "event_type": "Splice",
        "gene": "GENE1",
        "peptide": "AAAAAAAAA",
        "hla_allele": "HLA-A*02:01",
        "crosses_junction": "true",
        "rna_junction_reads": "12",
        "event_authenticity_state": "EVENT_CONFIRMED",
        "rna_support_state": "RNA_CONFIRMED",
        "presentation_consensus_state": "PRESENTATION_MODERATE",
        "mutant_specificity_state": "UNASSESSED",
        "evidence_grade": "R3-REVIEW",
    }]
    out = tmp_path / "patient_track_intro.html"
    make_patient_report(out, bundle, event_top_n=1)
    text = out.read_text(encoding="utf-8")
    section = text.split("<h3>Splice Top 1</h3>", 1)[1].split("</table>", 1)[0]
    intro, table = section.split("<table>", 1)
    assert "本赛道证据口径" in intro
    assert "普通点突变式DNA VAF" in intro
    assert "传统MT/WT配对" in intro
    assert "不再逐行重复不适用项" in intro
    assert "MT/WT=传统点突变式配对不适用" not in table
    assert "DNA证据：点突变VCF口径不适用" not in table


def test_patient_tool_table_discovers_source_tools_not_listed_in_provenance(tmp_path):
    bundle = _bundle()
    bundle.peptides[0]["source_tools"] = "EasyFuse;JAFFAL;SNAF;SpliceMutr"
    out = tmp_path / "patient_source_tools.html"
    make_patient_report(out, bundle)
    text = out.read_text(encoding="utf-8")
    for tool in ("EasyFuse", "JAFFAL", "SNAF", "SpliceMutr"):
        assert f"<td>{tool}</td>" in text
    assert "综合证据表已载入结果值（1/1 个适用组合）" in text


def test_patient_tool_coverage_uses_applicable_denominator_and_reports_nonapplicable():
    bundle = _bundle()
    bundle.peptides.extend([
        {"peptide_id": "LONG", "peptide": "AAAAAAAAAAAAAAA", "hla_allele": "HLA-A*02:01"},
        {"peptide_id": "NOHLA", "peptide": "AAAAAAAAA", "hla_allele": ""},
    ])
    rows = {row["流程/工具"]: row for row in _patient_tool_rows(bundle)}
    assert "1/1 个适用组合" in rows["NetMHCpan"]["状态"]
    assert "另 2 个" in rows["NetMHCpan"]["状态"]


def test_patient_tool_table_discovers_standard_runtime_outputs(tmp_path):
    (tmp_path / "rna" / "star").mkdir(parents=True)
    (tmp_path / "rna" / "star" / "Log.final.out").write_text("finished\n", encoding="utf-8")
    (tmp_path / "hla_loh" / "spechla").mkdir(parents=True)
    (tmp_path / "hla_loh" / "spechla" / "hla_loh.tsv").write_text(
        "hla_allele\tloh_status\nHLA-A*02:01\tno\n", encoding="utf-8"
    )
    base = _bundle()
    bundle = load_report_bundle(
        profile=base.profile, events=base.events, peptides=base.peptides, outdir=tmp_path,
    )
    out = tmp_path / "patient_runtime_tools.html"
    make_patient_report(out, bundle)
    text = out.read_text(encoding="utf-8")
    assert "<td>STAR</td>" in text
    assert "<td>SpecHLA</td>" in text
    assert "短读长RNA比对与junction提取" in text
    assert "HLA分型与HLA-LOH证据" in text


def test_patient_report_writes_machine_readable_release_audit(tmp_path):
    out = tmp_path / "patient.html"
    make_patient_report(out, _bundle())
    audit = out.with_suffix(".release_audit.json")
    assert audit.is_file()
    payload = json.loads(audit.read_text(encoding="utf-8"))
    assert len(payload["checks"]) == 10
    assert {row["编号"] for row in payload["checks"]} == {str(i) for i in range(1, 11)}


def test_patient_limitation_keeps_hard_fail_and_missing_dimensions():
    limitation = _patient_limitation({
        "hard_failure_codes": "FAIL_EVENT",
        "safety_state": "SAFETY_PARTIAL",
        "safety_missing_layers": "normal_hspc;normal_junction",
    })
    assert "FAIL_EVENT" in limitation
    assert "RNA证据未评估" in limitation
    assert "HSPC正常造血参考已接入但该候选未匹配到可判定记录" in limitation
    assert "正常融合/剪接连接背景未完成正式评估" in limitation
    assert "安全性证据不完整" not in limitation
    assert "NetMHCstabpan未评估" in limitation


def test_patient_key_gaps_describes_safety_review_signals():
    row = {
        "event_type": "Fusion",
        "safety_state": "SAFETY_REVIEW",
        "event_safety_reason": "critical_tissue_expression;normal_HSPC_expression",
        "critical_tissue_name": "Brain_Cerebellar_Hemisphere",
        "critical_tissue_max_tpm": "201.559",
        "normal_hspc_tpm": "160.7",
        "normal_hspc_unit": "HPA_nCPM",
    }
    gaps = _patient_key_gaps(row, _bundle())
    text = "；".join(gaps)
    assert "连接事件的自身相似性与正常组织风险需按精确新生序列复核" in text
    assert "伙伴基因正常组织/HSPC表达仅作背景" in text
    assert "201.5590 TPM" not in text
    assert "160.7000 HPA_nCPM" not in text
    assert "安全性证据不完整" not in text


def test_patient_limitation_describes_patient_specific_tool_conflict():
    limitation = _patient_limitation({
        "evidence_conflict_fields": "presentation_evidence_score,rna_alt_reads",
        "evidence_conflict_details": json.dumps([
            {
                "field": "presentation_evidence_score",
                "selected_source": "presentation_evidence",
                "selected_value": "0.82",
                "other_source": "ranked_peptides",
                "other_value": "0.31",
            },
            {
                "field": "rna_alt_reads",
                "selected_source": "rna_junction_evidence",
                "selected_value": "9",
                "other_source": "raw_events",
                "other_value": "0",
            },
        ]),
    })
    assert "具体证据冲突" in limitation
    assert "旧主排序副本" not in limitation
    assert "RNA突变等位基因reads：RNA位点/连接证据=9，原始事件表=0" in limitation


def test_patient_limitation_does_not_treat_stale_ranking_copy_as_tool_conflict():
    limitation = _patient_limitation({
        "evidence_conflict_fields": "presentation_evidence_score",
        "evidence_conflict_details": json.dumps([{
            "field": "presentation_evidence_score",
            "selected_source": "presentation_evidence",
            "selected_value": "0.82",
            "other_source": "ranked_peptides",
            "other_value": "0.31",
        }]),
    })
    assert "具体证据冲突" not in limitation
    assert "旧主排序副本" not in limitation


def test_patient_splice_change_uses_exon_path_and_coordinate():
    change = _patient_event_change({
        "event_type": "Splice",
        "peptide_consequence": "splice_junction",
        "source_event_id": "ENSG00000008300:E25.1_48645470-E30.3_48642042",
        "canonical_junction_id": "SJ|GRCh38|chr3|48642043|48645469|.",
        "event_name": "chr3:48642043",
    })
    assert "E25.1→E30.3" in change
    assert "chr3:48642043-48645469" in change
    assert "ORF/蛋白影响待确认" in change


def test_patient_splice_change_reports_formal_origin_chain():
    change = _patient_event_change({
        "event_type": "Splice",
        "canonical_junction_id": "SJ|GRCh38|chr1|101|200|+",
        "transcript_hypothesis_id": "STH|abc",
        "orf_id": "ORF|abc",
        "origin_peptide_id": "POR|abc",
        "peptide": "ABCDEFGHI",
    })
    assert "已完成局部转录本、ORF及跨junction肽段来源精确回链" in change
    assert "全长转录本真实性仍待独立验证" in change


def test_patient_key_gaps_translates_hard_failure_reason():
    bundle = ReportBundle(profile={}, events=[], peptides=[], appm_summary={}, validation_rows=[])
    gaps = _patient_key_gaps(
        {"hard_failure_codes": "HARD_REFERENCE_PROTEOME_MATCH", "event_type": "Splice"},
        bundle,
    )
    assert "阻断原因：候选肽与正常参考蛋白组存在精确匹配" in gaps


def test_new_patient_splice_gene_and_change_are_enriched_from_gtf(tmp_path, monkeypatch):
    gtf = tmp_path / "gencode.gtf"
    gtf.write_text(
        'chr7\ttest\tgene\t100\t500\t.\t+\t.\tgene_id "ENSG00999999999"; gene_name "GENE_NEW";\n',
        encoding="utf-8",
    )
    monkeypatch.setenv("GENCODE_GTF", str(gtf))
    event_id = "SJ|GRCh38|chr7|200|300|+"
    source_event_id = "ENSG00999999999:E2.1_199-E4.1_301"
    bundle = load_report_bundle(
        profile={},
        events=[{
            "event_id": event_id,
            "event_type": "Splice",
            "canonical_junction_id": event_id,
            "source_event_id": source_event_id,
        }],
        peptides=[{
            "peptide_id": "P_NEW",
            "event_id": event_id,
            "event_type": "Splice",
            "canonical_junction_id": event_id,
            "source_event_id": source_event_id,
            "peptide": "AAAAAAAAA",
            "hla_allele": "HLA-A*02:01",
        }],
        outdir=tmp_path,
        sample_id="NEW_PATIENT",
    )
    assert bundle.events[0]["gene"] == "GENE_NEW"
    assert bundle.peptides[0]["gene"] == "GENE_NEW"
    assert "E2.1→E4.1" in _patient_event_change(bundle.peptides[0])


def test_patient_splice_gene_symbol_is_discovered_from_recorded_asset_root(tmp_path):
    data_root = tmp_path / "portable_assets" / "data"
    gtf = data_root / "rna" / "gencode_v49" / "gencode.v49.annotation.gtf.gz"
    gtf.parent.mkdir(parents=True)
    with gzip.open(gtf, "wt", encoding="utf-8") as handle:
        handle.write(
            'chr7\ttest\tgene\t100\t500\t.\t+\t.\tgene_id "ENSG00999999998.1"; '
            'gene_type "protein_coding"; gene_name "PORTABLE_GENE";\n'
        )
    proteome = data_root / "normal" / "proteome" / "proteome.fa"
    proteome.parent.mkdir(parents=True)
    proteome.write_text(">P\nAAAA\n", encoding="utf-8")
    event_id = "SJ|GRCh38|chr7|200|300|+"
    bundle = load_report_bundle(
        profile={},
        events=[{
            "event_id": event_id,
            "event_type": "Splice",
            "gene": "ENSG00999999998",
            "canonical_junction_id": event_id,
        }],
        peptides=[{
            "peptide_id": "P_PORTABLE",
            "event_id": event_id,
            "event_type": "Splice",
            "gene": "ENSG00999999998",
            "canonical_junction_id": event_id,
            "peptide": "AAAAAAAAA",
            "hla_allele": "HLA-A*02:01",
        }],
        outdir=tmp_path,
        provenance={"input_files": {"reference_proteome": str(proteome)}},
        sample_id="PORTABLE_PATIENT",
    )
    assert bundle.events[0]["gene"] == "PORTABLE_GENE"
    assert bundle.peptides[0]["gene"] == "PORTABLE_GENE"


def test_consensus_report_streams_source_labels_for_unresolved_new_patient_gene(tmp_path, monkeypatch):
    gtf = tmp_path / "gencode.gtf"
    gtf.write_text(
        'chr7\ttest\tgene\t100\t500\t.\t+\t.\tgene_id "ENSG00999999999"; gene_name "GENE_NEW";\n',
        encoding="utf-8",
    )
    monkeypatch.setenv("GENCODE_GTF", str(gtf))
    event_id = "SJ|GRCh38|chr7|200|300|+"
    source_event_id = "ENSG00999999999:E2.1_199-E4.1_301"
    scoring = tmp_path / "scoring"
    scoring.mkdir()
    write_tsv(scoring / "all_tool_results.tsv", [{
        "event_id": event_id,
        "peptide_id": "P_NEW",
        "hla_allele": "HLA-A*02:01",
        "gene": "chr7:200",
        "canonical_junction_id": event_id,
        "source_record_id": source_event_id,
    }])
    bundle = load_report_bundle(
        profile={},
        events=[{
            "event_id": event_id,
            "event_type": "Splice",
            "gene": "chr7:200",
            "evidence_grade": "R3",
            "representative_1_peptide_id": "P_NEW",
        }],
        peptides=[{
            "peptide_id": "P_NEW",
            "event_id": event_id,
            "event_type": "Splice",
            "gene": "chr7:200",
            "peptide": "AAAAAAAAA",
            "hla_allele": "HLA-A*02:01",
            "evidence_grade": "R3",
        }],
        outdir=tmp_path,
        sample_id="NEW_PATIENT",
    )
    assert bundle.events[0]["gene"] == "GENE_NEW"
    assert bundle.peptides[0]["gene"] == "GENE_NEW"
    assert "E2.1→E4.1" in _patient_event_change(bundle.peptides[0])


def test_consensus_report_backfills_specific_protein_change_for_new_patient(tmp_path):
    scoring = tmp_path / "scoring"
    scoring.mkdir()
    write_tsv(scoring / "all_tool_results.tsv", [{
        "event_id": "chr1:100:A:T",
        "peptide_id": "P_SNV",
        "hla_allele": "HLA-A*02:01",
        "gene": "GENE1",
        "combined_protein_change": "ENST0001:p.Gly12Asp",
    }])
    bundle = load_report_bundle(
        profile={},
        events=[{
            "event_id": "chr1:100:A:T",
            "event_type": "SNV",
            "gene": "GENE1",
            "evidence_grade": "R3",
            "representative_1_peptide_id": "P_SNV",
        }],
        peptides=[{
            "peptide_id": "P_SNV",
            "event_id": "chr1:100:A:T",
            "event_type": "SNV",
            "gene": "GENE1",
            "peptide": "AAAAAAAAA",
            "hla_allele": "HLA-A*02:01",
            "evidence_grade": "R3",
        }],
        outdir=tmp_path,
        sample_id="NEW_PATIENT",
    )
    assert _patient_event_change(bundle.events[0]) == "p.Gly12Asp"
    assert _patient_event_change(bundle.peptides[0]) == "p.Gly12Asp"


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
    assert "Top 1证据维度可用性" in text
    assert "可作为当前分层证据" in text
    assert "尚不能作为可靠证据" in text
    assert "NetMHCstabpan" in text


def test_patient_evidence_audit_hides_candidate_level_ccf():
    rows = [
        {"ccf_estimate": "0.82", "ccf_confidence_state": "CCF_HIGH_CONFIDENCE"},
        {"ccf_estimate": "0.21", "ccf_confidence_state": "CCF_LOW_CONFIDENCE"},
        {"ccf_estimate": "", "ccf_status": "RNA_ONLY_UNRESOLVED"},
    ]
    assert all(row["证据维度"] != "克隆性/CCF" for row in _patient_evidence_audit_rows(rows))


def test_patient_report_uses_single_sample_level_clonality_boundary(tmp_path):
    bundle = _bundle()
    bundle.events = [
        {"event_id": "E1", "event_type": "SNV", "ccf_estimate": "0.82", "ccf_confidence_state": "CCF_HIGH_CONFIDENCE"},
        {"event_id": "E2", "event_type": "SNV", "ccf_estimate": "0.21", "ccf_confidence_state": "CCF_LOW_CONFIDENCE"},
    ]
    bundle.peptides[0].update({
        "ccf_estimate": "0.82", "ccf_ci_low": "0.61", "ccf_ci_high": "0.97",
        "ccf_confidence_state": "CCF_HIGH_CONFIDENCE", "total_cn": "3",
    })
    out = tmp_path / "patient_clonality_boundary.html"
    make_patient_report(out, bundle)
    text = out.read_text(encoding="utf-8")
    assert text.count("<b>克隆性证据边界：</b>") == 1
    assert "当前缺口是事件级CCF覆盖不足，而不是缺少样本纯度估计" not in text
    assert "CCF覆盖度与缺失影响" not in text
    assert "CCF=0.82" not in text
    assert "95%区间 0.61-0.97" not in text


def test_patient_clonality_boundary_separates_low_purity_from_missing_ccf():
    bundle = _bundle()
    bundle.purity_consensus = {
        "recommended_purity": "0.12", "recommended_ploidy": "2.1",
        "status": "LOW_PURITY_REVIEW",
    }
    bundle.events = [{"event_id": "E1", "event_type": "SNV", "ccf_status": "UNASSESSED"}]
    text = _patient_clonality_boundary(bundle)
    assert "低纯度审阅" in text
    assert "系统性降低CCF置信度" in text
    assert "不解释为候选不存在" in text


def test_patient_clonality_boundary_high_purity_names_event_level_gap():
    bundle = _bundle()
    bundle.purity_consensus = {
        "recommended_purity": "0.92", "recommended_ploidy": "2.0",
        "status": "CONCORDANT",
    }
    bundle.purity_tools = [
        {"tool": "FACETS", "purity": "0.91", "status": "FOUND"},
        {"tool": "PURPLE", "purity": "0.93", "status": "FOUND"},
        {"tool": "Sequenza", "purity": "", "status": "MISSING"},
        {"tool": "ASCAT", "purity": "", "status": "MISSING"},
    ]
    bundle.events = [{"event_id": "E1", "event_type": "SNV", "ccf_status": "UNASSESSED"}]
    text = _patient_clonality_boundary(bundle)
    assert "多工具一致或相互支持" in text
    assert "事件级CCF覆盖不足，而不是缺少样本纯度估计" in text


def test_patient_clonality_boundary_preserves_outlier_context():
    bundle = _bundle()
    bundle.purity_consensus = {
        "recommended_purity": "0.48", "recommended_ploidy": "2.4",
        "status": "MULTI_TOOL_REVIEW",
    }
    bundle.purity_tools = [
        {"tool": "FACETS", "purity": "0.47", "status": "FOUND"},
        {"tool": "Sequenza", "purity": "0.49", "status": "FOUND"},
        {"tool": "PURPLE", "purity": "0.81", "status": "OUTLIER_EXCLUDED", "note": "excluded outlier"},
    ]
    bundle.events = [{"event_id": "E1", "event_type": "SNV", "ccf_status": "UNASSESSED"}]
    text = _patient_clonality_boundary(bundle)
    assert "PURPLE" in text
    assert "冲突、离群或被排除" in text
    assert "并列保留全部工具结果" in text


def test_coding_variant_rna_coverage_is_event_level_and_excludes_junction_tracks():
    bundle = _bundle()
    bundle.events = [
        {"event_id": "S1", "event_type": "SNV", "rna_depth": "0", "rna_alt_reads": "0"},
        {"event_id": "S2", "event_type": "SNV", "rna_depth": "5", "rna_alt_reads": "0"},
        {"event_id": "S3", "event_type": "SNV", "rna_depth": "120", "rna_alt_reads": "0"},
        {"event_id": "S4", "event_type": "SNV", "rna_depth": "30", "rna_alt_reads": "4"},
        {"event_id": "I1", "event_type": "InDel"},
        {"event_id": "F1", "event_type": "Fusion", "rna_junction_reads": "100"},
    ]
    bundle.peptides = [
        {"event_id": "S4", "event_type": "SNV", "rna_depth": "30", "rna_alt_reads": "4"},
    ]
    rows = _patient_coding_variant_rna_rows(bundle)
    by_track = {row["编码变异赛道"]: row for row in rows}
    assert by_track["SNV"]["DNA检出事件"] == "4"
    assert by_track["SNV"]["RNA ALT≥3"] == "1"
    assert by_track["SNV"]["RNA无覆盖"] == "1"
    assert by_track["SNV"]["RNA低覆盖且ALT=0"] == "1"
    assert by_track["SNV"]["RNA充分覆盖且ALT=0"] == "1"
    assert by_track["InDel"]["位点RNA未计算"] == "1"
    assert by_track["SNV+InDel合计"]["DNA检出事件"] == "5"
    assert all("Fusion" not in row["编码变异赛道"] for row in rows)


def test_patient_report_explains_coding_variant_rna_gap_without_fusion_callers(tmp_path):
    bundle = _bundle()
    bundle.events = [
        {"event_id": "S1", "event_type": "SNV", "rna_depth": "100", "rna_alt_reads": "0"},
        {"event_id": "I1", "event_type": "InDel"},
    ]
    out = tmp_path / "coding_rna_gap.html"
    make_patient_report(out, bundle)
    text = out.read_text(encoding="utf-8")
    assert "编码变异的DNA–RNA证据覆盖" in text
    assert "不混入Fusion/Splice" in text
    assert "不把增加融合caller当作编码变异RNA补证" in text
    assert "DNA有、RNA无直接支持" in text


def test_experiment_entry_gates_separate_presentation_caps_and_safety_hard_failures():
    bundle = _bundle()
    bundle.appm_summary = {"appm_evidence_completeness": "LOW"}
    bundle.peptides = [
        {
            "event_id": "E1", "peptide": "AAAAAAAAA", "hla_allele": "HLA-A*02:01", "event_type": "SNV",
            "presentation_consensus_state": "PRESENTATION_DISCORDANT", "safety_state": "SAFETY_PARTIAL",
            "normal_proteome_exact_match_status": "NOT_DETECTED",
        },
        {
            "event_id": "E2", "peptide": "BBBBBBBBB", "hla_allele": "HLA-A*02:01", "event_type": "Fusion",
            "presentation_consensus_state": "PRESENTATION_CONSISTENT_STRONG", "safety_state": "SAFETY_REJECT",
            "reference_proteome_exact_match": "true", "hard_failure_codes": "HARD_REFERENCE_PROTEOME_MATCH",
            "normal_junction_assessment_status": "UNASSESSED",
        },
    ]
    rows = {row["实验入口门控"]: row for row in _patient_experiment_entry_gate_rows(bundle)}
    assert rows["NetMHCpan/MHCflurry核心呈递共识"]["谨慎/封顶"] == "1"
    assert rows["NetMHCpan/MHCflurry核心呈递共识"]["支持进入"] == "1"
    assert rows["正常蛋白组精确匹配"]["明确阻断"] == "1"
    assert rows["HLA/APPM证据完整度"]["谨慎/封顶"] == "2"
    assert rows["Fusion精确断点级正常背景"]["未评估"] == "1"


def test_patient_report_states_event_truth_does_not_equal_first_batch_candidate(tmp_path):
    bundle = _bundle()
    bundle.peptides[0].update({
        "presentation_consensus_state": "PRESENTATION_DISCORDANT",
        "normal_proteome_exact_match_status": "NOT_DETECTED",
        "safety_state": "SAFETY_PARTIAL",
    })
    out = tmp_path / "experiment_entry_gate.html"
    make_patient_report(out, bundle)
    text = out.read_text(encoding="utf-8")
    assert "首批实验入口门控" in text
    assert "事件真实存在不等于相应肽段已经适合进入首批实验" in text
    assert "不表示候选已被确认可注射、有效或安全" in text
    assert "NetMHCpan/MHCflurry核心呈递共识" in text


def test_lohhla_qc_gap_is_not_attributed_to_missing_purity_tools():
    bundle = _bundle()
    bundle.hla_loh_tool_results = [
        {"hla_allele": "HLA-A*02:01", "source_tool": "SpecHLA", "status": "RETAINED", "spechla_copy": "2"},
        {"hla_allele": "HLA-A*02:01", "source_tool": "LOHHLA", "status": "UNASSESSED", "lohhla_pval": "NA"},
    ]
    rows, _ = _patient_hla_loh_consensus(bundle)
    assert "CopyNumLoc/纯度倍性参数" in rows[0]["说明"]
    assert "不应归因于Sequenza或ASCAT未运行" in rows[0]["说明"]


def test_patient_evidence_audit_counts_later_tool_field_after_unassessed_placeholder():
    rows = [{
        "event_type": "Splice",
        "netmhcstabpan_rank": "UNASSESSED",
        "netchop_processing_status": "ASSESSED",
        "mutant_specificity_status": "UNASSESSED",
        "mutant_specificity_state": "NOVEL_SEQUENCE",
    }]
    audit = {row["证据维度"]: row for row in _patient_evidence_audit_rows(rows)}
    assert audit["加工/稳定性"]["可作为当前分层证据"].startswith("1（")
    assert "NetChop 1" in audit["加工/稳定性"]["判定口径"]
    assert audit["免疫原性"]["可作为当前分层证据"].startswith("0（")
    assert audit["MT/WT特异性"]["可作为当前分层证据"] == "1（传统MT/WT 0；异常连接新序列 1）"
    assert "正常连接或正常异构体肽对照" in audit["MT/WT特异性"]["判定口径"]


def test_patient_evidence_audit_counts_raw_immunogenicity_tools():
    rows = [
        {"prime_score": "0.82"},
        {"bigmhc_im_score": "0.71", "deepimmuno_score": "0.63"},
        {"immunogenicity_score": "UNASSESSED", "iedb_immunogenicity_score": "0.12"},
        {},
    ]
    audit = {row["证据维度"]: row for row in _patient_evidence_audit_rows(rows)}
    immunogenicity = audit["免疫原性"]
    assert immunogenicity["可作为当前分层证据"] == "3（至少一项候选级结果）"
    assert immunogenicity["尚不能作为可靠证据"] == "1（尚无候选级免疫原性值）"
    assert "PRIME 1" in immunogenicity["判定口径"]
    assert "BigMHC 1" in immunogenicity["判定口径"]
    assert "DeepImmuno 1" in immunogenicity["判定口径"]
    assert "IEDB 1" in immunogenicity["判定口径"]


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


def test_rerank_recovers_purity_tools_from_recorded_recommendation(tmp_path):
    base = _bundle()
    purity_dir = tmp_path / "upstream" / "evidence" / "purity_cnv"
    purity_dir.mkdir(parents=True)
    recommendation = purity_dir / "recommended_purity.tsv"
    write_tsv(recommendation, [{
        "purity": "0.6932", "ploidy": "1.96", "evidence_tool": "multi_tool_median",
        "consensus_status": "CONCORDANT",
    }])
    write_tsv(purity_dir / "purity_cnv_tool_summary.tsv", [
        {"tool": "FACETS", "status": "FOUND", "purity": "0.6864", "ploidy": "1.9648", "parse_method": "table columns"},
        {"tool": "PURPLE", "status": "FOUND", "purity": "0.7000", "ploidy": "1.9600", "parse_method": "table columns"},
        {"tool": "Sequenza", "status": "MISSING", "purity": "", "ploidy": "", "notes": "no parsed result found"},
    ])
    write_tsv(purity_dir / "purity_cnv_consensus.tsv", [{
        "status": "CONCORDANT", "recommended_purity": "0.6932",
        "interpretation": "FACETS and PURPLE are concordant.",
    }])
    peptides = [dict(base.peptides[0],
        purity="0.6932",
        purity_consensus_status="CONCORDANT",
        purity_tool_values='{"FACETS":0.6864,"PURPLE":0.7}',
        purity_recommendation_file=str(recommendation),
    )]
    final_dir = tmp_path / "rerank" / "final"
    final_dir.mkdir(parents=True)
    bundle = load_report_bundle(
        profile=base.profile, events=base.events, peptides=peptides,
        outdir=final_dir,
        provenance={
            "purity_cnv_tools": [{
                "tool": "FACETS", "purity": "", "ploidy": "",
                "status": "ASSESSED", "note": "stale run manifest",
            }],
            "purity_cnv_consensus": {
                "recommended_purity": "", "recommended_ploidy": "",
                "selected_tool": "FACETS", "status": "SINGLE_TOOL_NO_CROSSCHECK",
            },
        },
    )
    assert [(row["tool"], row["purity"], row["ploidy"]) for row in bundle.purity_tools] == [
        ("FACETS", "0.6864", "1.9648"),
        ("PURPLE", "0.7000", "1.9600"),
        ("Sequenza", "", ""),
    ]
    assert bundle.purity_consensus["recommended_purity"] == "0.6932"
    assert bundle.purity_consensus["recommended_ploidy"] == "1.96"
    assert bundle.purity_consensus["status"] == "CONCORDANT"


def test_explicit_patient_inputs_override_provenance_inputs():
    base = _bundle()
    bundle = load_report_bundle(
        profile=base.profile, events=base.events, peptides=base.peptides,
        provenance={"input_files": {"somatic_vcf": "/old/old.vcf.gz"}},
        patient_inputs={"somatic_vcf": "/new/new.vcf.gz"},
    )
    assert bundle.provenance["input_files"]["somatic_vcf"] == "/new/new.vcf.gz"


def test_patient_report_derives_normal_depth_from_paired_vcf(tmp_path):
    raw_events = tmp_path / "raw_events.tsv"
    write_tsv(raw_events, [
        {"event_id": "E1", "chrom": "chr1", "pos": "101", "ref": "A", "alt": "G"},
        {"event_id": "E2", "chrom": "1", "pos": "202", "ref": "C", "alt": "T"},
    ])
    vcf = tmp_path / "somatic.vcf.gz"
    with gzip.open(vcf, "wt", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tcase_blood\tcase_tumor\n")
        handle.write("chr1\t101\t.\tA\tG\t.\tPASS\t.\tGT:AD:DP\t0/0:20,0:20\t0/1:12,8:20\n")
        handle.write("chr1\t202\t.\tC\tT\t.\tPASS\t.\tGT:AD:DP\t0/0:40,0:40\t0/1:20,20:40\n")
    base = _bundle()
    bundle = load_report_bundle(
        profile=base.profile,
        events=base.events,
        peptides=base.peptides,
        provenance={"input_files": {"raw_events": str(raw_events), "somatic_vcf": str(vcf)}},
        outdir=tmp_path,
    )
    assert bundle.provenance["normal_dna_depth"] == (
        "候选事件位点中位正常样本有效深度 40x（n=2；VCF case_blood DP）"
    )


def test_patient_report_finds_paired_vcf_through_upstream_provenance(tmp_path):
    production = tmp_path / "case" / "pipeline" / "production"
    expression = production / "rna" / "expression" / "gene_tpm.tsv"
    expression.parent.mkdir(parents=True)
    expression.write_text("gene_id\ttpm\nGENE1\t1\n", encoding="utf-8")
    raw_events = production / "candidates" / "events.tsv"
    write_tsv(raw_events, [{"event_id": "E1", "chrom": "chr1", "pos": "101", "ref": "A", "alt": "G"}])
    vcf = production / "variants" / "somatic.vcf.gz"
    vcf.parent.mkdir(parents=True)
    with gzip.open(vcf, "wt", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tcase_blood\tcase_tumor\n")
        handle.write("chr1\t101\t.\tA\tG\t.\tPASS\t.\tGT:AD:DP\t0/0:31,0:31\t0/1:20,11:31\n")
    upstream = production / "final_original" / "final"
    upstream.mkdir(parents=True)
    (upstream / "provenance.json").write_text(json.dumps({
        "input_files": {"raw_events": str(raw_events), "somatic_vcf": str(vcf)},
    }), encoding="utf-8")
    base = _bundle()
    bundle = load_report_bundle(
        profile=base.profile,
        events=base.events,
        peptides=base.peptides,
        provenance={"input_files": {"expression": str(expression), "raw_events": str(raw_events)}},
        outdir=tmp_path / "derived",
    )
    assert bundle.provenance["normal_dna_depth"] == (
        "候选事件位点中位正常样本有效深度 31x（n=1；VCF case_blood DP）"
    )


def test_find_normal_bam_excludes_rna_and_hla_derived_bams(tmp_path):
    production = tmp_path / "case" / "pipeline" / "production"
    expression = production / "rna" / "expression" / "gene_tpm.tsv"
    expression.parent.mkdir(parents=True)
    expression.write_text("gene_id\ttpm\n", encoding="utf-8")
    normal = tmp_path / "inputs" / "case_blood_wgs.align.bam"
    rna = tmp_path / "inputs" / "case_blood_rna.bam"
    hla = tmp_path / "inputs" / "case_blood_lohhla_region.bam"
    normal.parent.mkdir()
    normal.write_bytes(b"normal")
    rna.write_bytes(b"rna" * 100)
    hla.write_bytes(b"hla" * 100)
    upstream = production / "final_original"
    upstream.mkdir(parents=True)
    (upstream / "provenance.json").write_text(json.dumps({
        "normal_dna_bam": str(normal),
        "normal_rna_bam": str(rna),
        "normal_lohhla_bam": str(hla),
    }), encoding="utf-8")
    provenance = {"input_files": {"expression": str(expression)}}
    assert _find_bam_input(provenance, "normal", tmp_path / "derived") == str(normal)


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
    assert "<td>未提供</td><td>未提示LOH</td>" in text
    assert "SpecHLA未提示相应限制性HLA-I等位基因LOH" in text
    assert "当前仅有单工具支持，不足以确认该等位基因在肿瘤中完整保留" in text
    assert "仅SpecHLA报告未提示LOH，证据有限" in text
    assert (
        "限制性HLA HLA-A*02:01 仅有单工具LOH判断："
        "LOHHLA=未提供；SpecHLA=未提示LOH；不足以确认该等位基因在肿瘤中完整保留或丢失"
    ) in text
    assert "拟进入高成本实验前建议补充" in text
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
    assert _replace_gene_ids("ENSG00000232325", {}) == "AC093627.1"


def test_scoring_all_tool_results_counts_as_canonical(tmp_path):
    scoring = tmp_path / "scoring"
    scoring.mkdir()
    evidence = scoring / "all_tool_results.tsv"
    write_tsv(evidence, [{"peptide_id": "P1", "hla_allele": "HLA-A*02:01", "netmhcpan_el_rank": "0.2"}])
    digest = hashlib.sha256(evidence.read_bytes()).hexdigest()
    (scoring / "all_tool_results.manifest.json").write_text(
        json.dumps({"output": {"sha256": digest}, "canonical": True}), encoding="utf-8"
    )
    base = _bundle()
    bundle = load_report_bundle(profile=base.profile, events=base.events, peptides=base.peptides, outdir=tmp_path)
    assert bundle.evidence_source_status == "CANONICAL_ALL_TOOL_RESULTS"
    assert bundle.evidence_integrity["status"] == "PASS"


def test_dotted_hla_loh_sidecars_are_loaded_for_both_tools(tmp_path):
    loh_dir = tmp_path / "hla_loh_consensus"
    loh_dir.mkdir()
    write_tsv(loh_dir / "lohhla.hla_loh.tsv", [
        {"hla_allele": "HLA-A*02:01", "loh_status": "loh", "evidence_tool": "lohhla"},
    ])
    write_tsv(loh_dir / "spechla.hla_loh.tsv", [
        {"hla_allele": "HLA-A*02:01", "loh_status": "no", "evidence_tool": "spechla"},
    ])
    base = _bundle()
    base.peptides[0]["hla_allele"] = "HLA-A*02:01"
    bundle = load_report_bundle(
        profile=base.profile, events=base.events, peptides=base.peptides, outdir=tmp_path,
    )
    assert {row.get("_report_tool") for row in bundle.hla_loh_tool_results} == {"LOHHLA", "SpecHLA"}
    rows, overall = _patient_hla_loh_consensus(bundle)
    assert rows[0]["综合判断"] == "工具结果冲突，暂不判定"
    assert "HLA-I LOH工具结果冲突" in overall


def test_hla_loh_consensus_tsv_is_expanded_into_tool_rows(tmp_path):
    consensus_dir = tmp_path / "hla_loh_consensus"
    consensus_dir.mkdir()
    write_tsv(consensus_dir / "hla_loh_consensus.tsv", [
        {"hla_allele": "HLA-A*02:01", "lohhla_status": "loh", "spechla_status": "no"},
        {"hla_allele": "HLA-B*07:02", "lohhla_status": "no", "spechla_status": "no"},
    ])
    base = _bundle()
    base.peptides[0]["hla_allele"] = "HLA-A*02:01"
    bundle = load_report_bundle(profile=base.profile, events=base.events, peptides=base.peptides, outdir=tmp_path)
    assert any(row.get("_report_tool") == "LOHHLA" for row in bundle.hla_loh_tool_results)
    assert any(row.get("_report_tool") == "SpecHLA" for row in bundle.hla_loh_tool_results)
    out = tmp_path / "patient_loh.html"
    make_patient_report(out, bundle)
    text = out.read_text(encoding="utf-8")
    assert "未提供逐等位基因HLA LOH结果" not in text
    assert "检出LOH" in text or "未见LOH" in text or "冲突" in text


def test_release_metadata_reads_parallel_ranking_rules_version():
    from neoag.reports_dual import _patient_release_metadata

    bundle = _bundle()
    bundle.provenance["parallel_rankings"] = {"rules_version": "2.1"}
    bundle.evidence_integrity = {"actual_sha256": "abc1234567890"}
    meta = _patient_release_metadata(bundle)
    assert meta["rules_version"] == "2.1"
    assert meta["run_id"].startswith("S1-")


def test_patient_assessed_accepts_a_valid_metric_after_unassessed_status():
    from neoag.reports_dual import _patient_assessed

    row = {"processing_status": "UNASSESSED", "processing_score": "0.82"}
    assert _patient_assessed(row, "processing_status", "processing_score")
    assert not _patient_assessed(row, "processing_status", "missing_score")


def test_junction_safety_uses_exact_sequence_not_partner_expression():
    gap = _patient_safety_gap({
        "event_type": "Fusion",
        "safety_status": "SAFETY_REVIEW",
        "safety_reason": "critical_tissue_expression;normal_hspc_expression",
        "critical_tissue_name": "Brain",
        "critical_tissue_max_tpm": "201.5",
        "normal_hspc_tpm": "160.7",
        "reference_proteome_status": "UNASSESSED",
        "normal_ligandome_status": "UNASSESSED",
        "normal_junction_assessment_status": "UNASSESSED",
    })
    assert "精确新生序列" in gap
    assert "伙伴基因正常组织/HSPC表达仅作背景" in gap
    assert "201.5" not in gap
    assert "160.7" not in gap

    partial_gap = _patient_safety_gap({
        "event_type": "Fusion", "safety_status": "SAFETY_PARTIAL",
        "normal_tissue_max_tpm": "213.994", "normal_tissue_max_tissue": "Nerve_Tibial",
        "normal_hspc_tpm": "18.0",
        "reference_proteome_status": "ASSESSED",
        "normal_ligandome_status": "ASSESSED",
        "normal_junction_assessment_status": "ASSESSED",
    })
    assert "精确新生序列" in partial_gap
    assert "213.994" not in partial_gap
    assert "18.0" not in partial_gap


def test_fusion_boundary_report_shows_both_gene_sides_and_flags_alternative_orfs():
    row = {
        "event_id": "F1", "event_type": "Fusion", "gene": "EWSR1::WT1",
        "peptide": "SSYGQQSEK", "crosses_junction": "yes",
        "junction_position_in_peptide_1based": "4",
        "fusion_left_gene": "EWSR1", "fusion_right_gene": "WT1",
        "fusion_left_peptide": "SSYG", "fusion_right_peptide": "QQSEK",
        "transcript_id": "FT1", "orf_id": "ORF1",
    }
    bundle = _bundle()
    bundle.peptides = [
        row,
        {**row, "peptide": "SSYGQQSVK", "transcript_id": "FT2", "orf_id": "ORF2"},
    ]
    text = _patient_fusion_boundary_evidence(row, bundle)
    assert "EWSR1来源 SSYG｜融合连接点（肽内位置 4）｜WT1来源 QQSEK" in text
    assert "SSYGQQSEK/SSYGQQSVK" in text
    assert "必须按患者精确融合转录本、阅读框和ORF逐一归因" in text


def test_presentation_quantitative_row_exposes_raw_values_and_allele_uncertainty():
    result = _patient_presentation_quantitative_row({
        "peptide": "SLYNTVATL", "hla_allele": "HLA-A*02:01",
        "mutation_positions_in_peptide": "4", "netmhcpan_mt_rank_el": "0.25",
        "netmhcpan_wt_rank_el": "2.5", "netmhcpan_mt_rank_ba": "0.7",
        "netmhcpan_wt_rank_ba": "3.5", "netmhcpan_mt_ic50": "48",
        "netmhcpan_wt_ic50": "810", "mhcflurry_mt_affinity_percentile": "0.8",
        "mhcflurry_wt_affinity_percentile": "4.2", "mhcflurry_mt_presentation_score": "0.91",
        "mhcflurry_wt_presentation_score": "0.12", "netmhcstabpan_rank": "0.6",
        "netmhcstabpan_score": "1.9", "netmhcstabpan_wt_rank": "3.1",
        "netmhcstabpan_wt_score": "0.4", "wildtype_peptide": "SLYATVATL", "prime_score": "0.72",
        "bigmhc_im_score": "0.81", "deepimmuno_score": "0.66",
        "mutation_position_role": "PUTATIVE_TCR_FACING",
        "wt_self_reactivity_risk_status": "WT_LOW_PREDICTED_BINDING",
    }, 1)
    assert "EL rank MT=0.25%, WT=2.5%" in result["NetMHCpan原始值"]
    assert "IC50 MT=48 nM, WT=810 nM" in result["NetMHCpan原始值"]
    assert "EL WT/MT=10" in result["MT/WT定量比较"]
    assert "presentation score MT=0.91, WT=0.12" in result["MHCflurry原始值"]
    assert "NetMHCstabpan rank=0.6%" in result["稳定性"]
    assert "WT rank=3.1%" in result["稳定性"]
    assert "PRIME MT=0.72" in result["免疫原性辅助模型"]
    assert "训练覆盖/外推状态未记录" in result["HLA模型覆盖"]
    assert result["肽长/变异位置"] == "9 aa；位置=4"
    assert "推定TCR暴露位置" in result["突变位置结构解释"]
    assert "不能据此排除耐受或交叉反应" in result["WT自身反应/耐受风险"]


def test_presentation_quantitative_row_flags_strong_wt_binding_without_precomputed_status():
    result = _patient_presentation_quantitative_row({
        "peptide": "SLYNTVATL", "wildtype_peptide": "SLYATVATL",
        "hla_allele": "HLA-A*02:01", "mutation_positions_in_peptide": "4",
        "mutation_tcr_facing": "yes", "netmhcpan_mt_rank_el": "0.4",
        "netmhcpan_wt_rank_el": "0.8", "netmhcpan_wt_ic50": "42",
        "mutant_specificity_status": "MARGINAL_MT_ADVANTAGE",
    }, 1)
    assert "WT仍预测为强结合" in result["WT自身反应/耐受风险"]
    assert "序列位置推断" in result["突变位置结构解释"]


def test_fusion_stability_does_not_render_missing_wt_as_percentage():
    result = _patient_presentation_quantitative_row({
        "event_type": "Fusion", "peptide": "DLDLGEKPY", "hla_allele": "HLA-A*01:01",
        "netmhcstabpan_rank": "0.8", "netmhcstabpan_score": "1.2",
    }, 1)
    assert "WT rank=未提供%" not in result["稳定性"]
    assert "传统点突变式WT稳定性不适用" in result["稳定性"]


def test_fusion_safety_dimensions_separate_direct_evidence_from_partner_expression():
    text = _patient_safety_dimensions({
        "event_type": "Fusion",
        "normal_proteome_exact_match_status": "NOT_DETECTED",
        "normal_transcript_junction_match_status": "NOT_DETECTED",
        "normal_immunopeptidome_match_status": "NOT_DETECTED",
        "similar_peptide_cross_reactivity_status": "UNASSESSED",
        "normal_tissue_max_tpm": "201.559",
        "normal_tissue_max_tissue": "Brain Cerebellar Hemisphere",
        "critical_tissue_max_tpm": "201.559",
        "critical_tissue_name": "Brain Cerebellar Hemisphere",
        "normal_hspc_tpm": "160.700",
        "normal_hspc_unit": "nCPM",
        "final_safety_conclusion": "PARTIAL_DIRECT_SAFETY_EVIDENCE",
        "safety_evidence_completeness": "0.75",
    })
    assert "完整肽正常蛋白组精确匹配=未检出" in text
    assert "相似肽交叉反应风险=未评估" in text
    assert "Brain Cerebellar Hemisphere最高 201.5590 TPM（伙伴基因辅助背景，不代表跨连接肽存在）" in text
    assert "160.7000 nCPM（伙伴基因辅助背景）" in text
    assert "数据库风险筛查证据不完整，暂不能评价实验性脱靶或交叉反应风险（直接证据完整度 75%）" in text


def test_patient_risk_screen_wording_does_not_claim_clinical_safety():
    text = _patient_safety_dimensions({
        "event_type": "SNV",
        "normal_proteome_exact_match_status": "NOT_DETECTED",
        "normal_transcript_junction_match_status": "NOT_DETECTED",
        "normal_immunopeptidome_match_status": "NOT_DETECTED",
        "similar_peptide_cross_reactivity_status": "LOW_RISK",
        "final_safety_conclusion": "NO_DIRECT_SAFETY_SIGNAL_DETECTED",
        "safety_evidence_completeness": "1.0",
    })
    assert "数据库风险筛查结论=当前数据库筛查未发现明确的同序列排除证据" in text
    assert "尚未完成实验性脱靶和交叉反应验证" in text
    assert "未见直接安全信号" not in text
    assert "已证明安全" not in text


def test_generic_fusion_artifact_review_flags_same_gene_and_complex_labels():
    same_gene = _patient_fusion_artifact_review({
        "event_type": "Fusion", "gene": "PPP1R9B::PPP1R9B",
    })
    complex_label = _patient_fusion_artifact_review({
        "event_type": "Fusion", "gene": "UNCX::ENSG00000233082(11139),MICALL2(148059)",
    })
    assert "同基因连接" in same_gene
    assert "复杂/多伙伴标注" in complex_label


def test_hla_loh_qc_only_record_does_not_count_as_assessed_tool():
    bundle = _bundle()
    bundle.hla_loh_tool_results = [
        {
            "_report_tool": "LOHHLA", "hla_allele": "HLA-A*02:01",
            "loh_status": "UNASSESSED", "lohhla_copy_number": "1.8",
        },
        {
            "_report_tool": "SpecHLA", "hla_allele": "HLA-A*02:01",
            "loh_status": "no",
        },
    ]
    rows, overall = _patient_hla_loh_consensus(bundle)
    assert overall == (
        "SpecHLA未提示相应限制性HLA-I等位基因LOH；另一工具因QC或输入对接不足未形成有效判断，"
        "因此当前仅有单工具支持，不足以确认该等位基因在肿瘤中完整保留。"
    )
    assert "多工具一致" not in overall
    assert rows[0]["LOHHLA"] == "未形成判断（QC不足）"
    audit = {
        row["证据维度"]: row
        for row in _patient_evidence_audit_rows(bundle.peptides[:1], bundle)
    }["限制性HLA状态"]
    assert audit["可作为当前分层证据"] == "0（逐等位基因多工具一致）"
    assert audit["尚不能作为可靠证据"] == "1（其中单工具或QC不足 1）"
    assert "仅单工具阴性或另一工具QC不足不计为完整保留证据" in audit["判定口径"]


def test_depth_and_rna_qc_are_summarized_by_unique_event():
    base = _bundle()
    duplicate = dict(base.peptides[0])
    duplicate["peptide_id"] = "P1_ALT"
    duplicate["hla_allele"] = "HLA-B*07:02"
    second = dict(base.peptides[0])
    second.update({
        "peptide_id": "P2", "event_id": "E2", "tumor_depth": "80",
        "rna_depth": "0", "rna_alt_reads": "0",
    })
    base.peptides[0]["tumor_depth"] = "100"
    bundle = load_report_bundle(
        profile=base.profile,
        events=[],
        peptides=[base.peptides[0], duplicate, second],
        provenance={"rna_qc_status": "已载入表达；位点级RNA pileup未评估"},
    )
    assert "n=2" in bundle.provenance["tumor_dna_depth"]
    assert "SNV/InDel位点级RNA已评估 2/2 个独立事件" in bundle.provenance["rna_qc_status"]
    assert "ALT reads≥1/3/5：1/1/0" in bundle.provenance["rna_qc_status"]


def test_rna_qc_denominator_keeps_dna_events_without_rna_metrics():
    base = _bundle()
    assessed = dict(base.peptides[0])
    missing = dict(base.peptides[0])
    missing.update({"peptide_id": "P2", "event_id": "E2"})
    for field in ("rna_depth", "rna_alt_reads", "rna_vaf"):
        missing.pop(field, None)
    bundle = load_report_bundle(
        profile=base.profile,
        events=[],
        peptides=[assessed, missing],
        provenance={"rna_qc_status": "位点级RNA pileup未评估"},
    )
    assert "SNV/InDel位点级RNA已评估 1/2 个独立事件" in bundle.provenance["rna_qc_status"]
    assert "位点RNA未计算：1" in bundle.provenance["rna_qc_status"]


def test_disease_knowledge_prioritizes_display_without_changing_r_grade(tmp_path):
    bundle = _bundle()
    common = {
        "event_type": "Fusion", "source_chain_confidence_tier": "C2",
        "netmhcpan_el_rank": "0.5", "evidence_grade": "R3",
        "best_evidence_grade": "R3", "safety_status": "PASS",
    }
    bundle.events = [
        {**common, "event_id": "F1", "gene": "PPP1R9B::PPP1R9B"},
        {**common, "event_id": "F2", "gene": "EWSR1::WT1"},
    ]
    bundle.peptides = [
        {**common, "peptide_id": "P1", "event_id": "F1", "gene": "PPP1R9B::PPP1R9B", "peptide": "AAAAAAAAA", "hla_allele": "HLA-A*11:01"},
        {**common, "peptide_id": "P2", "event_id": "F2", "gene": "EWSR1::WT1", "peptide": "SSYGQQSEK", "hla_allele": "HLA-A*11:01"},
    ]
    bundle.disease_knowledge = {
        "anchors": [{
            "event": "EWSR1::WT1", "aliases": ["EWSR1--WT1"],
            "label": "DSRCT disease anchor",
            "molecular_significance": "EWSR1::WT1是DSRCT的关键分子特征和核心致病机制事件",
            "patient_interpretation": (
                "检测到EWSR1::WT1融合。该融合是DSRCT的关键分子特征；"
                "如本事件经独立方法确认，应结合病理及临床资料进行诊断整合。"
                "本报告本身不替代病理诊断。"
            ),
            "public_neoantigens": [{
                "sequence": "SSYGQQSEK", "hla_prefixes": ["HLA-A*11"],
                "evidence": "HLA-IP/MS",
            }],
        }],
    }
    out = tmp_path / "patient_anchor.html"
    make_patient_report(out, bundle, event_top_n=2, candidate_top_n=2)
    text = out.read_text(encoding="utf-8")
    section = text.split("当前展示2个去重候选组合", 1)[1].split("</table>", 1)[0]
    assert section.index("EWSR1::WT1") < section.index("PPP1R9B::PPP1R9B")
    assert "仅优先展示，不自动提升R等级" in text
    assert "核心分子发现：EWSR1::WT1" in text
    assert "如本事件经独立方法确认，应结合病理及临床资料进行诊断整合" in text
    assert "本报告本身不替代病理诊断" in text
    assert "结构化临床诊断" in text
    assert "分子知识库锚定" in text
    assert "<td>R3-" in section
    assert "<td>R1</td>" not in section
    assert "<td>R2</td>" not in section


def test_disease_knowledge_is_auto_discovered_from_structured_disease():
    base = _bundle()
    bundle = load_report_bundle(
        profile=base.profile,
        events=base.events,
        peptides=base.peptides,
        provenance={"clinical_context": {"disease": "促结缔组织增生性小圆细胞肿瘤（DSRCT）"}},
    )
    assert bundle.disease_knowledge["status"] == "LOADED"
    assert bundle.disease_knowledge["disease_id"] == "DSRCT"
    assert bundle.disease_knowledge["anchors"][0]["event"] == "EWSR1::WT1"
