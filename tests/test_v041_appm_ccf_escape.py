from pathlib import Path

from neoag_v03.appm_v2 import build_appm_2
from neoag_v03.ccf_v2 import build_ccf_2
from neoag_v03.immune_escape import build_immune_escape_evidence
from neoag_v03.utils import read_tsv, write_tsv
from neoag_v03.config import load_profile


def test_appm_2_b2m_biallelic_loss_caps_mhc_i_peptide(tmp_path):
    vep = tmp_path / "vep.tsv"
    expr = tmp_path / "expr.tsv"
    cnv = tmp_path / "cnv.tsv"
    peps = tmp_path / "peptides.tsv"
    write_tsv(vep, [{"gene": "B2M", "consequence": "frameshift_variant"}])
    write_tsv(expr, [{"gene": "B2M", "TPM": "0.05"}, {"gene": "HLA-A", "TPM": "8"}])
    write_tsv(cnv, [{"gene": "B2M", "copy_number_status": "copy_neutral_loh", "total_cn": "1", "minor_cn": "0"}])
    write_tsv(peps, [{"peptide_id": "P1", "event_id": "E1", "peptide": "AAAAAAAAA", "hla_allele": "HLA-A*02:01", "mhc_class": "I"}])
    paths = build_appm_2(sample_id="S1", outdir=tmp_path / "appm", vep_tsv=vep, expression_tsv=expr, cnv_tsv=cnv, raw_peptides=peps, profile=load_profile("default"))
    gene = {r["gene"]: r for r in read_tsv(paths["appm_gene_status"])}
    flags = read_tsv(paths["peptide_appm_flags"])
    assert gene["B2M"]["biallelic_status"] == "BIALLELIC_LOSS"
    assert flags[0]["priority_cap"] == "D"
    assert flags[0]["appm_multiplier"] == "0.0000"


def test_ccf_2_reports_multiplicity_range_and_confidence(tmp_path):
    events = tmp_path / "events.tsv"
    purity = tmp_path / "purity.tsv"
    cnv = tmp_path / "cnv.tsv"
    write_tsv(events, [{"event_id": "E1", "sample_id": "S1", "mutation_source": "SNV", "chrom": "chr1", "pos": "150", "tumor_vaf": "0.35", "tumor_depth": "100", "tumor_alt_count": "35"}])
    write_tsv(purity, [{"purity": "0.70", "ploidy": "2.4"}])
    write_tsv(cnv, [{"chrom": "chr1", "start": "1", "end": "1000", "total_cn": "3", "major_cn": "2", "minor_cn": "1"}])
    rows = build_ccf_2(events, purity, cnv, load_profile("default"), tmp_path / "ccf.tsv")
    assert rows[0]["ccf_method"] == "SNV_INDEL_COPY_NUMBER_AWARE"
    assert rows[0]["multiplicity_best"]
    assert rows[0]["ccf_min"] and rows[0]["ccf_max"]
    assert rows[0]["ccf_confidence"] in {"medium", "high", "low"}


def test_immune_escape_2_context_retains_lost_hla_for_immunomonitoring(tmp_path):
    peps = tmp_path / "peptides.tsv"
    hla = tmp_path / "hla_loh.tsv"
    out = tmp_path / "escape"
    write_tsv(peps, [{"peptide_id": "P1", "event_id": "E1", "peptide": "AAAAAAAAA", "hla_allele": "HLA-A*02:01", "mhc_class": "I"}])
    write_tsv(hla, [{"hla_allele": "HLA-A*02:01", "loh_status": "loh"}])
    paths = build_immune_escape_evidence(sample_id="S1", raw_peptides=peps, outdir=out, hla_loh_tsv=hla, therapy_context="immunomonitoring")
    flags = read_tsv(paths["peptide_escape_flags"])
    assert flags[0]["restricting_hla_lost"] == "yes"
    assert flags[0]["priority_cap"] == "C_CAUTION"
    assert flags[0]["escape_multiplier"] != "0.0000"


def test_appm_2_writes_native_sidecars_and_avoids_hla_loh_double_penalty(tmp_path):
    hla = tmp_path / "hla_loh.tsv"
    peps = tmp_path / "peptides.tsv"
    write_tsv(hla, [{"hla_allele": "HLA-A*02:01", "loh_status": "loh"}])
    write_tsv(peps, [{"peptide_id": "P1", "event_id": "E1", "peptide": "AAAAAAAAA", "hla_allele": "HLA-A*02:01", "mhc_class": "I"}])
    paths = build_appm_2(sample_id="S1", outdir=tmp_path / "appm", hla_loh_tsv=hla, raw_peptides=peps, profile=load_profile("default"))
    for key in ["appm_module_scores", "appm_evidence_completeness", "appm_conflicts", "appm_peptide_modifiers"]:
        assert key in paths
        assert Path(paths[key]).exists()
    modifiers = read_tsv(paths["appm_peptide_modifiers"])
    assert modifiers[0]["restricting_locus_loh"] == "yes"
    assert modifiers[0]["appm_multiplier"] == "1.0000"
    assert "review_in_immune_escape" in modifiers[0]["appm_multiplier_reason"]
    completeness = read_tsv(paths["appm_evidence_completeness"])[0]
    assert completeness["hla_loh_assessed"] == "yes"




def test_appm_2_expression_missing_gene_is_not_called_low_expression(tmp_path):
    expr = tmp_path / "expr.tsv"
    peps = tmp_path / "peptides.tsv"
    write_tsv(expr, [{"gene": "HLA-A", "TPM": "10"}])
    write_tsv(peps, [{"peptide_id": "P1", "event_id": "E1", "peptide": "AAAAAAAAA", "hla_allele": "HLA-A*02:01", "mhc_class": "I"}])
    paths = build_appm_2(sample_id="S1", outdir=tmp_path / "appm", expression_tsv=expr, raw_peptides=peps, profile=load_profile("default"))
    gene = {r["gene"]: r for r in read_tsv(paths["appm_gene_status"])}
    assert gene["B2M"]["expression_status"] == "missing_from_expression_matrix"
    assert gene["B2M"]["expression_input_status"] == "gene_missing"
    assert gene["B2M"]["biallelic_status"] == "NO_EVIDENCE"
    assert gene["B2M"]["functional_validation_status"] == "computational_proxy"
    summary = read_tsv(paths["appm_summary"])[0]
    assert summary["rna_input_status"] == "gene_level_tpm"
    assert summary["functional_validation_status"] == "computational_proxy"


def test_appm_2_segment_only_cnv_is_flagged_not_treated_as_gene_level(tmp_path):
    cnv = tmp_path / "segments.tsv"
    write_tsv(cnv, [{"chrom": "chr15", "start": "1", "end": "1000", "total_cn": "1", "minor_cn": "0"}])
    paths = build_appm_2(sample_id="S1", outdir=tmp_path / "appm", cnv_tsv=cnv, profile=load_profile("default"))
    completeness = read_tsv(paths["appm_evidence_completeness"])[0]
    assert completeness["cnv_input_status"] == "segment_level_unmapped"
    assert completeness["cnv_assessed"] == "no"
    assert "cnv_gene_mapping" in completeness["missing_evidence"]
    summary = read_tsv(paths["appm_summary"])[0]
    assert summary["cnv_input_status"] == "segment_level_unmapped"



def test_appm_2_gene_status_compatible_mechanism_schema(tmp_path):
    vep = tmp_path / "vep.tsv"
    expr = tmp_path / "expr.tsv"
    cnv = tmp_path / "cnv.tsv"
    write_tsv(vep, [{"gene": "B2M", "consequence": "frameshift_variant"}])
    write_tsv(expr, [{"gene": "B2M", "TPM": "0.05"}, {"gene": "HLA-A", "TPM": "10"}])
    write_tsv(cnv, [{"gene": "B2M", "copy_number_status": "copy_neutral_loh", "total_cn": "1", "minor_cn": "0"}])
    paths = build_appm_2(sample_id="S1", outdir=tmp_path / "appm", vep_tsv=vep, expression_tsv=expr, cnv_tsv=cnv, profile=load_profile("default"))
    rows = read_tsv(paths["appm_gene_status"])
    required = {
        "module", "mutation_status", "expression_percentile", "protein_status",
        "ligandome_support", "gene_integrity_status", "gene_integrity_score",
        "evidence_completeness", "risk_reason",
    }
    assert required <= set(rows[0])
    gene = {r["gene"]: r for r in rows}
    assert gene["B2M"]["module"] == "MHC-I antigen presentation"
    assert gene["B2M"]["mutation_status"] == "damaging_variant"
    assert gene["B2M"]["gene_integrity_status"] in {"biallelic_loss", "conflicting"}
    assert float(gene["B2M"]["gene_integrity_score"]) <= 0.5
    assert gene["B2M"]["protein_status"] == "not_assessed"
    assert gene["B2M"]["ligandome_support"] == "not_assessed"
    assert gene["B2M"]["evidence_completeness"] == "PARTIAL"
    assert gene["HLA-A"]["gene_integrity_status"] == "intact"
    assert gene["HLA-A"]["expression_percentile"]


def test_appm_2_writes_input_evidence_status_layer(tmp_path):
    vep = tmp_path / "vep.tsv"
    expr = tmp_path / "expr.tsv"
    cnv = tmp_path / "cnv.tsv"
    hla_loh = tmp_path / "hla_loh.tsv"
    purity = tmp_path / "purity.tsv"
    peps = tmp_path / "peptides.tsv"
    write_tsv(vep, [{"gene": "B2M", "consequence": "frameshift_variant"}])
    write_tsv(expr, [{"gene": "B2M", "TPM": "0.5"}, {"gene": "HLA-A", "TPM": "12"}])
    write_tsv(cnv, [{"gene": "B2M", "copy_number_status": "loss", "loh_status": "loh"}])
    write_tsv(hla_loh, [{"hla_allele": "HLA-A*02:01", "loh_status": "loss"}])
    write_tsv(purity, [{"purity": "0.71", "ploidy": "2.1"}])
    write_tsv(peps, [{"peptide_id": "P1", "event_id": "E1", "peptide": "AAAAAAAAA", "hla_allele": "HLA-A*02:01", "mhc_class": "I"}])

    paths = build_appm_2(
        sample_id="S1",
        outdir=tmp_path / "appm",
        vep_tsv=vep,
        expression_tsv=expr,
        cnv_tsv=cnv,
        hla_loh_tsv=hla_loh,
        tumor_purity_tsv=purity,
        raw_peptides=peps,
        profile=load_profile("default"),
    )
    assert "appm_input_status" in paths
    rows = {r["input_type"]: r for r in read_tsv(paths["appm_input_status"])}
    expected = {
        "mutation_vep", "cnv_loh", "rna_expression", "hla_typing", "hla_loh",
        "tumor_purity", "proteomics", "phosphoproteomics", "hla_ligandome",
    }
    assert set(rows) == expected
    assert rows["mutation_vep"]["input_status"] == "assessed"
    assert rows["cnv_loh"]["input_status"] == "assessed"
    assert rows["rna_expression"]["input_status"] == "assessed"
    assert rows["hla_typing"]["input_status"] == "assessed"
    assert rows["hla_typing"]["assay_scope"] == "derived_from_raw_peptides"
    assert rows["hla_loh"]["input_status"] == "assessed"
    assert rows["tumor_purity"]["input_status"] == "assessed"
    assert rows["proteomics"]["input_status"] == "not_provided"
    assert rows["phosphoproteomics"]["input_status"] == "not_provided"
    assert rows["hla_ligandome"]["input_status"] == "not_provided"


def test_appm_2_input_evidence_status_distinguishes_empty_and_failed_parse(tmp_path):
    expr = tmp_path / "empty_expr.tsv"
    expr.write_text("gene\tTPM\n", encoding="utf-8")
    missing_proteomics = tmp_path / "missing_proteomics.tsv"
    paths = build_appm_2(
        sample_id="S1",
        outdir=tmp_path / "appm",
        expression_tsv=expr,
        proteomics_tsv=missing_proteomics,
        profile=load_profile("default"),
    )
    rows = {r["input_type"]: r for r in read_tsv(paths["appm_input_status"])}
    assert rows["rna_expression"]["input_status"] == "provided_empty"
    assert rows["rna_expression"]["parse_status"] == "parsed_empty"
    assert rows["proteomics"]["input_status"] == "failed_parse"
    assert rows["proteomics"]["parse_status"] == "missing_file"
    assert rows["cnv_loh"]["input_status"] == "not_provided"


def test_appm_2_fixed_gene_sets_and_immune_context_annotation_only(tmp_path):
    expr = tmp_path / "expr.tsv"
    write_tsv(expr, [
        {"gene": "HLA-A", "TPM": "12"},
        {"gene": "B2M", "TPM": "9"},
        {"gene": "STAT2", "TPM": "8"},
        {"gene": "CXCL9", "TPM": "50"},
        {"gene": "GZMB", "TPM": "20"},
    ])
    paths = build_appm_2(sample_id="S1", outdir=tmp_path / "appm", expression_tsv=expr, profile=load_profile("default"))
    rows = {r["gene"]: r for r in read_tsv(paths["appm_gene_status"])}
    assert rows["HLA-A"]["gene_set"] == "mhc_i_core"
    assert rows["B2M"]["gene_set"] == "mhc_i_core;mhc_i_processing"
    assert rows["STAT2"]["gene_set"] == "ifng_response"
    assert rows["CXCL9"]["gene_set"] == "optional_immune_context"
    assert rows["CXCL9"]["appm_integrity_role"] == "context_annotation_only"
    assert rows["CXCL9"]["functional_status"] == "context_annotation"
    assert rows["CXCL9"]["gene_integrity_status"] == "not_assessed"
    assert rows["CXCL9"]["gene_integrity_score"] == "1.0000"
    modules = {r["module"]: r for r in read_tsv(paths["appm_module_scores"])}
    assert modules["MHC-I"]["status"] == "MHC_I_INTACT"
    assert modules["IFNG-JAK-STAT"]["status"] == "IFNG_RESPONSE_INTACT"
    context = {r["gene"]: r for r in read_tsv(paths["appm_immune_context"])}
    assert context["CXCL9"]["context_marker_class"] == "interferon_inflamed_context"
    assert context["GZMB"]["context_marker_class"] == "cytotoxic_t_cell_context"
    assert context["CXCL9"]["context_interpretation"] == "background_annotation_not_appm_integrity"

