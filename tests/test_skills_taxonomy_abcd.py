from __future__ import annotations

import csv
import json
from pathlib import Path

from neoag.skill_taxonomy.registry import SKILLS_BY_NAME, registry_dict
from neoag.skill_taxonomy.runner import run_skill, validate_skill_dirs


def test_abcd_registry_complete():
    reg = registry_dict()
    assert set(reg["categories"]) == {"A", "B", "C", "D"}
    assert "neoag-vcf" in SKILLS_BY_NAME
    assert "neoag-ranking" in SKILLS_BY_NAME
    assert "open-neo-run" in SKILLS_BY_NAME
    assert "neoag-release-qc" in SKILLS_BY_NAME
    assert len(reg["skills"]) >= 20


def test_skill_dirs_exist():
    res = validate_skill_dirs(Path.cwd())
    assert res["status"] == "PASS", res["missing"]


def test_peptide_csv_skill(tmp_path: Path):
    inp = tmp_path / "peptides.tsv"
    inp.write_text("peptide_id\tpeptide\thla_allele\tgene\tel_rank\nP1\tSYFPEITHI\tHLA-A*02:01\tGENE1\t0.2\n", encoding="utf-8")
    outdir = tmp_path / "out"
    res = run_skill("neoag-peptide-csv", {"peptide_csv": str(inp), "outdir": str(outdir)}, dry_run=False)
    assert res["status"] == "PASS"
    assert (outdir / "raw_peptides.tsv").exists()
    assert (outdir / "presentation_evidence.tsv").exists()


def test_vcf_skill_minimal(tmp_path: Path):
    vcf = tmp_path / "sample.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t100\t.\tA\tT\t.\tPASS\t.\n", encoding="utf-8")
    outdir = tmp_path / "vcf_out"
    res = run_skill("neoag-vcf", {"vcf": str(vcf), "outdir": str(outdir), "sample_id": "S1"}, dry_run=False)
    assert res["status"] == "PASS"
    rows = list(csv.DictReader((outdir / "raw_events.tsv").open(), delimiter="\t"))
    assert rows and rows[0]["event_type"] == "SNV"


def test_experiment_design_skill(tmp_path: Path):
    ranked = tmp_path / "ranked.tsv"
    ranked.write_text("peptide_id\tpeptide\thla_allele\tgene\tevent_type\tfinal_priority\tefficacy_score\nP1\tSYFPEITHI\tHLA-A*02:01\tG1\tSNV\tB\t0.6\nP2\tABCDEFGH\tHLA-A*02:01\tG2\tfusion\tC\t0.5\n", encoding="utf-8")
    outdir = tmp_path / "design"
    res = run_skill("neoag-experiment-design", {"ranked_peptides": str(ranked), "outdir": str(outdir)}, dry_run=False)
    assert res["status"] == "PASS"
    assert (outdir / "short_peptide_pool.tsv").exists()
    assert (outdir / "minigene_design.tsv").exists()


def test_skill2_ranking_delegates_to_production_evidence_rank(tmp_path: Path):
    comprehensive = tmp_path / "comprehensive_peptide_evidence.tsv"
    comprehensive.write_text(
        "peptide_id\tevent_id\tevent_type\tpeptide\thla_allele\t"
        "cross_platform_status\trna_alt_reads\trna_vaf\t"
        "netmhcpan_mt_rank_el\tmhcflurry_presentation_score\t"
        "mutant_specificity_status\tccf_estimate\tccf_confidence\t"
        "appm_multiplier\thla_loh_status\trestricting_hla_lost\t"
        "safety_status\tsafety_evidence_completeness\n"
        "P1\tE1\tSNV\tSYFPEITHI\tHLA-A*02:01\t"
        "CROSS_PLATFORM_PASS_CONCORDANT\t8\t0.20\t0.2\t0.8\t"
        "MT_SPECIFIC\t0.9\thigh\t1.0\tRETAINED\tfalse\tPASS\tCOMPLETE\n",
        encoding="utf-8",
    )
    weighted = tmp_path / "ranked_peptides.tsv"
    weighted.write_text("peptide_id\tefficacy_score\tfinal_priority\nP1\t0.8\tB\n", encoding="utf-8")
    for skill_name in ("neoag-ranking", "open-neo-run"):
        outdir = tmp_path / skill_name
        result = run_skill(skill_name, {
            "comprehensive_evidence": str(comprehensive),
            "weighted_baseline": str(weighted),
            "outdir": str(outdir),
        })
        assert result["status"] == "PASS"
        assert result["production_command"] == "neoag evidence-rank"
        assert result["algorithm_owner"] == "src/neoag/evidence_consensus.py"
        for filename in (
            "all_tool_results.tsv",
            "ranked_peptides.weighted_baseline.tsv",
            "ranked_peptides.evidence_consensus.tsv",
            "ranked_events.evidence_consensus.tsv",
            "ranking_compare_weighted_vs_consensus.md",
        ):
            assert (outdir / filename).is_file()


def test_experiment_design_prefers_event_representatives(tmp_path: Path):
    ranked_events = tmp_path / "ranked_events.evidence_consensus.tsv"
    ranked_events.write_text(
        "event_evidence_rank\tevent_group_id\tevent_id\tgene\tevent_type\tbest_evidence_grade\t"
        "representative_1_peptide_id\trepresentative_1_peptide\trepresentative_1_hla_allele\t"
        "representative_2_peptide_id\trepresentative_2_peptide\trepresentative_2_hla_allele\n"
        "1\tEVENT:E1\tE1\tTBR1\tSNV\tR2\tP1\tSYFPEITHI\tHLA-A*02:01\tP2\tYFPEITHIA\tHLA-B*07:02\n"
        "2\tEVENT:E2\tE2\tGENE2\tfusion\tR3\tP3\tABCDEFGHI\tHLA-A*02:01\t\t\t\n",
        encoding="utf-8",
    )
    ranked_peptides = tmp_path / "ranked_peptides.evidence_consensus.tsv"
    ranked_peptides.write_text(
        "peptide_id\tpeptide\thla_allele\tgene\tevent_type\n"
        "WRONG\tAAAAAAAAA\tHLA-A*01:01\tDUPLICATE\tSNV\n",
        encoding="utf-8",
    )
    outdir = tmp_path / "design_events"
    res = run_skill("neoag-experiment-design", {
        "ranked_events": str(ranked_events), "ranked_peptides": str(ranked_peptides),
        "top_n": 3, "outdir": str(outdir),
    }, dry_run=False)
    assert res["status"] == "PASS"
    assert res["input_source"] == "ranked_events"
    rows = list(csv.DictReader((outdir / "experiment_candidates.tsv").open(), delimiter="\t"))
    assert [row["peptide_id"] for row in rows] == ["P1", "P3", "P2"]
    assert all(row["peptide_id"] != "WRONG" for row in rows)
    assert rows[0]["event_group_id"] == "EVENT:E1"



def test_vcf_skill_gz_multialt_csq(tmp_path: Path):
    import gzip
    vcf = tmp_path / "sample.vcf.gz"
    text = "##fileformat=VCFv4.2\n##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|SYMBOL|Feature|HGVSp\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t100\t.\tA\tT,G\t.\tPASS\tCSQ=T|missense_variant|GENE1|TX1|p.K1N,G|frameshift_variant|GENE2|TX2|p.K2fs\n"
    with gzip.open(vcf, "wt", encoding="utf-8") as fh:
        fh.write(text)
    outdir = tmp_path / "vcf_gz_out"
    res = run_skill("neoag-vcf", {"vcf": str(vcf), "outdir": str(outdir), "sample_id": "S1"}, dry_run=False)
    assert res["status"] == "PASS"
    rows = list(csv.DictReader((outdir / "raw_events.tsv").open(), delimiter="\t"))
    assert len(rows) == 2
    assert rows[0]["gene"] == "GENE1"
    assert rows[1]["gene"] == "GENE2"
    assert rows[1]["event_type"] == "SNV"


def test_fusion_splice_and_sv_entry_skills(tmp_path: Path):
    fusion = tmp_path / "fusion.tsv"
    fusion.write_text("FusionName\tJunctionReadCount\tFrame\tpeptide\thla\nGENE1::GENE2\t12\tin-frame\tSYFPEITHI\tA0201\n", encoding="utf-8")
    normal = tmp_path / "normal.tsv"
    normal.write_text("gene5\tgene3\nGENE1\tGENE2\n", encoding="utf-8")
    fout = tmp_path / "fusion_out"
    fres = run_skill("neoag-fusion", {"fusion": str(fusion), "normal_readthrough_db": str(normal), "outdir": str(fout), "sample_id": "S1"}, dry_run=False)
    assert fres["status"] == "PASS"
    frows = list(csv.DictReader((fout / "fusion_events.tsv").open(), delimiter="\t"))
    assert frows[0]["normal_background_hit"] == "true"
    prows = list(csv.DictReader((fout / "raw_peptides.tsv").open(), delimiter="\t"))
    assert prows and prows[0]["hla_allele"] == "HLA-A*02:01"

    splice = tmp_path / "splice.tsv"
    splice.write_text("chrom\tstart\tend\tgene\tread_count\tjunction_peptide\nchr1\t10\t20\tGENE3\t5\tAAAAAAAAA\n", encoding="utf-8")
    sout = tmp_path / "splice_out"
    sres = run_skill("neoag-splice", {"junctions": str(splice), "outdir": str(sout)}, dry_run=False)
    assert sres["status"] == "PASS"
    srows = list(csv.DictReader((sout / "splice_events.tsv").open(), delimiter="\t"))
    assert srows[0]["event_type"] == "Splice"

    sv = tmp_path / "sv.vcf"
    sv.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t100\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL\n2\t200\t.\tN\tN[3:300[\t.\tPASS\t.\n", encoding="utf-8")
    wgs_out = tmp_path / "sv_wgs"
    wgs = run_skill("neoag-sv-wgs", {"sv_vcf": str(sv), "outdir": str(wgs_out)}, dry_run=False)
    assert wgs["status"] == "PASS"
    svrows = list(csv.DictReader((wgs_out / "sv_events.tsv").open(), delimiter="\t"))
    assert {r["event_subtype"] for r in svrows} == {"DEL", "BND"}

    wes_fail = run_skill("neoag-sv-wes", {"sv_vcf": str(sv), "outdir": str(tmp_path / "sv_wes_fail"), "capture_bed": str(tmp_path / "missing.bed")}, dry_run=False)
    assert wes_fail["status"] == "FAIL"
    assert wes_fail["failure_reason"] == "CAPTURE_BED_NOT_FOUND"
    bed = tmp_path / "capture.bed"
    bed.write_text("chr1\t50\t150\n", encoding="utf-8")
    wes_out = tmp_path / "sv_wes"
    wes = run_skill("neoag-sv-wes", {"sv_vcf": str(sv), "capture_bed": str(bed), "outdir": str(wes_out)}, dry_run=False)
    assert wes["status"] == "PASS"
    wes_rows = list(csv.DictReader((wes_out / "sv_events.tsv").open(), delimiter="\t"))
    assert wes_rows[0]["capture_limited"] == "true"
    assert wes_rows[0]["capture_overlap"] == "true"
    assert wes_rows[1]["capture_overlap"] == "false"
