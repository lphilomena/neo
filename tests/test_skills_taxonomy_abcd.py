from __future__ import annotations

import csv
import json
from pathlib import Path

from neoag_v03.skill_taxonomy.registry import SKILLS_BY_NAME, registry_dict
from neoag_v03.skill_taxonomy.runner import run_skill, validate_skill_dirs


def test_abcd_registry_complete():
    reg = registry_dict()
    assert set(reg["categories"]) == {"A", "B", "C", "D"}
    assert "neoag-vcf" in SKILLS_BY_NAME
    assert "neoag-ranking" in SKILLS_BY_NAME
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
