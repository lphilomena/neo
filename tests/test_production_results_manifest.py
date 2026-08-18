from __future__ import annotations

import subprocess
import sys
import tomllib
from pathlib import Path

from neoag.utils import read_tsv


ROOT = Path(__file__).resolve().parents[1]


def write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def test_generator_builds_all_three_upstream_consensus_stages(tmp_path):
    optitype = write(tmp_path / "optitype/result.tsv", "A1\tA2\tB1\tB2\tC1\tC2\nA*02:01\tA*11:01\tB*15:01\tB*40:01\tC*03:04\tC*08:01\n").parent
    spechla = write(tmp_path / "spechla/hla.result.txt", "A\tA*02:01\tA*11:01\nB\tB*15:01\tB*40:01\nC\tC*03:04\tC*08:01\n").parent
    facets = write(tmp_path / "facets/facets_purity.txt", "purity\t0.60\n").parent
    sequenza = write(tmp_path / "sequenza/result.tsv", "purity\tploidy\n0.62\t2.1\n").parent
    purple = write(tmp_path / "purple/purity.tsv", "purity\tploidy\n0.61\t2.0\n").parent
    lohhla = write(tmp_path / "lohhla/hla_loh.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    spechla_loh = write(tmp_path / "spechla_loh/hla_loh.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    star = write(tmp_path / "star/star-fusion.fusion_predictions.tsv", "#FusionName\tJunctionReadCount\nEWSR1--WT1\t20\n")
    output = tmp_path / "production.toml"
    command = [sys.executable, str(ROOT / "scripts/generate_production_from_results_manifest.py"), "--project-root", str(ROOT), "--sample-id", "S1", "--outdir", str(tmp_path / "run"), "--output", str(output), "--optitype", str(optitype), "--spechla-typing", str(spechla), "--facets", str(facets), "--sequenza", str(sequenza), "--purple", str(purple), "--lohhla", str(lohhla), "--spechla-loh", str(spechla_loh), "--star-fusion", str(star)]
    subprocess.run(command, check=True)
    manifest = tomllib.loads(output.read_text(encoding="utf-8"))
    assert manifest["run"]["hla_file"].endswith("evidence/hla_typing/recommended_hla.txt")
    assert {"hla_typing_consensus", "purity_cnv_consensus", "hla_loh_consensus", "fusion_candidates"} <= set(manifest["stages"])
    assert manifest["evidence"]["purity"].endswith("evidence/purity_cnv/recommended_purity.tsv")
    assert manifest["evidence"]["hla_loh"].endswith("evidence/hla_loh/hla_loh_consensus.tsv")
    assert manifest["run"]["profile"].endswith("profiles/sarcoma_rna_supported_v2_provisional.toml")
    assert manifest["evidence"]["evidence_consensus_rules"].endswith(
        "configs/ranking/sarcoma_evidence_consensus_v3_source_chain.toml"
    )


def test_generator_filters_splice_candidates_before_production_scoring(tmp_path):
    hla = write(tmp_path / "hla.txt", "HLA-A*02:01\n")
    facets = write(tmp_path / "facets/facets_purity.txt", "purity\t0.60\n").parent
    ascat = write(tmp_path / "ascat/ascat_summary.tsv", "sample_id\tpurity\tploidy\nS1\t0.64\t2.1\n").parent
    lohhla = write(tmp_path / "lohhla/hla_loh.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    spechla_loh = write(tmp_path / "spechla_loh/hla_loh.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    junctions = write(tmp_path / "junctions.tsv", "chrom\tstart\tend\tstrand\tjunction_reads\nchr1\t10\t20\t+\t12\n")
    snaf = write(tmp_path / "snaf.tsv", "event_id\tpeptide\thla_allele\tbinding_rank\nE1\tACDEFGHIK\tHLA-A*02:01\t0.8\n")
    output = tmp_path / "production.toml"
    subprocess.run([sys.executable, str(ROOT / "scripts/generate_production_from_results_manifest.py"), "--project-root", str(ROOT), "--sample-id", "S1", "--outdir", str(tmp_path / "run"), "--output", str(output), "--hla-file", str(hla), "--facets", str(facets), "--ascat", str(ascat), "--lohhla", str(lohhla), "--spechla-loh", str(spechla_loh), "--junctions", str(junctions), "--snaf", str(snaf)], check=True)
    manifest = tomllib.loads(output.read_text(encoding="utf-8"))
    splice = manifest["stages"]["splice_candidates"]
    assert "filter_splice_production_candidates.py" in splice["command"]
    assert splice["outputs"]["raw_peptides"].endswith("production_selected/raw_peptides.tsv")


def test_generator_accepts_facets_and_ascat_when_other_purity_tools_failed(tmp_path):
    hla = write(tmp_path / "hla.txt", "HLA-A*02:01\n")
    fasta = write(tmp_path / "GRCh38.fa", ">chr1\nACGT\n")
    facets = write(tmp_path / "facets/facets_purity.txt", "purity\t0.60\n").parent
    ascat = write(tmp_path / "ascat/ascat_summary.tsv", "sample_id\tpurity\tploidy\nS1\t0.64\t2.1\n").parent
    lohhla = write(tmp_path / "lohhla/hla_loh.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    spechla_loh = write(tmp_path / "spechla_loh/hla_loh.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    star = write(tmp_path / "star/star-fusion.fusion_predictions.tsv", "#FusionName\tJunctionReadCount\nEWSR1--WT1\t20\n")
    output = tmp_path / "production.toml"
    command = [sys.executable, str(ROOT / "scripts/generate_production_from_results_manifest.py"), "--project-root", str(ROOT), "--sample-id", "S1", "--outdir", str(tmp_path / "run"), "--output", str(output), "--hla-file", str(hla), "--reference-fasta", str(fasta), "--facets", str(facets), "--ascat", str(ascat), "--lohhla", str(lohhla), "--spechla-loh", str(spechla_loh), "--somatic-vcf", str(star), "--star-fusion", str(star)]
    subprocess.run(command, check=True)
    manifest = tomllib.loads(output.read_text(encoding="utf-8"))
    assert manifest["run"]["required_tool_groups"]["purity_cnv"]["tools"] == ["facets", "ascat"]
    assert "purity_facets" in manifest["stages"]
    assert "purity_ascat" in manifest["stages"]
    assert "purity_sequenza" not in manifest["stages"]
    assert "purity_purple" not in manifest["stages"]
    assert manifest["run"]["profile"].endswith("profiles/sarcoma_rna_supported_v2_provisional.toml")
    assert manifest["evidence"]["evidence_consensus_rules"].endswith(
        "configs/ranking/sarcoma_evidence_consensus_v3_source_chain.toml"
    )
    assert f'NEOAG_REFERENCE_FASTA="{fasta}"' in manifest["stages"]["snv_indel_candidates"]["command"]


def test_generator_builds_star_and_rna_allele_count_stages_from_fastq(tmp_path):
    hla = write(tmp_path / "hla.txt", "HLA-A*02:01\n")
    purity = write(tmp_path / "purity.tsv", "sample_id\tpurity\tploidy\nS1\t0.60\t2.0\n")
    cnv = write(tmp_path / "cnv.tsv", "chrom\tstart\tend\ttotal_cn\nchr1\t1\t1000\t2\n")
    lohhla = write(tmp_path / "lohhla.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    spechla_loh = write(tmp_path / "spechla_loh.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    vcf = write(tmp_path / "somatic.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\tT\t.\tPASS\t.\n")
    fastq1 = write(tmp_path / "rna_R1.fastq", "@r1\nACGT\n+\nFFFF\n")
    fastq2 = write(tmp_path / "rna_R2.fastq", "@r1\nTGCA\n+\nFFFF\n")
    star_index = tmp_path / "star_index"
    star_index.mkdir()
    gtf = write(tmp_path / "gencode.gtf", 'chr1\ttest\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG1";\n')
    output = tmp_path / "production.toml"
    subprocess.run([
        sys.executable, str(ROOT / "scripts/generate_production_from_results_manifest.py"),
        "--project-root", str(ROOT), "--sample-id", "S1",
        "--outdir", str(tmp_path / "run"), "--output", str(output),
        "--hla-file", str(hla), "--purity", str(purity), "--cnv", str(cnv),
        "--lohhla", str(lohhla), "--spechla-loh", str(spechla_loh),
        "--somatic-vcf", str(vcf), "--rna-fastq1", str(fastq1),
        "--rna-fastq2", str(fastq2), "--star-index", str(star_index),
        "--gencode-gtf", str(gtf), "--rna-threads", "8",
    ], check=True)
    manifest = tomllib.loads(output.read_text(encoding="utf-8"))
    assert "rna_star_alignment" in manifest["stages"]
    assert "run_star_rna_fastq.sh" in manifest["stages"]["rna_star_alignment"]["command"]
    assert manifest["stages"]["rna_alt_vaf"]["depends_on"] == ["rna_star_alignment"]
    assert "rna_allele_counts_pysam.py" in manifest["stages"]["rna_alt_vaf"]["command"]
    assert manifest["evidence"]["rna_vaf"] == "{outdir}/rna/rna_alt_vaf.tsv"


def test_generator_reuses_existing_rna_vaf_without_star(tmp_path):
    hla = write(tmp_path / "hla.txt", "HLA-A*02:01\n")
    purity = write(tmp_path / "purity.tsv", "sample_id\tpurity\tploidy\nS1\t0.60\t2.0\n")
    cnv = write(tmp_path / "cnv.tsv", "chrom\tstart\tend\ttotal_cn\nchr1\t1\t1000\t2\n")
    lohhla = write(tmp_path / "lohhla.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    spechla_loh = write(tmp_path / "spechla_loh.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    vcf = write(tmp_path / "somatic.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    rna_vaf = write(tmp_path / "rna_alt_vaf.tsv", "chrom\tpos\tref\talt\trna_depth\trna_alt_reads\trna_vaf\nchr1\t10\tA\tT\t20\t5\t0.25\n")
    output = tmp_path / "production.toml"
    subprocess.run([
        sys.executable, str(ROOT / "scripts/generate_production_from_results_manifest.py"),
        "--project-root", str(ROOT), "--sample-id", "S1",
        "--outdir", str(tmp_path / "run"), "--output", str(output),
        "--hla-file", str(hla), "--purity", str(purity), "--cnv", str(cnv),
        "--lohhla", str(lohhla), "--spechla-loh", str(spechla_loh),
        "--somatic-vcf", str(vcf), "--rna-vaf", str(rna_vaf),
    ], check=True)
    manifest = tomllib.loads(output.read_text(encoding="utf-8"))
    assert "rna_alt_vaf_input" in manifest["stages"]
    assert "rna_star_alignment" not in manifest["stages"]
    assert manifest["evidence"]["rna_vaf"] == str(rna_vaf.resolve())


def test_fusion_union_keeps_event_only_callers_and_provided_peptides(tmp_path):
    hla = write(tmp_path / "hla.txt", "HLA-A*02:01\n")
    star = write(tmp_path / "star.tsv", "#FusionName\tLeftBreakpoint\tRightBreakpoint\tJunctionReadCount\nEWSR1--WT1\tchr22:1:+\tchr11:2:-\t20\n")
    arriba = write(tmp_path / "arriba.tsv", "gene1\tgene2\tbreakpoint1\tbreakpoint2\tjunction_peptide\tsplit_reads1\nCPSF6\tRARG\tchr12:3:+\tchr12:4:-\tQYYSTPWTF\t12\n")
    outdir = tmp_path / "union"
    subprocess.run([sys.executable, str(ROOT / "scripts/build_fusion_caller_union.py"), "--sample-id", "S1", "--hla-file", str(hla), "--star-fusion", str(star), "--arriba", str(arriba), "--outdir", str(outdir)], cwd=ROOT, check=True)
    events = read_tsv(outdir / "raw_events.tsv")
    peptides = read_tsv(outdir / "raw_peptides.tsv")
    audit = read_tsv(outdir / "fusion_caller_union.tsv")
    assert {row["gene"] for row in events} == {"EWSR1::WT1", "CPSF6::RARG"}
    assert any(row["peptide"] == "QYYSTPWTF" and row["hla_allele"] == "HLA-A*02:01" for row in peptides)
    assert any(row["gene_pair"] == "EWSR1::WT1" and row["peptide_status"] == "ORF_PEPTIDE_UNAVAILABLE_REVIEW_ONLY" for row in audit)


def test_existing_upstream_results_produce_hla_purity_and_loh_consensus(tmp_path):
    optitype = write(tmp_path / "optitype/result.tsv", "A1\tA2\tB1\tB2\tC1\tC2\nA*02:01\tA*11:01\tB*15:01\tB*40:01\tC*03:04\tC*08:01\n").parent
    spechla = write(tmp_path / "spechla/hla.result.txt", "A\tA*02:01\tA*11:01\nB\tB*15:01\tB*40:01\nC\tC*03:04\tC*08:01\n").parent
    hla_out = tmp_path / "hla_consensus"
    subprocess.run([sys.executable, "-m", "neoag.agent_skills.hla_typing_compare", "--result-dir", str(optitype), "--result-dir", str(spechla), "--sample-id", "S1", "--outdir", str(hla_out)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    subprocess.run([sys.executable, str(ROOT / "scripts/hla_consensus_to_file.py"), "--consensus", str(hla_out / "hla_typing_consensus.tsv"), "--output", str(hla_out / "recommended_hla.txt")], cwd=ROOT, check=True)
    assert "HLA-A*02:01" in (hla_out / "recommended_hla.txt").read_text(encoding="utf-8")

    facets = write(tmp_path / "facets/facets_purity.txt", "purity\t0.60\n").parent
    sequenza = write(tmp_path / "sequenza/result.tsv", "purity\tploidy\n0.62\t2.1\n").parent
    purple = write(tmp_path / "purple/purity.tsv", "purity\tploidy\n0.61\t2.0\n").parent
    purity_out = tmp_path / "purity_consensus"
    subprocess.run([sys.executable, "-m", "neoag.agent_skills.purity_cnv_review", "--result-dir", str(facets), "--result-dir", str(sequenza), "--result-dir", str(purple), "--sample-id", "S1", "--outdir", str(purity_out)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    purity_rows = read_tsv(purity_out / "purity_cnv_consensus.tsv")
    assert purity_rows and int(purity_rows[0]["n_tools"]) >= 2

    lohhla = write(tmp_path / "lohhla.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    spechla_loh = write(tmp_path / "spechla_loh.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    loh_out = tmp_path / "loh_consensus"
    subprocess.run([sys.executable, str(ROOT / "scripts/build_hla_loh_consensus.py"), "--sample-id", "S1", "--lohhla", str(lohhla), "--spechla", str(spechla_loh), "--outdir", str(loh_out)], cwd=ROOT, check=True)
    assert read_tsv(loh_out / "hla_loh_consensus.tsv")[0]["consensus_status"] == "CONSENSUS_RETAINED"
